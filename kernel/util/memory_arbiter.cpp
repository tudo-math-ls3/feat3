#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/os_windows.hpp>

#include <cstdint>
#include <mutex>
#include <unordered_set>

namespace FEAT
{
  namespace Memory
  {
    /// \cond internal
    /**
     * \brief main memory manager mutex object
     *
     * This mutex is used by the Manager class whenever a new block is created during a call of
     * Manager::create() or an existing one is destroyed during a call of Manager::release(),
     * because this may involve a resize of the block deque and a change of the free handle deque.
     *
     * This mutex is not locked by any action of the Block class, as that class used the arbiter_mutex
     * array instead.
     */
    static std::mutex manager_mutex;

    static constexpr std::size_t arbiter_mutex_shift = 6;

    /// arbiter mutex index mask
    static constexpr std::size_t arbiter_mutex_mask = 0x7;

    /**
     * \brief memory block arbiter object arrays
     */
    static std::mutex arbiter_mutex[arbiter_mutex_mask + 1u];

    typedef std::unique_lock<std::mutex> Lock;

    /**
     * \brief Set of all currently allocated memory blocks
     */
    static std::unordered_set<Block*> manager_block_set;

    /// \endcond

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * \brief Memory Arbiter memory block class
     *
     * This class is used internally by the Memory::Arbiter class to store all relevant data that is managed by the
     * arbiter. This class is only defined internally and cannot be used directly by anyone else than the arbiter.
     *
     * \author Peter Zajac
     */
    class Block
    {
    public:
      /// full size
      static constexpr std::size_t full_size = ~std::size_t(0);

      /**
       * \brief total number of references to this memory block
       */
      std::size_t ref_count;

      /**
       * \brief total size of memory block in bytes
       */
      std::size_t bytes;

      /**
       * \brief memory block pointer to main (host) memory
       */
      void* mem_main;

      /**
       * \brief memory block pointer to cuda device memory
       */
      void* mem_cuda;

      /**
       * \brief combination of all currently up-to-date memory locations
       *
       * This enum holds the combination of all memory locations which are currently up-to-date.
       */
      Location valid_locations;

      /**
       * \brief location of current write view
       *
       * This holds the location of the currently open write view(s).
       */
      Location write_location;

      /**
       * \brief number of currently active write views
       */
      std::size_t write_views;

      /**
      * \brief number of currently active main memory views with read access
      */
      std::size_t read_views_main;

      /**
      * \brief number of currently active cuda memory views with read access
      */
      std::size_t read_views_cuda;

      /**
       * \brief specifies whether the write view(s) can be overlapped
       */
      bool write_overlap;

    public:
      Block() :
        ref_count(0u),
        bytes(0u),
        mem_main(nullptr),
        mem_cuda(nullptr),
        valid_locations(Location::none),
        write_location(Location::none),
        write_views(0u),
        read_views_main(0u),
        read_views_cuda(0u),
        write_overlap(false)
      {
      }

      ~Block()
      {
        XASSERT(ref_count == 0u);
        XASSERT(mem_main == nullptr);
        XASSERT(mem_cuda == nullptr);
      }

      // no move constructor
      Block(Block&& other) = delete;
      // no move-assign operator
      Block& operator=(Block&& other) = delete;
      /// no copy constructor
      Block(const Block&) = delete;
      /// no copy-assign operator
      Block& operator=(const Block&) = delete;

      /// locks the arbiter mutex assigned to this memory block
      Lock lock_mutex() const
      {
        return Lock(arbiter_mutex[(reinterpret_cast<std::size_t>(this) >> arbiter_mutex_shift) & arbiter_mutex_mask]);
      }

      /// checks wether another block is assigned to the same mutex as this block
      bool uses_same_lock(const Block& other) const
      {
        return
          ((reinterpret_cast<std::size_t>(this)   >> arbiter_mutex_shift) & arbiter_mutex_mask) ==
          ((reinterpret_cast<std::size_t>(&other) >> arbiter_mutex_shift) & arbiter_mutex_mask);
      }

      void create(std::size_t bytes_, Location loc_, Init all_)
      {
        std::unique_lock<std::mutex> lock(lock_mutex());

        // make sure that this block is empty
        if(this->bytes > 0u)
          throw MemoryException("Memory::Block::create: memory block is already allocated");

        // save byte count and try to allocate
        this->bytes = bytes_;
        if(loc_ != Location::none)
          this->_alloc(loc_, all_);
        ++this->ref_count;

        // finally, register this block with the manager
        {
          Lock manager_lock(manager_mutex);
          manager_block_set.insert(this);
        }
      }

      static Block* attach(Block* block)
      {
        if(block == nullptr)
          return nullptr;

        Lock lock(block->lock_mutex());
        //XASSERTM(block != nullptr, "cannot attach to null block");
        ++block->ref_count;
        return block;
      }

      static bool detach(Block* block)
      {
        if(block == nullptr)
          return false;

        Lock lock(block->lock_mutex());

        //XASSERTM(block != nullptr, "cannot detach from null block");
        if(block->ref_count == 0u)
          throw MemoryException("Memory::Block::detach: reference count is already zero");

        if(0u < --block->ref_count)
          return true;

        // unregister this block with the manager
        {
          Lock manager_lock(manager_mutex);
          manager_block_set.erase(block);
        }

        // all references have been detached, so let's free the block
        block->_free_main();
        block->_free_cuda();

        // finally, free the block itself
        delete block;

        return false;
      }

      void alloc(Location mem_loc, Init mem_ini)
      {
        std::unique_lock<std::mutex> lock(lock_mutex());
        _alloc(mem_loc, mem_ini);
      }

      void free(Location mem_loc)
      {
        std::unique_lock<std::mutex> lock(lock_mutex());
        if(*(mem_loc & Location::main))
          _free_main();
        if(*(mem_loc & Location::cuda))
          _free_cuda();
      }

      void free_unmapped(Location mem_loc)
      {
        std::unique_lock<std::mutex> lock(lock_mutex());

        // ensure that no view is currently mapped, otherwise scram
        if(_is_mapped_any())
          return;

        // ensure that a valid copy remains in main memory
        _copy_main_from_uptodate();

        // free requested memory locations (but never free main memory)
        if(*(mem_loc & Location::cuda))
          _free_cuda();
      }

      void copy_from(const Block& source, Location mem_loc, std::size_t this_offset, std::size_t src_offset, std::size_t size)
      {
        // first lock the mutex for this memory block
        std::unique_lock<std::mutex> lock(lock_mutex());

        // make sure that this block is not mapped
        if(this->_is_mapped_any())
          throw AccessException("Memory::Block::copy_from: cannot copy to mapped memory block");

        // now lock the mutex for the source block, but only if it actually uses another mutex than ours
        std::unique_lock<std::mutex> source_lock;
        if(!this->uses_same_lock(source))
          source_lock = source.lock_mutex();

        // ensure that the source block is not mapped for writing
        if(source._is_mapped_any_write())
          throw AccessException("Memory::Block::copy_from: cannot copy from write-mapped memory block");

        // if no target location is specified, then simply copy the existing source location
        if(mem_loc == Location::none)
          mem_loc = source.valid_locations;

        // ensure the ranges are valid
        if(size == full_size)
          throw RangeException("Memory::Block::copy_from: memory size must be given explicitly");
        if(this->bytes < this_offset + size)
          throw RangeException("Memory::Block::copy_from: invalid destination memory range");
        if(source.bytes < src_offset + size)
          throw RangeException("Memory::Block::copy_from: invalid source memory range");

        // perform the actual copy
        if(size > std::size_t(0))
        {
          this->_copy_from(source, mem_loc, this_offset, src_offset, size);
        }

        // store our up-to-date locations
        this->valid_locations = mem_loc;
      }

      void copy_from(const void* src_ptr, Location src_loc, Location mem_loc, std::size_t this_offset, std::size_t size)
      {
        // first lock the mutex for this memory block
        std::unique_lock<std::mutex> lock(lock_mutex());

        // make sure that this block is not mapped
        if(this->_is_mapped_any())
          throw AccessException("Memory::Block::copy_from: cannot copy to mapped memory block");

        // if no target location is specified, then simply copy the existing source location
        if(mem_loc == Location::none)
          mem_loc = src_loc;

        // ensure the ranges are valid
        if(this->bytes < this_offset)
          throw RangeException("Memory::Block::copy_from: invalid destination memory offset");
        if(size == full_size)
          size = this->bytes - this_offset;
        if(this->bytes < this_offset + size)
          throw RangeException("Memory::Block::copy_from: invalid destination memory range");

        // perform the actual copy
        if(size > std::size_t(0))
        {
          this->_copy_from(src_ptr, src_loc, mem_loc, this_offset, size);
        }

        // store our up-to-date locations
        this->valid_locations = mem_loc;
      }

      void copy_to(void* dst_ptr, Location dst_loc, std::size_t this_offset, std::size_t size) const
      {
        // first lock the mutex for this memory block
        std::unique_lock<std::mutex> lock(lock_mutex());

        // make sure that this block is not mapped
        if(this->_is_mapped_any_write())
          throw AccessException("Memory::Block::copy_to: cannot copy to mapped memory block");

        // ensure the ranges are valid
        if(this->bytes < this_offset)
          throw RangeException("Memory::Block::copy_from: invalid destination memory offset");
        if(size == full_size)
          size = this->bytes - this_offset;
        if(this->bytes < this_offset + size)
          throw RangeException("Memory::Block::copy_from: invalid destination memory range");

        // perform the actual copy
        if(size > std::size_t(0))
        {
          this->_copy_to(dst_ptr, dst_loc, this_offset, size);
        }
      }

      void format(Location mem_loc, std::size_t this_offset, std::size_t size)
      {
        std::unique_lock<std::mutex> lock(lock_mutex());
        if(this->_is_mapped_any())
          throw AccessException("Memory::Block::format: cannot format mapped memory block");

        // ensure the ranges are valid
        if(this->bytes < this_offset)
          throw RangeException("Memory::Block::format: invalid destination memory offset");
        if(size == full_size)
          size = this->bytes - this_offset;
        if(this->bytes < this_offset + size)
          throw RangeException("Memory::Block::format: invalid destination memory range");

        // if no target location is specified, then simply format the existing location
        if(mem_loc == Location::none)
          mem_loc = this->valid_locations;

        // if there was no existing location, then let's continue with main memory
        if(mem_loc == Location::none)
          mem_loc = Location::main;

        // reset memory location
        this->valid_locations = Location::none;

        // format main memory?
        if(*(mem_loc & Location::main))
        {
          if(this->mem_main != nullptr)
            Memory::memset_main(reinterpret_cast<char*>(this->mem_main) + this_offset, 0, size);
          else
            this->_alloc_main(Init::format_to_zero);
          this->valid_locations = this->valid_locations | Location::main;
        }

        // format cuda memory?
        if(*(mem_loc & Location::cuda))
        {
          if(this->mem_cuda != nullptr)
            Memory::memset_cuda(reinterpret_cast<char*>(this->mem_cuda) + this_offset, 0, size);
          else
            this->_alloc_cuda(Init::format_to_zero);
          this->valid_locations = this->valid_locations | Location::cuda;
        }
      }

      void format_random(Location mem_loc, std::size_t this_offset, std::size_t size)
      {
        std::unique_lock<std::mutex> lock(lock_mutex());
        if(this->_is_mapped_any())
          throw AccessException("Memory::Block::format_random: cannot format mapped memory block");

        // ensure the ranges are valid
        if(this->bytes < this_offset)
          throw RangeException("Memory::Block::format_random: invalid destination memory offset");
        if(size == full_size)
          size = this->bytes - this_offset;
        if(this->bytes < this_offset + size)
          throw RangeException("Memory::Block::format: invalid destination memory range");

        // ensure that main memory location is allocated
        this->_alloc(Location::main, Init::uninitialized);

        // format main memory block to random
        if(size > 0u)
          memset_random_main(reinterpret_cast<char*>(this->mem_main) + this_offset, size);

        // remember that we wrote to main memory last and all other memory locations are out-of-date now
        this->valid_locations = Location::main;

        // copy from main to whatever the user asked for
        if(*(mem_loc & Location::cuda))
        {
          _copy_cuda_from_main(reinterpret_cast<const char*>(this->mem_main) + this_offset, this_offset, size);
        }
      }

      View map_view(Location mem_loc, Access access, std::size_t offset, std::size_t size)
      {
        std::unique_lock<std::mutex> lock(lock_mutex());

        // do we need to create a view for a sub-region?
        if((offset != std::size_t(0)) || (size != Arbiter::full_size))
        {
          // We are trying to map a write-accessible view to a sub-region of the memory block, so we need to add the
          // read access flag independently of whether the user asked for it, because this will ensure that the
          // entire memory block is up-to-date at the target memory location. If the read access flag was missing and
          // the memory block was not already up-to-date on the target memory location, then the user's write
          // operations would only overwrite the sub-region memory, but the region outside of the mapped view would
          // remain in an undefined state.
          if(*(access & Access::write))
            access = access | Access::read;
        }

        // ensure that the offset is not out-of-bounds
        if(this->bytes < offset)
          throw RangeException("Memory::Block::map_view: memory offset is out of bounds");

        // set size to entire (rest of) memory block if desired (which may be zero, which is fine)
        if(size == full_size)
          size = this->bytes - offset;
        else if(this->bytes < offset + size)
          throw RangeException("Memory::Block::map_view: memory range is out of bounds");

        // do we need write access?
        if(*(access & Access::write))
        {
          // do we allow overlapping write access?
          if(*(access & Access::overlap))
          {
            // yes, overlap is allowed, but we must not have any mapped read views
            if(this->_is_mapped_any_read())
            {
              throw AccessException("Memory::Block::map_view: cannot map read-mapped memory block for overlapping write access");
            }
            else if(this->_is_mapped_any_write())
            {
              // we can live with overlapping access, but some other write view is already mapped to this block,
              // so ensure that a) it also allows overlapped access and b) it is mapped to the same location
              if(!this->write_overlap)
                throw AccessException("Memory::Block::map_view: cannot map exclusively write-mapped memory block for overlapping write access");
              if(this->write_location != mem_loc)
                throw AccessException("Memory::Block::map_view: memory block is write-mapped onto a different memory location");
            }
          }
          else if(this->_is_mapped_any())
          {
            // we require exclusive write access, but there is some other view mapped to this block
            throw AccessException("Memory::Block::map_view: cannot map locked memory block for exclusive write access");
          }
        }

        // do we need read access? if so, then make sure that there is no exclusive write lock on the memory block
        if(*(access & Access::read) && !this->write_overlap)
        {
          if(this->_is_mapped_any_write())
            throw AccessException("Memory::Block::map_view: cannot map write-locked memory block for read access");
        }

        // zero size? create an empty view in this case
        if(size == std::size_t(0))
          return View(this, mem_loc, access, size, nullptr);

        // ensure that the target memory location is allocated
        this->_alloc(mem_loc, Init::dont_care);

        // read access requested?
        if(*(access & Access::read))
        {
          // do we have an up-to-date memory location?
          if(this->valid_locations == Location::none)
          {
            // emit a warning
            std::cerr << "\n>>> MEMORY ARBITER: map_view: trying to map uninitialized memory for reading!\n";
            //XABORTM("\n>>> MEMORY ARBITER: map_view: trying to map uninitialized memory for reading!\n");
#ifdef _WIN32
            Windows::dump_call_stack_to_stderr();
#endif
          }
          // check whether the location that we want to map is out-of-date
          else if(!*(this->valid_locations & mem_loc))
          {
            // the desired target memory location is out of date, so first make sure that
            // main memory is up-to-date and then copy from main to target location
            _copy_main_from_uptodate();
            _copy_main_to(mem_loc);
          }

          // increase read lock counters
          if(*(mem_loc & Location::main))
            ++this->read_views_main;
          if(*(mem_loc & Location::cuda))
            ++this->read_views_cuda;
        }

        // write access requested?
        if(*(access & Access::write))
        {
          if(this->write_views == std::size_t(0))
          {
            // remember whether we asked for overlapping or exclusive access
            this->write_overlap = *(access & Access::overlap);

            // remember memory write location
            this->write_location = mem_loc;
          }

          // increase number of write views
          ++this->write_views;
        }

        // create view for the requested memory location
        if(mem_loc == Location::main)
        {
          return View(this, Location::main, access, size, reinterpret_cast<char*>(this->mem_main) + offset);
        }
        if(mem_loc == Location::cuda)
        {
          return View(this, Location::cuda, access, size, reinterpret_cast<char*>(this->mem_cuda) + offset);
        }

        // we should never arrive here
        XABORTM("Memory::Block::map_view: invalid memory location");
      }

      void release_view(const View& view)
      {
        std::unique_lock<std::mutex> lock(lock_mutex());

        // if the view had a write lock, then remove it
        if(*(view.access() & Access::write))
        {
          if(--this->write_views == std::size_t(0))
          {
            // ok, this was the last open write view
            this->valid_locations = this->write_location;
            this->write_location = Memory::Location::none;
            this->write_overlap = false;
          }
        }

        // if the view had a read lock, then decrease the number of read views by one
        if(*(view.access() & Access::read))
        {
          if(view.location() == Location::main)
          {
            XASSERT(this->read_views_main > 0u);
            --this->read_views_main;
          }
          if(view.location() == Location::cuda)
          {
            XASSERT(this->read_views_cuda > 0u);
            --this->read_views_cuda;
          }
        }
      }

    private:
      // ATTENTION: none of the private functions are thread-safe, so they have to be called from
      // within a mutex-locked scope!
      void _reset()
      {
        this->bytes = 0u;
        this->ref_count = 0u;
        this->mem_main = nullptr;
        this->mem_cuda = nullptr;
        this->valid_locations = Location::none;
        this->write_location = Location::none;
        this->write_views = 0u;
        this->read_views_main = 0u;
        this->read_views_cuda = 0u;
        this->write_overlap = false;
      }

      bool _is_mapped_any() const
      {
        return _is_mapped_any_read() || _is_mapped_any_write();
      }

      bool _is_mapped_any_read() const
      {
        return this->read_views_main + this->read_views_cuda > 0u;
      }

      bool _is_mapped_any_write() const
      {
        return this->write_views > 0u;
      }

      bool _is_mapped_cuda() const
      {
        return _is_mapped_cuda_read() || _is_mapped_cuda_write();
      }

      bool _is_mapped_cuda_read() const
      {
        return read_views_cuda > 0u;
      }

      bool _is_mapped_cuda_write() const
      {
        return (write_views > 0u) && (write_location == Memory::Location::cuda);
      }

      bool _is_mapped_main() const
      {
        return _is_mapped_main_read() || _is_mapped_main_write();
      }

      bool _is_mapped_main_read() const
      {
        return read_views_main > 0u;
      }

      bool _is_mapped_main_write() const
      {
        return (write_views > 0u) && (write_location == Memory::Location::main);
      }

      bool _is_write_overlapped() const
      {
        return write_overlap;
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _alloc_main(Init mem_ini)
      {
        // try to allocate memory
        if((this->mem_main = Memory::alloc_main(this->bytes)) == nullptr)
          throw AllocException(this->bytes, "main");

        // pre-format memory to ones if required
        switch(mem_ini)
        {
        case Init::dont_care:
#ifdef DEBUG
          Memory::memset_main(this->mem_main, 0xFF, this->bytes);
#endif
          break;

        case Init::uninitialized:
          // nothing to do here
          break;

        case Init::format_to_zero:
          Memory::memset_main(this->mem_main, 0, this->bytes);
          this->valid_locations = this->valid_locations | Location::main;
          break;

        case Init::format_to_one:
          Memory::memset_main(this->mem_main, 0xFF, this->bytes);
          this->valid_locations = this->valid_locations | Location::main;
          break;
        }
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _alloc_cuda(Init mem_ini)
      {
        if((this->mem_cuda = Memory::alloc_cuda(this->bytes)) == nullptr)
          throw AllocException(this->bytes, "cuda");

        // pre-format memory to ones if required
        switch(mem_ini)
        {
        case Init::dont_care:
#ifdef DEBUG
          Memory::memset_cuda(this->mem_cuda, 0xFF, this->bytes);
#endif
          break;

        case Init::uninitialized:
          // nothing to do here
          break;

        case Init::format_to_zero:
          Memory::memset_cuda(this->mem_cuda, 0, this->bytes);
          this->valid_locations = this->valid_locations | Location::cuda;
          break;

        case Init::format_to_one:
          Memory::memset_cuda(this->mem_cuda, 0xFF, this->bytes);
          this->valid_locations = this->valid_locations | Location::cuda;
          break;
        }
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _alloc(Location mem_loc, Init mem_ini)
      {
        // nothing to allocate?
        if(this->bytes <= 0u)
          return;

        // allocate main memory?
        if(*(mem_loc & Location::main) && (this->mem_main == nullptr))
          this->_alloc_main(mem_ini);

        // allocate cuda memory?
        if(*(mem_loc & Location::cuda) && (this->mem_cuda == nullptr))
          this->_alloc_cuda(mem_ini);
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _free_main()
      {
        if(this->mem_main)
        {
          XASSERTM(this->read_views_main == std::size_t(0), "Memory::Block::_free_main: Cannot free main memory block that has a read view attached!");
          XASSERTM(this->write_location != Memory::Location::main, "Memory::Block::_free_main: Cannot free main memory block that has a write view attached!");
          Memory::free_main(this->mem_main);
          this->mem_main = nullptr;
          this->valid_locations = this->valid_locations & ~Memory::Location::main;
        }
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _free_cuda()
      {
        if(this->mem_cuda)
        {
          XASSERTM(this->read_views_cuda == std::size_t(0), "Memory::Block::_free_cuda: Cannot free cuda memory block that has a read view attached!");
          XASSERTM(this->write_location != Memory::Location::cuda, "Memory::Block::_free_cuda: Cannot free cuda memory block that has a write view attached!");
          Memory::free_cuda(this->mem_cuda);
          this->mem_cuda = nullptr;
          this->valid_locations = this->valid_locations & ~Memory::Location::cuda;
        }
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_main_from_main(const void* mem_main_src, std::size_t this_offset, std::size_t size)
      {
        //XASSERT(!this->_is_mapped_main());
        XASSERT(mem_main_src != nullptr);
        XASSERT(this_offset + size <= this->bytes);
        if(this->mem_main == nullptr)
          this->_alloc_main(Init::uninitialized);
        Memory::memcopy_main(reinterpret_cast<char*>(this->mem_main) + this_offset, mem_main_src, size);

        // remember that main is an up-to-date location now
        this->valid_locations = this->valid_locations | Location::main;
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_cuda_from_main(const void* mem_main_src, std::size_t this_offset, std::size_t size)
      {
        //XASSERT(!this->_is_mapped_cuda());
        XASSERT(mem_main_src != nullptr);
        XASSERT(this_offset + size <= this->bytes);
        if(this->mem_cuda == nullptr)
          this->_alloc_cuda(Init::uninitialized);
        Memory::memcopy_main_to_cuda(reinterpret_cast<char*>(this->mem_cuda) + this_offset, mem_main_src, size);

        // remember that cuda is an up-to-date location now
        this->valid_locations = this->valid_locations | Location::cuda;
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_main_from_cuda(const void* mem_cuda_src, std::size_t this_offset, std::size_t size)
      {
        //XASSERT(!this->_is_mapped_cuda_write());
        //XASSERT(!this->_is_mapped_main());
        XASSERT(mem_cuda_src != nullptr);
        XASSERT(this_offset + size <= this->bytes);
        if(this->mem_main == nullptr)
          this->_alloc_main(Init::uninitialized);
        Memory::memcopy_cuda_to_main(reinterpret_cast<char*>(this->mem_main) + this_offset, mem_cuda_src, size);

        // remember that main is an up-to-date location now
        this->valid_locations = this->valid_locations | Location::main;
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_cuda_from_cuda(const void* mem_cuda_src, std::size_t this_offset, std::size_t size)
      {
        //XASSERT(!this->_is_mapped_cuda());
        XASSERT(mem_cuda_src != nullptr);
        XASSERT(this_offset + size <= this->bytes);
        if(this->mem_cuda == nullptr)
          this->_alloc_cuda(Init::uninitialized);
        Memory::memcopy_cuda(reinterpret_cast<char*>(this->mem_cuda) + this_offset, mem_cuda_src, size);

        // remember that cuda is an up-to-date location now
        this->valid_locations = this->valid_locations | Location::cuda;
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_main_from_uptodate()
      {
        // is main memory already up to date?
        if(*(this->valid_locations & Location::main))
        {
          // there's nothing left to do then
          return;
        }

        // is cuda memory up-to-date?
        if(*(this->valid_locations & Location::cuda))
        {
          // copy cuda to main then
          XASSERT(!_is_mapped_cuda_write());
          _copy_main_from_cuda(this->mem_cuda, std::size_t(0), this->bytes);
          return;
        }

        // if we come out here, then there is no up-to-date memory to copy
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_main_to(Location mem_loc)
      {
        if(*(mem_loc & Location::cuda))
        {
          XASSERT(!_is_mapped_main_write());
          _copy_cuda_from_main(this->mem_main, std::size_t(0), this->bytes);
        }
      }

      //void _copy(const Block& other, Location mem_loc, std::size_t offset)
      void _copy_from(const Block& source, Location mem_loc, std::size_t this_offset, std::size_t src_offset, std::size_t size)
      {
        if(*(mem_loc & Location::main))
        {
          if(*(source.valid_locations & Location::main))
          {
            XASSERT(!source._is_mapped_main_write());
            _copy_main_from_main(reinterpret_cast<char*>(source.mem_main) + src_offset, this_offset, size);
          }
          else if(*(source.valid_locations & Location::cuda))
          {
            XASSERT(!source._is_mapped_cuda_write());
            _copy_main_from_cuda(reinterpret_cast<char*>(source.mem_cuda) + src_offset, this_offset, size);
          }
        }

        if(*(mem_loc & Location::cuda))
        {
          if(*(source.valid_locations & Location::cuda))
          {
            XASSERT(!source._is_mapped_cuda_write());
            _copy_cuda_from_cuda(reinterpret_cast<char*>(source.mem_cuda) + src_offset, this_offset, size);
          }
          else if(*(source.valid_locations & Location::main))
          {
            XASSERT(!source._is_mapped_main_write());
            _copy_cuda_from_main(reinterpret_cast<char*>(source.mem_main) + src_offset, this_offset, size);
          }
        }
      }

      void _copy_from(const void* src_ptr, Location src_loc, Location mem_loc, std::size_t this_offset, std::size_t size)
      {
        if(*(mem_loc & Location::main))
        {
          switch(src_loc)
          {
          case Location::main:
            _copy_main_from_main(src_ptr, this_offset, size);
            break;

          case Location::cuda:
            _copy_main_from_cuda(src_ptr, this_offset, size);
            break;

          default:
            XABORTM("Memory::Block::_copy_from: invalid source memory location");
            break;
          }
        }

        if(*(mem_loc & Location::cuda))
        {
          switch(src_loc)
          {
          case Location::cuda:
            _copy_cuda_from_cuda(src_ptr, this_offset, size);
            break;

          case Location::main:
            _copy_cuda_from_main(src_ptr, this_offset, size);
            break;

          default:
            XABORTM("Memory::Block::_copy_from: invalid source memory location");
            break;
          }
        }
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_main_to_main(void* mem_main_dst, std::size_t this_offset, std::size_t size) const
      {
        XASSERT(!this->_is_mapped_main_write());
        XASSERT(mem_main_dst != nullptr);
        XASSERT(this->mem_main != nullptr);
        XASSERT(this_offset + size <= this->bytes);
        Memory::memcopy_main(mem_main_dst, reinterpret_cast<char*>(this->mem_main) + this_offset, size);
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_cuda_to_main(void* mem_main_dst, std::size_t this_offset, std::size_t size) const
      {
        XASSERT(!this->_is_mapped_cuda_write());
        XASSERT(mem_main_dst != nullptr);
        XASSERT(this->mem_cuda != nullptr);
        XASSERT(this_offset + size <= this->bytes);
        Memory::memcopy_cuda_to_main(mem_main_dst, reinterpret_cast<char*>(this->mem_cuda) + this_offset, size);
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_main_to_cuda(void* mem_cuda_dst, std::size_t this_offset, std::size_t size) const
      {
        XASSERT(!this->_is_mapped_cuda_write());
        XASSERT(mem_cuda_dst != nullptr);
        XASSERT(this->mem_main != nullptr);
        XASSERT(this_offset + size <= this->bytes);
        Memory::memcopy_main_to_cuda(mem_cuda_dst, reinterpret_cast<char*>(this->mem_main) + this_offset, size);
      }

      // not thread-safe; must be called in a mutex-locked region!
      void _copy_cuda_to_cuda(void* mem_cuda_dst, std::size_t this_offset, std::size_t size) const
      {
        XASSERT(!this->_is_mapped_cuda_write());
        XASSERT(mem_cuda_dst != nullptr);
        XASSERT(this->mem_cuda != nullptr);
        XASSERT(this_offset + size <= this->bytes);
        Memory::memcopy_cuda(mem_cuda_dst, reinterpret_cast<char*>(this->mem_cuda) + this_offset, size);
      }

      void _copy_to(void* dst_ptr, Location dst_loc, std::size_t this_offset, std::size_t size) const
      {
        switch(dst_loc)
        {
        case Location::main:
          if(*(this->valid_locations & Location::main))
            _copy_main_to_main(dst_ptr, this_offset, size);
          else if(*(this->valid_locations & Location::cuda))
            _copy_cuda_to_main(dst_ptr, this_offset, size);
          else
            XABORTM("Memory::Block::_copy_to: invalid destination memory location");
          break;

        case Location::cuda:
          if(*(this->valid_locations & Location::cuda))
            _copy_cuda_to_cuda(dst_ptr, this_offset, size);
          else if(*(this->valid_locations & Location::main))
            _copy_cuda_to_main(dst_ptr, this_offset, size);
          else
            XABORTM("Memory::Block::_copy_to: invalid destination memory location");
          break;

        default:
          XABORTM("Memory::Block::_copy_to: no up-to-date source memory location");
          break;
        }
      }
    }; // class Block

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    View::View() :
      _block(nullptr),
      _location(Location::none),
      _access(Access::none),
      _bytes(0u),
      _ptr(nullptr)
    {
    }

    View::View(Block* blk, Location l, Access a, std::size_t b, void* p) :
      _block(blk),
      _location(l),
      _access(a),
      _bytes(b),
      _ptr(p)
    {
    }

    View::~View()
    {
      if((this->_block != nullptr) && (this->_ptr != nullptr))
      {
        this->_block->release_view(*this);
      }
    }

    /// move-constructor
    View::View(View&& other) :
      _block(other._block),
      _location(other._location),
      _access(other._access),
      _bytes(other._bytes),
      _ptr(other._ptr)
    {
      other._reset();
    }

    /// move-assign operator
    View& View::operator=(View&& other)
    {
      if(this != &other)
      {
        this->release();
        this->_block = other._block;
        this->_location = other._location;
        this->_access = other._access;
        this->_bytes = other._bytes;
        this->_ptr = other._ptr;
        other._reset();
      }
      return *this;
    }

    void View::release() const
    {
      if(this->_block)
        this->_block->release_view(*this);
      const_cast<View*>(this)->_reset();
    }

    void View::_reset()
    {
      this->_block = nullptr;
      this->_location = Location::none;
      this->_access = Access::none;
      this->_bytes = 0u;
      this->_ptr = nullptr;
    }

    void View::format_to_zero()
    {
      if(!this->writeable())
        throw AccessException("Memory::View::format_to_zero(): write access required for format");
      if(this->empty())
        return;

      if(this->_location == Location::main)
        Memory::memset_main(this->_ptr, 0, this->_bytes);
      else if(this->_location == Location::cuda)
        Memory::memset_cuda(this->_ptr, 0, this->_bytes);
    }

    void View::format_to_one()
    {
      if(!this->writeable())
        throw AccessException("Memory::View::format_to_zero(): write access required for format");
      if(this->empty())
        return;

      if(this->_location == Location::main)
        Memory::memset_main(this->_ptr, 0xFF, this->_bytes);
      else if(this->_location == Location::cuda)
        Memory::memset_cuda(this->_ptr, 0xFF, this->_bytes);
    }

    void View::copy_raw_from(const void* src_ptr, Location src_loc)
    {
      XASSERT(src_ptr != nullptr);
      if(!this->writeable())
        throw AccessException("Memory::View::format_to_zero(): write access required for copy");

      // where do we copy to?
      switch(this->_location)
      {
      case Location::main:
        switch(src_loc)
        {
        case Location::main:
          Memory::memcopy_main(this->_ptr, src_ptr, this->_bytes);
          break;

        case Location::cuda:
          Memory::memcopy_cuda_to_main(this->_ptr, src_ptr, this->_bytes);
          break;

        default:
          XABORTM("Memory::View::copy_raw_from: unknown source memory view location");
          break;
        }
        break;

      case Location::cuda:
        switch(src_loc)
        {
        case Location::main:
          Memory::memcopy_main_to_cuda(this->_ptr, src_ptr, this->_bytes);
          break;

        case Location::cuda:
          Memory::memcopy_cuda(this->_ptr, src_ptr, this->_bytes);
          break;

        default:
          XABORTM("Memory::View::copy_raw_from: unknown source memory view location");
          break;
        }
        break;

      default:
        XABORTM("Memory::View::copy_raw_from: unknown destination memory view location");
        break;
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Arbiter::Arbiter(Block* block, std::size_t offset, std::size_t size) :
      _block(block),
      _offset(offset),
      _size(size)
    {
      if(block != nullptr)
      {
        if(block->bytes < offset)
          throw RangeException("Memory::Arbiter: invalid block memory offset");
        if(size != full_size)
        {
          if(block->bytes < offset + size)
            throw RangeException("Memory::Arbiter: invalid block memory range");
        }
      }
      else
      {
        XASSERT(offset == std::size_t(0));
        XASSERT(size == std::size_t(0));
      }
    }

    Arbiter::Arbiter(std::size_t bytes, Location mem_loc, Init mem_ini) :
      _block(new Block()),
      _offset(0),
      _size(full_size)
    {
      _block->create(bytes, mem_loc, mem_ini);
    }

    Arbiter::Arbiter(Arbiter&& other) :
      Arbiter()
    {
      this->_move(std::forward<Arbiter>(other));
    }

    Arbiter& Arbiter::operator=(Arbiter&& other)
    {
      if(this->_block != other._block)
      {
        Block::detach(this->_block);
        this->_block = nullptr;
        this->_move(std::forward<Arbiter>(other));
      }
      return *this;
    }

    Arbiter::~Arbiter()
    {
      Block::detach(this->_block);
    }

    void Arbiter::_move(Arbiter&& other)
    {
      XASSERT(this->_block == nullptr);
      this->_block = other._block;
      this->_offset = other._offset;
      this->_size = other._size;
      other._block = nullptr;
      other._offset = 0u;
      other._size = 0u;
    }

    Arbiter Arbiter::attach(std::size_t offset, std::size_t size) const
    {
      if(this->_block == nullptr)
      {
        if(offset != std::size_t(0))
          throw RangeException("Memory::Arbiter::attach(): invalid block memory offset");
        if((size != full_size) && (size != std::size_t(0)))
          throw RangeException("Memory::Arbiter::attach(): invalid block memory size");
        return Arbiter();
      }

      if(this->_size != full_size)
      {
        if(this->_size < this->_offset + offset)
          throw RangeException("Memory::Arbiter::attach(): invalid block memory offset");
        if(size == full_size)
          return Arbiter(Block::attach(this->_block), this->_offset + offset, this->_size);
        if(this->_size < size)
          throw RangeException("Memory::Arbiter::attach(): invalid block memory range");
        return Arbiter(Block::attach(this->_block), this->_offset + offset, size);
      }
      else // this->_size == full_size
      {
        if(this->_block->bytes < this->_offset + offset)
          throw RangeException("Memory::Arbiter::attach(): invalid block memory offset");
        if(size == full_size)
          return Arbiter(Block::attach(this->_block), this->_offset + offset, full_size);
        if(this->_block->bytes < this->_offset + offset + size)
          throw RangeException("Memory::Arbiter::attach(): invalid block memory range");
        return Arbiter(Block::attach(this->_block), this->_offset + offset, size);
      }
    }

    void Arbiter::attach(const Arbiter& other, std::size_t offset, std::size_t size)
    {
      // self-attach?
      if(this->_block == other._block)
      {
        if(this->_size != full_size)
        {
          if(this->_size < this->_offset + offset)
            throw RangeException("Memory::Arbiter::attach(): invalid block memory offset");
          if(size != full_size)
          {
            if(this->_size < size)
              throw RangeException("Memory::Arbiter::attach(): invalid block memory range");
            this->_size = size;
          }
        }
        else // this->_size == full_size
        {
          if(this->_block->bytes < this->_offset + offset)
            throw RangeException("Memory::Arbiter::attach(): invalid block memory offset");
          if(size != full_size)
          {
            if(this->_block->bytes < this->_offset + offset + size)
              throw RangeException("Memory::Arbiter::attach(): invalid block memory range");
            this->_size = size;
          }
        }
        this->_offset += offset;
        return;
      }

      // detach old block
      Block::detach(this->_block);

      // attach new block
      this->_move(other.attach(offset, size));
    }

    Arbiter Arbiter::clone(std::size_t offset, std::size_t size, bool clone_block) const
    {
      // no block assigned?
      if(this->_block == nullptr)
      {
        XASSERTM(offset == std::size_t(0), "arbiter clone range offset is out-of-bounds");
        XASSERTM((size == std::size_t(0)) || (size == full_size), "arbiter clone range size is out-of-bounds");
        return Arbiter();
      }

      // validate range offset and size
      if(this->_block->bytes < offset)
        throw RangeException("Memory::Arbiter::clone(): invalid block memory offset");
      if((size != full_size) && (this->_block->bytes < offset + size))
        throw RangeException("Memory::Arbiter::clone(): invalid block memory range");

      // create new block and copy content
      if(clone_block)
      {
        // clone entire block
        Arbiter arb(this->_block->bytes);
        arb._block->copy_from(*this->_block, Memory::Location::none, std::size_t(0), std::size_t(0), this->_block->bytes);
        arb._offset = this->_offset + offset;
        arb._size = size;
        return arb;
      }
      else
      {
        // clone selected memory range
        if(size == full_size)
        {
          if(this->_size == full_size)
            size = this->_block->bytes - offset;
          else
            size = this->_size - offset;
        }
        Arbiter arb(size);
        arb._block->copy_from(*this->_block, Memory::Location::none, std::size_t(0), this->_offset + offset, size);
        return arb;
      }
    }

    void Arbiter::clone(const Arbiter& source, std::size_t offset, std::size_t size, bool clone_block)
    {
      // self-clone?
      if((this->_block == source._block) && clone_block)
      {
        // in this case, we only need to update the offset and size, if necessary
        if(this->_size != full_size)
        {
          if(this->_size < this->_offset + offset)
            throw RangeException("Memory::Arbiter::attach(): invalid block memory offset");
          if(size != full_size)
          {
            if(this->_size < size)
              throw RangeException("Memory::Arbiter::attach(): invalid block memory range");
            this->_size = size;
          }
        }
        else // this->_size == full_size
        {
          if(this->_block->bytes < this->_offset + offset)
            throw RangeException("Memory::Arbiter::attach(): invalid block memory offset");
          if(size != full_size)
          {
            if(this->_block->bytes < this->_offset + offset + size)
              throw RangeException("Memory::Arbiter::attach(): invalid block memory range");
            this->_size = size;
          }
        }
        this->_offset += offset;
        return;
      }

      // clone to temporary, detach and move temporary to properly handle self-clone
      Arbiter temp(source.clone(offset, size, clone_block));
      Block::detach(this->_block);
      this->_block = nullptr;
      this->_move(std::move(temp));
    }

    void Arbiter::release()
    {
      Block::detach(this->_block);
      this->_block = nullptr;
    }

    std::size_t Arbiter::bytes() const
    {
      if(_block == nullptr)
        return std::size_t(0u);
      else if(this->_size == full_size)
        return _block->bytes - this->_offset;
      else
        return this->_size;
    }

    std::size_t Arbiter::block_bytes() const
    {
      return (_block != nullptr ? _block->bytes : std::size_t(0u));
    }

    Location Arbiter::location() const
    {
      return (this->_block ? this->_block->valid_locations : Location::none);
    }

    void Arbiter::alloc(Location mem_loc, Init mem_ini)
    {
      if(this->_block)
        this->_block->alloc(mem_loc, mem_ini);
    }

    void Arbiter::free(Location mem_loc)
    {
      if(this->_block)
        this->_block->free(mem_loc);
    }

    void Arbiter::copy(const Arbiter& source, Location mem_loc)
    {
      if(this->_block == source._block)
        return;

      XASSERT(this->_block != nullptr);
      XASSERT(source._block != nullptr);

      std::size_t dst_size = (this->_size == full_size ? this->_block->bytes : this->_size);
      std::size_t src_size = (source._size == full_size ? source._block->bytes : source._size);
      XASSERTM(dst_size == src_size, "source and destination arbiters have different size");

      this->_block->copy_from(*source._block, mem_loc, this->_offset, source._offset, dst_size);
    }

    void Arbiter::copy(const void* src_ptr, Location src_loc, Location mem_loc)
    {
      XASSERT(this->_block != nullptr);
      this->_block->copy_from(src_ptr, src_loc, mem_loc, this->_offset, this->_size);
    }

    void Arbiter::copy_to(void* dst_ptr, Location dst_loc) const
    {
      XASSERT(this->_block != nullptr);
      this->_block->copy_to(dst_ptr, dst_loc, this->_offset, this->_size);
    }

    void Arbiter::format(Location mem_loc)
    {
      if(this->_block)
        this->_block->format(mem_loc, this->_offset, this->_size);
    }

    void Arbiter::format_random(Location mem_loc)
    {
      if(this->_block)
        this->_block->format_random(mem_loc, this->_offset, this->_size);
    }

    View Arbiter::view(Location mem_loc, Access access)
    {
      XASSERT(this->_block != nullptr);
      return this->_block->map_view(mem_loc, access, this->_offset, this->_size);
    }

    View Arbiter::view(Location mem_loc, Access access) const
    {
      XASSERT(this->_block != nullptr);
      XASSERTM(!*(access & Memory::Access::write), "Memory::Arbiter::view(): const overload does not offer write access");
      return this->_block->map_view(mem_loc, access, this->_offset, this->_size);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void Manager::initialize(int, char**)
    {
      // nothing to do here (for now)
    }

    void Manager::finalize()
    {
      // all memory blocks should have been freed by now
      if(manager_block_set.empty())
        return;

      // print a warning
      std::cout << "\nWARNING: Memory::Manager: " << manager_block_set.size()
        << " memory arbiter blocks are still allocated upon program termination\n\n";
    }

    void Manager::free_all_unmapped(Location mem_loc)
    {
      std::unique_lock<std::mutex> lock(manager_mutex);

      for(auto& mbi : manager_block_set)
      {
        mbi->free_unmapped(mem_loc);
      }
    }

    std::size_t Manager::bytes_reserved()
    {
      std::unique_lock<std::mutex> lock(manager_mutex);

      std::size_t r = 0u;
      for(const auto& mbi : manager_block_set)
      {
        r += mbi->bytes;
      }

      return r;
    }

    std::size_t Manager::bytes_allocated_main()
    {
      std::unique_lock<std::mutex> lock(manager_mutex);

      std::size_t r = 0u;
      for(const auto& mbi : manager_block_set)
      {
        if(mbi->mem_main != nullptr)
          r += mbi->bytes;
      }

      return r;
    }

    std::size_t Manager::bytes_allocated_cuda()
    {
      std::unique_lock<std::mutex> lock(manager_mutex);

      std::size_t r = 0u;
      for(const auto& mbi : manager_block_set)
      {
        if(mbi->mem_cuda != nullptr)
          r += mbi->bytes;
      }

      return r;
    }
  } // namespace Memory
} // namespace FEAT
