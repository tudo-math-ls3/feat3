// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/backend.hpp>

#include <memory>
#include <deque>

namespace FEAT
{
  namespace Memory
  {
    // forward declarations
    class Arbiter;
    class Block;
    class Manager;
    class View;

    /**
     * \brief Memory location enumeration
     *
     * This enum is used by the Memory::Arbiter class and its related classes to specify the
     * locations of a memory block and views that are mapped to it.
     *
     * See \ref sec_memarb_location for more information.
     */
    enum class Location : std::uint32_t
    {
      /// no or unknown or unspecified location
      none = 0,
      /// main memory
      main = 0x0001,
      /// cuda memory
      cuda = 0x0002
    }; // enum class Location

    /// \cond internal
    inline Location operator&(Location a, Location b)
    {
      return static_cast<Location>(static_cast<std::uint32_t>(a) & static_cast<std::uint32_t>(b));
    }

    inline Location operator|(Location a, Location b)
    {
      return static_cast<Location>(static_cast<std::uint32_t>(a) | static_cast<std::uint32_t>(b));
    }

    inline Location& operator|=(Location& a, Location b)
    {
      return a = (a | b);
    }

    inline Location operator~(Location a)
    {
      return static_cast<Location>(~static_cast<std::uint32_t>(a));
    }

    inline bool operator*(Location a)
    {
      // only check the relevant bits
      return (static_cast<std::uint32_t>(a) & 0x3u) != 0u;
    }

    /// \endcond

    /**
     * \brief Memory access enumeration
     *
     * This enum is used by the Memory::Arbiter class and its related classes to specify the access
     * rights requested/granted to a view of a memory block.
     *
     * See \ref sec_memarb_access for more details.
     */
    enum class Access : std::uint32_t
    {
      /// no access
      none       = 0,
      /// read-only access
      read       = 0x1,
      /// write-only access (see attention note above)
      write      = 0x2,
      /// read-write access
      read_write = 0x3,
      /// explicitly allow overlapping writes (see \ref sec_memarb_overlap)
      overlap    = 0x4
    }; // enum class Access

    /// \cond internal
    inline Access operator&(Access a, Access b)
    {
      return static_cast<Access>(static_cast<std::uint32_t>(a) & static_cast<std::uint32_t>(b));
    }

    inline Access operator|(Access a, Access b)
    {
      return static_cast<Access>(static_cast<std::uint32_t>(a) | static_cast<std::uint32_t>(b));
    }

    inline bool operator*(Access a)
    {
      // only check the relevant bits
      return (static_cast<std::uint32_t>(a) & 0x7u) != 0u;
    }
    /// \endcond

    /**
     * \brief Memory initialization enumeration
     *
     * This enum is used by the Memory::Arbiter class and its related classes to specify how newly
     * allocated memory arrays are to be initialized.
     *
     * See \ref memarb_arbiter_alloc_init for more information.
     */
    enum class Init : std::uint32_t
    {
      /// don't care
      dont_care = 0,
      /// leave memory uninitialized
      uninitialized = 1,
      /// allocate and format all bits to zero
      format_to_zero = 2,
      /// allocate and format all bits to one
      format_to_one = 3
    }; // enum class Init

    /**
     * \brief Memory related exception base class
     */
    class MemoryException :
      public Exception
    {
    public:
      explicit MemoryException(const String& msg) :
        Exception(msg)
      {
      }
    }; // class MemoryException

    /**
     * \brief Memory allocation exception
     *
     * An exception of this class is thrown when a memory allocation fails.
     */
    class AllocException :
      public MemoryException
    {
    public:
      explicit AllocException(std::size_t bytes, const String& loc) :
        MemoryException(String("Failed to allocate ") + stringify(bytes) + " in " + loc + " memory")
      {
      }
    }; // class AllocException

    /**
     * \brief Memory access exception
     *
     * An exception of this class is thrown when the required memory access cannot be granted or
     * when an operation misses the required access rights.
     */
    class AccessException :
      public MemoryException
    {
    public:
      explicit AccessException(const String& msg) :
        MemoryException(msg)
      {
      }
    }; // class AccessException

    /**
     * \brief Memory range exception
     *
     * An exception of this class is thrown when an invalid memory range (offset and size) are
     * specified.
     */
    class RangeException :
      public MemoryException
    {
    public:
      explicit RangeException(const String& msg) :
        MemoryException(msg)
      {
      }
    }; // class RangeException

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * \brief Raw memory view class
     *
     * This class represents a raw view for a (region of a bigger) memory block that is managed by the Memory::Arbiter
     * class. This class stores a void pointer and the size of the memory block in bytes and it therefore does not
     * offer any functionality related to array elements etc. as this class has no knowledge of how the data is to
     * be interpreted.
     *
     * If you need a view for a memory block that is interpreted as an (array of) object(s) of a certain data type or
     * class, then consider the documentation of the Memory::TypedView class template, which is derived from this class.
     *
     * See \ref memarb_views for more information.
     *
     * \author Peter Zajac
     */
    class View
    {
      // I've got friends!
      friend class Arbiter;
      friend class Block;

    protected:
      /// the block that this view belongs to
      Block* _block;
      /// the current memory location
      Location _location;
      /// the granted memory access
      Access _access;
      /// the size of the memory block in bytes
      std::size_t _bytes;
      /// the actual memory pointer (may be a main or cuda or whatever pointer)
      void* _ptr;

      /// protected constructor; only the Arbiter is allowed to create views via this constructor
      explicit View(Block* blk, Location l, Access a, std::size_t b, void* p);

      /// reset all attributes to null
      void _reset();

    public:
      /// default constructor
      View();

      /// virtual destructor
      virtual ~View();

      /// no copy allowed
      View(const View&) = delete;
      /// no copy-assign allowed
      View& operator=(const View&) = delete;

      /// move-constructor
      View(View&& other);
      /// move-assign operator
      View& operator=(View&& other);

      /**
       * \brief Checks whether the memory view has a size of 0 bytes.
       */
      bool empty() const
      {
        return _bytes <= std::size_t(0u);
      }

      /**
       * \brief Checks whether the memory pointer is a nullptr.
       */
      bool is_nullptr() const
      {
        return _ptr == nullptr;
      }

      /// \returns The size of the memory view in bytes
      std::size_t bytes() const
      {
        return _bytes;
      }

      /// Releases this view
      void release() const;

      /**
       * \brief Returns a raw pointer to the memory block that grants write (and read) access
       *
       * \attention
       * This function will throw an exception if you call this function on a memory view that
       * was created without write access rights!
       */
      void* raw_w()
      {
        if(!this->writeable())
          throw AccessException("Memory::View::raw_w(): write access required");
        return _ptr;
      }

      /**
       * \brief Returns a const raw pointer to the memory block that grants read access
       */
      const void* raw_r() const
      {
        //XASSERT(this->readable());
        return _ptr;
      }

      /// \returns The location of the memory that this view maps to; is always a single location
      Location location() const
      {
        return _location;
      }

      /// \returns The access rights that this view was granted to the memory region
      Access access() const
      {
        return _access;
      }

      /// \returns Checks whether this view was granted read access rights to the memory region
      bool readable() const
      {
        return *(this->_access & Access::read);
      }

      /// \returns Checks whether this view was granted write access rights to the memory region
      bool writeable() const
      {
        return *(this->_access & Access::write);
      }

      /// Formats all bits of the memory view to zero
      void format_to_zero();

      /// Formats all bits of the memory view to one (yields NaN for floating point types)
      void format_to_one();

      /**
       * \brief Copies the contents of a given memory region into the memory of this memory view
       *
       * \param[in] src_ptr
       * A \transient pointer to a memory block that is to be copied into this view's memory
       *
       * \param[in] src_loc
       * The location of the source memory pointer src_ptr.
       */
      void copy_raw_from(const void* src_ptr, Location src_loc);

      /**
       * \brief Copies the contents of another view's memory region into the memory of this memory view
       *
       * In other words: this function performs a deep copy of the source view.
       *
       * \note This function has no effect if &src == this.
       *
       * \param[in] src
       * A \transient reference to the view whose memory is to be copied into this view's memory.
       * The source view may be larger than this view, but it must not be smaller.
       */
      void copy_from(const View& src)
      {
        if(src.bytes() < this->bytes())
          throw RangeException("Memory::View::copy_from(): source memory view is smaller than the destination memory view");
        if(this != &src)
          this->copy_raw_from(src.raw_r(), src.location());
      }
    }; // class View

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * \brief Typed memory view class
     *
     * This class represents a typed view for a (sub-) array of elements of a certain data type. This class derives
     * from the Memory::View class and extends its functionality by granting typed pointers, array element access
     * as well as array element conversion functions.
     *
     * See \ref memarb_views for more information.
     *
     * \author Peter Zajac
     */
    template<typename T_>
    class TypedView :
      public View
    {
    public:
      /// default constructor
      TypedView() = default;

      /// default destructor
      virtual ~TypedView() = default;

      /// no copy
      TypedView(const TypedView&) = delete;
      /// no copy-assign
      TypedView& operator=(const TypedView&) = delete;

      /**
       * \brief View move constructor
       */
      explicit TypedView(View&& other) :
        View(std::forward<View>(other))
      {
      }

      /**
       * \brief View move-assign operator
       */
      TypedView& operator=(View&& other)
      {
        static_cast<View&>(*this) = std::forward<View>(other);
        return *this;
      }

      /**
       * \brief View move constructor
       */
      explicit TypedView(TypedView&& other) :
        View(std::forward<View>(other))
      {
      }

      /**
       * \brief View move-assign operator
       */
      TypedView& operator=(TypedView&& other)
      {
        static_cast<View&>(*this) = std::forward<View>(other);
        return *this;
      }

      /**
       * \brief Deleted type-conversion move constructor
       *
       * This constructor is explicitly deleted to prevent accidental conversions between typed
       * views of different underlying data types, as this is usually unintended behavior.
       *
       * If you really intend to convert a typed view into another typed view with another underlying
       * data type, then you need to manually wrap the source object in a <c>std::forward<Memory::View&&></c>.
       *
       * \see \ref sec_memarb_view_move
       */
      template<typename Tother_>
      TypedView(TypedView<Tother_>&&) = delete;

      /**
       * \brief Deleted type-conversion move-assign operator
       *
       * This operator is explicitly deleted to prevent accidental conversions between typed
       * views of different underlying data types, as this is usually unintended behavior.
       *
       * If you really intend to convert a typed view into another typed view with another underlying
       * data type, then you need to manually wrap the source object in a <c>std::forward<Memory::View&&></c>.
       *
       * \see \ref sec_memarb_view_move
       */
      template<typename Tother_>
      TypedView& operator=(TypedView<Tother_>&&) = delete;

      /**
       * \brief Returns the number of array elements mapped to this view
       */
      Index size() const
      {
        return Index(this->bytes() / sizeof(T_));
      }

      /**
       * \brief Returns a typed pointer to the array that grants write (and read) access
       *
       * \attention
       * This function will fire an assertion if you call this function on a memory view that
       * was created without write access rights!
       */
      T_* get_w()
      {
        if(!this->writeable())
          throw AccessException("Memory::TypedView::get_w(): write access required");
        return reinterpret_cast<T_*>(this->_ptr);
      }

      /**
       * \brief Returns a const typed pointer to the array that grants read access
       */
      const T_* get_r() const
      {
        //XASSERT(this->readable());
        return reinterpret_cast<const T_*>(this->_ptr);
      }

      /**
       * \brief Returns a mutable reference to an array element that grants write (and read) access
       *
       * \param[in] i
       * The index of the array element that is to be accessed
       *
       * \returns
       * A mutable reference to the i-th array element
       */
      T_& operator[](Index i)
      {
        ASSERT(i < size());
        return get_w()[i];
      }

      /**
       * \brief Returns a const reference to an array elements that grants read access
       *
       * \param[in] i
       * The index of the array element that is to be accessed
       *
       * \returns
       * A const reference to the i-th array element
       */
      const T_& operator[](Index i) const
      {
        ASSERT(i < size());
        return get_r()[i];
      }

      /**
       * \brief Returns a const reference to an array elements that grants read access
       *
       * \param[in] i
       * The index of the array element that is to be accessed
       *
       * \returns
       * A const reference to the i-th array element
       */
      const T_& operator()(Index i) const
      {
        ASSERT(i < size());
        return get_r()[i];
      }

      /// overload of convert_from for the identical type
      void convert_from(const T_* ptr_src)
      {
        this->copy_raw_from(ptr_src, this->location());
      }

      /**
       * \brief Converts the array elements of another type to this view's array elements
       *
       * \param[in] ptr_src
       * A \transient pointer to the source array whose elements are to be converted into this view's array.
       * Must be a pointer to memory in the same memory location as this view, otherwise the behavior is
       * undefined and will very likely result in a crash due to a segmentation fault.
       */
      template<typename Tsrc_>
      void convert_from(const Tsrc_* ptr_src)
      {
        switch(this->location())
        {
        case Location::main:
          Memory::convert_main(this->get_w(), ptr_src, this->size());
          break;

        case Location::cuda:
          Memory::convert_cuda(this->get_w(), ptr_src, this->size());
          break;

        default:
          XABORTM("TypedView::convert_from: unknown memory location");
          break;
        }
      }

      /**
       * \brief Converts the array elements of another view to this view's array elements
       *
       * \param[in] view_src
       * A \transient reference to a view whose elements are to be converted into this view's array.
       * Must reside in the same memory location as this view, otherwise an assertion will fire.
       */
      template<typename Tsrc_>
      void convert_from(const TypedView<Tsrc_>& view_src)
      {
        if((const void*)this != (const void*)&view_src)
        {
          XASSERTM(this->location() == view_src.location(), "cannot convert between different memory locations!");
          this->convert_from(view_src.get_r());
        }
      }

      /**
       * \brief Formats the array elements to a given value, i.e. sets all elements to the same value
       *
       * \param[in] value
       * A \transient reference to the value that this view's elements are to be set to
       */
      void format(const T_& value)
      {
        switch(this->location())
        {
        case Location::main:
          Memory::format_main(this->get_w(), value, this->size());
          break;

        case Location::cuda:
          Memory::format_cuda(this->get_w(), value, this->size());
          break;

        default:
          XABORTM("TypedView::format: unknown memory location");
          break;
        }
      }
    }; // class TypedView<T_>

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * \brief Memory arbiter class
     *
     * This class implements a memory-location-aware reference-counted smart pointer, which is capable of automatically
     * transferring memory between host and device memory when needed without explicit synchronization by the user.
     *
     * See \ref memory_arbiter for more information.
     *
     * \author Peter Zajac
     */
    class Arbiter
    {
    public:
      /// full size constant for size parameters
      static constexpr std::size_t full_size = ~std::size_t(0);

    private:
      /// a pointer to the memory block associated with this arbiter
      mutable Block* _block;

      /// the offset of this arbiter's memory region within the block
      std::size_t _offset;

      /// size of this arbiter's memory region within the block
      std::size_t _size;

      /// private constructor
      explicit Arbiter(Block* block, std::size_t offset, std::size_t size);

      /// moves another arbiter into this one
      void _move(Arbiter&& other);

    public:
      /// standard constructor
      Arbiter() :
        _block(nullptr),
        _offset(0u),
        _size(0u)
      {
      }

      /**
       * \brief Creates a new memory arbiter block of a given size
       *
       * \note
       * Although there is little point in doing that, it is perfectly legal to create an arbiter
       * with bytes = 0, in which case a new but empty memory block is created. It is also possible
       * to create views for this memory block and the pointer stored in these views will be a
       * valid non-null pointer in the requested memory location, but the mapped memory block will
       * have a size of 0 bytes (as requested), so no reasonable operations can be performed on
       * that memory block. In consequence, it is also necessary to release the memory block of
       * size 0, as the memory manager will report a memory leak otherwise.
       *
       * \param[in] bytes
       * The size of the memory block to be created in bytes.
       *
       * \param[in] mem_loc
       * A combination of memory locations on which the memory block has to be allocated at time
       * of creation, i.e. during this function call. If set to Location::none, then no memory is
       * allocated initially and the allocation is performed when the first view to this memory
       * block is mapped.
       *
       * \param[in] mem_ini
       * Specifies how the memory of this block is to be initialized if memory is to be allocated
       * during this function call, i.e. if mem_loc != Location::none.
       */
      explicit Arbiter(std::size_t bytes, Location mem_loc = Location::none, Init mem_ini = Init::dont_care);

      /// deleted copy constructor
      Arbiter(const Arbiter& other) = delete;
      /// deleted copy-assign operator
      Arbiter& operator=(const Arbiter& other) = delete;
      /// move constructor
      Arbiter(Arbiter&& other);
      /// move-assign operator
      Arbiter& operator=(Arbiter&& other);

      /// destructor
      ~Arbiter();

      /**
       * \brief Creates a new arbiter that is attached to the same memory block as this arbiter
       *
       * \note
       * This function effectively performs a <em>shallow copy</em> of the memory block, i.e. both
       * this arbiter as well as the returned arbiter will point to the same memory block.
       * If you indent to perform a <em>deep copy</em> of the memory block instead, use the clone()
       * function instead.
       *
       * \returns
       * A new arbiter that is attached to this arbiters memory block.
       */
      Arbiter attach(std::size_t offset = 0u, std::size_t size = full_size) const;

      /**
       * \brief Attaches this arbiter to another arbiter's memory block
       *
       * This function is effectively equivalent to
       * \code{.cpp}
         this->release();
         (*this) = other.attach(offset, size);
       * \endcode
       *
       * \note
       * This function effectively performs a <em>shallow copy</em> of the memory block, i.e. both
       * the input arbiter \p other as well as this arbiter will point to the same memory block.
       * If you indent to perform a <em>deep copy</em> instead, use the clone() function instead.
       *
       * \param[in] other
       * A \transient reference to another arbiter to whose block this arbiter should attach to.
       */
      void attach(const Arbiter& other, std::size_t offset = 0u, std::size_t size = full_size);

      /**
       * \brief Creates a new arbiter that contains a deep copy to (a sub-region of) this arbiters memory block
       *
       * \note
       * This function performs a <em>deep copy</em> of (a sub-region of) the memory block, i.e.
       * this arbiter and the returned arbiter will point to independent memory blocks. If you
       * indent to perform a <em>shallow copy</em> of the memory block instead, use the attach()
       * function instead.
       *
       * \param[in] offset
       * Specifies the offset of the memory sub-region in bytes that is to be copied.
       *
       * \param[in] size
       * Specifies the size of the memory sub-region in bytes that is to be copied.
       * Set to Arbiter::full_size if the memory sub-region should extend to the end of the input
       * memory block (sub-region).
       *
       * \param[in] clone_block
       * Specifies if the entire memory block is to be cloned or only the sub-region of the block
       * that is defined by the range [offset, offset+size).
       *
       * \returns
       * A new arbiter that contains a deep copy of (a sub-region of) this arbiters memory block.
       */
      Arbiter clone(std::size_t offset = 0u, std::size_t size = full_size, bool clone_block = false) const;

      /**
       * \brief Sets this arbiter to a deep copy of (a sub-region of) another arbiter's memory block
       *
       * This function is effectively equivalent to
       * \code{.cpp}
         this->release();
         (*this) = other.clone(offset, size);
       * \endcode
       *
       * \note
       * This function performs a <em>deep copy</em> of the memory block, i.e. the input arbiter
       * \p other and this arbiter will point to independent memory blocks. If you indent to perform
       * a <em>shallow copy</em> instead, use the attach() function instead.
       *
       * \param[in] source
       * A \transient reference to another arbiter to whose block this arbiter should clone from.
       *
       * \param[in] offset
       * Specifies the offset of the memory sub-region in bytes that is to be copied.
       *
       * \param[in] size
       * Specifies the size of the memory sub-region in bytes that is to be copied.
       * Set to Arbiter::full_size if the memory sub-region should extend to the end of the input
       * memory block (sub-region).
       *
       * \param[in] clone_block
       * Specifies if the entire memory block is to be cloned or only the sub-region of the block
       * that is defined by the range [offset, offset+size).
       */
      void clone(const Arbiter& source, std::size_t offset = 0u, std::size_t size = full_size, bool clone_block = false);

      /**
       * \brief Releases this arbiter from its associated memory block and frees it if necessary
       *
       * This function decrements the reference counter of this arbiters memory block by one and,
       * if the reference counter reaches zero, frees all memory allocated within the memory block
       * as well as the block itself.
       *
       * This function has no effect if the arbiter is not attached to a memory block.
       */
      void release();

      /**
       * \brief Returns the size of the arbiter's range in the attached memory block in bytes
       */
      std::size_t bytes() const;

      /**
      * \brief Returns the size of the entire attached memory block in bytes
      */
      std::size_t block_bytes() const;

      /**
       * \brief Returns a combination of all up-to-date memory locations
       */
      Location location() const;

      /**
       * \brief Allocates memory for a given memory block
       *
       * \param[in] mem_loc
       * A combination of memory locations on which the memory block has to be allocated
       *
       * \param[in] mem_ini
       * Specifies how the memory of this block is to be initialized if memory is to be allocated
       */
      void alloc(Location mem_loc, Init mem_ini = Init::dont_care);

      /**
       * \brief Frees memory for a given memory block
       *
       * \param[in] mem_loc
       * A combination of memory locations on which the memory block has to be freed
       */
      void free(Location mem_loc);

      /**
       * \brief Formats the memory block to zero
       *
       * \param[in] mem_loc
       * A combination of memory locations which have to be formatted. If set to Location::none,
       * then all memory locations that are currently allocated will be formatted.
       */
      void format(Location mem_loc = Location::none);

      /**
       * \brief Formats the memory block to a given value
       *
       * This function interprets the memory block as an array of elements of type \p T_ and sets
       * all elements of the array to the given value
       *
       * \tparam T_
       * Specifies the type of the array elements that the memory block's contents are to be
       * interpreted as
       *
       * \param[in] value_
       * The value that all array elements are to be set to.
       *
       * \param[in] mem_loc
       * A combination of memory locations which have to be formatted. If set to Location::none,
       * then all memory locations that are currently allocated will be formatted.
       */
      template<typename T_>
      void format(const T_& value, Location mem_loc)
      {
        if(mem_loc == Location::none)
          mem_loc = this->location();

        if(*(mem_loc & Location::main))
          this->template typed_view<T_>(Location::main, Access::write).format(value);
        if(*(mem_loc & Location::cuda))
          this->template typed_view<T_>(Location::cuda, Access::write).format(value);
      }

      /**
       * \brief Formats the memory block to random values
       *
       * \param[in] mem_loc
       * A combination of memory locations which have to be formatted. If set to Location::none,
       * then all memory locations that are currently allocated will be formatted.
       */
      void format_random(Location mem_loc = Location::none);

      /**
       * \brief Copies of the contents of one memory block into another one
       *
       * This function performs a deep copy of the currently up-to-date memory location(s) of the
       * source memory block into the memory locations of the destination memory block that are
       * specified by mem_loc.
       *
       * \param[in] source
       * The handle of the source memory block that is to be copied.
       *
       * \param[in] mem_loc
       * Specifies which memory location of the destination memory block should receive the copy.
       * If set to Location::none, then all up-to-date locations of the source memory block are
       * copied into the destination memory block.
       */
      void copy(const Arbiter& source, Location mem_loc = Location::none);

      /**
       * \brief Copies of the contents of a raw memory region into this memory block
       *
       * This function performs a deep copy of the given source memory region into the memory
       * locations of the destination memory block that are specified by mem_loc.
       *
       * \param[in] src_ptr
       * A \transient pointer to the source memory region that is to be deep copied.
       *
       * \param[in] src_loc
       * The memory location of the source memory region.
       * Must be a single memory location and must not be Location::none.
       *
       * \param[in] mem_loc
       * Specifies which memory location of the destination memory block should receive the copy.
       * If set to Location::none, then the destination memory location will be equal to the source
       * memory location \p src_loc
       */
      void copy(const void* src_ptr, Location src_loc, Location mem_loc = Location::none);

      /**
       * \brief Copies the contents of this memory block into a raw memory region
       *
       * This function performs a deep copy of this memory block's contents into the given memory
       * region. This function automatically chooses a valid source memory location in this block.
       *
       * \param[out] dst_ptr
       * A \transient pointer to the destination memory region that is to be deep copied to.
       *
       * \param[in] dst_loc
       * The memory location of the source memory region.
       * Must be a single memory location and must not be Location::none.
       */
      void copy_to(void* dst_ptr, Location dst_loc) const;

      /**
       * \brief Converts the contents of one memory block into another one by casting the underlying data types
       *
       * This function interprets the contents of the destination memory block as an array of
       * elements of type \p Tdst_ and the contents of the source memory block as an array of elements
       * of type \p Tsrc_ and it converts each element of the source array into its corresponding
       * element of the destination array by performing a static cast to the destination data type.
       *
       * \tparam Tdst_
       * Specifies the underlying type of the destination array.
       *
       * \tparam Tsrc_
       * Specifies the underlying type of the source array; must be convertible to Tdst_.
       *
       * \param[in] source
       * The arbiter of the source memory block that is to be converted.
       *
       * \param[in] mem_loc
       * Specifies which memory location of the destination memory block should receive the copy.
       * Must be a single memory location and must not be Location::none.
       */
      template<typename Tdst_, typename Tsrc_>
      void convert(const Arbiter& source, Location mem_loc)
      {
        TypedView<Tdst_> view_dst = this->template typed_view<Tdst_>(mem_loc, Access::write);
        TypedView<Tsrc_> view_src = source.template typed_view<Tsrc_>(mem_loc, Access::read);
        view_dst.convert_from(view_src);
      }

      /**
       * \brief Converts the contents of one memory block into another one by casting the underlying data types
       *
       * This function interprets the contents of the destination memory block as an array of
       * elements of type Tdst_ and the contents of the source memory block as an array of elements
       * of type Tsrc_ and it converts each element of the source array into its corresponding
       * element of the destination array by performing a static cast to the destination data type.
       *
       * \tparam Tdst_
       * Specifies the underlying type of the destination.
       *
       * \tparam Tsrc_
       * Specifies the underlying type of the source array; must be convertible to Tdst_.
       *
       * \param[in] source
       * The arbiter of the source memory block that is to be converted.
       */
      template<typename Tdst_, typename Tsrc_>
      void convert(const Arbiter& source)
      {
        Location mem_loc = this->location() | source.location();
        if(*(mem_loc & Location::main))
        {
          TypedView<Tdst_> view_dst = this->template typed_view<Tdst_>(Location::main, Access::write);
          TypedView<Tsrc_> view_src = source.template typed_view<Tsrc_>(Location::main, Access::read);
          view_dst.convert_from(view_src);
        }
        if(*(mem_loc & Location::cuda))
        {
          TypedView<Tdst_> view_dst = this->template typed_view<Tdst_>(Location::cuda, Access::write);
          TypedView<Tsrc_> view_src = source.template typed_view<Tsrc_>(Location::cuda, Access::read);
          view_dst.convert_from(view_src);
        }
      }

      /**
       * \brief Maps a raw view to a memory block
       *
       * \param[in] mem_loc
       * The memory location that the view is to be mapped for. Must be a single location and must
       * not be Location::none.
       *
       * \param[in] access
       * The access rights that are to be granted for the view.
       *
       * \returns
       * A raw view for the (sub-region of the) memory block.
       */
      View view(Location mem_loc, Access access);

      /**
       * \brief Maps a raw read-only view to the arbiters memory block
       *
       * \attention
       * For this const-overload the access must be equal to Memory::access::read, otherwise an
       * assertion will fire.
       *
       * \param[in] mem_loc
       * The memory location that the view is to be mapped for.
       * Must be a single location and must not be Location::none.
       *
       * \param[in] access
       * The access rights that are to be granted for the view. Must be equal to Memory::Access:read
       * for this overload.
       *
       * \returns
       * A raw view for the memory block.
       */
      View view(Location mem_loc, Access access = Memory::Access::read) const;

      /**
       * \brief Maps a typed view to a memory block
       *
       * \tparam T_
       * The type that is to be used for the typed view.
       *
       * \param[in] mem_loc
       * The memory location that the view is to be mapped for. Must be a single location and must
       * not be Location::none.
       *
       * \param[in] access
       * The access rights that are to be granted for the view.
       *
       * \returns
       * A typed view for the (sub-region of the) memory block of the requested type T_.
       */
      template<typename T_>
      TypedView<T_> typed_view(Location mem_loc, Access access)
      {
        return TypedView<T_>(this->view(mem_loc, access));
      }

      /**
       * \brief Maps a read-only typed view to a memory block
       *
       * \tparam T_
       * The type that is to be used for the typed view.
       *
       * \param[in] mem_loc
       * The memory location that the view is to be mapped for. Must be a single location and must
       * not be Location::none.
       *
       * \param[in] access
       * The access rights that are to be granted for the view. Must be equal to Memory::Access:read
       * for this overload.
       *
       * \returns
       * A typed view for the (sub-region of the) memory block of the requested type T_.
       */
      template<typename T_>
      TypedView<T_> typed_view(Location mem_loc, Access access = Access::read) const
      {
        return TypedView<T_>(this->view(mem_loc, access));
      }

      /**
       * \brief Checks whether this arbiter is not the null arbiter
       */
      operator bool() const
      {
        return this->_block != nullptr;
      }

      /// Checks whether two arbiters refer to the same memory region
      bool operator==(const Arbiter& other) const
      {
        return (this->_block == other._block) && (this->_offset == other._offset) && (this->_size == other._size);
      }

      /// Checks whether two arbiters refer to the different memory regions
      bool operator!=(const Arbiter& other) const
      {
        return (this->_block != other._block) || (this->_offset != other._offset) || (this->_size != other._size);
      }

      /// Checks whether two arbiters refer to the same memory block ignoring offset and size
      bool same_block(const Arbiter& other) const
      {
        return (this->_block == other._block);
      }
    }; // class Arbiter

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * \brief Memory management interface class
     *
     * This class acts as an interface for all kinds on memory-management settings.
     *
     * Please note that there is no need to directly interact with this class except for debugging and profiling
     * purposes.
     *
     * \author Peter Zajac
     */
    class Manager
    {
    public:
      /**
       * \brief Initializes FEAT's memory management
       *
       * This function is automatically called by Runtime::initialize(), so there is no point in calling it by yourself.
       */
      static void initialize(int argc, char** argv);

      /**
       * \brief Finalizes FEAT's memory management
       *
       * This function is automatically called by Runtime::finalize(), so there is no point in calling it by yourself.
       */
      static void finalize();

      /**
       * \brief Frees all memory blocks for a given set of memory location that are currently not mapped to a view
       *
       * \param[in] mem_loc
       * A combination of memory locations on which the currently unmapped memory blocks have to be freed
       *
       * \note
       * This function never frees main memory, even if requested.
       */
      static void free_all_unmapped(Location mem_loc);

      /**
       * \brief Returns the total number of bytes reserved over all memory blocks
       */
      static std::size_t bytes_reserved();

      /**
       * \brief Returns the total number of bytes currently allocated in main memory over all blocks
       */
      static std::size_t bytes_allocated_main();

      /**
       * \brief Returns the total number of bytes currently allocated in cuda memory over all blocks
       */
      static std::size_t bytes_allocated_cuda();
    }; // class Manager
  } // namespace Memory
} // namespace FEAT
