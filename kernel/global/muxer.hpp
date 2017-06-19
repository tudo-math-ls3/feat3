#pragma once
#ifndef KERNEL_GLOBAL_MUXER_HPP
#define KERNEL_GLOBAL_MUXER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/lafem/dense_vector.hpp>

#include <vector>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Global multiplexer implementation
     *
     * \tparam LocalVector_
     * The local vector type that the muxer is to be applied onto.
     *
     * \tparam Mirror_
     * The mirror type that is to be used.
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename Mirror_>
    class Muxer
    {
    public:
      /// the memory type
      typedef typename LocalVector_::MemType MemType;
      /// the data type
      typedef typename LocalVector_::DataType DataType;
      /// the index type
      typedef typename LocalVector_::IndexType IndexType;
      /// the internal buffer vector type (possibly in device memory)
      typedef LAFEM::DenseVector<MemType, DataType, IndexType> BufferType;
      /// the internal buffer vector type in main memory
      typedef LAFEM::DenseVector<Mem::Main, DataType, IndexType> BufferMain;

      /// Our 'base' class type
      template <typename LocalVector2_, typename Mirror2_>
      using MuxerType = class Muxer<LocalVector2_, Mirror2_>;

      /// this typedef lets you create a gate container with new Memory, Data and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using MuxerTypeByMDI = class Muxer<typename LocalVector_::template ContainerType<Mem2_, DataType2_, IndexType2_>, typename Mirror_::template MirrorType<Mem2_, DataType2_, IndexType2_> >;

    public:
      /// the sibling communicator
      const Dist::Comm* _sibling_comm;

      /// the rank of the parent
      int _parent_rank;

      /// buffer size
      Index _buffer_size;

      /// parent mirror (on all children)
      Mirror_ _parent_mirror;

      /// child mirrors (only on parent)
      std::vector<Mirror_> _child_mirrors;

    public:
      /// standard constructor
      explicit Muxer() :
        _sibling_comm(nullptr),
        _parent_rank(-1),
        _buffer_size(0),
        _parent_mirror(),
        _child_mirrors()
      {
      }

      /// move-constructor
      Muxer(Muxer&& other) :
        _sibling_comm(other._sibling_comm),
        _parent_rank(other._parent_rank),
        _buffer_size(other._buffer_size),
        _parent_mirror(std::forward<Mirror_>(other._parent_mirror)),
        _child_mirrors(std::forward<std::vector<Mirror_>>(other._child_mirrors))
      {
      }

      /// virtual destructor
      virtual ~Muxer()
      {
      }


      /// move-assignment operator
      Muxer& operator=(Muxer&& other)
      {
        if(this == &other)
          return *this;

        _sibling_comm = other._sibling_comm;
        _parent_rank = other._parent_rank;
        _buffer_size = other._buffer_size;
        _parent_mirror = std::forward<Mirror_>(other._parent_mirror);
        _child_mirrors = std::forward<std::vector<Mirror_>>(other._child_mirrors);

        return *this;
      }

      /// conversion function
      template<typename LVT2_, typename MT2_>
      void convert(const Muxer<LVT2_, MT2_>& other)
      {
        if((void*)this == (void*)&other)
          return;

        this->_sibling_comm = other._sibling_comm;
        this->_parent_rank = other._parent_rank;
        this->_buffer_size = other._buffer_size;
        this->_parent_mirror.convert(other._parent_mirror);

        for(std::size_t i(0); i < other._child_mirrors.size(); ++i)
        {
          this->_child_mirrors.at(i).convert(other._child_mirrors.at(i));
        }
      }

      /// Returns the internal data size in bytes.
      std::size_t bytes() const
      {
        std::size_t b = _parent_mirror.bytes();
        for(const auto& cm : _child_mirrors)
          b += cm.bytes();
        return b;
      }

      /**
       * \brief Specifies whether this process represents a child in the muxer.
       *
       * \returns
       * \c true, if this process is a child, otherwise \c false.
       */
      bool is_child() const
      {
        return (_sibling_comm != nullptr) && (_sibling_comm->size() > 0);
      }

      /**
       * \brief Specifies whether this process represents a parent in the muxer.
       *
       * \returns
       * \c true, if this process is a parent, otherwise \c false.
       */
      bool is_parent() const
      {
        return (_sibling_comm != nullptr) && (_sibling_comm->rank() == _parent_rank);
      }

      /**
       * \brief Specifies whether this process is a ghost in the muxer.
       *
       * \note
       * A ghost is a child process, which is not a parent process.
       *
       * \returns
       * \c true, if this process is a ghost, otherwise \c false.
       */
      bool is_ghost() const
      {
        return is_child() && (!is_parent());
      }

      /**
       * \brief Returns the sibling communicator.
       *
       * \returns The sibling communicator.
       */
      const Dist::Comm* get_sibling_comm() const
      {
        return _sibling_comm;
      }

      /**
       * \brief Sets the sibling communicator.
       *
       * \param[in] sibling_comm_
       * The sibling communicator.
       *
       * \param[in] parent_rank
       * The rank of the parent process in the sibling comm.
       *
       * \param[in] parent_mirror
       * The parent mirror.
       */
      void set_parent(const Dist::Comm* sibling_comm_, int parent_rank, Mirror_&& parent_mirror)
      {
        _sibling_comm = sibling_comm_;
        XASSERT(_sibling_comm != nullptr);
        XASSERT(parent_rank >= 0);
        XASSERT(parent_rank < _sibling_comm->size());
        _parent_rank = parent_rank;
        _parent_mirror = std::forward<Mirror_>(parent_mirror);
      }

      /**
       * \brief Adds a child rank and mirror for a parent process.
       *
       * \param[in] child_rank
       * The rank of the child patch of this parent process.
       *
       * \param[in] child_mirror
       * The mirror of the child patch of this parent process.
       */
      void push_child(Mirror_&& child_mirror)
      {
        _child_mirrors.push_back(std::forward<Mirror_>(child_mirror));
      }

      /**
       * \brief Compiles the muxer.
       *
       * \param[in] vec_tmp_
       * A temporary vector for internal use.
       */
      void compile(const LocalVector_& vec_tmp_)
      {
        if(!is_child())
          return; // nothing to

        XASSERT(_sibling_comm != nullptr);
        XASSERT(_parent_rank >= 0);
        XASSERT(_parent_rank < _sibling_comm->size());

        Index bufsize(0);

        // Is this the parent?
        if(_sibling_comm->rank() == _parent_rank)
        {
          // loop over all child buffers and compute the largest buffer size
          for(const auto& cm : _child_mirrors)
          {
            bufsize = Math::max(bufsize, cm.buffer_size(vec_tmp_));
          }
        }

        // broadcast buffer size to all siblings
        _sibling_comm->bcast(&bufsize, std::size_t(1), _parent_rank);

        // verify against length of parent buffer
        XASSERT(bufsize >= _parent_mirror.buffer_size(vec_tmp_));

        // store buffer size
        this->_buffer_size = bufsize;
      }

      /**
       * \brief Sends a join operation to the parent process.
       *
       * \param[in] vec_src
       * The type-0 child vector that is to be joined on the parent process.
       */
      void join_send(const LocalVector_& vec_src) const
      {
        XASSERT(_sibling_comm != nullptr);
        XASSERT(_sibling_comm->size() > 1);
        XASSERT(_sibling_comm->rank() != _parent_rank); // parent must call join() instead

        TimeStamp ts_execute;

        // create parent buffer in device memory
        BufferType parent_buffer(_buffer_size);

        // gather source to parent buffer
        _parent_mirror.gather(parent_buffer, vec_src);

        // convert buffer to main memory
        BufferMain parent_buffer_main;
        parent_buffer_main.convert(parent_buffer);

        Statistics::add_time_mpi_execute(ts_execute.elapsed_now());

        // gather to parent sibling
        TimeStamp ts_collective;
        DataType* dummy = nullptr;
        _sibling_comm->gather(parent_buffer_main.elements(), _buffer_size, dummy, std::size_t(0), _parent_rank);
        Statistics::add_time_mpi_wait_collective(ts_collective.elapsed_now());
      }

      /**
       * \brief Performs a join operation on the parent process.
       *
       * \param[in] vec_src
       * The type-0 child vector that is to be joined on the parent process.
       *
       * \param[in] vec_trg
       * The joined type-0 parent vector.
       */
      void join(const LocalVector_& vec_src, LocalVector_& vec_trg) const
      {
        // if this muxer is not a child, then this operation is a simple copy
        if((_sibling_comm == nullptr) || (_sibling_comm->size() <= 1))
        {
          vec_trg.copy(vec_src);
          return;
        }

        XASSERT(_sibling_comm != nullptr);
        XASSERT(_sibling_comm->size() > 0);
        XASSERT(_sibling_comm->rank() == _parent_rank);

        TimeStamp ts_execute;

        const Index num_children = Index(_child_mirrors.size());

        // create parent buffer in device memory
        BufferType parent_buffer(_buffer_size);

        // gather source to parent buffer
        _parent_mirror.gather(parent_buffer, vec_src);

        // convert buffer to main memory
        BufferMain parent_buffer_main;
        parent_buffer_main.convert(parent_buffer);

        // allocate child buffers
        BufferMain child_buffers_main(_buffer_size * num_children);

        Statistics::add_time_mpi_execute(ts_execute.elapsed_now());

        // gather from siblings
        TimeStamp ts_collective;
        _sibling_comm->gather(parent_buffer_main.elements(), _buffer_size, child_buffers_main.elements(), _buffer_size, _parent_rank);
        Statistics::add_time_mpi_wait_collective(ts_collective.elapsed_now());

        ts_execute.stamp();

        // download child buffers to device memory
        BufferType child_buffers;
        child_buffers.convert(child_buffers_main);

        // format target vector
        vec_trg.format();

        // scatter child buffers into target vector
        for(Index i(0); i < num_children; ++i)
        {
          // scatter buffers
          _child_mirrors.at(i).scatter_axpy(vec_trg, child_buffers, DataType(1), i*_buffer_size);
        }
        Statistics::add_time_mpi_execute(ts_execute.elapsed_now());
      }

      /**
       * \brief Receives a split operation from the parent process.
       *
       * \param[out] vec_trg
       * The split type-1 child vector.
       */
      void split_recv(LocalVector_& vec_trg) const
      {
        XASSERT(_sibling_comm != nullptr);
        XASSERT(_sibling_comm->size() > 1);
        XASSERT(_sibling_comm->rank() != _parent_rank); // parent must call split() instead

        // allocate parent buffer
        BufferMain parent_buffer_main(_buffer_size);

        // receive scatter from parent sibling
        TimeStamp ts_collective;
        DataType* dummy = nullptr;
        _sibling_comm->scatter(dummy, std::size_t(0), parent_buffer_main.elements(), _buffer_size, _parent_rank);
        Statistics::add_time_mpi_wait_collective(ts_collective.elapsed_now());

        TimeStamp ts_execute;

        // convert to device memory
        BufferType parent_buffer;
        parent_buffer.convert(parent_buffer_main);

        // scatter into target vector
        vec_trg.format();
        _parent_mirror.scatter_axpy(vec_trg, parent_buffer);
        Statistics::add_time_mpi_execute(ts_execute.elapsed_now());
      }

      /**
       * \brief Performs a split operation on the parent process.
       *
       * \param[out] vec_trg
       * The split type-1 child vector.
       *
       * \param[in] vec_src
       * The type-1 parent vector that is to be split
       */
      void split(LocalVector_& vec_trg, const LocalVector_& vec_src) const
      {
        // if this muxer is not a child, then this operation is a simple copy
        if((_sibling_comm == nullptr) || (_sibling_comm->size() <= 1))
        {
          vec_trg.copy(vec_src);
          return;
        }

        XASSERT(_sibling_comm != nullptr);
        XASSERT(_sibling_comm->size() > 0);
        XASSERT(_sibling_comm->rank() == _parent_rank);

        TimeStamp ts_execute;

        const Index num_children = Index(_child_mirrors.size());

        // allocate child buffers
        BufferType child_buffers(_buffer_size * num_children);

        // gather child buffers from target vector
        for(Index i(0); i < num_children; ++i)
        {
          _child_mirrors.at(i).gather(child_buffers, vec_src, i*_buffer_size);
        }

        // convert child buffers to main memory
        BufferMain child_buffers_main;
        child_buffers_main.convert(child_buffers);

        // create parent buffer
        BufferMain parent_buffer_main(_buffer_size);

        Statistics::add_time_mpi_execute(ts_execute.elapsed_now());

        // scatter to siblings
        TimeStamp ts_collective;
        _sibling_comm->scatter(child_buffers_main.elements(), _buffer_size, parent_buffer_main.elements(), _buffer_size, _parent_rank);
        Statistics::add_time_mpi_wait_collective(ts_collective.elapsed_now());

        ts_execute.stamp();

        // convert to device memory
        BufferType parent_buffer;
        parent_buffer.convert(parent_buffer_main);

        // scatter into target vector
        vec_trg.format();
        _parent_mirror.scatter_axpy(vec_trg, parent_buffer);

        Statistics::add_time_mpi_execute(ts_execute.elapsed_now());
      }
    }; // class Muxer<...>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_MUXER_HPP
