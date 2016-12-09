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
      /// the internal buffer vector type
      typedef LAFEM::DenseVector<Mem::Main, DataType, IndexType> BufferVectorType;

      /// Our 'base' class type
      template <typename LocalVector2_, typename Mirror2_>
      using MuxerType = class Muxer<LocalVector2_, Mirror2_>;

      /// this typedef lets you create a gate container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using MuxerTypeByMDI = class Muxer<typename LocalVector_::template ContainerType<Mem2_, DataType2_, IndexType2_>, typename Mirror_::template MirrorType<Mem2_, DataType2_, IndexType2_> >;

    public:
      /// the child communicator
      const Dist::Comm* _comm;

      /// parent rank (on all children)
      int _parent_rank;
      /// parent mirror (on all children)
      Mirror_ _parent_mirror;
      /// parent buffer (on all children)
      mutable BufferVectorType _parent_buffer;

      /// child ranks (only on parent)
      std::vector<int> _child_ranks;
      /// child mirrors (only on parent)
      std::vector<Mirror_> _child_mirrors;
      /// child buffers (only on parent)
      mutable std::vector<BufferVectorType> _child_buffers;
      /// child requests (only on parent)
      mutable Dist::RequestVector _child_reqs;

    public:
      /// standard constructor
      explicit Muxer() :
        _comm(nullptr),
        _parent_rank(-1)
      {
      }

      /// move-constructor
      Muxer(Muxer&& other) :
        _comm(other._comm),
        _parent_rank(other._parent_rank),
        _parent_mirror(std::forward<Mirror_>(other._parent_mirror)),
        _parent_buffer(std::forward<BufferVectorType>(other._parent_buffer)),
        _child_ranks(std::forward<std::vector<int>>(other._child_ranks)),
        _child_mirrors(std::forward<std::vector<Mirror_>>(other._child_mirrors)),
        _child_buffers(std::forward<std::vector<BufferVectorType>>(other._child_buffers)),
        _child_reqs(std::forward<Dist::RequestVector>>(other._child_reqs))
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

        _comm = other._comm;
        _parent_rank = other._parent_rank;
        _parent_mirror = std::forward<Mirror_>(other._parent_mirror);
        _parent_buffer = std::forward<BufferVectorType>(other._parent_buffer);
        _child_ranks = std::forward<std::vector<int>>(other._child_ranks);
        _child_mirrors = std::forward<std::vector<Mirror_>>(other._child_mirrors);
        _child_buffers = std::forward<std::vector<BufferVectorType>>(other._child_buffers);
        _child_reqs = std::forward<Dist::RequestVector>>(other._child_reqs);
        return *this;
      }

      /// conversion function
      template<typename LVT2_, typename MT2_>
      void convert(const Muxer<LVT2_, MT2_>& other)
      {
        if((void*)this == (void*)&other)
          return;

        this->_comm = other._comm;

        this->_parent_rank = other._parent_rank;
        this->_parent_mirror.convert(other._parent_mirror);
        this->_parent_buffer.convert(other._parent_buffer);

        this->_child_ranks = other._child_ranks;
        this->_child_mirrors.resize(other._child_mirrors.size());
        this->_child_buffers.resize(other._child_buffers.size());
        this->_child_reqs.resize(other._child_reqs.size());

        for(std::size_t i(0); i < other._child_ranks.size(); ++i)
        {
          this->_child_mirrors.at(i).convert(other._child_mirrors.at(i));
          this->_child_buffers.at(i).convert(other._child_buffers.at(i));
        }
      }

      /// Returns the internal data size in bytes.
      std::size_t bytes() const
      {
        std::size_t b = _parent_mirror.bytes() + _parent_buffer.bytes();
        b += _child_ranks.size() * sizeof(int);
        for(const auto& cm : _child_mirrors)
          b += cm.bytes();
        for(const auto& cb : _child_buffers)
          b += cb.bytes();
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
        return (_parent_rank >= 0);
      }

      /**
       * \brief Specifies whether this process represents a parent in the muxer.
       *
       * \returns
       * \c true, if this process is a parent, otherwise \c false.
       */
      bool is_parent() const
      {
        return !_child_ranks.empty();
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
       * \brief Returns the child communicator.
       *
       * \returns The child communicator.
       */
      const Dist::Comm* get_comm() const
      {
        return _comm;
      }

      /**
       * \brief Sets the child communicator.
       *
       * \param[in] comm_
       * The child communicator.
       */
      void set_comm(const Dist::Comm* comm_)
      {
        _comm = comm_;
      }

      /**
       * \brief Sets the parent rank and mirror for a child process.
       *
       * \param[in] parent_rank
       * The rank of the parent of this child process.
       *
       * \param[in] parent_mirror
       * The mirror of the parent patch of this child process.
       */
      void set_parent(int parent_rank, Mirror_&& parent_mirror)
      {
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
      void push_child(int child_rank, Mirror_&& child_mirror)
      {
        _child_ranks.push_back(child_rank);
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
        if(is_child())
        {
          // create parent buffer
          _parent_buffer = _parent_mirror.create_buffer(vec_tmp_);
        }

        if(is_parent())
        {
          // create child buffers
          _child_reqs.resize(_child_ranks.size());
          _child_buffers.resize(_child_ranks.size());

          for(std::size_t i(0); i < _child_ranks.size(); ++i)
            _child_buffers.at(i) = _child_mirrors.at(i).create_buffer(vec_tmp_);
        }
      }

      /**
       * \brief Sends a join operation to the parent process.
       *
       * \param[in] vec_src
       * The type-0 child vector that is to be joined on the parent process.
       */
      bool join_send(const LocalVector_& vec_src) const
      {
        XASSERT(is_child());
        XASSERT(!is_parent());
        XASSERT(_comm != nullptr);

        // pack our source vector into the parent buffer
        _parent_mirror.gather(_parent_buffer, vec_src, Index(0));

        // and set it to the parent
        _comm->send(_parent_buffer.elements(), _parent_buffer.used_elements(), _parent_rank);

        return true;
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
      bool join(const LocalVector_& vec_src, LocalVector_& vec_trg) const
      {
        // if this muxer is not a child, then this operation is a simple copy
        if(!is_child() || (_comm == nullptr))
        {
          vec_trg.copy(vec_src);
          return true;
        }

        XASSERT(is_parent());
        XASSERT(_comm != nullptr);

        // post the receives for our children
        for(std::size_t i(0); i < _child_ranks.size(); ++i)
        {
          // post receive request from child
          _child_reqs[i] = _comm->irecv(
            _child_buffers.at(i).elements(),
            _child_buffers.at(i).used_elements(),
            _child_ranks.at(i));
        }

        // pack our source vector into the parent buffer
        _parent_mirror.gather(_parent_buffer, vec_src, Index(0));

        // and set it to the parent
        _comm->send(_parent_buffer.elements(), _parent_buffer.used_elements(), _parent_rank);

        // format target vector
        vec_trg.format();

        // process all receive requests
        for(std::size_t idx; _child_reqs.wait_any(idx); )
        {
          // scatter buffers
          _child_mirrors.at(idx).scatter_axpy(vec_trg, _child_buffers.at(idx));
        }

        // okay
        return true;
      }

      /**
       * \brief Receives a split operation from the parent process.
       *
       * \param[out] vec_trg
       * The split type-1 child vector.
       *
       * \returns
       * \c true, if the split operation was performed or \c false,
       * if the split operation was cancelled by the parent.
       */
      bool split_recv(LocalVector_& vec_trg) const
      {
        XASSERT(is_child());
        XASSERT(!is_parent());
        XASSERT(_comm != nullptr);

        // receive status from parent
        int parent_split_status(1);
        _comm->recv(&parent_split_status, std::size_t(1), _parent_rank);
        if(parent_split_status != 0)
          return false;

        // receive from parent
        _comm->recv(_parent_buffer.elements(), _parent_buffer.used_elements(), _parent_rank);

        // format correction
        vec_trg.format();

        // scatter buffer
        _parent_mirror.scatter_axpy(vec_trg, _parent_buffer);

        // okay
        return true;
      }

      /**
       * \brief Cancels a split operation.
       *
       * This function may be called by a parent process to notify the child processes
       * that the split operation is not going to be completed. In this case, the
       * return value of the #split_recv call on the child processes will return \c false.
       */
      void split_cancel() const
      {
        XASSERT(is_child());
        XASSERT(is_parent());
        XASSERT(_comm != nullptr);

        // post negative status to all children
        int child_split_status(1);
        for(std::size_t i(0); i < _child_ranks.size(); ++i)
        {
          _child_reqs[i] = _comm->isend(&child_split_status, std::size_t(1), _child_ranks.at(i));
        }

        // receive status from parent (must be != 0)
        int parent_split_status(1);
        _comm->recv(&parent_split_status, std::size_t(1), _parent_rank);
        XASSERT(parent_split_status != 0);

        // wait until all sends are finished
        _child_reqs.wait_all();
      }

      /**
       * \brief Performs a split operation on the parent process.
       *
       * \param[out] vec_trg
       * The split type-1 child vector.
       *
       * \param[in] vec_src
       * The type-1 parent vector that is to be split
       *
       * \returns
       * \c true, if the split operation was performed or \c false,
       * if the split operation was cancelled by the parent.
       */
      bool split(LocalVector_& vec_trg, const LocalVector_& vec_src) const
      {
        // if this muxer is not a child, then this operation is a simple copy
        if(!is_child() || (_comm == nullptr))
        {
          vec_trg.copy(vec_src);
          return true;
        }

        XASSERT(is_parent());
        XASSERT(_comm != nullptr);

        // post positive status to all children
        int child_split_status(0);
        for(std::size_t i(0); i < _child_ranks.size(); ++i)
        {
          _child_reqs[i] = _comm->isend(&child_split_status, std::size_t(1), _child_ranks.at(i));
        }

        // receive status from parent (must be 0)
        int parent_split_status(1);
        _comm->recv(&parent_split_status, std::size_t(1), _parent_rank);
        XASSERT(parent_split_status == 0);

        // wait until all sends are finished
        _child_reqs.wait_all();

        // post the sends for our children
        for(std::size_t i(0); i < _child_ranks.size(); ++i)
        {
          auto& child_buf = _child_buffers.at(i);

          // gather vector into child buffer
          _child_mirrors.at(i).gather(child_buf, vec_src, Index(0));

          // post send request
          _child_reqs[i] = _comm->isend(child_buf.elements(), child_buf.used_elements(), _child_ranks.at(i));
        }

        // receive from parent
        _comm->recv(_parent_buffer.elements(), _parent_buffer.used_elements(), _parent_rank);

        // format correction
        vec_trg.format();

        // scatter buffer
        _parent_mirror.scatter_axpy(vec_trg, _parent_buffer);

        // wait until all sends are finished
        _child_reqs.wait_all();

        // okay
        return true;
      }
    }; // class Muxer<...>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_MUXER_HPP
