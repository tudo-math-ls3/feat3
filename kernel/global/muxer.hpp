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

      /// parent ranks (on all children)
      std::vector<int> _parent_ranks;
      /// parent mirrors (on all children)
      std::vector<Mirror_> _parent_mirrors;
      /// parent buffers (on all children)
      mutable std::vector<BufferVectorType> _parent_buffers;
      /// parent requests (on all children)
      mutable Dist::RequestVector _parent_reqs;

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
        _comm(nullptr)
      {
      }

      /// move-constructor
      Muxer(Muxer&& other) :
        _comm(other._comm),
        _parent_ranks(std::forward<std::vector<int>>(other._parent_ranks)),
        _parent_mirrors(std::forward<std::vector<Mirror_>>(other._parent_mirrors)),
        _parent_buffers(std::forward<std::vector<BufferVectorType>>(other._parent_buffers)),
        _parent_reqs(std::forward<Dist::RequestVector>(other._parent_reqs)),
        _child_ranks(std::forward<std::vector<int>>(other._child_ranks)),
        _child_mirrors(std::forward<std::vector<Mirror_>>(other._child_mirrors)),
        _child_buffers(std::forward<std::vector<BufferVectorType>>(other._child_buffers)),
        _child_reqs(std::forward<Dist::RequestVector>(other._child_reqs))
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
        _parent_ranks = std::forward<std::vector<int>>(other._parent_ranks);
        _parent_mirrors = std::forward<std::vector<Mirror_>>(other._parent_mirrors);
        _parent_buffers = std::forward<std::vector<BufferVectorType>>(other._parent_buffers);
        _parent_reqs = std::forward<Dist::RequestVector>(other._parent_reqs);
        _child_ranks = std::forward<std::vector<int>>(other._child_ranks);
        _child_mirrors = std::forward<std::vector<Mirror_>>(other._child_mirrors);
        _child_buffers = std::forward<std::vector<BufferVectorType>>(other._child_buffers);
        _child_reqs = std::forward<Dist::RequestVector>(other._child_reqs);
        return *this;
      }

      /// conversion function
      template<typename LVT2_, typename MT2_>
      void convert(const Muxer<LVT2_, MT2_>& other)
      {
        if((void*)this == (void*)&other)
          return;

        this->_comm = other._comm;

        this->_parent_ranks = other._parent_ranks;
        this->_parent_mirrors.resize(other._parent_mirrors.size());
        this->_parent_buffers.resize(other._parent_buffers.size());

        for(std::size_t i(0); i < other._parent_ranks.size(); ++i)
        {
          this->_parent_mirrors.at(i).convert(other._parent_mirrors.at(i));
          this->_parent_buffers.at(i).convert(other._parent_buffers.at(i));
        }

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
        std::size_t b = (_parent_ranks.size() + _child_ranks.size()) * sizeof(int);
        for(const auto& pm : _parent_mirrors)
          b += pm.bytes();
        for(const auto& pb : _parent_buffers)
          b += pb.bytes();
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
        return !_parent_ranks.empty();
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
       * \brief Adds a parent rank and mirror for a child process.
       *
       * \param[in] parent_rank
       * The rank of the parent of this child process.
       *
       * \param[in] parent_mirror
       * The mirror of the parent patch of this child process.
       */
      void push_parent(int parent_rank, Mirror_&& parent_mirror)
      {
        _parent_ranks.push_back(parent_rank);
        _parent_mirrors.push_back(std::forward<Mirror_>(parent_mirror));
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
          // create parent buffers
          _parent_reqs.resize(_parent_ranks.size());
          _parent_buffers.resize(_parent_ranks.size());

          for(std::size_t i(0); i < _parent_ranks.size(); ++i)
            _parent_buffers.at(i) = _parent_mirrors.at(i).create_buffer(vec_tmp_);
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

        // post sends to all parents
        for(std::size_t i(0); i < _parent_ranks.size(); ++i)
        {
          // gather vector into buffer
          _parent_mirrors.at(i).gather(_parent_buffers.at(i), vec_src, Index(0));

          // post send
          _parent_reqs[i] = _comm->isend(
            _parent_buffers.at(i).elements(),
            _parent_buffers.at(i).used_elements(),
            _parent_ranks.at(i));
        }

        // wait for all sends to finish
        _parent_reqs.wait_all();

        // okay
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

        // post sends to all parents
        for(std::size_t i(0); i < _parent_ranks.size(); ++i)
        {
          // gather vector into buffer
          _parent_mirrors.at(i).gather(_parent_buffers.at(i), vec_src, Index(0));

          // post send
          _parent_reqs[i] = _comm->isend(
            _parent_buffers.at(i).elements(),
            _parent_buffers.at(i).used_elements(),
            _parent_ranks.at(i));
        }

        // format target vector
        vec_trg.format();

        // process all receive requests
        for(std::size_t idx; _child_reqs.wait_any(idx); )
        {
          // scatter buffers
          _child_mirrors.at(idx).scatter_axpy(vec_trg, _child_buffers.at(idx));
        }

        // wait for all sends to finish
        _parent_reqs.wait_all();

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

        // receive the statuses from all parents
        std::vector<int> parent_split_status(_parent_ranks.size(), 1);
        for(std::size_t i(0); i < _parent_ranks.size(); ++i)
        {
          _parent_reqs[i] = _comm->irecv(&parent_split_status.at(i), std::size_t(1), _parent_ranks.at(i));
        }

        // wait for all status receives to finish
        _parent_reqs.wait_all();

        // did any of our parents cancel?
        bool parent_split_cancelled = false;

        // receive all statuses, one by one
        for(std::size_t i(0); i < _parent_ranks.size(); ++i)
        {
          // check received status
          if(parent_split_status.at(i) == 0)
          {
            // okay, parent successful, so post another receive for the data
            _parent_reqs[i] = _comm->irecv(
              _parent_buffers.at(i).elements(),
              _parent_buffers.at(i).used_elements(),
              _parent_ranks.at(i));
          }
          else
          {
            // parent cancelled
            parent_split_cancelled = true;
          }
        }

        // If at least one of our parents cancelled, then we also cancel,
        // but we still need to wait for the receives from the successful
        // parents to finish:
        if(parent_split_cancelled)
        {
          // wait for successful parents
          _parent_reqs.wait_all();

          // cancelled by at least one parent
          return false;
        }

        // If we come out here, then all parents were successful, so
        // receive their data and scatter it into the target vector

        // format target vector
        vec_trg.format();

        // process all receive requests
        for(std::size_t idx; _parent_reqs.wait_any(idx); )
        {
          // scatter buffers
          _parent_mirrors.at(idx).scatter_axpy(vec_trg, _parent_buffers.at(idx));
        }

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

        // receive statuses from our parents
        std::vector<int> parent_split_status(_parent_ranks.size(), 1);
        for(std::size_t i(0); i < _parent_ranks.size(); ++i)
        {
          _parent_reqs[i] = _comm->irecv(&parent_split_status.at(i), std::size_t(1), _parent_ranks.at(i));
        }

        // post negative status to all our children
        int child_split_status(1);
        for(std::size_t i(0); i < _child_ranks.size(); ++i)
        {
          _child_reqs[i] = _comm->isend(&child_split_status, std::size_t(1), _child_ranks.at(i));
        }

        // wait for all status receives to finish
        _parent_reqs.wait_all();

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

        // receive statuses from our parents
        std::vector<int> parent_split_status(_parent_ranks.size(), 1);
        for(std::size_t i(0); i < _parent_ranks.size(); ++i)
        {
          _parent_reqs[i] = _comm->irecv(&parent_split_status.at(i), std::size_t(1), _parent_ranks.at(i));
        }

        // post positive status to all children
        int child_split_status(0);
        for(std::size_t i(0); i < _child_ranks.size(); ++i)
        {
          _child_reqs[i] = _comm->isend(&child_split_status, std::size_t(1), _child_ranks.at(i));
        }

        // wait for all status receives to finish
        _parent_reqs.wait_all();

        // did any of our parents cancel?
        bool parent_split_cancelled = false;

        // receive all statuses, one by one
        for(std::size_t i(0); i < _parent_ranks.size(); ++i)
        {
          // check received status
          if(parent_split_status.at(i) == 0)
          {
            // okay, parent successful, so post another receive for the data
            _parent_reqs[i] = _comm->irecv(
              _parent_buffers.at(i).elements(),
              _parent_buffers.at(i).used_elements(),
              _parent_ranks.at(i));
          }
          else
          {
            // parent cancelled
            parent_split_cancelled = true;
          }
        }

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


        // If at least one of our parents cancelled, then we also cancel,
        // but we still need to wait for the receives from the successful
        // parents to finish:
        if(parent_split_cancelled)
        {
          // wait for successful parents
          _parent_reqs.wait_all();

          // wait until all sends are finished
          _child_reqs.wait_all();

          // cancelled by at least one parent
          return false;
        }

        // If we come out here, then all parents were successful, so
        // receive their data and scatter it into the target vector

        // format target vector
        vec_trg.format();

        // process all receive requests
        for(std::size_t idx; _parent_reqs.wait_any(idx); )
        {
          // scatter buffers
          _parent_mirrors.at(idx).scatter_axpy(vec_trg, _parent_buffers.at(idx));
        }

        // wait until all sends are finished
        _child_reqs.wait_all();

        // okay
        return true;
      }
    }; // class Muxer<...>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_MUXER_HPP
