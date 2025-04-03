// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/lafem/dense_vector.hpp>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Ticket class for asynchronous global operations on vectors
     *
     * \todo statistics
     *
     * \author Dirk Ribbrock, Peter Zajac
     */
    template <typename VT_, typename VMT_>
    class SynchVectorTicket
    {
    public:
      /// the buffer vector type
      using BufferType = LAFEM::DenseVector<typename VT_::DataType, typename VT_::IndexType>;

    protected:
      /// signals, whether wait was already called
      bool _finished;

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// the vector to be synchronized
      VT_* _target;
      /// our communicator
      const Dist::Comm* _comm;
      /// the vector mirrors
      const std::vector<VMT_>* _mirrors;
      /// send and receive request vectors
      Dist::RequestVector _send_reqs, _recv_reqs;
      /// send and receive buffers
      std::vector<BufferType> _send_bufs, _recv_bufs;
#endif // FEAT_HAVE_MPI || DOXYGEN

    public:
      /// default constructor
      SynchVectorTicket() :
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
        _finished(true),
        _target(nullptr),
        _comm(nullptr),
        _mirrors(nullptr),
        _send_reqs(),
        _recv_reqs(),
        _send_bufs(),
        _recv_bufs()
#else
        _finished(true)
#endif // FEAT_HAVE_MPI || DOXYGEN
      {
      }

      /**
       * \brief Constructor
       *
       * \param[inout] target
       * The type-0 vector to be synchronized
       *
       * \param[in] comm
       * The communicator
       *
       * \param[in] ranks
       * The neighbor ranks within the communicator
       *
       * \param[in] mirrors
       * The vector mirrors to be used for synchronization
       */
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      SynchVectorTicket(VT_ & target, const Dist::Comm& comm, const std::vector<int>& ranks, const std::vector<VMT_> & mirrors) :
        _finished(false),
        _target(&target),
        _comm(&comm),
        _mirrors(&mirrors)
      {
        TimeStamp ts_start;
        const std::size_t n = ranks.size();

        XASSERTM(_mirrors->size() == n, "invalid vector mirror count");

        // post receives
        _recv_reqs.reserve(n);
        _recv_bufs.resize(n);
        for(std::size_t i(0); i < n; ++i)
        {
          // create buffer vector in main memory
          _recv_bufs.at(i) = BufferType(_mirrors->at(i).buffer_size(*_target));

          // post receive
          _recv_reqs.push_back(_comm->irecv(_recv_bufs.at(i).elements(), _recv_bufs.at(i).size(), ranks.at(i)));
        }

        // post sends
        _send_reqs.reserve(n);
        _send_bufs.resize(n);
        for(std::size_t i(0); i < n; ++i)
        {
          // create buffer in device memory
          _send_bufs.at(i) = BufferType(_mirrors->at(i).buffer_size(*_target));

          // gather from mirror
          _mirrors->at(i).gather(_send_bufs.at(i), *_target);

          // post send
          _send_reqs.push_back(_comm->isend(_send_bufs.at(i).elements(), _send_bufs.at(i).size(), ranks.at(i)));
        }

        Statistics::add_time_mpi_execute_blas2(ts_start.elapsed_now());
      }
#else // non-MPI version
      SynchVectorTicket(VT_ &, const Dist::Comm&, const std::vector<int>& ranks, const std::vector<VMT_> &) :
        _finished(false)
      {
        XASSERT(ranks.empty());
      }
#endif // FEAT_HAVE_MPI
      /// Unwanted copy constructor: Do not implement!
      SynchVectorTicket(const SynchVectorTicket &) = delete;
      /// Unwanted copy assignment operator: Do not implement!
      SynchVectorTicket & operator=(const SynchVectorTicket &) = delete;

      /// move constructor
      SynchVectorTicket(SynchVectorTicket&& other) :
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
        _finished(other._finished),
        _target(other._target),
        _comm(other._comm),
        _mirrors(other._mirrors),
        _send_reqs(std::forward<Dist::RequestVector>(other._send_reqs)),
        _recv_reqs(std::forward<Dist::RequestVector>(other._recv_reqs)),
        _send_bufs(std::forward<std::vector<BufferType>>(other._send_bufs)),
        _recv_bufs(std::forward<std::vector<BufferType>>(other._recv_bufs))
      {
        other._finished = true;
        other._comm = nullptr;
        other._target = nullptr;
        other._mirrors = nullptr;
      }
#else
        _finished(other._finished)
      {
        other->_finished = true;
      }
#endif // FEAT_HAVE_MPI

      /// move-assign operator
      SynchVectorTicket& operator=(SynchVectorTicket&& other)
      {
        if(this == &other)
          return *this;

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
        _finished = other._finished;
        _target = other._target;
        _comm = other._comm;
        _mirrors = other._mirrors;
        _send_reqs = std::forward<Dist::RequestVector>(other._send_reqs);
        _recv_reqs = std::forward<Dist::RequestVector>(other._recv_reqs);
        _send_bufs = std::forward<std::vector<BufferType>>(other._send_bufs);
        _recv_bufs = std::forward<std::vector<BufferType>>(other._recv_bufs);

        other->_finished = true;
        other->_comm = nullptr;
        other->_target = nullptr;
        other->_mirrors = nullptr;
#else
        _finished = other._finished;
        other._finished = true;
#endif // FEAT_HAVE_MPI

        return *this;
      }

      /**
       * \brief wait method
       *
       * wait for completion barrier
       * on return, the target vector has been updated
       */
      void wait()
      {
        XASSERTM(!_finished, "ticket was already completed by a wait call");

#ifdef FEAT_HAVE_MPI
        TimeStamp ts_start;

        // process all pending receives
        for(std::size_t idx(0u); _recv_reqs.wait_any(idx); )
        {
          // scatter the receive buffer
          _mirrors->at(idx).scatter_axpy(*_target, _recv_bufs.at(idx));
        }

        // wait for all sends to finish
        _send_reqs.wait_all();

        Statistics::add_time_mpi_wait_blas2(ts_start.elapsed_now());
#endif // FEAT_HAVE_MPI

        _finished = true;
      }

      /// Destructor
      ~SynchVectorTicket()
      {
        XASSERTM(_finished, "trying to destroy an unfinished SynchVectorTicket");
      }
    }; // class SynchVectorTicket
  } // namespace Global
} // namespace FEAT
