// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef GLOBAL_SYNCH_VEC_HPP
#define GLOBAL_SYNCH_VEC_HPP 1

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
      /// the buffer vector type (possibly in device memory)
      using BufferType = LAFEM::DenseVector<typename VT_::MemType, typename VT_::DataType, typename VT_::IndexType>;

      /// the buffer vector type in main memory
      using BufferMain = LAFEM::DenseVector<Mem::Main, typename VT_::DataType, typename VT_::IndexType>;

    protected:
      /// signals, whether wait was already called
      bool _finished;

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// the vector to be synchronized
      VT_& _target;
      /// our communicator
      const Dist::Comm& _comm;
      /// the vector mirrors
      const std::vector<VMT_>& _mirrors;
      /// send and receive request vectors
      Dist::RequestVector _send_reqs, _recv_reqs;
      /// send and receive buffers
      std::vector<BufferMain> _send_bufs, _recv_bufs;
#endif // FEAT_HAVE_MPI || DOXYGEN

    public:
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
        _target(target),
        _comm(comm),
        _mirrors(mirrors)
      {
        TimeStamp ts_start;
        const std::size_t n = ranks.size();

        XASSERTM(mirrors.size() == n, "invalid vector mirror count");

        // post receives
        _recv_reqs.reserve(n);
        _recv_bufs.resize(n);
        for(std::size_t i(0); i < n; ++i)
        {
          // create buffer vector in main memory
          _recv_bufs.at(i) = BufferMain(_mirrors.at(i).buffer_size(target), LAFEM::Pinning::disabled);

          // post receive
          _recv_reqs.push_back(_comm.irecv(_recv_bufs.at(i).elements(), _recv_bufs.at(i).size(), ranks.at(i)));
        }

        // post sends
        _send_reqs.reserve(n);
        _send_bufs.resize(n);
        for(std::size_t i(0); i < n; ++i)
        {
          // create buffer in device memory
          BufferType buffer(_mirrors.at(i).buffer_size(target), LAFEM::Pinning::disabled);

          // gather from mirror
          _mirrors.at(i).gather(buffer, _target);

          // convert buffer to main memory
          _send_bufs.at(i).convert(buffer);

          // post send
          _send_reqs.push_back(_comm.isend(_send_bufs.at(i).elements(), _send_bufs.at(i).size(), ranks.at(i)));
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
          // convert buffer to device memory
          BufferType buffer;
          buffer.convert(_recv_bufs.at(idx));

          // scatter the receive buffer
          _mirrors.at(idx).scatter_axpy(_target, buffer);
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
        XASSERT(_finished);
      }
    }; // class SynchVectorTicket

    /**
     * \brief Synchronizes a type-0 vector
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
    template<typename VT_, typename VMT_>
    void synch_vector(VT_& target, const Dist::Comm& comm, const std::vector<int>& ranks, const std::vector<VMT_>& mirrors)
    {
      SynchVectorTicket<VT_, VMT_> ticket(target, comm, ranks, mirrors);
      ticket.wait();
    }
  } // namespace Global
} // namespace FEAT

#endif // GLOBAL_SYNCH_VEC_HPP
