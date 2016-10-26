#pragma once
#ifndef GLOBAL_SYNCH_VEC_HPP
#define GLOBAL_SYNCH_VEC_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>

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
      using BufferVectorType = LAFEM::DenseVector<Mem::Main, typename VT_::DataType, typename VT_::IndexType>;

    protected:
      /// signals, whether wait was already called
      bool _finished;

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// the vector to be synchronised
      VT_& _target;
      /// our communicator
      const Dist::Comm& _comm;
      /// the vector mirrors
      const std::vector<VMT_>& _mirrors;
      /// send and receive request vectors
      Dist::RequestVector _send_reqs, _recv_reqs;
      /// send and receive buffers
      std::vector<BufferVectorType> _send_bufs, _recv_bufs;
#endif // FEAT_HAVE_MPI || DOXYGEN

    public:
      /**
       * \brief Constructor
       *
       * \param[inout] target
       * The type-0 vector to be synchronised
       *
       * \param[in] comm
       * The communicator
       *
       * \param[in] ranks
       * The neighbour ranks within the communicator
       *
       * \param[in] mirrors
       * The vector mirrors to be used for synchronisation
       */
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      SynchVectorTicket(VT_ & target, const Dist::Comm& comm, const std::vector<int>& ranks, const std::vector<VMT_> & mirrors) :
        _finished(false),
        _target(target),
        _comm(comm),
        _mirrors(mirrors)
      {
        const std::size_t n = ranks.size();

        XASSERTM(mirrors.size() == n, "invalid vector mirror count");

        _recv_reqs.reserve(n);
        _send_reqs.reserve(n);
        _recv_bufs.reserve(n);
        _send_bufs.reserve(n);

        // post receives
        for(std::size_t i(0); i < n; ++i)
        {
          // create receive buffer vector
          _recv_bufs.emplace_back(_mirrors.at(i).create_buffer_vector());
          BufferVectorType& buf = _recv_bufs.back();

          // post receive
          _recv_reqs.push_back(_comm.irecv(buf.elements(), buf.size(), ranks.at(i)));
        }

        // post sends
        for(std::size_t i(0); i < n; ++i)
        {
          // create receive buffer vector
          _send_bufs.emplace_back(_mirrors.at(i).create_buffer_vector());
          BufferVectorType& buf = _send_bufs.back();

          // gather from mirror
          _mirrors.at(i).gather_dual(buf, _target);

          // post send
          _send_reqs.push_back(_comm.isend(buf.elements(), buf.size(), ranks.at(i)));
        }
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
        for(std::size_t idx; _recv_reqs.wait_any(idx); )
        {
          // scatter the receive buffer
          _mirrors.at(idx).scatter_axpy_dual(_target, _recv_bufs.at(idx));
        }

        // wait for all sends to finish
        _send_reqs.wait_all();
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
     * \brief Synchronises a type-0 vector
     *
     * \param[inout] target
     * The type-0 vector to be synchronised
     *
     * \param[in] comm
     * The communicator
     *
     * \param[in] ranks
     * The neighbour ranks within the communicator
     *
     * \param[in] mirrors
     * The vector mirrors to be used for synchronisation
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
