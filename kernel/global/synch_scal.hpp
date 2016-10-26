#pragma once
#ifndef GLOBAL_SYNCH_SCAL_HPP
#define GLOBAL_SYNCH_SCAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Ticket class for asynchronous global operations on scalars
     *
     * \author Dirk Ribbrock, Peter Zajac
     */
    template <typename DT_>
    class SynchScalarTicket
    {
    protected:
      /// buffer containing the send data
      DT_ _x;
      /// buffer containing the received data
      DT_ _r;
      /// should we compute the sqrt of the result
      bool _sqrt;
      /// signals, whether wait was already called
      bool _finished;

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// our communicator
      const Dist::Comm& _comm;
      /// Our request for the corresponding iallreduce mpi call
      Dist::Request _req;
#endif // FEAT_HAVE_MPI || DOXYGEN

    public:
      /**
       * \brief Constructor
       *
       * \param[in] x
       * The value to be synchronised.
       *
       * \param[in] comm
       * The communicator to be used for synchronisation.
       *
       * \param[in] op
       * The reduction operation to be applied.
       *
       * \param[in] sqrt
       * Specifies whether to apply the square-root onto the reduction result.
       */
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      explicit SynchScalarTicket(DT_ x, const Dist::Comm& comm, const Dist::Operation& op, bool sqrt = false) :
        _x(x),
        _r(),
        _sqrt(sqrt),
        _finished(false),
        _comm(comm),
        _req()
      {
        TimeStamp ts_start;
        _req = _comm.iallreduce(&_x, &_r, std::size_t(1), op);
        Statistics::add_time_mpi_execute(ts_start.elapsed_now());
      }
#else // non-MPI version
      explicit SynchScalarTicket(const DT_ & in, const Dist::Comm&, const Dist::Operation&, bool sqrt = false) :
        _x(in),
        _r(in),
        _sqrt(sqrt),
        _finished(false)
      {
      }
#endif // FEAT_HAVE_MPI || DOXYGEN

      /// Unwanted copy constructor: Do not implement!
      SynchScalarTicket(const SynchScalarTicket &) = delete;
      /// Unwanted copy assignment operator: Do not implement!
      SynchScalarTicket & operator=(const SynchScalarTicket &) = delete;

      /**
       * \brief wait method
       *
       * wait for completion barrier
       *
       * \returns the accumulated data
       */
      DT_ wait()
      {
        XASSERTM(!_finished, "ticket was already completed by a wait call");

#ifdef FEAT_HAVE_MPI
        TimeStamp ts_start;

        _req.wait();

        Statistics::add_time_mpi_wait_reduction(ts_start.elapsed_now());
#endif // FEAT_HAVE_MPI

        _finished = true;

        return (_sqrt ?  Math::sqrt(_r) : _r);
      }

      /// Destructor
      ~SynchScalarTicket()
      {
        XASSERT(_finished);
      }
    }; // class SynchScalarTicket

    /**
     * \brief Synchronises a scalar value by applying a reduction operation
     *
     * \param[in] x
     * The value to be synchronised.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation.
     *
     * \param[in] op
     * The reduction operation to be applied.
     *
     * \param[in] sqrt
     * Specifies whether to apply the square-root onto the reduction result.
     *
     * \returns
     * The synchronised value.
     */
    template<typename DT_>
    DT_ synch_scalar(DT_ x, const Dist::Comm& comm, const Dist::Operation& op, bool sqrt = false)
    {
      SynchScalarTicket<DT_> ticket(x, comm, op, sqrt);
      return ticket.wait();
    }
  } // namespace Global
} // namespace FEAT

#endif // GLOBAL_SYNCH_SCAL_HPP
