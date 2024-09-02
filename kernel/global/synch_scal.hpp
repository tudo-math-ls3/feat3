// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef GLOBAL_SYNCH_SCAL_HPP
#define GLOBAL_SYNCH_SCAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>

#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Ticket class for asynchronous global operations on scalars
     *
     * Internally a cpp thread is used to process the mpi call and ensure proper processing by calling MPI_Wait ahead of the actual tickets wait call.
     *
     * \author Dirk Ribbrock, Peter Zajac
     */
    template <typename DT_>
    class SynchScalarTicket
    {
    protected:
      /// buffer containing the received data
      DT_ _r;
      /// buffer containing the send data
      DT_ _x;
      /// should we compute the sqrt of the result
      bool _sqrt;
      /// holds our mpi execution toe
      double _mpi_exec;
      /// holds our mpi reduction wait toe
      double _mpi_wait;
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
#ifdef FEAT_MPI_THREAD_MULTIPLE
      /// mutex for our cv
      std::mutex _mutex;
      /// did the thread call allreduce already?
      std::condition_variable _cv_allreduce_called;
      /// predicate for our cv
      bool _flag_allreduce_called;
      /// thread handle
      std::thread _thread;
#else // no MPI_THREAD_MULTIPLE
      /// Our request for the corresponding iallreduce mpi call
      Dist::Request _req;
#endif // MPI_THREAD_MULTIPLE
#endif //(defined(FEAT_HAVE_MPI) && defined(FEAT_MPI_THREAD_MULTIPLE)) || defined(DOXYGEN)
      /// signals, whether wait was already called
      bool _finished;

    public:
      /// standard constructor
      SynchScalarTicket() :
        _r(DT_(0)),
        _x(DT_(0)),
        _sqrt(false),
        _mpi_exec(0.0),
        _mpi_wait(0.0),
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
#ifdef FEAT_MPI_THREAD_MULTIPLE
        _mutex(),
        _cv_allreduce_called(),
        _flag_allreduce_called(false),
        _thread(),
#else // no MPI_THREAD_MULTIPLE
        _req(),
#endif // MPI_THREAD_MULTIPLE
#endif //(defined(FEAT_HAVE_MPI) && defined(FEAT_MPI_THREAD_MULTIPLE)) || defined(DOXYGEN)
        _finished(true)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] x
       * The value to be synchronized.
       *
       * \param[in] comm
       * The communicator to be used for synchronization.
       *
       * \param[in] op
       * The reduction operation to be applied.
       *
       * \param[in] sqrt
       * Specifies whether to apply the square-root onto the reduction result.
       */
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      explicit SynchScalarTicket(DT_ x, const Dist::Comm& comm, const Dist::Operation& op, bool sqrt = false) :
        _r(DT_(0)),
        _x(x),
        _sqrt(sqrt),
        _mpi_exec(double(0)),
        _mpi_wait(double(0)),
#ifdef FEAT_MPI_THREAD_MULTIPLE
        _flag_allreduce_called(false),
        _thread(_wait_function, std::cref(_x), std::cref(comm), std::cref(op), std::ref(_r), std::ref(_mpi_exec), std::ref(_mpi_wait),
            std::ref(_mutex), std::ref(_cv_allreduce_called), std::ref(_flag_allreduce_called)),
#else
        _req(),
#endif
        _finished(false)
      {
#ifdef FEAT_MPI_THREAD_MULTIPLE
        std::unique_lock<std::mutex> l(_mutex);
        _cv_allreduce_called.wait(l, [this]() {return _flag_allreduce_called == true; });
#else // no FEAT_MPI_THREAD_MULTIPLE
        TimeStamp ts_start;
        _req = comm.iallreduce(&_x, &_r, std::size_t(1), op);
        Statistics::add_time_mpi_execute_reduction(ts_start.elapsed_now());
#endif // FEAT_MPI_THREAD_MULTIPLE
      }
#else // non-MPI version
      explicit SynchScalarTicket(const DT_ & in, const Dist::Comm&, const Dist::Operation&, bool sqrt = false) :
        _r(in),
        _x(in),
        _sqrt(sqrt),
        _finished(false)
      {
      }
#endif // FEAT_HAVE_MPI || DOXYGEN

      /// Unwanted copy constructor: Do not implement!
      SynchScalarTicket(const SynchScalarTicket &) = delete;
      /// Unwanted copy assignment operator: Do not implement!
      SynchScalarTicket & operator=(const SynchScalarTicket &) = delete;

      /// move constructor
      SynchScalarTicket(SynchScalarTicket&& other) :
        _r(other._r),
        _x(other._x),
        _sqrt(other._sqrt),
        _mpi_exec(other._mpi_exec),
        _mpi_wait(other._mpi_wait),
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
#ifdef FEAT_MPI_THREAD_MULTIPLE
        _mutex(),
        _cv_allreduce_called(),
        _flag_allreduce_called(other._flag_allreduce_called),
        _thread(std::move(other._thread)),
#else // no MPI_THREAD_MULTIPLE
        _req(std::move(other._req)),
#endif // MPI_THREAD_MULTIPLE
#endif //(defined(FEAT_HAVE_MPI) && defined(FEAT_MPI_THREAD_MULTIPLE)) || defined(DOXYGEN)
        _finished(other._finished)
      {
        #ifdef FEAT_MPI_THREAD_MULTIPLE
        XASSERTM(other._mutex.try_lock(), "Other is still locked by other thread!");
        XASSERTM(other._flag_allreduce_called, "Allreduce not called!");
        #endif // FEAT_MPI_THREAD_MULTIPLE
        other._finished = true;
      }

      /// move-assign operator
      SynchScalarTicket& operator=(SynchScalarTicket&& other)
      {
        if(this == &other)
          return *this;

        _r = other._r;
        _x = other._x;
        _sqrt = other._sqrt;
        _mpi_exec = other._mpi_exec;
        _mpi_wait = other._mpi_wait;
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
#ifdef FEAT_MPI_THREAD_MULTIPLE
        XASSERTM(other._mutex.try_lock(), "Other is still locked by other thread!");
        XASSERTM(other._flag_allreduce_called, "Allreduce not called!");
        _flag_allreduce_called = other._flag_allreduce_called;
        _thread = std::move(other._thread);
#else // no MPI_THREAD_MULTIPLE
        _req = std::move(other._req);
#endif // MPI_THREAD_MULTIPLE
#endif //(defined(FEAT_HAVE_MPI) && defined(FEAT_MPI_THREAD_MULTIPLE)) || defined(DOXYGEN)
        _finished = other._finished;

        other._finished = true;
        return *this;
      }

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
#ifdef FEAT_MPI_THREAD_MULTIPLE
        _thread.join();
        Statistics::add_time_mpi_execute_reduction(_mpi_exec);
        Statistics::add_time_mpi_wait_reduction(_mpi_wait);
#else // no FEAT_MPI_THREAD_MULTIPLE
        TimeStamp ts_start;
        _req.wait();
        Statistics::add_time_mpi_wait_reduction(ts_start.elapsed_now());
#endif // FEAT_MPI_THREAD_MULTIPLE
#endif // FEAT_HAVE_MPI
        _finished = true;
        return (_sqrt ?  Math::sqrt(_r) : _r);
      }

      /// Destructor
      ~SynchScalarTicket()
      {
        XASSERT(_finished);
      }

    private:
#ifdef FEAT_HAVE_MPI
      static void _wait_function(const DT_ & x, const Dist::Comm& comm, const Dist::Operation & op, DT_ & r, double & mpi_exec, double & mpi_wait, std::mutex & mutex,
          std::condition_variable & cv_allreduce_called, bool & flag_allreduce_called)
      {
        TimeStamp ts_start;
        std::unique_lock<std::mutex> l(mutex);
        Dist::Request req = comm.iallreduce(&x, &r, std::size_t(1), op);
        flag_allreduce_called = true;
        l.unlock();
        cv_allreduce_called.notify_all();
        mpi_exec = ts_start.elapsed_now();
        ts_start.stamp();
        req.wait();
        mpi_wait = ts_start.elapsed_now();
      }
#endif // FEAT_HAVE_MPI
    }; // class SynchScalarTicket
  } // namespace Global
} // namespace FEAT

#endif // GLOBAL_SYNCH_SCAL_HPP
