#pragma once
#ifndef KERNEL_UITL_STOP_WATCH_HPP
#define KERNEL_UITL_STOP_WATCH_HPP 1

// includes, FEAST
#include <kernel/util/time_stamp.hpp>

namespace FEAST
{
  /**
   * \brief Stop-Watch class
   *
   * This class implements an incremental stop-watch.
   *
   * \author Peter Zajac
   */
  class StopWatch
  {
  private:
    TimeStamp _start_time;
    TimeStamp _stop_time;
    long long _micros;

  public:
    StopWatch() :
      _micros(0ll)
    {
    }

    /// Resets the elapsed time.
    void reset()
    {
      _micros = 0ll;
    }

    /// Starts the stop-watch.
    void start()
    {
      _start_time.stamp();
    }

    /// Stops the stop-watch and increments elapsed time.
    void stop()
    {
      _stop_time.stamp();

      // update elapsed time
      _micros += _stop_time.elapsed_micros(_start_time);
    }

    /// Returns the total elapsed time in seconds.
    double elapsed() const
    {
      return 1E-6 * double(_micros);
    }

    /// Returns the total elapsed time in micro-seconds.
    long long elapsed_micros() const
    {
      return _micros;
    }
  }; // class StopWatch
} // namespace FEAST

#endif // KERNEL_UITL_STOP_WATCH_HPP
