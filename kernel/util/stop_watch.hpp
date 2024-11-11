// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/util/time_stamp.hpp>

namespace FEAT
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
    bool _running;
    long long _micros;

  public:
    StopWatch() :
      _running(false),
      _micros(0ll)
    {
    }

    /// Resets the elapsed time.
    void reset()
    {
      _running = false;
      _micros = 0ll;
    }

    /// Starts the stop-watch.
    void start()
    {
      if(!_running)
        _start_time.stamp();
      _running = true;
    }

    /// Stops the stop-watch and increments elapsed time.
    void stop()
    {
      if(_running)
      {
        _stop_time.stamp();

        // update elapsed time
        _micros += _stop_time.elapsed_micros(_start_time);
      }
      _running = false;
    }

    /// Returns \c true if the stop-watch is currently running.
    bool running() const
    {
      return _running;
    }

    /// Returns the total elapsed time in seconds.
    double elapsed() const
    {
      return 1E-6 * double(elapsed_micros());
    }

    /// Returns the total elapsed time in micro-seconds.
    long long elapsed_micros() const
    {
      if(!_running)
        return _micros;
      else
      {
        // stop-watch is currently running, so compute current time
        return _micros + TimeStamp().elapsed_micros(_start_time);
      }
    }

    /**
     * \brief Return the time elapsed in the stop-watch.
     *
     * See TimeStamp::format_micros() for more information about the formatting options.
     *
     * \param[in] format
     * Specifies the output string format to be used.
     *
     * \returns
     * The time elapsed between in the stop-watch as a formatted string <code>h:mm:ss.nnn</code>.
     */
    String elapsed_string(TimeFormat format = TimeFormat::s_m) const
    {
      return TimeStamp::format_micros(elapsed_micros(), format);
    }

    /**
     * \brief Returns the formatted percentage that this stop watch elapsed time relative to another total runtime stop watch
     *
     * \param[in] total_watch
     * The stop watch that represents the total runtime
     *
     * \param[in] prefix
     * The prefix string
     *
     * \param[in] postfix
     * The postfix string
     *
     * \returns The formatted percentage of this stop watch elapsed time relative to watch_total
     */
    String percent_of(const StopWatch& total_watch, const String& prefix = " [", const String& postfix = "% ]") const
    {
      return prefix + stringify_fp_fix(100.0 * this->elapsed() / total_watch.elapsed(), 3, 7) + postfix;
    }
  }; // class StopWatch
} // namespace FEAT
