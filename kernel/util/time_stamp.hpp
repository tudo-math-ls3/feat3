#pragma once
#ifndef KERNEL_UTIL_TIMESTAMP_HPP
#define KERNEL_UTIL_TIMESTAMP_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#if defined(__linux) || defined(__unix__)
#  define FEAST_HAVE_GETTIMEOFDAY 1
#  include <sys/time.h>
#elif defined(_WIN32)
// Do not include <windows.h> to avoid namespace pollution -- define the necessary prototypes by hand instead.
extern "C" int __stdcall QueryPerformanceCounter(long long int*);
extern "C" int __stdcall QueryPerformanceFrequency(long long int*);
#else
#  include <ctime>
#endif

namespace FEAST
{
  /**
   * \brief Time stamp class
   *
   * This class is used to store time stamps and compute elapsed times. The implementation of this
   * class depends on the current platform:
   *  - For unix/linux systems, this class makes use of the \c gettimeofday() function.
   *  - For Windows systems, this class makes use of the \c QueryPerformanceCounter() function.
   *  - For other systems, this class makes use of the ANSI-C \c clock() function as a fallback implementation.
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  class TimeStamp
  {
  private:
    /// Our time-stamp.
#if defined(FEAST_HAVE_GETTIMEOFDAY)
    timeval _time;
#elif defined(_WIN32)
    long long int _counter;
#else
    clock_t _clock;
#endif

  public:
    /// Constructor.
    TimeStamp()
    {
      stamp();
    }

    /**
     * \brief Stamps the current time-stamp.
     *
     * This function updates the time-stamp to the current time.
     *
     * \returns \c *this
     */
    TimeStamp& stamp()
    {
#if defined(FEAST_HAVE_GETTIMEOFDAY)
      gettimeofday(&_time, 0);
#elif defined(_WIN32)
      QueryPerformanceCounter(&_counter);
#else
      _clock = ::clock();
#endif
      return *this;
    }

    /**
     * \brief Calculates the time elapsed between two time stamps.
     *
     * \param[in] before
     * A time stamp that represents a previous moment.
     *
     * \returns
     * The time elapsed between the time stamps \p before and \c this in seconds.
     */
    double elapsed(const TimeStamp& before) const
    {
#if defined(FEAST_HAVE_GETTIMEOFDAY)
      return double(_time.tv_sec - before._time.tv_sec) + 1E-6 * double(_time.tv_usec - before._time.tv_usec);
#elif defined(_WIN32)
      long long int freq = 0;
      QueryPerformanceFrequency(&freq);
      return (freq == 0) ? 0.0 : (double(_counter - before._counter) / double(freq));
#else
      return double(_clock - before._clock) / double(CLOCKS_PER_SEC);
#endif
    }

    /**
     * \brief Calculate the time elapsed between two time stamps in microseconds.
     *
     * \param[in] before
     * A time stamp that represents a previous moment.
     *
     * \returns
     * The time elapsed between the time stamps \p before and \c this in microseconds.
     */
    long long int elapsed_micros(const TimeStamp& before) const
    {
#if defined(FEAST_HAVE_GETTIMEOFDAY)
      return 1000000ll * long long int(_time.tv_sec - before._time.tv_sec)
        + long long int(_time.tv_usec - before._time.tv_usec);
#elif defined(_WIN32)
      long long int freq = 0ll;
      QueryPerformanceFrequency(&freq);
      return (freq == 0ll) ? 0ll : (1000000ll * (_counter - before._counter) / freq);
#else
      return 1000000ll * long long int(_clock - before._clock) / long long int(CLOCKS_PER_SEC);
#endif
    }

    /**
     * \brief Comparison operator.
     *
     * \param[in] other
     * Another time-stamp.
     *
     * \returns
     * \c true, if \c this time-stamp has been taken earlier than \p other, otherwise \c false.
     */
    bool operator< (const TimeStamp & other) const
    {
#if defined(FEAST_HAVE_GETTIMEOFDAY)
      return (_time.tv_sec < other._time.tv_sec) ||
        ((_time.tv_sec == other._time.tv_sec) && (_time.tv_usec < other._time.tv_usec));
#elif defined(_WIN32)
      return (_counter < other._counter);
#else
      return (_clock < other._clock);
#endif
    }
  }; // class TimeStamp
} // namespace FEAST

#endif // KERNEL_UTIL_TIMESTAMP_HPP
