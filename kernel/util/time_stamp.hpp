#pragma once
#ifndef KERNEL_UTIL_TIMESTAMP_HPP
#define KERNEL_UTIL_TIMESTAMP_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <sstream>
#include <iomanip>

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
      long long freq = 0ll;
      QueryPerformanceFrequency(&freq);
      return (freq == 0ll) ? 0.0 : (double(_counter - before._counter) / double(freq));
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
    long long elapsed_micros(const TimeStamp& before) const
    {
#if defined(FEAST_HAVE_GETTIMEOFDAY)
      return 1000000ll * (long long)(_time.tv_sec - before._time.tv_sec)
        + (long long)(_time.tv_usec - before._time.tv_usec);
#elif defined(_WIN32)
      long long freq = 0ll;
      QueryPerformanceFrequency(&freq);
      return (freq == 0ll) ? 0ll : (1000000ll * (_counter - before._counter)) / freq;
#else
      return 1000000ll * (long long)(_clock - before._clock) / (long long)CLOCKS_PER_SEC;
#endif
    }

    /**
     * \brief Formats an elapsed time in microseconds as a string.
     *
     * This function formats the given elapsed time as a string, where the format depends
     * on the \p bhours and \p bmillis parameters:
     * - If \p bhours = \c true, then the least granular unit is an hour and the string is
     *   of the format <c>h:mm:ss[.nnn]</c>.
     * - If \p bhours = \c false, then the least granular unit is a minute and the string
     *   is of the format <c>m:ss[.nnn]</c>.
     * - If \p bmillis = \c true, then the most granular unit is a millisecond and the string is
     *   of the format <c>[h:m]m:ss.nnn</c>.
     * - If \p bmillis = \c false, then the most granular unit is a second and the string is
     *   of the format <c>[h:m]m:ss</c>.
     *
     * \param[in] micros
     * The elapsed time to be formatted.
     *
     * \param[in] bhours
     * Specifies whether hours are to be printed.
     *
     * \param[in] bmillis
     * Specifies whether milliseconds are to be printed.
     *
     * \returns
     * The time elapsed as a formatted string <code>[h:m]m:ss[.nnn]</code>.
     */
    static String format_micros(long long micros, bool bhours = true, bool bmillis = true)
    {
      // check whether the time is negative
      bool neg(false);
      if(micros < 0)
      {
        neg = true;
        micros = -micros;
      }

      std::ostringstream oss;
      oss << std::setfill('0');
      if(bhours)
      {
        // write hours
        oss << (micros / 3600000000ll);
        // write minutes
        oss << ":" << std::setw(2) << ((micros / 60000000ll) % 60ll);
      }
      else
      {
        // write minutes
        oss << ":" << (micros / 60000000ll);
      }
      if(bmillis)
      {
        // write seconds
        oss << ":" << std::setw(2) << (((micros + 500) / 1000000ll) % 60ll);
        // write milli-seconds by rounding
        oss << "." << std::setw(3) << (((micros + 500) / 1000ll) % 1000ll);
      }
      else
      {
        // write seconds by rounding
        oss << ":" << std::setw(2) << (((micros + 500000) / 1000000ll) % 60ll);
      }

      // return formatted string
      String str(oss.str());
      if(neg)
        str.push_front('-');
      return str;
    }

    /**
     * \brief Return the time elapsed between two time stamps as a string.
     *
     * See TimeStamp::format_micros() for more information about the formatting options.
     *
     * \param[in] before
     * A time stamp that represents a previous moment.
     *
     * \param[in] bhours
     * Specifies whether hours are to be printed.
     *
     * \param[in] bmillis
     * Specifies whether milliseconds are to be printed.
     *
     * \returns
     * The time elapsed between the time stamps \p before and \c this as a formatted string
     * <code>h:mm:ss.nnn</code>.
     */
    String elapsed_string(const TimeStamp& before, bool bhours = true, bool bmillis = true) const
    {
      return format_micros(elapsed_micros(before), bhours, bmillis);
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
