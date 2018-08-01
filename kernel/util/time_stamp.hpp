#pragma once
#ifndef KERNEL_UTIL_TIMESTAMP_HPP
#define KERNEL_UTIL_TIMESTAMP_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/os_windows.hpp>

// includes, system
#include <sstream>
#include <iomanip>

#if defined(__linux) || defined(__unix__)
#  define FEAT_HAVE_GETTIMEOFDAY 1
#  include <sys/time.h>
#else
#  include <ctime>
#endif

namespace FEAT
{
  /**
   * Supported time string formatting.
   */
  enum class TimeFormat
  {
    h_m_s_m, /**< <c>h:mm:ss.mmm</c> */
    h_m_s, /**< <c>h:mm:ss</c> */
    m_s_m, /**< <c>m:ss.mmm</c> */
    m_s, /**< <c>m:ss</c> */
    s_m, /**< <c>s.mmm</c> */
  };

  /**
   * \brief Time stamp class
   *
   * This class is used to store time stamps and compute elapsed times. The implementation of this
   * class depends on the current platform:
   *  - For unix/linux systems, this class makes use of the \c gettimeofday() function.
   *  - For Windows systems, this class makes use of the \c QueryPerformanceCounter() function.
   *  - For other systems, this class makes use of the ANSI-C \c clock() function as a fallback implementation.
   *
   * \platformswitch Windows does not use gettimeofday
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  class TimeStamp
  {
  private:
    /// Our time-stamp.
#if defined(FEAT_HAVE_GETTIMEOFDAY)
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
#if defined(FEAT_HAVE_GETTIMEOFDAY)
      gettimeofday(&_time, 0);
#elif defined(_WIN32)
      _counter = Windows::query_performance_counter();
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
#if defined(FEAT_HAVE_GETTIMEOFDAY)
      return double(_time.tv_sec - before._time.tv_sec) + 1E-6 * double(_time.tv_usec - before._time.tv_usec);
#elif defined(_WIN32)
      long long freq = Windows::query_performance_frequency();
      return (freq == 0ll) ? 0.0 : (double(_counter - before._counter) / double(freq));
#else
      return double(_clock - before._clock) / double(CLOCKS_PER_SEC);
#endif
    }

    /**
     * \brief Calculates the time elapsed between the time stamp and now.
     *
     * \note
     * This function does \b not update (stamp) this time stamp.
     *
     * \returns
     * The time elapsed between this time stamp and now in seconds.
     */
    double elapsed_now() const
    {
      return TimeStamp().elapsed(*this);
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
#if defined(FEAT_HAVE_GETTIMEOFDAY)
      return 1000000ll * (long long)(_time.tv_sec - before._time.tv_sec)
        + (long long)(_time.tv_usec - before._time.tv_usec);
#elif defined(_WIN32)
      long long freq = Windows::query_performance_frequency();
      return (freq == 0ll) ? 0ll : (1000000ll * (_counter - before._counter)) / freq;
#else
      return 1000000ll * (long long)(_clock - before._clock) / (long long)CLOCKS_PER_SEC;
#endif
    }

    /**
     * \brief Calculates the time elapsed between the time stamp and now in microseconds.
     *
     * \note
     * This function does \b not update (stamp) this time stamp.
     *
     * \returns
     * The time elapsed between this time stamp and now in microseconds.
     */
    long long elapsed_micros_now() const
    {
      return TimeStamp().elapsed_micros(*this);
    }

    /**
     * \brief Formats an elapsed time in microseconds as a string.
     *
     * This function formats the given elapsed time as a string, where the format depends
     * on the TimeFormat parameter:
     * - For TimeFormat::h_m_s_m the string is
     *   of the format <c>h:mm:ss.mmm</c>.
     * - For TimeFormat::h_m_s the string is
     *   of the format <c>h:mm:ss</c>.
     * - For TimeFormat::m_s_m the string is
     *   of the format <c>m:ss.mmm</c>.
     * - For TimeFormat::m_s the string is
     *   of the format <c>m:ss</c>.
     * - For TimeFormat::s_m the string is
     *   of the format <c>s.mmm</c>.
     *
     * \param[in] micros
     * The elapsed time to be formatted.
     *
     * \param[in] format
     * Specifies string formatting to be used.
     *
     * \returns
     * The time elapsed as a formatted string, for example <code>h:m:ss.mmm</code>.
     */
    static String format_micros(long long micros, TimeFormat format = TimeFormat::s_m)
    {
      std::ostringstream oss;
      oss << std::setfill('0');

      // check whether the time is negative
      if(micros < 0)
      {
        oss << '-';
        micros = -micros;
      }

      switch (format)
      {
        case TimeFormat::h_m_s_m:
          oss << (micros / 3600000000ll);
          oss << ":" << std::setw(2) << ((micros / 60000000ll) % 60ll);
          oss << ":" << std::setw(2) << ((micros / 1000000ll) % 60ll);
          oss << "." << std::setw(3) << ((micros / 1000ll) % 1000ll);
          break;

        case TimeFormat::h_m_s:
          oss << (micros / 3600000000ll);
          oss << ":" << std::setw(2) << ((micros / 60000000ll) % 60ll);
          oss << ":" << std::setw(2) << ((micros / 1000000ll) % 60ll);
          break;

        case TimeFormat::m_s_m:
          oss << (micros / 60000000ll);
          oss << ":" << std::setw(2) << ((micros / 1000000ll) % 60ll);
          oss << "." << std::setw(3) << ((micros / 1000ll) % 1000ll);
          break;

        case TimeFormat::m_s:
          oss << (micros / 60000000ll);
          oss << ":" << std::setw(2) << ((micros / 1000000ll) % 60ll);
          break;

        case TimeFormat::s_m:
          oss << (micros / 1000000ll);
          oss << "." << std::setw(3) << ((micros / 1000ll) % 1000ll);
          break;

        default:
          throw InternalError("TimeFormat not supported!");
      }

      // return formatted string
      return oss.str();
    }

    /**
     * \brief Return the time elapsed between two time stamps as a string.
     *
     * See TimeStamp::format_micros() for more information about the formatting options.
     *
     * \param[in] before
     * A time stamp that represents a previous moment.
     *
     * \param[in] format
     * Specifies string formatting to be used.
     *
     * \returns
     * The time elapsed between the time stamps \p before and \c this as a formatted string,
     * for example <code>h:mm:ss.mm</code>.
     */
    String elapsed_string(const TimeStamp& before, TimeFormat format = TimeFormat::s_m) const
    {
      return format_micros(elapsed_micros(before), format);
    }

    /**
     * \brief Calculates the time elapsed between the time stamp and now as a string.
     *
     * \note
     * This function does \b not update (stamp) this time stamp.
     *
     * See TimeStamp::format_micros() for more information about the formatting options.
     *
     * \param[in] format
     * Specifies string formatting to be used.
     *
     * \returns
     * The time elapsed between this time stamp and now as a formatted string,
     * for example <code>h:mm:ss.mm</code>.
     */
    String elapsed_string_now(TimeFormat format = TimeFormat::s_m) const
    {
      return TimeStamp().elapsed_string(*this, format);
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
#if defined(FEAT_HAVE_GETTIMEOFDAY)
      return (_time.tv_sec < other._time.tv_sec) ||
        ((_time.tv_sec == other._time.tv_sec) && (_time.tv_usec < other._time.tv_usec));
#elif defined(_WIN32)
      return (_counter < other._counter);
#else
      return (_clock < other._clock);
#endif
    }
  }; // class TimeStamp
} // namespace FEAT

#endif // KERNEL_UTIL_TIMESTAMP_HPP
