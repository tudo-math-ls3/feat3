#pragma once
#ifndef KERNEL_UTIL_TIMESTAMP_HPP
#define KERNEL_UTIL_TIMESTAMP_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#include <limits>
#ifdef __linux
#  include <sys/time.h>
#else
#  error "TimeStamp is not yet supported on your platform!"
#endif


namespace FEAST
{
  class TimeStamp
  {
    private:
      /// Our time-stamp.
      timeval _time;

    public:
      /// Constructor.
      TimeStamp()
      {
        _time.tv_sec = std::numeric_limits<typeof(_time.tv_sec)>::max();
        _time.tv_usec = std::numeric_limits<typeof(_time.tv_usec)>::max();
      }

      /// Take a new time stamp.
      TimeStamp take()
      {
        gettimeofday(&_time, 0);

        return *this;
      }

      /// Returns our seconds.
      unsigned long sec() const
      {
        return _time.tv_sec;
      }

      /// Returns our useconds.
      long usec() const
      {
        return _time.tv_usec;
      }

      /// Return total time in seconds.
      double total() const
      {
        return _time.tv_sec + (_time.tv_usec / 1e6);
      }

      /**
       * Our comparison operator.
       *
       * Return true if our time-stamp has been taken earlier than the
       * other.
       *
       * \param other Another time-stamp.
       */
      bool operator< (const TimeStamp & other)
      {
        double this_time(_time.tv_sec + (_time.tv_usec / 1e6));
        double other_time(other._time.tv_sec + (other._time.tv_usec / 1e6));
        return this_time < other_time;
      }
  };
}

#endif
