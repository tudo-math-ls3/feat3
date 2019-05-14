// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_THREAD_HPP
#define KERNEL_UTIL_THREAD_HPP 1

#include <kernel/base_header.hpp>

#include <condition_variable>
#include <mutex>
#include <thread>

namespace FEAT
{
  /**
   * \brief Thread fence synchronization utility class
   *
   * This class implements a "fence" for manual thread synchronization,
   * which is characterized by the following properties:
   * - The fence can be either open or closed.
   * - If the fence is closed, any thread calling the 'wait' function
   *   will block until some other thread opens the fence.
   * - If the fence is open, any thread calling the 'wait' function
   *   will return immediately without blocking.
   * - Any thread can open the fence by calling the 'open' function,
   *   which also notifies and unblocks all waiting threads.
   *
   * In addition to the above properties, this class also stores an
   * additional internal boolean state variable, which is specified by
   * the open() function and is returned by the wait() function.
   *
   * \author Peter Zajac
   */
  class ThreadFence
  {
  private:
    /// the internal mutex
    std::mutex _mtx;
    /// the internal condition variable
    std::condition_variable _cvar;
    /// current state of the fence: open or closed
    bool _open;
    /// additional state variable
    bool _okay;

  public:
    /// constructor
    explicit ThreadFence() : _open(false), _okay(false) {}

    /// delete copy-constructor
    ThreadFence(const ThreadFence&) = delete;
    /// delete copy-assignment operator
    ThreadFence& operator=(const ThreadFence&) = delete;

    /**
     * \brief Wait for the fence to be opened
     *
     * This function causes the calling thread to wait for
     * the fence to be opened by another thread.
     * If the fence is already open at the time when this
     * function is called, then this function returns
     * immediately without blocking the calling thread.
     *
     * \returns
     * The internal state as set by the \c open() function.
     */
    bool wait()
    {
      std::unique_lock<std::mutex> lock(_mtx);
      // we must to use a while loop here to correctly handle spurious wake-ups
      while(!_open)
        _cvar.wait(lock);
      return _okay;
    }

    /**
     * \brief Open the fence and notify all waiting threads
     *
     * This function opens the fence, set the internal state
     * and notifies all waiting threads.
     *
     * \param[in] okay
     * The state that is to be returned by the wait() function.
     */
    void open(bool okay = true)
    {
      std::unique_lock<std::mutex> lock(_mtx);
      _open = true;
      _okay = okay;
      _cvar.notify_all();
    }

    /**
     * \brief Close the fence and reset internal state.
     */
    void close()
    {
      std::unique_lock<std::mutex> lock(_mtx);
      _open = _okay = false;
    }
  }; // class ThreadFence
} // namespace FEAT

#endif // KERNEL_UTIL_THREAD_HPP
