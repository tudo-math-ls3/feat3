// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>

namespace FEAT
{
  /**
   * \brief FEAT Runtime management class
   *
   * \author Peter Zajac, Dirk Ribbrock
   */
  class Runtime
  {
  private:
    /// signals, if initialize was called
    static bool _initialized;

    /// signals, if finalize was called
    static bool _finalized;

  public:
    /**
     * \brief Runtime scope guard class
     */
    class ScopeGuard
    {
    public:
      /**
       * \brief Runtime scope guard constructor
       *
       * This constructor effectively calls Runtime::initialize().
       *
       * \param[in] argc, argv
       * The argument parameters of the calling \p main function.
       */
      explicit ScopeGuard(int& argc, char**& argv)
      {
        Runtime::initialize(argc, argv);
      }

      /**
       * \brief Runtime scope guard destructor
       *
       * This destructor effectively calls Runtime::finalize().
       */
      ~ScopeGuard()
      {
        Runtime::finalize();
      }
    }; // class Runtime::Guard

    /**
     * \brief FEAT initialization
     *
     * This function performs the basic initialization of the FEAT library.
     *
     * \attention
     * This function should be the first functional called in an application's
     * \p main function.
     *
     * \param[in] argc, argv
     * The argument parameters of the calling \p main function.
     */
    static void initialize(int& argc, char**& argv);

    /**
     * \brief FEAT abortion
     *
     * This function terminates this process and, in a MPI-based run, also
     * all other processes belonging to this group.
     *
     * \param[in] dump_call_stack
     * Specifies whether to dump the call-stack to stderr prior to process termination.\n
     * Note that a call-stack dump may not be available on all platforms.
     */
    [[noreturn]] static void abort(bool dump_call_stack = true);

    /**
     * \brief FEAT finalization
     *
     * This function finalizes the FEAT library.
     *
     * \attention
     * This function should be the last function called in an application's
     * \p main function.
     *
     * \note
     * Internally this functions calls the MemoryPool::finalize function.
     * To get proper warnings one memory that is still in use and not freed correctly
     * one has to make sure, that any FEAT Container has been destructed when calling the finalize method.
     * This is usually achieved by keeping the C++ main function slim and kicking off all the fancy application stuff in a separate function/method.
     * Thus (at most) every FEAT related stuff is destructed when this separate function/method ends.
     *
     * \returns
     * An exit code (<c>EXIT_SUCCESS</c>) that can be returned by the \p main function.
     */
    static int finalize();
  }; // class Runtime
} // namespace FEAT
