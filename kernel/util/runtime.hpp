// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_RUNTIME_HPP
#define KERNEL_RUNTIME_HPP 1

#include <kernel/base_header.hpp>

namespace FEAT
{
  /// The class Runtime encapsulates various settings and functionality needed to run FEAT properly
  class Runtime
  {
  private:
    /// signals, if initialise was called
    static bool _initialised;

    /// signals, if finalise was called
    static bool _finished;

  public:
    /**
     * \brief FEAT initialisation
     *
     * This function performs the basic initialisation of the FEAT library.
     *
     * \attention
     * This function should be the first functional called in an application's
     * \p main function.
     *
     * \param[in] argc, argv
     * The argument parameters of the calling \p main function.
     */
    static void initialise(int& argc, char**& argv);

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
     * \brief FEAT finalisation
     *
     * This function finalises the FEAT library.
     *
     * \attention
     * This function should be the last function called in an application's
     * \p main function.
     *
     * \note
     * Internally this functions calls the MemoryPool::finalise function.
     * To get proper warnings one memory that is still in use and not freed correctly
     * one has to make sure, that any FEAT Container has been destructed when calling the finalise method.
     * This is usually achieved by keeping the C++ main function slim and kicking off all the fancy application stuff in a separate function/method.
     * Thus (at most) every FEAT related stuff is destructed when this separate function/method ends.
     *
     * \returns
     * An exit code (<c>EXIT_SUCCESS</c>) that can be returned by the \p main function.
     */
    static int finalise();
  };

} // namespace FEAT

#endif // KERNEL_RUNTIME_HPP
