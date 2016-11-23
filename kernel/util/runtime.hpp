#pragma once
#ifndef KERNEL_RUNTIME_HPP
#define KERNEL_RUNTIME_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/property_map.hpp>

namespace FEAT
{
  /// The class Runtime encapsulates various settings and functionality needed to run FEAT properly
  class Runtime
  {
    private:

      /// Feat's global property map, reads in feat.ini at startup
      static PropertyMap _global_property_map;

      /// signals, if initialise was called
      static bool _initialised;

      /// signals, if finalise was called
      static bool _finished;

    public:

    /// Returns global property map, containing initial FEAT configuration, read in from feat.ini file
    static PropertyMap * global_property();

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
    static void abort(bool dump_call_stack = true);

    /**
     * \brief FEAT finalisation
     *
     * This function finalises the FEAT library.
     *
     * \attention
     * This function should be the last function called in an application's
     * \p main function.
     *
     * \returns
     * An exit code (<c>EXIT_SUCCESS</c>) that can be returned by the \p main function.
     */
    static int finalise();
  };

} // namespace FEAT

#endif // KERNEL_RUNTIME_HPP
