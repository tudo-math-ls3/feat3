#pragma once
#ifndef VISUAL_STUDIO_FEAST_CONFIG_HPP
#define VISUAL_STUDIO_FEAST_CONFIG_HPP 1

/**
 * \file
 * \brief FEAST configuration file for Visual Studio projects.
 *
 * \author Peter Zajac
 */

// Hide the internals of this file from doxygen - it might conflict with other definitions.
/// \cond nodoxy

// disable the context stack - in VS we can use a *real* debugger for that...
#define FEAST_NO_CONTEXT 1

// set the root directory - as we do not want to use environment variables, we'll set it to "."
// this works as long as the application is launched in the FEAST root directory
#define FEAST_SRC_DIR "."

// the binary dir for windows applications is called "win"
#define FEAST_BINARY_DIR "./win"

/// \endcond

#endif // VISUAL_STUDIO_FEAST_CONFIG_HPP
