#pragma once
#ifndef KERNEL_UTIL_COMPILER_INTEL_HPP
#define KERNEL_UTIL_COMPILER_INTEL_HPP 1

/**
 * \file compiler_intel.hpp
 *
 * \brief Compiler detection header for Intel C++ compiler.
 *
 * \author Dominik Goeddeke
 */

#if !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

// define FEAST_COMPILER_INTEL macro
// Note that __ICC is already linear sortable
#  define FEAST_COMPILER_INTEL __INTEL_COMPILER

// map version to human-readable string
#  if(__INTEL_COMPILER >= 1210)
#    define FEAST_COMPILER "Intel C/C++ compiler 12.1 (or newer)"
#  elif(__INTEL_COMPILER >= 1200)
#    define FEAST_COMPILER "Intel C/C++ compiler 12.0"
#  elif(__INTEL_COMPILER >= 1110)
#    define FEAST_COMPILER "Intel C/C++ compiler 11.1"
#  elif(__INTEL_COMPILER >= 1100)
#    define FEAST_COMPILER "Intel C/C++ compiler 11.0"
#  elif(__INTEL_COMPILER >= 1010)
#    define FEAST_COMPILER "Intel C/C++ compiler 10.1"
#  elif(__INTEL_COMPILER >= 1000)
#    define FEAST_COMPILER "Intel C/C++ compiler 10.0"
#  else
// too old to have a chance to support FEAST anyway
#    define FEAST_COMPILER "Intel C/C++ compiler"
#  endif

// Note: Up to version 12.1, the ICC compiler has only "experimental" C++0x support, which needs
// to be enabled explicitly via "-std=c++0x", including the following features:
//   1. static_assert: since v11.0
//   2. nullptr: since v12.1
// At least to our knowledge there is no way to detect here whether this optional has been enabled
// or not and therefore we do not define HAVE_CPP11_NULLPTR or HAVE_CPP11_STATIC_ASSERT here.

// Fortunately, at least from v10.1 on, the ICC supports the __func__ variable.
#define HAVE_CPP11_FUNC 1

#endif // !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

#endif // KERNEL_UTIL_COMPILER_INTEL_HPP
