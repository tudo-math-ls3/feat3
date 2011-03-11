#pragma once
#ifndef UTIL_COMPILER_INTEL_HPP
/// Header guard
#define UTIL_COMPILER_INTEL_HPP

/**
 * \file compiler_gnu.hpp
 *
 * \brief Compiler detection header for Intel C++ compiler.
 *
 * \author Dominik Goeddeke
 */

#if !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

// define FEAST_COMPILER_INTEL macro
// Note that __ICC is already linear sortable
#  define FEAST_COMPILER_INTEL __ICC

// map version to human-readable string and add a few
// C++0x-specific details.
#  if (__ICC >= 1200)
//   no Intel compiler up to 12.0.1.107 supports "nullptr"
//   so do not define HAVE_CPP0X_NULLPTR 1
#    define FEAST_COMPILER "Intel C/C++ compiler 12.0"
#  elif (__ICC >= 1110)
#    define FEAST_COMPILER "Intel C/C++ compiler 11.1"
#  elif (__ICC >= 1100)
#    define FEAST_COMPILER "Intel C/C++ compiler 11.0"
#  elif (__ICC >= 1010)
#    define FEAST_COMPILER "Intel C/C++ compiler 10.1"
#  elif (__ICC >= 1000)
#    define FEAST_COMPILER "Intel C/C++ compiler 10.0"
#  else
// too old to have a chance to support FEAST anyway
#    define FEAST_COMPILER "Intel C/C++ compiler"
#  endif

#endif // !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

#endif // UTIL_COMPILER_INTEL_HPP
