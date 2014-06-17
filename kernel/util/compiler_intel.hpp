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

#  define FEAST_IVDEP _Pragma("ivdep")

// Note: The ICC 14.0.x compilers have a bug which causes the compiler to choke
// on _Pragma statements in preprocessed include files; see
// https://software.intel.com/en-us/forums/topic/515154?language=en
// Therefore, we skip implementing the FEAST_DISABLE/RESTORE_WARNINGS macros
// for these versions.
#if (__INTEL_COMPILER != 1400)

#define FEAST_DISABLE_WARNINGS _Pragma("warning(push,0)") \
  _Pragma("warning(disable:177)") \
  _Pragma("warning(disable:858)")

#define FEAST_RESTORE_WARNINGS _Pragma("warning(pop)")

#endif //  (__INTEL_COMPILER != 1400)

#endif // !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

#endif // KERNEL_UTIL_COMPILER_INTEL_HPP
