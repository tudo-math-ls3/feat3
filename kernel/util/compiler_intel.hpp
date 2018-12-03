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

#if !defined(FEAT_COMPILER) && defined(__INTEL_COMPILER)

// define FEAT_COMPILER_INTEL macro
// Note that __ICC is already linear sortable
#  define FEAT_COMPILER_INTEL __INTEL_COMPILER

// map version to human-readable string
#  if(__INTEL_COMPILER >= 2000)
#    define FEAT_COMPILER "Intel C/C++ compiler 20.x (or newer)"
#  elif(__INTEL_COMPILER >= 1900)
#    define FEAT_COMPILER "Intel C/C++ compiler 19.x"
#  elif(__INTEL_COMPILER >= 1800)
#    define FEAT_COMPILER "Intel C/C++ compiler 18.x"
#  elif(__INTEL_COMPILER >= 1700)
#    define FEAT_COMPILER "Intel C/C++ compiler 17.x"
#  elif(__INTEL_COMPILER >= 1600)
#    define FEAT_COMPILER "Intel C/C++ compiler 16.x"
#  elif(__INTEL_COMPILER >= 1500)
#    define FEAT_COMPILER "Intel C/C++ compiler 15.x"
#  elif(__INTEL_COMPILER >= 1400)
#    define FEAT_COMPILER "Intel C/C++ compiler 14.0"
#  elif(__INTEL_COMPILER >= 1310)
#    define FEAT_COMPILER "Intel C/C++ compiler 13.1"
#  elif(__INTEL_COMPILER >= 1300)
#    define FEAT_COMPILER "Intel C/C++ compiler 13.0"
#  elif(__INTEL_COMPILER >= 1210)
#    define FEAT_COMPILER "Intel C/C++ compiler 12.1"
#  elif(__INTEL_COMPILER >= 1200)
#    define FEAT_COMPILER "Intel C/C++ compiler 12.0"
#  elif(__INTEL_COMPILER >= 1110)
#    define FEAT_COMPILER "Intel C/C++ compiler 11.1"
#  elif(__INTEL_COMPILER >= 1100)
#    define FEAT_COMPILER "Intel C/C++ compiler 11.0"
#  elif(__INTEL_COMPILER >= 1010)
#    define FEAT_COMPILER "Intel C/C++ compiler 10.1"
#  elif(__INTEL_COMPILER >= 1000)
#    define FEAT_COMPILER "Intel C/C++ compiler 10.0"
#  else
// too old to have a chance to support FEAT anyway
#    define FEAT_COMPILER "Intel C/C++ compiler"
#  endif

#  define FEAT_DISABLE_WARNINGS _Pragma("warning(push,0)") \
    _Pragma("warning(disable:177)") \
    _Pragma("warning(disable:2259)") \
    _Pragma("warning(disable:1478)") \
    _Pragma("warning(disable:1599)") \
    _Pragma("warning(disable:1944)") \
    _Pragma("warning(disable:3280)") \
    _Pragma("warning(disable:858)")

#  define FEAT_RESTORE_WARNINGS _Pragma("warning(pop)")

#  define FEAT_IVDEP _Pragma("ivdep")


// disable warning #2196 (routine is both "inline" and "noinline") unconditionally
/// \todo evaluate, if the icc finally inlines or not the corresponding routines, marked by NOINLINE
_Pragma("warning(disable:2196)")

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

#define FORCE_INLINE inline __forceinline

#endif // !defined(FEAT_COMPILER) && defined(__INTEL_COMPILER)

#endif // KERNEL_UTIL_COMPILER_INTEL_HPP
