// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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

#if !defined(FEAT_COMPILER) && defined(__INTEL_LLVM_COMPILER)

// define FEAT_COMPILER_INTEL macro
// Note that __ICC is already linear sortable
#  define FEAT_COMPILER_INTEL __INTEL_LLVM_COMPILER

// map version to human-readable string
#  if(__INTEL_LLVM_COMPILER >= 20250000)
#    define FEAT_COMPILER "Intel OneAPI C/C++ compiler 2025.x (or newer)"
#  elif(__INTEL_LLVM_COMPILER >= 20240000)
#    define FEAT_COMPILER "Intel OneAPI C/C++ compiler 2024.x"
#  elif(__INTEL_LLVM_COMPILER >= 20230000)
#    define FEAT_COMPILER "Intel OneAPI C/C++ compiler 2023.x"
#  elif(__INTEL_LLVM_COMPILER >= 20220000)
#    define FEAT_COMPILER "Intel OneAPI C/C++ compiler 2022."
#  elif(__INTEL_LLVM_COMPILER >= 20210000)
#    define FEAT_COMPILER "Intel OneAPI C/C++ compiler 2021."
#  elif(__INTEL_LLVM_COMPILER >= 20200000)
#    define FEAT_COMPILER "Intel OneAPI C/C++ compiler 2020."
#  else
// too old to have a chance to support FEAT anyway
#    define FEAT_COMPILER "Intel OneAPI C/C++ compiler"
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

#  define FEAT_PRAGMA_IVDEP _Pragma("ivdep")


// disable warning #2196 (routine is both "inline" and "noinline") unconditionally
/// \todo evaluate, if the icc finally inlines or not the corresponding routines, marked by NOINLINE
_Pragma("warning(disable:2196)")

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

#define FORCE_INLINE inline __forceinline

#endif // !defined(FEAT_COMPILER) && defined(__INTEL_LLVM_COMPILER)

#endif // KERNEL_UTIL_COMPILER_INTEL_HPP
