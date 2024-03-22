// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_BASE_HEADER_HPP
#define KERNEL_BASE_HEADER_HPP 1

/**
 * \file
 * \brief FEAT Kernel base header.
 *
 * This file is the base header for the FEAT kernel, which is included by all other FEAT header and source files.
 * It defines macros and data types which are frequently used in other files.
 */

// Include FEAT configuration header.
#ifndef FEAT_NO_CONFIG
#  include <feat_config.hpp>
#endif

// Make sure the DOXYGEN macro is not defined at compile-time;
// it is reserved for doxygen's preprocessor.
#ifdef DOXYGEN
#  error The DOXYGEN macro must not be defined at compile-time
#else
#  define DOXY(x)
#endif

/// \cond nodoxy
// Activate DEBUG macro if the build system tells us to do so.
#if defined(FEAT_DEBUG_MODE) && !defined(DEBUG)
#  define DEBUG 1
#endif
/// \endcond

// include compiler detection headers
#include <kernel/util/compiler_pgi.hpp>        // PGI C/C++ compiler
#include <kernel/util/compiler_intel_oneapi.hpp>   // Intel(R) OneAPI C/C++ compiler // Needs to be included before classic intel and llvm and gnu compilers
#include <kernel/util/compiler_intel.hpp>      // Intel(R) C/C++ compiler
#include <kernel/util/compiler_clang.hpp>      // Clang/LLVM Compiler.
#include <kernel/util/compiler_microsoft.hpp>  // Microsoft(R) (Visual) C/C++ compiler
#include <kernel/util/compiler_oracle.hpp>     // SunStudio/OracleStudio C/C++ compiler
// The GNU compiler must be the last one in this list, because other compilers (e.g. Intel and Clang)
// also define the __GNUC__ macro used to identify the GNU C/C++ compiler, thus leading to incorrect
// compiler detection.
#include <kernel/util/compiler_gnu.hpp>        // GNU C/C++ compiler

// hide the following block from doxygen
/// \cond nodoxy

// If the compiler does not support a 'noinline' specifier, we'll define it as an empty macro.
#ifndef NOINLINE
#define NOINLINE
#endif

// If the compiler does not support a 'force-inline' specifier, we'll define it as a simple inline.
#ifndef FORCE_INLINE
#define FORCE_INLINE inline
#endif

// If the compiler does not support a loop vectorization specifier, we'll define it as an empty macro.
#ifndef FEAT_PRAGMA_IVDEP
#define FEAT_PRAGMA_IVDEP
#endif

// If the compiler does not support disabling/restoring warnings, we'll define the corresponding
// macros as empty.
#ifndef FEAT_DISABLE_WARNINGS
#define FEAT_DISABLE_WARNINGS
#endif
#ifndef FEAT_RESTORE_WARNINGS
#define FEAT_RESTORE_WARNINGS
#endif

// Unless the compiler detection header explicitly defined 'FEAT_PRAGMA_OMP', we define it as an OpenMP
// pragma here, if the FEAT_HAVE_OMP define is set, otherwise we will define it as an empty macro.
#ifndef FEAT_PRAGMA_OMP
#ifdef FEAT_HAVE_OMP
#define FEAT_PRAGMA_OMP_HELPER(x) _Pragma(#x)
#define FEAT_PRAGMA_OMP(x) FEAT_PRAGMA_OMP_HELPER(omp x)
#else
#define FEAT_PRAGMA_OMP(x)
#endif
#endif

// If we use the quadmath library, we can use this macro to define 128-bit floating point constants.
#ifdef FEAT_HAVE_QUADMATH
#  define FEAT_F128C(x) x##Q
#else
#  define FEAT_F128C(x) x
#endif

///\endcond
// end of block hidden from doxygen

// include integer definition header from C/C++ standard library
#include <cstdint>

/**
 * \brief FEAT namespace
 */
namespace FEAT
{
  /// FEAT major version number
  static constexpr int version_major = 3;
  /// FEAT minor version number
  static constexpr int version_minor = 0;
  /// FEAT patch version number
  static constexpr int version_patch = 0;

  /**
   * \brief Index data type.
   */
#ifdef FEAT_INDEX_U32
  typedef std::uint32_t Index;
#else
  typedef std::uint64_t Index;
#endif

  /**
   * \brief Real data type.
   */
  typedef double Real;
} // namespace FEAT

#endif // KERNEL_BASE_HEADER_HPP
