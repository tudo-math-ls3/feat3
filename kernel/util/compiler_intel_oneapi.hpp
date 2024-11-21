// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

/**
 * \file compiler_intel.hpp
 *
 * \brief Compiler detection header for Intel OneAPI C++ compiler.
 *
 * \author Dirk Ribbrock, Peter Zajac
 */

#if !defined(FEAT_COMPILER) && defined(__INTEL_LLVM_COMPILER)

// define FEAT_COMPILER_INTEL_ONEAPI macro
#define FEAT_COMPILER_INTEL_ONEAPI __INTEL_LLVM_COMPILER

// define compiler string
#define FEAT_COMPILER __VERSION__

// oneAPI compilers are based on clang, so use clang-style diagnostics

#define FEAT_DISABLE_WARNINGS _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wall\"") \
  _Pragma("clang diagnostic ignored \"-Wunknown-pragmas\"") \
  _Pragma("clang diagnostic ignored \"-Wshadow\"") \
  _Pragma("clang diagnostic ignored \"-Wunused-parameter\"") \
  _Pragma("clang diagnostic ignored \"-Wdeprecated-builtins\"") \
  _Pragma("clang diagnostic ignored \"-Wdeprecated-copy-with-user-provided-dtor\"")

#define FEAT_RESTORE_WARNINGS _Pragma("clang diagnostic pop")

#define FEAT_PRAGMA_IVDEP _Pragma("ivdep")

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

#define FORCE_INLINE inline __forceinline

#endif // !defined(FEAT_COMPILER) && defined(__INTEL_LLVM_COMPILER)
