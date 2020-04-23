// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_COMPILER_CLANG_HPP
#define KERNEL_UTIL_COMPILER_CLANG_HPP 1

/**
 * \file compiler_clang.hpp
 *
 * \brief Compiler detection header for Clang C++ compiler.
 *
 * \author Dirk Ribbrock
 */

#if !defined(FEAT_COMPILER) && defined(__clang__)

// calc linear sortable clang version
#define _CLANG_VER (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)

#define FEAT_COMPILER_CLANG _CLANG_VER

#define FEAT_COMPILER "Clang Compiler" // __clang_version__ contains details

#if(__clang_major__ > 8)
#define FEAT_DISABLE_WARNINGS _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wunused-variable\"") \
  _Pragma("clang diagnostic ignored \"-Wunused-parameter\"") \
  _Pragma("clang diagnostic ignored \"-Wsign-compare\"") \
  _Pragma("clang diagnostic ignored \"-Wconversion\"") \
  _Pragma("clang diagnostic ignored \"-Wmismatched-tags\"") \
  _Pragma("clang diagnostic ignored \"-Wignored-qualifiers\"") \
  _Pragma("clang diagnostic ignored \"-Wcast-qual\"") \
  _Pragma("clang diagnostic ignored \"-Wdeprecated-declarations\"") \
  _Pragma("clang diagnostic ignored \"-Wshadow\"") \
  _Pragma("clang diagnostic ignored \"-Wundef\"") \
  _Pragma("clang diagnostic ignored \"-Wimplicit-fallthrough\"") \
  _Pragma("clang diagnostic ignored \"-Wcomma\"") \
  _Pragma("clang diagnostic ignored \"-Wextra-semi\"") \
  _Pragma("clang diagnostic ignored \"-Wextra-semi-stmt\"") \
  _Pragma("clang diagnostic ignored \"-Wc++98-compat-extra-semi\"")
#endif
#if(__clang_major__ >= 4) && (__clang_major__ < 8)
#define FEAT_DISABLE_WARNINGS _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wunused-variable\"") \
  _Pragma("clang diagnostic ignored \"-Wunused-parameter\"") \
  _Pragma("clang diagnostic ignored \"-Wsign-compare\"") \
  _Pragma("clang diagnostic ignored \"-Wconversion\"") \
  _Pragma("clang diagnostic ignored \"-Wmismatched-tags\"") \
  _Pragma("clang diagnostic ignored \"-Wignored-qualifiers\"") \
  _Pragma("clang diagnostic ignored \"-Wcast-qual\"") \
  _Pragma("clang diagnostic ignored \"-Wdeprecated-declarations\"") \
  _Pragma("clang diagnostic ignored \"-Wshadow\"") \
  _Pragma("clang diagnostic ignored \"-Wundef\"") \
  _Pragma("clang diagnostic ignored \"-Wimplicit-fallthrough\"") \
  _Pragma("clang diagnostic ignored \"-Wcomma\"")
#endif
#if(__clang_major__ < 4)
#define FEAT_DISABLE_WARNINGS _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wunused-variable\"") \
  _Pragma("clang diagnostic ignored \"-Wconversion\"") \
  _Pragma("clang diagnostic ignored \"-Wmismatched-tags\"") \
  _Pragma("clang diagnostic ignored \"-Wignored-qualifiers\"")
#endif

#define FEAT_RESTORE_WARNINGS _Pragma("clang diagnostic pop")

#define FEAT_IVDEP _Pragma("omp simd")

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

#define FORCE_INLINE inline __attribute__((always_inline))

#endif // !defined(FEAT_COMPILER) && defined(__clang__)

#endif // KERNEL_UTIL_COMPILER_CLANG_HPP
