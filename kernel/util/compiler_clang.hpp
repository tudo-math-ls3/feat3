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

#if(__clang_major__ > 3) || (clang_minor > 8)
#define FEAT_DISABLE_WARNINGS _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wunused-variable\"") \
  _Pragma("clang diagnostic ignored \"-Wconversion\"") \
  _Pragma("clang diagnostic ignored \"-Wmismatched-tags\"") \
  _Pragma("clang diagnostic ignored \"-Wignored-qualifiers\"") \
  _Pragma("clang diagnostic ignored \"-Wcomma\"")
#else
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
