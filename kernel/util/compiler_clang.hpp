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

#if !defined(FEAST_COMPILER) && defined(__clang__)

// calc linear sortable clang version
#define _CLANG_VER (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)

#define FEAST_COMPILER_CLANG _CLANG_VER

#define FEAST_COMPILER "Clang Compiler" // __clang_version__ contains details

#define FEAST_DISABLE_WARNINGS _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wunused-variable\"") \
  _Pragma("clang diagnostic ignored \"-Wignored-qualifiers\"")

#define FEAST_RESTORE_WARNINGS _Pragma("clang diagnostic pop")

#  define FEAST_IVDEP _Pragma("omp simd")

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

#endif // !defined(FEAST_COMPILER) && defined(__clang__)

#endif // KERNEL_UTIL_COMPILER_CLANG_HPP
