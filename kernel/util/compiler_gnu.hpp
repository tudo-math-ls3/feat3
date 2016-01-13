#pragma once
#ifndef KERNEL_UTIL_COMPILER_GNU_HPP
#define KERNEL_UTIL_COMPILER_GNU_HPP 1

/**
 * \file compiler_gnu.hpp
 *
 * \brief Compiler detection header for GNU C++ compiler.
 *
 * \author Dirk Ribbrock
 * \author Dominik Goeddeke
 */

#if !defined(FEAST_COMPILER) && defined(__GNUC__)

// calc linear sortable gcc version
#  define _GCC_VER (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

// define FEAST_COMPILER_GCC
#  define FEAST_COMPILER_GNU _GCC_VER

#if(__GNUC__ >= 5)
#  define FEAST_COMPILER "GNU C++ compiler 5.x.x (or newer)"
#elif(__GNUC__ >= 4)
#  define FEAST_COMPILER "GNU C++ compiler 4.x.x"
#else
// too old to compile FEAST anyway...
#  define FEAST_COMPILER "GNU C++ compiler"
#endif

#if(_GCC_VER >= 40900)
#  define FEAST_IVDEP _Pragma("GCC ivdep")
#else
#  define FEAST_IVDEP
#endif

#define FEAST_DISABLE_WARNINGS _Pragma("GCC diagnostic push") \
  _Pragma("GCC diagnostic ignored \"-Wunused-variable\"") \
  _Pragma("GCC diagnostic ignored \"-Wignored-qualifiers\"")

#define FEAST_RESTORE_WARNINGS _Pragma("GCC diagnostic pop")

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

#define FORCE_INLINE inline __attribute__((always_inline))

#endif // !defined(FEAST_COMPILER) && defined(__GNUC__)

#endif // KERNEL_UTIL_COMPILER_GNU_HPP
