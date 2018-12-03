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

#if !defined(FEAT_COMPILER) && defined(__GNUC__)

// calc linear sortable gcc version
#  define _GCC_VER (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

// define FEAT_COMPILER_GCC
#  define FEAT_COMPILER_GNU _GCC_VER

#if(__GNUC__ >= 9)
#  define FEAT_COMPILER "GNU C++ compiler 9.x.x (or newer)"
#elif(__GNUC__ >= 8)
#  define FEAT_COMPILER "GNU C++ compiler 8.x.x"
#elif(__GNUC__ >= 7)
#  define FEAT_COMPILER "GNU C++ compiler 7.x.x"
#elif(__GNUC__ >= 6)
#  define FEAT_COMPILER "GNU C++ compiler 6.x.x"
#elif(__GNUC__ >= 5)
#  define FEAT_COMPILER "GNU C++ compiler 5.x.x"
#elif(__GNUC__ >= 4)
#  define FEAT_COMPILER "GNU C++ compiler 4.x.x"
#else
// too old to compile FEAT anyway...
#  define FEAT_COMPILER "GNU C++ compiler"
#endif

#if(_GCC_VER >= 40900)
#  define FEAT_IVDEP _Pragma("GCC ivdep")
#else
#  define FEAT_IVDEP
#endif


#if(_GCC_VER >= 50000)
#define FEAT_DISABLE_WARNINGS _Pragma("GCC diagnostic push") \
  _Pragma("GCC diagnostic ignored \"-Wunused-variable\"") \
  _Pragma("GCC diagnostic ignored \"-Wundef\"") \
  _Pragma("GCC diagnostic ignored \"-Wparentheses\"") \
  _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"") \
  _Pragma("GCC diagnostic ignored \"-Wshadow\"") \
  _Pragma("GCC diagnostic ignored \"-Wsuggest-override\"") \
  _Pragma("GCC diagnostic ignored \"-Wdouble-promotion\"") \
  _Pragma("GCC diagnostic ignored \"-Wpedantic\"") \
  _Pragma("GCC diagnostic ignored \"-Wignored-qualifiers\"")
#else
#define FEAT_DISABLE_WARNINGS _Pragma("GCC diagnostic push") \
  _Pragma("GCC diagnostic ignored \"-Wunused-variable\"") \
  _Pragma("GCC diagnostic ignored \"-Wundef\"") \
  _Pragma("GCC diagnostic ignored \"-Wparentheses\"") \
  _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"") \
  _Pragma("GCC diagnostic ignored \"-Wshadow\"") \
  _Pragma("GCC diagnostic ignored \"-Wdouble-promotion\"") \
  _Pragma("GCC diagnostic ignored \"-Wpedantic\"") \
  _Pragma("GCC diagnostic ignored \"-Wignored-qualifiers\"")
#endif

#define FEAT_RESTORE_WARNINGS _Pragma("GCC diagnostic pop")

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

#define FORCE_INLINE inline __attribute__((always_inline))

#endif // !defined(FEAT_COMPILER) && defined(__GNUC__)

#endif // KERNEL_UTIL_COMPILER_GNU_HPP
