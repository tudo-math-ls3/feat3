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

#if(__GNUC__ >= 4)
#  define FEAST_COMPILER "GNU C++ compiler 4.x.x (or newer)"
#elif(__GNUC__ >= 3)
#  define FEAST_COMPILER "GNU C++ compiler 3.x.x"
#else
// too old to compile FEAST anyway...
#  define FEAST_COMPILER "GNU C++ compiler 2.x.x (or older)"
#endif

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

#endif // !defined(FEAST_COMPILER) && defined(__GNUC__)

#endif // KERNEL_UTIL_COMPILER_GNU_HPP
