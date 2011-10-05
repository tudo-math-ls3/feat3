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

// Note: Open64 is a fork of GCC, so defines __GNUC__ as well.
// Consequence: Also check for __OPEN64__ here to be able to
// distinguish them.
#if !defined(FEAST_COMPILER) && defined(__GNUC__) && !defined(__OPEN64__)

// calc linear sortable gcc version
#  define _GCC_VER (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

// define FEAST_COMPILER_GCC
#  define FEAST_COMPILER_GNU _GCC_VER

#  define FEAST_COMPILER "GNU C++ compiler"

// Now let's see what C++0x features the compiler offers
#ifdef __GXX_EXPERIMENTAL_CXX0X__
#  if (_GCC_VER > 40600)
#    define HAVE_CPP0X_NULLPTR 1
#  endif
#  if (_GCC_VER > 40300)
#    define HAVE_CPP0X_STATIC_ASSERT 1
#  endif
#  if (_GCC_VER > 40100)
#    define HAVE_CPP0X_EXTERN_TEMPLATE 1
#  endif
#endif // __GXX_EXPERIMENTAL_CXX0X__

#endif // !defined(FEAST_COMPILER) && defined(__GNUC__)

#endif // KERNEL_UTIL_COMPILER_GNU_HPP
