#pragma once
#ifndef UTIL_COMPILER_GNU_HPP
/// Header guard
#define UTIL_COMPILER_GNU_HPP 1

/**
 * \file compiler_gnu.hpp
 *
 * \brief Compiler detection header for GNU C++ compiler.
 *
 * \author Dirk Ribbrock
 */

#if !defined(FEAST_COMPILER) && defined(__GNUC__)

// calc linear sortable gcc version
#  define _GCC_VER (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

// define FEAST_COMPILER_GCC
#  define FEAST_COMPILER_GNU _GCC_VER

#  define FEAST_COMPILER "GNU C++ compiler"

// Now let's see what C++0x features the compiler offers
#  if (_GCC_VER > 40600)
#    define HAVE_CPP0X_NULLPTR 1
#  endif
#  if (_GCC_VER > 40300)
#    define HAVE_CPP0X_STATIC_ASSERT 1
#  endif
#  if (_GCC_VER > 40100)
#    define HAVE_CPP0X_EXTERN_TEMPLATE 1
#  endif

#endif // !defined(FEAST_COMPILER) && defined(__GNUC__)

#endif // UTIL_COMPILER_GNU_HPP
