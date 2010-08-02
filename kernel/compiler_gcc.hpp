#pragma once
#ifndef KERNEL_COMPILER_GCC_HPP
#define KERNEL_COMPILER_GCC_HPP 1

// Compiler detection header for GNU C++ compiler.

/// \todo detect implement compiler detection for GCC

#if !defined(FEAST_COMPILER) && defined(__GNUC__)

// calc linear sortable gcc version
#  define _GCC_VER (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

// define FEAST_COMPILER_GCC
#  define FEAST_COMPILER_MSC _GCC_VER

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

#endif // KERNEL_COMPILER_GCC_HPP
