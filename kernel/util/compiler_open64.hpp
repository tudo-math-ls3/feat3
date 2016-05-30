#pragma once
#ifndef KERNEL_UTIL_COMPILER_OPEN64_HPP
#define KERNEL_UTIL_COMPILER_OPEN64_HPP 1

/**
 * \file compiler_open64.hpp
 *
 * \brief Compiler detection header for Open64 C++ compiler.
 *
 * \author Dominik Goeddeke
 */

#if !defined(FEAT_COMPILER) && defined(__OPEN64__)

// calc linear sortable open64 version
#  define _OPEN64_VER (__OPENCC__ * 10000 + __OPENCC_MINOR__ * 100 + __OPENCC_PATCHLEVEL__)

// define FEAT_COMPILER_OPEN64 macro
#  define FEAT_COMPILER_OPEN64 _OPEN64_VER

// map version to human-readable string and add a few
// C++0x-specific details.
#  if (_OPEN64_VER >= 40204)
//   Open64 does not support "nullptr" so do not define HAVE_CPP0X_NULLPTR 1
#    define FEAT_COMPILER "Open64 C/C++ compiler 4.2.4 (or greater)"
#  else
#    define FEAT_COMPILER "Open64 C/C++ compiler (unknown version)"
#  endif

#endif // !defined(FEAT_COMPILER) && defined(__OPEN64__)

#endif // KERNEL_UTIL_COMPILER_OPEN64_HPP
