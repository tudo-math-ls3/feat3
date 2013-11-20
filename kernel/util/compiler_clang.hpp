#pragma once
#ifndef KERNEL_UTIL_COMPILER_CLANG_HPP
#define KERNEL_UTIL_COMPILER_CLANG_HPP 1

/**
 * \file compiler_gnu.hpp
 *
 * \brief Compiler detection header for Clang C++ compiler.
 *
 * \author Dirk Ribbrock
 */

#if !defined(FEAST_COMPILER) && defined(__clang__)

// calc linear sortable clang version
#  define _CLANG_VER (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)

#  define FEAST_COMPILER_CLANG _CLANG_VER

#  define FEAST_COMPILER "Clang Compiler" // __clang_version__ contains details

// define the noinline specifier
#define NOINLINE __attribute__((noinline))

// Now claim, that all C++11 features are offered
#define HAVE_CPP11_NULLPTR 1
#define HAVE_CPP11_STATIC_ASSERT 1
#define HAVE_CPP11_EXTERN_TEMPLATE 1
#define HAVE_CPP11_SMART_POINTER 1
#define HAVE_CPP11_FUNC 1

#endif // !defined(FEAST_COMPILER) && defined(__clang__)

#endif // KERNEL_UTIL_COMPILER_CLANG_HPP
