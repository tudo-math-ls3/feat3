#pragma once
#ifndef KERNEL_BASE_HEADER_HPP
#define KERNEL_BASE_HEADER_HPP 1

/**
 * \file base_header.hpp
 * \brief Feast Kernel base header.
 * \details
 * This file is the base header for the Feast kernel, which is included by all other
 * Feast header and source files.
 * This file defines macros and data types which are frequently used in other files.
 */

// Make sure the DOXYGEN macro is not defined at compile-time;
// it is reserved for doxygen's preprocessor.
#ifdef DOXYGEN
#  error The DOXYGEN macro must not be defined at compile-time
#endif // DOXYGEN

// The DEBUG and NDEBUG macros are mutually exclusive
#if defined(DEBUG) && defined(NDEBUG)
#  error The DEBUG and NDEBUG macros must not be defined at the same time.
#endif // defined(DEBUG) && defined(NDEBUG)

// include compiler detection headers
#include <kernel/compiler_gcc.hpp>  // GNU C/C++ compiler
#include <kernel/compiler_icc.hpp>  // Intel(R) C/C++ compiler
#include <kernel/compiler_msc.hpp>  // Microsoft(R) (Visual) C/C++ compiler

// If the compiler doesn't support the C++0x nullptr, we have to define it via pre-processor.
#if !defined(HAVE_CPP0X_NULLPTR) && !defined(DOXYGEN)
#  define nullptr 0
#endif

// If the compiler doesn't support the C++0x static_assert statement, we will define an
// empty pre-processor macro for it.
#if !defined(HAVE_CPP0X_STATIC_ASSERT) && !defined(DOXYGEN)
#  define static_assert(const_expr, err_msg)
#endif

// If the compiler doesn't support the 'extern template' specifier or its usage is
// disabled by defining the FEAST_NO_EXPL_TEMPL_INST macro, both the DEF_EXPL_TEMPL_INST
// and DECL_EXPL_TEMPL_INST macros are defined as empty macros.
#if !defined(HAVE_CPP0X_EXTERN_TEMPLATE) || defined(FEAST_NO_EXPL_TEMPL_INST)
#  define DEF_EXPL_TEMPL_INST(expr)
#  define DECL_EXPL_TEMPL_INST(expr)
#else
#  define DEF_EXPL_TEMPL_INST(expr) template<> expr
#  define DECL_EXPL_TEMPL_INST(expr) extern template<> expr
#endif

// define WINDOWS macro on windows platform
#ifdef _WIN32
#  define WINDOWS 1
#endif

// Include the C header assert.h, so that the assert() function will always be defined.
#include <assert.h>

// If the NDEBUG macro is defined, ASSERT(x) will be redefined as an empty macro, otherwise
// the ASSERT(x) definition of the assert.h header will be used.
#ifdef NDEBUG
#  undef ASSERT
#  define ASSERT(x)
#endif // NDEBUG

// include STL string header
#include <string>

/**
 * \brief Feast namespace.
 */
namespace Feast
{
  // Feast version
  enum
  {
    /// Feast major version number
    version_major = 1,
    /// Feast minor version number
    version_minor = 0,
    /// Feast patch version number
    version_patch = 0
  };

  /**
   * \brief Feast String type.
   */
  typedef std::string String;

  /// signed 8-bit integer
  typedef signed char int8;

  /// unsigned 8-bit integer
  typedef unsigned char uint8;

  /// signed 16-bit integer
  typedef signed short int16;

  /// unsigned 16-bit integer
  typedef unsigned short uint16;

  /// signed 32-bit integer
  typedef signed int int32;

  /// unsigned 32-bit integer
  typedef unsigned int uint32;

  /// signed 64-bit integer
  typedef signed long long int64;

  /// unsigned 64-bit integer
  typedef unsigned long long uint64;

  /**
   * \brief Nil class definition.
   * \details
   * This is an empty tag class which may be used for templates with optional parameters.
   */
  class Nil
  {
  }; // class Nil

} // namespace Feast

#endif // KERNEL_BASE_HEADER_HPP
