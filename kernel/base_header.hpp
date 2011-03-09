#pragma once
#ifndef KERNEL_BASE_HEADER_HPP
/// Header guard
#define KERNEL_BASE_HEADER_HPP 1

#include <feast_config.hpp>

/**
* \brief FEAST Kernel base header.
*
* This file is the base header for the FEAST kernel, which is included by all other FEAST header and source files.
* It defines macros and data types which are frequently used in other files.
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

// Assure that DEBUG or NDEBUG is defined
// In consequence, if DEBUG and NDEBUG are note defined, NDEBUG is defined.
#ifndef DEBUG
#  define NDEBUG 1
#endif
#if !defined (NDEBUG) && !defined(DEBUG)
#  define DEBUG 1
#endif

// include compiler detection headers
#include <kernel/util/compiler_gcc.hpp>  // GNU C/C++ compiler
#include <kernel/util/compiler_icc.hpp>  // Intel(R) C/C++ compiler
#include <kernel/util/compiler_msc.hpp>  // Microsoft(R) (Visual) C/C++ compiler

// If the compiler doesn't support the C++0x nullptr, we have to define it via pre-processor.
#if !defined(HAVE_CPP0X_NULLPTR) && !defined(DOXYGEN)
#  define nullptr 0
#endif

// define __PRETTY_FUNCTION if not defined
#ifndef __PRETTY_FUNCTION__
#  ifndef DOXYGEN
#    define __PRETTY_FUNCTION__ __FUNCTION__
#  endif
#endif

// If the compiler doesn't support the 'extern template' specifier or its usage is
// disabled by defining the FEAST_NO_EXPL_TEMPL_INST macro, both the DEF_EXPL_TEMPL_INST
// and DECL_EXPL_TEMPL_INST macros are defined as empty macros.
#if !defined(DOXYGEN)
#  if !defined(HAVE_CPP0X_EXTERN_TEMPLATE) || defined(FEAST_NO_EXPL_TEMPL_INST)
#    define DEF_EXPL_TEMPL_INST(expr)
#    define DECL_EXPL_TEMPL_INST(expr)
#  else
#    define DEF_EXPL_TEMPL_INST(expr) template<> expr
#    define DECL_EXPL_TEMPL_INST(expr) extern template<> expr
#  endif
#endif

// define WINDOWS macro on windows platform
#ifdef _WIN32
#  define WINDOWS 1
#endif

// let Doxygen ignore the following block
//\cond
// In order to increase performance one can ignore the status object in some MPI functions
// (see Section 3.2.6 in the MPI 2.2 standard). In debug mode, we use a real status object, otherwise
// the special MPI status MPI_STATUS_IGNORE. Requires that the MPI_Status object *always* has the name 'status'.
// How to use this feature?
// Whenever a status object is required, define it like this:
//   #ifndef NDEBUG
//     MPI_Status status;
//   #endif
// Within the corresponding MPI routines, then use "MPI_STATUS_MACRO" instead of "&status", e.g.
//   MPI_Recv(..., some_source, some_tag, some_comm, MPI_STATUS_MACRO);
#ifdef NDEBUG
#  define MPI_STATUS_MACRO MPI_STATUS_IGNORE
#else
#  define MPI_STATUS_MACRO &status
#endif
//\endcond

/// FEAST namespace
namespace FEAST
{
  /// FEAST version
  enum
  {
    /// FEAST major version number
    version_major = 1,
    /// FEAST minor version number
    version_minor = 0,
    /// FEAST patch version number
    version_patch = 0
  };

  /**
  * \brief index data type
  *
  * This type is used for indexing entities on a single matrix patch like e.g. vertice indices, degrees of freedom
  * or column indices in a sparse matrix implementation.
  */
  typedef unsigned long index_t;

  /**
  * \brief global index data type
  *
  * This type is used for indexing entities in a global parallel simulation.
  */
  typedef unsigned long index_glob_t;

  /**
  * \brief Nil class definition.
  *
  * This is an empty tag class which may be used for templates with optional parameters.\n
  * Some template implementations might recognise the usage of a \c Nil parameter as <em>parameter not given</em>.
  */
  class Nil
  {
  }; // class Nil
} // namespace FEAST

#endif // KERNEL_BASE_HEADER_HPP
