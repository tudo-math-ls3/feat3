#pragma once
#ifndef KERNEL_UTIL_ASSERTION_HPP
#define KERNEL_UTIL_ASSERTION_HPP 1

// The following line is necessary - otherwise doxygen won't document the #define's in this file.
/** \file */

// includes, FEAST
#include <kernel/util/exception.hpp>

// includes, system

// check for standard C assert usage
#ifdef FEAST_STDC_ASSERT
// ensure that NDEBUG is defined unless DEBUG is defined
#  if !defined(DEBUG) && !defined(NDEBUG)
#    define NDEBUG
#  endif
#  include <cassert>
#endif // FEAST_STDC_ASSERT

namespace FEAST
{
  /**
  * \brief defining assertion
  *
  * An assertion is thrown when a critical condition is not fulfilled. Together with the macro defined below, it
  * replaces the standard C assert(...).
  *
  * \author Dirk Ribbrock
  */
  class Assertion
    : public Exception
  {

  public:

    /**
    * \brief CTOR
    *
    * \param[in] function
    * name of the function in which the assertion failed
    *
    * \param[in] file
    * name of the source file that contains the failed assertion
    *
    * \param[in] line
    * line number of the failed assertion
    *
    * \param[in] message_in
    * message that shall be displayed
    */
    Assertion(
      const char * const function,
      const char * const file,
      const long line,
      const String & message_in)
      : Exception(stringify(file) + ":" + stringify(line) + ": in " + stringify(function) + ": " + message_in)
    {
#ifndef FEAST_NO_CONTEXT
      std::cout << backtrace("\n") << this->message() << std::endl;
#else
      std::cout << this->message() << std::endl;
#endif
    }
  };

/**
 * \def ASSERT
 * \brief Convenience definition that provides a way to throw Assertion exceptions.
 *
 * The thrown Assertion will be automatically provided with the correct filename,
 * line number and function name.
 *
 * \param expr Boolean expression that shall be asserted.
 * \param msg Error message that will be display in case that expr evaluates to false.
 *
 * \note This macro will only be compiled in debug mode; it is an empty macro in no-debug mode.
 */
/**
 * \def ASSERT_
 * \brief Convenience definition that provides a way to throw Assertion exceptions.
 *
 * The thrown Assertion will be automatically provided with the correct filename, line number and function name.\n
 * In contrast to the #ASSERT macro, this macro has only one parameter, whereas the error message is a stringified
 * version of the expression to be asserted.
 *
 * \param expr Boolean expression that shall be asserted.
 *
 * \note This macro will only be compiled in debug mode; it is an empty macro in no-debug mode.
 */
#if defined (FEAST_STDC_ASSERT)
// define as standard C assert
/// \compilerhack microsoft compiler complains about its own flaws
#  if defined(FEAST_COMPILER_MICROSOFT)
//   Note: The MSC assert macro has a "bug": the __LINE__ macro, which is internally used by the
//         "assert" macro, is of type 'long', whereas the function receiving the value of
//         __LINE__ interprets it as 'unsigned int', which triggers a warning C4365.
//         The following is a hack to make the MSC compiler shut up about its own flaws.
#    define ASSERT(expr, msg) __pragma(warning(push)) __pragma(warning(disable:4365)) assert(expr) __pragma(warning(pop))
#    define ASSERT_(expr) __pragma(warning(push)) __pragma(warning(disable:4365)) assert(expr) __pragma(warning(pop))
#  else
#    define ASSERT(expr, msg) assert(expr)
#    define ASSERT_(expr) assert(expr)
#  endif
#elif defined (DEBUG)
// use FEAST::Assertion exception
#  define ASSERT(expr, msg) \
    do { \
        if (! (expr)) \
            throw FEAST::Assertion(__func__, __FILE__, __LINE__, msg); \
    } while (false)
#  define ASSERT_(expr) ASSERT(expr, #expr)
#else
#  define ASSERT(expr, msg)
#  define ASSERT_(expr)
#endif

} // namespace FEAST

#endif // KERNEL_UTIL_ASSERTION_HPP
