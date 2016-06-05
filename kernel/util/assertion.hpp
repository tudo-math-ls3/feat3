#pragma once
#ifndef KERNEL_UTIL_ASSERTION_HPP
#define KERNEL_UTIL_ASSERTION_HPP 1

// The following line is necessary - otherwise doxygen won't document the #define's in this file.
/** \file */

// includes, FEAT
#include <kernel/base_header.hpp>

namespace FEAT
{
  /**
   * \brief Assertion function
   *
   * This function implements the actual assertion that is called by the ASSERT, ASSERTM,
   * XASSERT and XASSERTM macros.
   *
   * The behaviour of this function is as follows:
   * - If \p expr evaluates to \c true, then this function does nothing.
   * - If \p expr evaluates to \c false, then this function prints an informative error message
   *   to stderr and calls Runtime::abort to terminate the process.
   *
   * \param[in] expr
   * The expression that is asserted.
   *
   * \param[in] expr_str
   * The stringified expression, usually <c>#expr</c>.
   *
   * \param[in] func
   * The name of the function that contains the assertion, usually <c>__func__</c>.
   *
   * \param[in] file
   * The name of the source/header file that contains the assertion, usually <c>__FILE__</c>.
   *
   * \param[in] line
   * The line number of the assertion in the source/header file, usually <c>__LINE__</c>.
   *
   * \param[in] msg
   * A custom error message to be displayed in addition to the standard information.
   * May be \c nullptr if no additional information is available for the assertion.
   */
  void assertion(
    bool expr,
    const char * const expr_str,
    const char * const func,
    const char * const file,
    const int line,
    const char * const msg = nullptr);

  /**
   * \def ASSERT
   * \brief Debug-Assertion macro definition
   *
   * This macro defines a debug-mode assertion that will abort program execution if
   * the asserted expression evaluates to \c false.
   *
   * \param expr
   * Boolean expression that shall be asserted.
   *
   * \note
   * This macro will only be compiled in debug mode; it is an empty macro in non-debug mode.
   * Use the XASSERT macro if you want to use the assertion in both debug and non-debug modes.
   */
  /**
   * \def ASSERTM
   * \brief Debug-Assertion macro definition with custom message
   *
   * This macro defines a debug-mode assertion that will abort program execution if
   * the asserted expression evaluates to \c false.
   *
   * \param expr
   * Boolean expression that shall be asserted.
   * \param msg
   * An error message that is to be displayed if the assertion fails.
   *
   * \note This macro will only be compiled in debug mode; it is an empty macro in non-debug mode.
   * Use the XASSERTM macro if you want to use the assertion in both debug and non-debug modes.
   */
#if defined(DEBUG)
#  define ASSERT(expr) assertion(expr, #expr, __func__, __FILE__, __LINE__)
#  define ASSERTM(expr, msg)  assertion(expr, #expr, __func__, __FILE__, __LINE__, msg)
#else
#  define ASSERT(expr)
#  define ASSERTM(expr, msg)
#endif

  /**
   * \def XASSERT
   * \brief Assertion macro definition
   *
   * This macro defines an assertion that will abort program execution if
   * the asserted expression evaluates to \c false.
   *
   * \param expr
   * Boolean expression that shall be asserted.
   *
   * \note
   * This macro will be compiled in both debug and non-debug mode builds.
   * Use the ASSERT macro if you want to use the assertion ony in debug builds.
   */
  /**
   * \def XASSERTM
   * \brief Assertion macro definition with custom message
   *
   * This macro defines a n assertion that will abort program execution if
   * the asserted expression evaluates to \c false.
   *
   * \param expr
   * Boolean expression that shall be asserted.
   * \param msg
   * An error message that is to be displayed if the assertion fails.
   *
   * This macro will be compiled in both debug and non-debug mode builds.
   * Use the ASSERTM macro if you want to use the assertion ony in debug builds.
   */
#define XASSERT(expr) assertion(expr, #expr, __func__, __FILE__, __LINE__)
#define XASSERTM(expr, msg)  assertion(expr, #expr, __func__, __FILE__, __LINE__, msg)

} // namespace FEAT

#endif // KERNEL_UTIL_ASSERTION_HPP
