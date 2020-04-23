// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_ASSERTION_HPP
#define KERNEL_UTIL_ASSERTION_HPP 1

// The following line is necessary - otherwise doxygen won't document the #define's in this file.
/** \file */

// includes, FEAT
#include <kernel/base_header.hpp>

// includes, system
#include <string>

namespace FEAT
{
  /**
   * \brief Abortion function
   *
   * This function implements the actual abortion that is called by the XABORTM macro, which prints
   * an informative error message to stderr and calls Runtime::abort to terminate the process.
   *
   * \param[in] func
   * The name of the function that contains the abortion, usually <c>__func__</c>.
   *
   * \param[in] file
   * The name of the source/header file that contains the abortion, usually <c>__FILE__</c>.
   *
   * \param[in] line
   * The line number of the abortion in the source/header file, usually <c>__LINE__</c>.
   *
   * \param[in] msg
   * A custom error message to be displayed in addition to the standard information.
   */
  [[noreturn]] void abortion(const char* const func, const char* const file, const int line, const char* const msg);

  [[noreturn]] inline void abortion(const char* const func, const char* const file, const int line, const std::string& msg)
  {
    abortion(func, file, line, msg.c_str());
  }

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
   * The stringified expression, usually <c>\#expr</c>.
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

  inline void assertion(
    bool expr,
    const char* const expr_str,
    const char* const func,
    const char* const file,
    const int line,
    const std::string& msg)
  {
    assertion(expr, expr_str, func, file, line, msg.c_str());
  }

  /**
   * \def XABORTM
   * \brief Abortion macro definition with custom message
   *
   * This macro prints an errors message and aborts program execution.
   *
   * \param msg
   * An error message that is to be displayed.
   *
   * This macro will be compiled in both debug and non-debug mode builds.
   *
   * \compilerhack Intel C++ compiler is too dumb for [[noreturn]] attribute,
   * so attach an additional call to std::abort() to silence warning 1011
   * that complains about missing return statements in non-void functions.
    */
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 2000)
#  define XABORTM(msg) do {FEAT::abortion(__func__, __FILE__, __LINE__, msg); std::abort();} while(false)
#else // any smart C++ compiler
#  define XABORTM(msg) FEAT::abortion(__func__, __FILE__, __LINE__, msg)
#endif


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
#  define ASSERT(expr) FEAT::assertion(expr, #expr, __func__, __FILE__, __LINE__)
#  define ASSERTM(expr, msg) FEAT::assertion(expr, #expr, __func__, __FILE__, __LINE__, msg)
#else
#  define ASSERT(expr) void(0)
#  define ASSERTM(expr, msg) void(0)
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
#define XASSERT(expr) FEAT::assertion(expr, #expr, __func__, __FILE__, __LINE__)
#define XASSERTM(expr, msg) FEAT::assertion(expr, #expr, __func__, __FILE__, __LINE__, msg)

} // namespace FEAT

#endif // KERNEL_UTIL_ASSERTION_HPP
