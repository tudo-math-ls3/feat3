#pragma once
#ifndef UTIL_ASSERTION_HH
/// Header guard
#define UTIL_ASSERTION_HH 1

#include <kernel/util/exception.hpp>
#include <string>
#include <iostream>

namespace Feast
{
    /**
     * Assertion is thrown when a critical condition is not fulfilled.
     */
    class Assertion :
        public Exception
    {
        public:
            /**
             * Constructor.
             *
             * \param function Name of the function in which the assertion failed.
             * \param file Name of the source file that contains the failed assertion.
             * \param line Line number of the failed assertion.
             * \param message Message that shall be displayed.
             */
            Assertion(const char * const function, const char * const file,
                const long line, const std::string & message) :
              Exception(stringify(file) + ":" + stringify(line) + ": in " + stringify(function) + ": " + message)
            {
              std::cout << backtrace("\n") << this->message() << std::endl;
            }
    };


/**
 * \def ASSERT
 *
 * \brief Convenience definition that provides a way to throw Assertion exceptions.
 *
 * The thrown Assertion will be automatically provided with the correct filename,
 * line number and function name.
 *
 * \param expr Boolean expression that shall be asserted.
 * \param msg Error message that will be display in case that expr evaluates to false.
 *
 * \warning Will only be compiled in when debug support is enabled.
 */
#if defined (DEBUG)
#define ASSERT(expr, msg) \
    do { \
        if (! (expr)) \
            throw Assertion(__PRETTY_FUNCTION__, __FILE__, __LINE__, msg); \
    } while (false)
#else
#define ASSERT(expr, msg)
#endif

/**
 * \def STATIC_ASSERT
 *
 * \brief Convenience definition that provides a way to throw compile-time Assertion exceptions.
 *
 * The thrown Assertion will be automatically provided with the correct filename,
 * line number and function name.
 *
 * \param expr Boolean expression that shall be asserted.
 * \param msg Error message that will be display in case that expr evaluates to false.
 *
 * \warning Will only be compiled in when debug support is enabled.
 */
#if defined (DEBUG)
#define STATIC_ASSERT(const_expr, msg) \
    do { \
        if (! (const_expr)) \
            throw Assertion(__PRETTY_FUNCTION__, __FILE__, __LINE__, msg); \
    } while (false)
#else
#define STATIC_ASSERT(const_expr, msg)
#endif
}

#endif //UTIL_ASSERTION_HPP
