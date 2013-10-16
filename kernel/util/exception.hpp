#pragma once
#ifndef KERNEL_UTIL_EXCEPTION_HPP
#define KERNEL_UTIL_EXCEPTION_HPP 1

// The following line is necessary - otherwise doxygen won't document the #define's in this file.
/** \file */

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <cstdlib>
#include <iostream>
#include <list>

namespace FEAST
{
#ifndef FEAST_NO_CONTEXT

  struct ContextData;
#endif // FEAST_NO_CONTEXT

  /**
  * \brief Base exception class.
  *
  * \author Dirk Ribbrock
  */
  class Exception :
    public std::exception
  {
  private:
#ifndef FEAST_NO_CONTEXT
    /// Our (local) context data.
    ContextData * const _context_data;
#endif // FEAST_NO_CONTEXT

    /// descriptive error message
    const String _message;

    /// Our what string (for std::exception).
    mutable String _what_str;

  protected:
    /**
    * \brief CTOR
    *
    * \param message
    * the exception's message
    */
    Exception(const String & message);

    /// copy CTOR
    Exception(const Exception & other);

  public:
    /// DTOR
    virtual ~Exception() throw();

    /// returns error message
    const String message() const;

#ifndef FEAST_NO_CONTEXT
    /// returns backtrace
    String backtrace(const String & delimiter) const;
#endif // FEAST_NO_CONTEXT

    /// returns true if the backtrace is empty
    bool empty() const;

    /// return descriptive exception name
    const char * what() const throw();
  };


  /**
  * \brief exception that is thrown if something that is never supposed to happen happens
  *
  * It simply prefixes the exception message with "Internal error: ", otherwise it does not differ from its
  * parent class Exception.
  *
  * \author Dirk Ribbrock
  */
  class InternalError :
    public Exception
  {
  public:
    /**
    * \brief Constructor.
    *
    * \param message
    * A short error message.
    */
    InternalError(const String & message) :
      Exception("Internal error: " + message)
    {
    }
  };

  /**
   * \brief Syntax Error exception class
   *
   * This class derives from FEAST::Exception and is thrown by the parser classes when a syntax error
   * is detected.
   *
   * \author Constantin Christof
   * \author Peter Zajac
   */
  class SyntaxError :
    public Exception
  {
  protected:
    /// name of the file containing the syntax error
    String _filename;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] message
     * A description of the syntax error.
     *
     * \param[in] filename
     * The name of the file in which the syntax error has been detected.
     */
     explicit SyntaxError(
     String message,
      String filename = ""):
      Exception(message + (filename.empty() ? "": " in file " + filename)),
      _filename(filename)
    {
    }

    /// virtual destructor
    virtual ~SyntaxError() throw()
    {
    }

    /// returns the filename
    String get_filename() const
    {
      return _filename;
    }
  }; // class SyntaxError


#if !defined(FEAST_NO_CONTEXT) || defined(FEAST_TRACE_CONTEXT)
  /**
  * \brief Backtrace class context.
  *
  * \author Dirk Ribbrock
  * \author Peter Zajac
  */
  class Context
    : public InstantiationPolicy<Context, NonCopyable>
  {
#  ifdef FEAST_TRACE_CONTEXT
  private:
    /// context livetrace prefix string
    static String _prefix;

  public:
    /// specifies whether context live-tracing is enabled
    static bool live_trace;
#  endif // defined(FEAST_TRACE_CONTEXT)

  public:
    /// specifies whether to print file names and line numbers
    static bool file_line;

  public:
    /**
    * \brief Constructor.
    *
    * \param file
    * name of the source file that contains the context
    * \param line
    * line number of the context
    * \param context
    * description of the context
    */
    Context(const char * const file, const long line, const String & context);

    /// DTOR
    ~Context();

#  ifndef FEAST_NO_CONTEXT
    /**
    * \brief Current context
    *
    * \param[in] delimiter
    * A delimiter added between to context strings
    */
    static String backtrace(const String & delimiter);
#endif // !defined(FEAST_NO_CONTEXT)
  }; // class Context
#endif // !defined(FEAST_NO_CONTEXT) || defined(FEAST_TRACE_CONTEXT)

/**
 * \def CONTEXT
 * \brief Convenience definition that provides a way to declare uniquely-named instances of class Context.
 *
 * \param msg
 * Context message that is to be displayed by an exception-triggered backtrace or during livetracing.
 */
/**
 * \def CONTEXT_FILE_LINE
 * \brief Enables or disables decoration of context strings with file names and line numbers.
 *
 * \param enable
 * Specifies whether to enable (\c true) or disable (\c false) file name and line number decoration.
 */
#if (defined (DEBUG) && !defined(FEAST_NO_CONTEXT)) || defined(FEAST_TRACE_CONTEXT)
  // C preprocessor abomination following...
#  define CONTEXT_NAME_(x) ctx_##x
#  define CONTEXT_NAME(x) CONTEXT_NAME_(x)
#  define CONTEXT(msg) Context CONTEXT_NAME(__LINE__)(__FILE__, __LINE__, (msg))
#  define CONTEXT_FILE_LINE(enable) Context::file_line = (enable)
#else
#  define CONTEXT(msg)
#  define CONTEXT_FILE_LINE(enable)
#endif // (defined (DEBUG) && !defined(FEAST_NO_CONTEXT)) || defined(FEAST_TRACE_CONTEXT)

/**
 * \def CONTEXT_LIVE_TRACE
 * \brief Enables or disables context livetracing.
 *
 * \param enable
 * Specifies whether to enable (\c true) or disable (\c false) context livetracing.
 */
#ifdef FEAST_TRACE_CONTEXT
#  define CONTEXT_LIVE_TRACE(enable) Context::live_trace = (enable)
#else
#  define CONTEXT_LIVE_TRACE(enable)
#endif // defined(FEAST_TRACE_CONTEXT)

} // namespace FEAST

#endif // KERNEL_UTIL_EXCEPTION_HPP
