#pragma once
#ifndef KERNEL_UTIL_EXCEPTION_HPP
#define KERNEL_UTIL_EXCEPTION_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/string_utils.hpp>

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

#ifndef FEAST_NO_CONTEXT
  /**
  * \brief Backtrace class context.
  *
  * \author Dirk Ribbrock
  */
  class Context
    : public InstantiationPolicy<Context, NonCopyable>
  {
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

    /**
    * \brief Current context
    *
    * \param[in] delimiter
    * A delimiter added between to context strings
    */
    static String backtrace(const String & delimiter);
  };
#endif // FEAST_NO_CONTEXT

  /**
  * \def CONTEXT
  *
  * \brief Convenience definition that provides a way to declare uniquely-named instances of class Context.
  *
  * The created Context will be automatically provided with the correct filename and line number.
  *
  * \param s
  * Context message that can be display by an exception-triggered backtrace.
  *
  * \warning Will only be compiled in when debug support is enabled.
  */
#if defined (DEBUG) && !defined(FEAST_NO_CONTEXT)
  // C preprocessor abomination following...
#define CONTEXT_NAME_(x) ctx_##x
#define CONTEXT_NAME(x) CONTEXT_NAME_(x)
#define CONTEXT(s) \
  Context CONTEXT_NAME(__LINE__)(__FILE__, __LINE__, (s))
#else
#define CONTEXT(s)
#endif

} // namespace FEAST

#endif // KERNEL_UTIL_EXCEPTION_HPP
