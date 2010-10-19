#pragma once
#ifndef UTIL_EXCEPTION_HHP
/// Header guard
#define UTIL_EXCEPTION_HHP 1

// includes, system
#include <string>
#include <cstdlib>
#include <iostream>
#include <list>

// includes, Feast
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/base_header.hpp>

/// FEAST namespace
namespace FEAST
{
  /// The global context stack.
  /// \todo Ist der stack global oder compile-unit lokal?
  std::list<std::string> * context_stack = 0;

  /**
  * \brief This structs holds the actual context history
  *
  * \author Dirk Ribbrock
  */
  struct ContextData
  {
    /// The local context stack
    std::list<std::string> local_context_stack;

    /// CTOR
    ContextData()
    {
      if (context_stack)
      {
        local_context_stack.assign(context_stack->begin(), context_stack->end());
      }
    }

    /// returns a full context stack aka backtrace
    std::string backtrace(const std::string& delimiter) const
    {
      return StringUtils::join(local_context_stack.begin(), local_context_stack.end(), delimiter);
    }
  };


  /**
  * \brief Base exception class.
  *
  * \author Dirk Ribbrock
  */
  class Exception :
    public std::exception
  {
  private:

    /// Our (local) context data.
    ContextData * const _context_data;

    /// descriptive error message
    const std::string _message;

    /// Our what string (for std::exception).
    mutable std::string _what_str;

  protected:
    /**
    * \brief CTOR
    *
    * \param message
    * the exception's message
    */
    Exception(const std::string & message) throw () :
      _context_data(new ContextData),
      _message(message)
    {
    }

    /// copy CTOR
    Exception(const Exception & other) :
      std::exception(other),
      _context_data(new ContextData(*other._context_data)),
      _message(other._message)
    {
    }

  public:
    /// DTOR
    virtual ~Exception() throw ()
    {
      delete _context_data;
    }

    /// returns error message
    const std::string message() const throw ()
    {
      return _message;
    }

    /// returns backtrace
    std::string backtrace(const std::string & delimiter) const
    {
      return _context_data->backtrace(delimiter);
    }

    /// returns true if the backtrace is empty
    bool empty() const;

    /// return descriptive exception name
    const char * what() const throw ()
    {
      /// \todo Add a working win32 alternative (see http://www.int0x80.gr/papers/name_mangling.pdf)
      /*if (_what_str.empty())
      {
        int status(0);
        char * const name(abi::__cxa_demangle(
              ("_Z" + StringUtils::stringify(std::exception::what())).c_str(), 0, 0, &status));
        if (0 == status)
        {
          _what_str = name;
          _what_str += " (" + message() + ")";
          std::free(name);
        }
      }*/
      if (_what_str.empty())
      {
        _what_str = StringUtils::stringify(std::exception::what());
        _what_str += " (" + message() + ")";
      }
      return _what_str.c_str();
    }
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
    InternalError(const std::string & message) throw () :
      Exception("Internal error: " + message)
    {
    }
  };


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
    Context(const char * const file, const long line, const std::string & context)
    {
      if (! context_stack)
      {
        context_stack = new std::list<std::string>;
      }

      context_stack->push_back(context + " (" + StringUtils::stringify(file) + ":"
                               + StringUtils::stringify(line) +")");
    }

    /// DTOR
    ~Context()
    {
      if (! context_stack)
        throw InternalError("no context!");

      context_stack->pop_back();

      if (context_stack->empty())
      {
        delete context_stack;
        context_stack = 0;
      }
    }

    /**
    * \brief Current context
    *
    * \param[in] delimiter
    * A delimiter added between to context strings
    */
    static std::string backtrace(const std::string & delimiter)
    {
      if (! context_stack)
        return "";

      return StringUtils::join(context_stack->begin(), context_stack->end(), delimiter);
    }
  };


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
#if defined (DEBUG)
  // C preprocessor abomination following...
#define CONTEXT_NAME_(x) ctx_##x
#define CONTEXT_NAME(x) CONTEXT_NAME_(x)
#define CONTEXT(s) \
  Context CONTEXT_NAME(__LINE__)(__FILE__, __LINE__, (s))
#else
#define CONTEXT(s)
#endif

} // namespace FEAST

#endif //UTIL_EXCEPTION_HPP
