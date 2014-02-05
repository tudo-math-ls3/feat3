// includes, FEAST
#include <kernel/util/exception.hpp>

namespace FEAST
{
#ifndef FEAST_NO_CONTEXT
  /// The global context stack.
  std::list<String> * context_stack = 0;

  /**
   * \brief This structs holds the actual context history
   *
   * \author Dirk Ribbrock
   */
  struct ContextData
  {
    /// The local context stack
    std::list<String> local_context_stack;

    /// CTOR
    ContextData()
    {
      if (context_stack)
      {
        local_context_stack.assign(context_stack->begin(), context_stack->end());
      }
    }

    /// returns a full context stack aka backtrace
    String backtrace(const String& delimiter) const
    {
      if (context_stack)
      {
        String str;
        return str.join(local_context_stack.begin(), local_context_stack.end(), delimiter);
      }
      else return "";
    }
  };
#endif // FEAST_NO_CONTEXT

  Exception::Exception(const String & message) :
#ifndef FEAST_NO_CONTEXT
    _context_data(new ContextData),
#endif // FEAST_NO_CONTEXT
    _message(message)
  {
  }

  Exception::Exception(
      const char* const function,
      const char* const file,
      const long line,
      const String & message) :
#ifndef FEAST_NO_CONTEXT
    _context_data(new ContextData),
#endif // FEAST_NO_CONTEXT
    _message(stringify(file) + ":" + stringify(line) + ": in " + stringify(function) + ": " + message )
  {
  }

  Exception::Exception(const Exception & other) :
    std::exception(other),
#ifndef FEAST_NO_CONTEXT
    _context_data(new ContextData(*other._context_data)),
#endif // FEAST_NO_CONTEXT
    _message(other._message)
  {
  }

  Exception::~Exception() throw()
  {
#ifndef FEAST_NO_CONTEXT
    delete _context_data;
#endif // FEAST_NO_CONTEXT
  }

  const String Exception::message() const
  {
    return _message;
  }

#ifndef FEAST_NO_CONTEXT
  String Exception::backtrace(const String & delimiter) const
  {
    return _context_data->backtrace(delimiter);
  }
#endif // FEAST_NO_CONTEXT

  const char * Exception::what() const throw()
  {
    /// \todo Add a working win32 alternative (see http://www.int0x80.gr/papers/name_mangling.pdf)
    /*if (_what_str.empty())
      {
      int status(0);
      char * const name(abi::__cxa_demangle(("_Z" + stringify(std::exception::what())).c_str(), 0, 0, &status));
      if (0 == status)
      {
      _what_str = name;
      _what_str += " (" + message() + ")";
      std::free(name);
      }
      }*/
    if (_what_str.empty())
    {
      _what_str = stringify(std::exception::what());
      _what_str += " (" + message() + ")";
    }
    return _what_str.c_str();
  }

#if !defined(FEAST_NO_CONTEXT) || defined(FEAST_TRACE_CONTEXT)

#  ifdef FEAST_TRACE_CONTEXT
  // the prefix is initially empty
  String Context::_prefix("");

  // by default, disable context livetracting
  bool Context::live_trace(false);
#  endif // defined(FEAST_TRACE_CONTEXT)

  // by default, disable file names and line numbers
  bool Context::file_line(false);

  Context::Context(const char * const file, const long line, const String & context)
  {
    // format the trace message
    String msg;
    if(file_line)
    {
      msg = String(file) + "[" + stringify(line) + "]: " + context;
    }
    else
    {
      msg = context;
    }

#  ifndef FEAST_NO_CONTEXT
    {
      // create a context stack if we don't already have one
      if(context_stack == nullptr)
      {
        context_stack = new std::list<String>;
      }

      // push the context onto the stack
      context_stack->push_back(msg);
    }
#  endif // !defined(FEAST_NO_CONTEXT)

#  ifdef FEAST_TRACE_CONTEXT
    {
      // increment prefix
      _prefix.push_back('*');

      // trace context
      if(live_trace)
      {
        std::clog << _prefix << " " << msg << std::endl;
      }
    }
#  endif // defined(FEAST_TRACE_CONTEXT)
  } // Context::Context()

  Context::~Context()
  {
#  ifndef FEAST_NO_CONTEXT
    {
      // ensure that we have a context stack
      if (context_stack == nullptr)
        throw InternalError("no context!");

      // remove the last context from the stack
      context_stack->pop_back();

      // and delete the stack if its empty
      if (context_stack->empty())
      {
        delete context_stack;
        context_stack = nullptr;
      }
    }
#  endif // !defined(FEAST_NO_CONTEXT)

#  ifdef FEAST_TRACE_CONTEXT
    {
      // decrement prefix
      // Note: In C++03, the std::string has a push_back function but misses a pop_back function,
      // so we'll use erase() here.
      if(!_prefix.empty())
      {
        _prefix.erase(_prefix.size() - 1);
      }
    }
#  endif // defined(FEAST_TRACE_CONTEXT)
  } // Context::~Context()


#  ifndef FEAST_NO_CONTEXT
  String Context::backtrace(const String & delimiter)
  {
    if (! context_stack)
      return "";

    String str;
    return str.join(context_stack->begin(), context_stack->end(), delimiter);
  }
#  endif // !defined(FEAST_NO_CONTEXT)
#endif // !defined(FEAST_NO_CONTEXT) || defined(FEAST_TRACE_CONTEXT)

} // namespace FEAST
