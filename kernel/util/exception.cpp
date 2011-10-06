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


#ifndef FEAST_NO_CONTEXT
  Context::Context(const char * const file, const long line, const String & context)
  {
    if (! context_stack)
    {
      context_stack = new std::list<String>;
    }

    context_stack->push_back(context + " (" + stringify(file) + ":" + stringify(line) +")");
  }

  Context::~Context()
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

  String Context::backtrace(const String & delimiter)
  {
    if (! context_stack)
      return "";

    String str;
    return str.join(context_stack->begin(), context_stack->end(), delimiter);
  }
#endif // FEAST_NO_CONTEXT
} // namespace FEAST
