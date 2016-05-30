// includes, FEAT
#include <kernel/util/exception.hpp>

namespace FEAT
{
  Exception::Exception(const String & message_in) :
    _message(message_in)
  {
  }

  Exception::Exception(
      const char* const function,
      const char* const file,
      const long line,
      const String & message_in) :
    _message(stringify(file) + ":" + stringify(line) + ": in " + stringify(function) + ": " + message_in)
  {
  }

  Exception::Exception(const Exception & other) :
    std::exception(other),
    _message(other._message)
  {
  }

  Exception::~Exception() throw()
  {
  }

  const String Exception::message() const
  {
    return _message;
  }

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
} // namespace FEAT
