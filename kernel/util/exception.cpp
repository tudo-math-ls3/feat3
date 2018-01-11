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
    return _message.empty() ? std::exception::what() : _message.c_str();
  }
} // namespace FEAT
