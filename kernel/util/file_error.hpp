#pragma once
#ifndef KERNEL_UTIL_FILE_ERROR_HPP
#define KERNEL_UTIL_FILE_ERROR_HPP 1

// includes, FEAST
#include <kernel/util/exception.hpp>

namespace FEAST
{
  /**
   * \brief Base class for file related errors
   * \author Peter Zajac
   */
  class FileError :
    public Exception
  {
  public:
    /**
     * \brief Constructor
     * \param[in] message_in
     * The error message.
     */
    explicit FileError(const String& message_in) :
      Exception(message_in)
    {
    }
    virtual ~FileError() throw()
    {
    }
  }; // class FileError

  /**
   * \brief File-Not-Found exception
   *
   * This exception is thrown when a file was not found.
   *
   * \author Peter Zajac
   */
  class FileNotFound :
    public FileError
  {
  public:
    /**
     * \brief Constructor
     * \param[in] filename
     * The name (and path) of the file that was not found.
     */
    explicit FileNotFound(const String& filename) :
      FileError("File not found: '" + filename + "'")
    {
    }
    virtual ~FileNotFound() throw()
    {
    }
  }; // class FileNotFound
} // namespace FEAST

#endif // KERNEL_UTIL_FILE_ERROR_HPP
