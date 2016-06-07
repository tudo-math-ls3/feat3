#pragma once
#ifndef KERNEL_UTIL_EXCEPTION_HPP
#define KERNEL_UTIL_EXCEPTION_HPP 1

// The following line is necessary - otherwise doxygen won't document the #define's in this file.
/** \file */

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <cstdlib>
#include <iostream>
#include <list>

namespace FEAT
{
  /**
  * \brief Base exception class.
  *
  * \author Dirk Ribbrock
  */
  class Exception :
    public std::exception
  {
  private:
    /// descriptive error message
    const String _message;

    /// Our what string (for std::exception).
    mutable String _what_str;

  protected:
    /**
     * \brief CTOR
     *
     * \param message
     * the exception's message.
     */
    explicit Exception(const String & message);

    /**
     * \brief CTOR
     *
     * \param function the current function name.
     * \param file the current file name.
     * \param line the current line number.
     * \param message the exception's message.
     */
    Exception(
        const char* const function,
        const char* const file,
        const long line,
        const String & message);

    /// copy CTOR
    Exception(const Exception & other);

  public:
    /// DTOR
    virtual ~Exception() throw();

    /// returns error message
    const String message() const;

    /// return descriptive exception name
    virtual const char * what() const throw() override;
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
    * \param message_in
    * A short error message.
    */
    InternalError(const String & message_in) :
      Exception("Internal error: " + message_in)
    {
    }

    /**
    * \brief Constructor.
    *
    * \param function the current function name.
    * \param file the current file name.
    * \param line the current line number.
    * \param message_in A short error message.
    */
    InternalError(
        const char* const function,
        const char* const file,
        const long line,
        const String & message_in) :
      Exception(function, file, line, "Internal error: " + message_in)
    {
    }
  };

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

  /**
   * \brief Syntax Error exception class
   *
   * This class derives from FEAT::Exception and is thrown by the parser classes when a syntax error
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
     String message_in,
      String filename = ""):
      Exception(message_in + (filename.empty() ? "": " in file " + filename)),
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
} // namespace FEAT

#endif // KERNEL_UTIL_EXCEPTION_HPP
