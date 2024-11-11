// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// The following line is necessary - otherwise doxygen won't document the #define's in this file.
/** \file */

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <iostream>

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

  protected:
    /**
     * \brief CTOR
     *
     * \param message_in
     * the exception's message.
     */
    explicit Exception(const String & message_in) :
      _message(message_in)
    {
    }

    /**
     * \brief CTOR
     *
     * \param function the current function name.
     * \param file the current file name.
     * \param line the current line number.
     * \param message_in the exception's message.
     */
    Exception(
        const char* const function,
        const char* const file,
        const long line,
        const String & message_in) :
      _message(stringify(file) + ":" + stringify(line) + ": in " + stringify(function) + ": " + message_in)
    {
    }

    /// copy CTOR
    Exception(const Exception &)
    {
    }

  public:
    /// DTOR
    virtual ~Exception() noexcept
    {
    }

    /// returns error message
    const String message() const
    {
      return _message;
    }

    /// return descriptive exception name
    virtual const char * what() const noexcept override
    {
      return _message.empty() ? std::exception::what() : _message.c_str();
    }
  }; // class Exception


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
    explicit InternalError(const String & message_in) :
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
  }; // class InternalError

  /**
   * \brief Class for parser related errors
   */
  class ParseError :
    public Exception
  {
  public:
    /**
     * \brief Constructor.
     *
     * \param[in] message_in
     * A short error message.
     */
    explicit ParseError(const String& message_in) :
      Exception("ParseError: " + message_in)
    {
    }

    /**
     * \brief Constructor
     *
     * \param[in] name
     * The name of the entry that was to be parsed.
     *
     * \param[in] got
     * The string that was meant to be parsed.
     *
     * \param[in] expect
     * What was expected instead, e.g. "a non-negative integer"
     */
    explicit ParseError(const String& name, const String& got, const String& expect) :
      Exception("ParseError: " + name + ": expected " + expect + " but got '" + got + "'")
    {
    }

    virtual ~ParseError() noexcept
    {
    }
  }; // class ParseError

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
    virtual ~FileError() noexcept
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
    virtual ~FileNotFound() noexcept
    {
    }
  }; // class FileNotFound

  /**
   * \brief File-Not-Created exception
   *
   * This exception is thrown when a file could not be created.
   *
   * \author Peter Zajac
   */
  class FileNotCreated :
    public FileError
  {
  public:
    /**
     * \brief Constructor
     * \param[in] filename
     * The name (and path) of the file that could not be created.
     */
    explicit FileNotCreated(const String& filename) :
      FileError("Could not create file: '" + filename + "'")
    {
    }
    virtual ~FileNotCreated() noexcept
    {
    }
  }; // class FileNotCreated

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
     * \param[in] message_in
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
    virtual ~SyntaxError() noexcept
    {
    }

    /// returns the filename
    String get_filename() const
    {
      return _filename;
    }
  }; // class SyntaxError
} // namespace FEAT
