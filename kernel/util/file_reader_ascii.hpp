/* GENERAL_REMARK_BY_HILMAR:
 * Rudimentary class for reading ascii files. Should be improved/replaced when the 'real' ascii/xml file reader is
 * implemented.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef UTIL_FILE_READER_ASCII_HPP
#define UTIL_FILE_READER_ASCII_HPP 1

// includes, system
#include <iostream>
#include <fstream>

// includes, FEAST
#include <kernel/util/exception.hpp>

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief wrapper class for ASCII (line-based) input
  *
  * TODO: Das hier ist eher sehr unhuebsch, und noch nicht wirklich typsicher. Kommt spaeter. Dito templatisierung.
  *
  * \author Dominik Goeddeke
  */
  class FileReaderASCII
  {

  private:
    /// line end delimiter
    /// TODO: std::endl for cross-platform compatibility? // DG Oct 19, 2010
    static const char ENDL = '\n';

    /// filename associated with this reader
    const std::string& _filename;

    /// stream associated with this reader
    std::ifstream myfile;

    /// leading comment char indicating lines (and remainder of lines) to be ignored
    char _comment_char;

    /// skip empty lines?
    bool _skip_empty_lines;

    /// main util function (returns false if end of file is reached, true otherwise)
    inline bool get_next_line(std::string& sbuf)
    {
      CONTEXT("FileReaderASCII::get_next_line()");
      // skip leading spaces
      while (myfile.peek() == ' ')
      {
        myfile.get();
      }
      // skip comment lines
      if(myfile.peek() == _comment_char)
      {
        getline(myfile, sbuf);
        return get_next_line(sbuf);
      }
      // skip empty lines
      else if(_skip_empty_lines && myfile.peek() == ENDL)
      {
        getline(myfile, sbuf);
        return get_next_line(sbuf);
      }
      // stop at EOF
      else if(!myfile.eof())
      {
        getline(myfile, sbuf);
        // trim line of trailing comments
        size_t pos = sbuf.find(_comment_char);
        sbuf = sbuf.substr(0,pos-1);
        return true;
      }
      else
      {
        return false;
      }
    }


  public:

    /**
    * \brief CTOR
    *
    * \param[in] filename
    * name of the file to be connected to this file reader
    *
    * \param[in] comment_char
    * character to be interpreted as comment character
    *
    * \param[in] skip_empty_lines
    * flag whether to skip empty lines
    */
    FileReaderASCII(
      const std::string& filename,
      const char comment_char,
      const bool skip_empty_lines)
      : _filename(filename),
        _comment_char(comment_char),
        _skip_empty_lines(skip_empty_lines)
    {
      CONTEXT("FileReaderASCII::FileReaderASCII()");
      myfile.open(filename.c_str(), std::ios::in);
      if(myfile.fail())
      {
        throw InternalError("Problem occured when opening file " + _filename + ".");
      }
    }


    /// DTOR
    ~FileReaderASCII()
    {
      CONTEXT("FileReaderASCII::~FileReaderASCII()");
      myfile.close();
    }


    /// skips one line
    inline void skip_line()
    {
      CONTEXT("FileReaderASCII::skip_line()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
    }


    /**
    * \brief skips empty lines, and checks if first non-empty line starts with given keyword
    *
    * \param[in] keyword
    * keyword to be checked
    */
    inline void read(std::string const& keyword)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      if(line.find(keyword) == std::string::npos)
      {
        throw InternalError("Keyword <" + keyword + "> not found");
      }
    }


    /**
    * \brief stores first non-empty non-comment line in given string
    *
    * \param[out] value
    * stores the line
    */
    inline void read(std::string& value)
    {
      CONTEXT("FileReaderASCII::read()");
      if(!get_next_line(value))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
    }


    /**
    * \brief extracts two strings from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first string value
    *
    * \param[out] value2
    * stores second string value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      std::string& value1,
      std::string& value2)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }

      // since "std::istringstream(line) >> value1 >> value2" does not work for strings, we have to extract the two
      // key-value strings manually

      // remove leading and trailing white spaces from line
      line = trim(line);

      if(line.length() == 0)
      {
        throw InternalError("Line is emty!");
      }

      // find position of the first blank after the key
      size_t pos = line.find_first_of(' ');
      if(pos == std::string::npos)
      {
        throw InternalError("No white space found after key in line '" + line + "'.");
      }
      // set first value to first token in the line
      value1 = line.substr(0, pos);

      // set second value to the trimmed remaining line
      value2 = trim(line.substr(pos+1));
    }


    /**
    * \brief extracts one integer from first non-empty non-comment line
    *
    * \param[out] value
    * stores first integer value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(int& value)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value;
    }


    /**
    * \brief extracts two integers from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first integer value
    *
    * \param[out] value2
    * stores second integer value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      int& value1,
      int& value2)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2;
    }


    /**
    * \brief extracts three integers from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first integer value
    *
    * \param[out] value2
    * stores second integer value
    *
    * \param[out] value3
    * stores first integer value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      int& value1,
      int& value2,
      int& value3)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3;
    }


    /**
    * \brief extracts four integers from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first integer value
    *
    * \param[out] value2
    * stores second integer value
    *
    * \param[out] value3
    * stores first integer value
    *
    * \param[out] value4
    * stores second integer value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      int& value1,
      int& value2,
      int& value3,
      int& value4)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3 >> value4;
    }


    /**
    * \brief extracts one unsigned integer from first non-empty non-comment line
    *
    * \param[out] value
    * stores first integer value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(unsigned int& value)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value;
    }


    /**
    * \brief extracts two unsigned integers from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first unsigned integer value
    *
    * \param[out] value2
    * stores second unsigned integer value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      unsigned int& value1,
      unsigned int& value2)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2;
    }


    /**
    * \brief extracts four unsigned integers from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first unsigned integer value
    *
    * \param[out] value2
    * stores second unsigned integer value
    *
    * \param[out] value3
    * stores first unsigned integer value
    *
    * \param[out] value4
    * stores second unsigned integer value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      unsigned int& value1,
      unsigned int& value2,
      unsigned int& value3,
      unsigned int& value4)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3 >> value4;
    }


    /**
    * \brief extracts two doubles and one unsigned integer from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first double value
    *
    * \param[out] value2
    * stores second double value
    *
    * \param[out] value3
    * stores first unsigned integer value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      double& value1,
      double& value2,
      unsigned int& value3)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3;
    }


    /**
    * \brief extracts two doubles from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first double value
    *
    * \param[out] value2
    * stores second double value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      double& value1,
      double& value2)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2;
    }


    /**
    * \brief extracts two doubles from first non-empty non-comment line and stores them in the given array
    *
    * \param[out] values
    * array of double values
    *
    * \note Assumes that the line contains nothing else except possible trailing comments. Assumes the array to have
    * at least length 2.
    */
    inline void read2(double values[])
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> values[0] >> values[1];
    }

    /**
    * \brief extracts two integers and three doubles from first non-empty non-comment line
    *
    * \param[out] value1
    * stores first integer value
    *
    * \param[out] value2
    * stores second integer value
    *
    * \param[out] value3
    * stores first double value
    *
    * \param[out] value4
    * stores second double value
    *
    * \param[out] value5
    * stores third double value
    *
    * \note Assumes that the line contains nothing else except possible trailing comments.
    */
    inline void read(
      int& value1,
      int& value2,
      double& value3,
      double& value4,
      double& value5)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3 >> value4 >> value5;
    }
  }; // class FileReaderASCII
} // namespace FEAST

#endif // guard UTIL_FILE_READER_ASCII_HPP
