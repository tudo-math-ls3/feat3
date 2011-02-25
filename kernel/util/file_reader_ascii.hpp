#pragma once
#ifndef UTIL_FILE_READER_ASCII_HPP
#define UTIL_FILE_READER_ASCII_HPP 1

// includes, system
#include <iostream>
#include <fstream>

// includes, Feast
#include <kernel/util/exception.hpp>

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

    /// CTOR
    FileReaderASCII(const std::string& filename, const char comment_char, const bool skip_empty_lines)
      : _filename(filename), _comment_char(comment_char), _skip_empty_lines(skip_empty_lines)
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


    /// skips empty lines, and checks if first non-empty line starts with given keyword
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


    /// stores first non-empty non-comment line in given string
    inline void read(std::string& value)
    {
      CONTEXT("FileReaderASCII::read()");
      if(!get_next_line(value))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
    }


    /// extracts one integer from first non-empty non-comment line (assumes that the line contains nothing else except
    /// possible trailing comments)
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


    /// extracts two integers from first non-empty non-comment line (assumes that the line contains nothing else except
    /// possible trailing comments)
    inline void read(int& value1, int& value2)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2;
    }


    /// extracts three integers from first non-empty non-comment line (assumes that the line contains nothing else except
    /// possible trailing comments)
    inline void read(int& value1, int& value2, int& value3)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3;
    }


    /// extracts four integers from first non-empty non-comment line (assumes that the line contains nothing else except
    /// possible trailing comments)
    inline void read(int& value1, int& value2, int& value3, int& value4)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3 >> value4;
    }


    /// extracts one integer from first non-empty non-comment line (assumes that the line contains nothing else except
    /// possible trailing comments)
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


    /// extracts two integers from first non-empty non-comment line (assumes that the line contains nothing else except
    /// possible trailing comments)
    inline void read(unsigned int& value1, unsigned int& value2)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2;
    }


    /// extracts four integers from first non-empty non-comment line (assumes that the line contains nothing else except
    /// possible trailing comments)
    inline void read(unsigned int& value1, unsigned int& value2, unsigned int& value3, unsigned int& value4)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3 >> value4;
    }


    /// extracts two doubles and one int from first non-empty non-comment line
    /// (assumes that the line contains nothing else except possible trailing comments)
    inline void read(double& value1, double& value2, unsigned int& value3)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2 >> value3;
    }


    /// extracts two doubles from first non-empty non-comment line
    /// (assumes that the line contains nothing else except possible trailing comments)
    inline void read(double& value1, double& value2)
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> value1 >> value2;
    }


    /// extracts two doubles from first non-empty non-comment line, assumes the array values to have length 2
    /// (assumes that the line contains nothing else except possible trailing comments)
    inline void read(double values[])
    {
      CONTEXT("FileReaderASCII::read()");
      std::string line;
      if(!get_next_line(line))
      {
        throw InternalError("Unexpected end of file (in " + _filename + ").");
      }
      std::istringstream(line) >> values[0] >> values[1];
    }


    /// extracts two ints and three doubles from first non-empty non-comment line
    /// (assumes that the line contains nothing else except possible trailing comments)
    inline void read(int& value1, int& value2, double& value3, double& value4, double& value5)
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
