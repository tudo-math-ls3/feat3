// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_XML_SCANNER_HPP
#define KERNEL_UTIL_XML_SCANNER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <memory>
#include <iostream>
#include <fstream>
#include <deque>
#include <map>

namespace FEAT
{
  /**
   * \brief Xml namespace
   */
  namespace Xml
  {
    /**
     * \brief Xml Error base class
     *
     * This class acts as a common base class for all errors related to the
     * parsing of XML files.
     *
     * \author Peter Zajac
     */
    class Error :
      public FEAT::Exception
    {
    private:
      /// erroneous line number
      int my_iline;
      /// erroneous line string
      String my_sline;
      /// error message
      String my_msg;
      /// what string
      String my_what;

    public:
      Error(int iline, const String& sline, const String& msg, const String& err_type) :
        FEAT::Exception("XML Error"),
        my_iline(iline),
        my_sline(sline),
        my_msg(msg),
        my_what(err_type + ": " + msg + " in line " + stringify(iline) + ": '" + sline + "'")
      {
      }

      /// \returns the erroneous line number
      int get_line_number() const
      {
        return my_iline;
      }

      /// \returns the erroneous line string
      const String& get_line_string() const
      {
        return my_sline;
      }

      /// \returns the error message
      const String& get_message() const
      {
        return my_msg;
      }

      /// \returns the what message
      virtual const char * what() const throw() override
      {
        return my_what.c_str();
      }
    }; // class Error

    /**
     * \brief Xml syntax error class
     *
     * This error class is thrown if the XML scanner encounters a syntax error.
     *
     * A "syntax error" is any type of error that violates the syntax of a valid XML file
     * under the limitations of the XML scanner class.
     *
     * \author Peter Zajac
     */
    class SyntaxError
      : public Error
    {
    public:
      SyntaxError(int iline, const String& sline, const String& msg) :
        Error(iline, sline, msg, "XML Syntax Error")
      {
      }
    }; // class SyntaxError

    /**
     * \brief Xml grammar error class
     *
     * This error class is thrown if the XML scanner encounters a grammar error.
     *
     * A "grammar error" is any type of error that does not violate the XML syntax,
     * but that cannot be interpreted by the parser.
     *
     * \author Peter Zajac
     */
    class GrammarError
      : public Error
    {
    public:
      GrammarError(int iline, const String& sline, const String& msg) :
        Error(iline, sline, msg, "XML Grammar Error")
      {
      }
    }; // class GrammarError

    /**
     * \brief Xml content error class
     *
     * This error class is thrown if the XML parser encounters a content error.
     *
     * A "content error" is thrown when a parser node cannot interpret a content line.
     *
     * \author Peter Zajac
     */
    class ContentError
      : public Error
    {
    public:
      ContentError(int iline, const String& sline, const String& msg) :
        Error(iline, sline, msg, "XML Content Error")
      {
      }
    }; // class ContentError

    /**
     * \brief XML Markup Parser interface
     *
     * This class acts as an interface for all XML markup parsers.
     *
     * \author Peter Zajac
     */
    class MarkupParser
    {
    public:
      /// virtual destructor
      virtual ~MarkupParser() {}

      /**
       * \brief Specifies the mandatory and optional attributes
       *
       * \param[out] attrs
       * A map of all supported attribute key names. The second component specifies
       * whether the attribute is mandatory (\c true) or optional (\c false).
       *
       * \returns
       * \c true, if the scanner should check for valid arguments, otherwise \c false.
       */
      virtual bool attribs(std::map<String,bool>& attrs) const = 0;

      /**
       * \brief Creates this markup parser node
       *
       * \param[in] iline
       * The line number of the markup in the XML file.
       *
       * \param[in] sline
       * The line string of the markup in the XML file.
       *
       * \param[in] name
       * The name of the markup.
       *
       * \param[in] attrs
       * A map of all attributes of the markup.
       *
       * \param[in] closed
       * Specifies whether the markup is closed.
       */
      virtual void create(
        int iline,
        const String& sline,
        const String& name,
        const std::map<String, String>& attrs,
        bool closed) = 0;

      /**
       * \brief Closes this markup parser node
       *
       * \param[in] iline
       * The line number of the markup terminator in the XML file.
       *
       * \param[in] sline
       * The line string of the markup in the XML file.
       */
      virtual void close(int iline, const String& sline) = 0;

      /**
       * \brief Called to process a child markup node
       *
       * This function is called by the XML scanner when a child markup
       * inside this markup parser node is detected.
       *
       * \param[in] iline
       * The line number of the child markup in the XML file.
       *
       * \param[in] sline
       * The line string of the child markup in the XML file.
       *
       * \param[in] name
       * The markup name of the child markup.
       *
       * \returns
       * A new MarkupParser object if the markup is a valid child or \c nullptr
       * if a child markup with this name is not supported.
       *
       * \note
       * Returning \c nullptr will cause the scanner to throw a Xml::GammarError.
       */
      virtual std::shared_ptr<MarkupParser> markup(int iline, const String& sline, const String& name) = 0;

      /**
       * \brief Called to process a content line
       *
       * This function is called by the XML scanner when a content line inside
       * this markup parser node is detected.
       *
       * \param[in] iline
       * The line number of the content line in the XML file.
       *
       * \param[in] sline
       * The line string of the content line in the XML file.
       *
       * \returns
       * \c true, if the content line was parsed successfully, otherwise \c false.
       *
       * \note
       * Returning \c false will cause the scanner to throw a Xml::GammarError.
       */
      virtual bool content(int iline, const String& sline) = 0;
    }; // class MarkupParser

    /**
     * \brief Dummy XML markup parser class
     *
     * This class implements the MarkupParser interface an implements
     * a "dummy" parser which accepts all kinds of markups and content,
     * but does not do anything with it - think of it as <c>/dev/null</c>
     *
     * This class may be somewhat useful for development purposes.
     *
     * \author Peter Zajac
     */
    class DummyParser :
      public MarkupParser
    {
    public:
      explicit DummyParser();
      virtual bool attribs(std::map<String,bool>&) const override;
      virtual void create(int, const String&, const String&, const std::map<String, String>&, bool) override;
      virtual void close(int, const String&) override;
      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override;
      virtual bool content(int, const String&) override;
    }; // class DummyParser

    /**
     * \brief Dump XML markup parser class
     *
     * This class implements the MarkupParser interface and simply 'dumps'
     * the parsed XML data to std::cout.
     * This class may be somewhat useful for debugging purposes.
     *
     * \author Peter Zajac
     */
    class DumpParser :
      public MarkupParser
    {
    protected:
      /// the name of this markup parser node
      String _name;
      /// specifies whether to write line numbers
      bool _lines;
      /// indentation for this parser node
      std::size_t _indent;

    public:
      explicit DumpParser(const String& name, bool lines = true, std::size_t indent = 0);
      virtual bool attribs(std::map<String,bool>&) const override;
      virtual void create(int iline, const String&, const String&, const std::map<String, String>& attrs, bool closed) override;
      virtual void close(int iline, const String&) override;
      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String& name) override;
      virtual bool content(int iline, const String& sline) override;
    }; // class DumpParser

    /**
     * \brief XML Scanner class
     *
     * \todo document
     *
     * \author Peter Zajac
     */
    class Scanner
    {
    private:
      /**
       * \brief Markup info class
       *
       * This class is used internally for the management of the markup parser stack.
       */
      class MarkupInfo
      {
      private:
        int my_line;
        String my_name;
        std::shared_ptr<MarkupParser> my_parser;

      public:
        explicit MarkupInfo(int _line, const String& _name, std::shared_ptr<MarkupParser> _parser) :
          my_line(_line), my_name(_name), my_parser(_parser)
        {
        }

        int line() const
        {
          return my_line;
        }

        const String& name() const
        {
          return my_name;
        }

        std::shared_ptr<MarkupParser> parser()
        {
          return my_parser;
        }
      }; // class MarkupInfo

    private:
      /// the input stream
      std::istream& _stream;
      /// the markup stack
      std::vector<MarkupInfo> _markups;
      /// specifies whether the root has been read
      bool _have_read_root;
      /// current line number
      int _cur_iline;
      /// current line string
      String _cur_sline;
      /// specifies whether this line is a markup
      bool _cur_is_markup;
      /// specifies whether the markup is closed
      bool _cur_is_closed;
      /// specifies whether the markup is a terminator
      bool _cur_is_termin;
      /// specifies the markup name of the current line
      String _cur_name;
      /// specifies the markup attributes of the current line
      std::map<String,String> _cur_attribs;

    public:
      /**
       * \brief Creates a XML scanner for an input stream
       *
       * \param[in] stream
       * The input stream that is to be read from.
       *
       * \param[in] lines_read
       * The number of lines already read from the stream.
       */
      explicit Scanner(std::istream& stream, int lines_read = 0) :
        _stream(stream),
        _markups(),
        _have_read_root(false),
        _cur_iline(lines_read),
        _cur_sline(),
        _cur_is_markup(false),
        _cur_is_closed(false),
        _cur_is_termin(false),
        _cur_name(),
        _cur_attribs()
      {
        _markups.reserve(32);
        _cur_sline.reserve(1024);
      }

      /// virtual destructor
      virtual ~Scanner()
      {
      }

      /// \returns the current line number
      int get_cur_iline() const
      {
        return _cur_iline;
      }

      /// \returns the current line string
      const String& get_cur_sline() const
      {
        return _cur_sline;
      }

      /// \returns \c true, if the current line is a markup line, otherwise \c false
      bool is_cur_markup() const
      {
        return _cur_is_markup;
      }

      /// \returns \c true, if the current markup line is a closed markup
      bool is_cur_closed() const
      {
        return _cur_is_closed;
      }

      /// \returns \c true, if the current markup line is a terminator
      bool is_cur_termin() const
      {
        return _cur_is_termin;
      }

      /// \returns the name of the current markup line
      const String& get_cur_name() const
      {
        return _cur_name;
      }

      /// \returns the attribute map of the current markup line
      const std::map<String,String>& get_cur_attribs() const
      {
        return _cur_attribs;
      }

      /**
       * \brief Tries to read the XML root markup.
       *
       * This function reads the first non-empty line from the input stream
       * and tries to interpret it as a XML markup representing the root markup.
       *
       * If this function returns without throwing an Xml::SyntaxError, the
       * information about the XML root markup can be read from this object by
       * using the getter functions.
       */
      void read_root();

      /**
       * \brief Sets the root markup parser node
       *
       * This function sets the parser for the root markup.
       *
       * \attention
       * This function can only be called after the #read_root function has
       * been called.
       *
       * \param[in] parser
       * The MarkupParser object responsible for parsing the root node.
       */
      void set_root_parser(std::shared_ptr<MarkupParser> parser);

      /**
       * \brief Scans the stream.
       *
       * This function scans the streams and processes all lines until the
       * root markup terminator is found.
       *
       * \attention
       * This function can only be called after a root parser has been set
       * by calling the #set_root_parser function, otherwise an InternalError
       * will be thrown.
       */
      void scan();

      /**
       * \brief Scans the stream with a given root parser.
       *
       * This function effectively just calls the three functions
       * - #read_root()
       * - #set_root_parser()
       * - #scan()
       */
      void scan(std::shared_ptr<MarkupParser> root_parser);

      /**
       * \brief Throws a SyntaxError for the current line
       *
       * \param[in] msg
       * The error message.
       */
      void throw_syntax(const String& msg) const
      {
        throw SyntaxError(_cur_iline, _cur_sline, msg);
      }

      /**
       * \brief Throws a GrammarError for the current line
       *
       * \param[in] msg
       * The error message.
       */
      void throw_grammar(const String& msg) const
      {
        throw GrammarError(_cur_iline, _cur_sline, msg);
      }

      /**
       * \brief Throws a ContentError for the current line
       *
       * \param[in] msg
       * The error message.
       */
      void throw_content(const String& msg) const
      {
        throw ContentError(_cur_iline, _cur_sline, msg);
      }

    protected:
      /**
       * \brief Reads the next non-empty line from the stream.
       *
       * \returns
       * \c true, if another non-empty line was read or
       * \c false, if the end of the stream was reached
       */
      bool read_next_line();

      /**
       * \brief Tries to interpret the current line as a XML markup line.
       *
       * This function checks whether the current line seems to be a XML markup line
       * and, if so, tries to decompose it into its parts.
       *
       * \returns
       * \c true, if the current line is a valid XML markup line or
       * \c false, if the current line seems to be a content line
       *
       * \note
       * This function throws a Xml::SyntaxError if the line seems to be
       * a markup line but violates the syntax rules.
       */
      bool scan_markup();

      /**
       * \brief Tries to process the current XML markup line.
       *
       * This function passes the current markup line to the top parser
       * and checks whether a new top parser node has to be allocated.
       *
       * This function throws a Xml::GrammarError if the current markup
       * is not accepted by the parser.
       */
      void process_markup();

      /**
       * \brief Creates the top parser node.
       *
       * This function checks whether the attributes match the top parser's
       * requirements and, if so, calls the create() method of the parser.
       */
      void create_top_parser();

      /**
       * \brief Tries to process the current XML content line.
       *
       * This function passes the current content line to the top parser.
       */
      void process_content();
    }; // class Scanner
  } // namespace Xml
} // namespace FEAT

#endif // KERNEL_UTIL_XML_SCANNER_HPP
