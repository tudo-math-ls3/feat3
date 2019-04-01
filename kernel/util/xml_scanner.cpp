// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/xml_scanner.hpp>
#include <cctype>

namespace FEAT
{
  namespace Xml
  {
    DummyParser::DummyParser()
    {
    }

    bool DummyParser::attribs(std::map<String,bool>&) const
    {
      return false;
    }

    void DummyParser::create(int, const String&, const String&, const std::map<String, String>&, bool)
    {
    }

    void DummyParser::close(int, const String&)
    {
    }

    std::shared_ptr<MarkupParser> DummyParser::markup(int, const String&, const String&)
    {
      // return a new DummyParser node
      return std::make_shared<DummyParser>();
    }

    bool DummyParser::content(int, const String&)
    {
      return true;
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    DumpParser::DumpParser(const String& name, bool lines, std::size_t indent) :
      _name(name), _lines(lines), _indent(indent)
    {
    }

    bool DumpParser::attribs(std::map<String,bool>&) const
    {
      return false;
    }

    void DumpParser::create(int iline, const String&, const String& name, const std::map<String, String>& attrs, bool closed)
    {
      if(_lines) std::cout << std::setw(5) << iline << ":";
      std::cout << std::string(_indent, ' ');
      std::cout << "<" << name;
      for(auto it = attrs.begin(); it != attrs.end(); ++it)
      {
        std::cout << " " << it->first << "=\"" << it->second << "\"";
      }
      std::cout << (closed ? " />" : ">") << std::endl;
    }

    void DumpParser::close(int iline, const String&)
    {
      if(_lines) std::cout << std::setw(5) << iline << ":";
      std::cout << std::string(_indent, ' ');
      std::cout << "</" << _name << ">" << std::endl;
    }

    std::shared_ptr<MarkupParser> DumpParser::markup(int, const String&, const String& name)
    {
      return std::make_shared<DumpParser>(name, _lines, _indent+2);
    }

    bool DumpParser::content(int iline, const String& sline)
    {
      if(_lines) std::cout << std::setw(5) << iline << ":";
      std::cout << std::string(_indent+2, ' ');
      std::cout << sline << std::endl;
      return true;
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    void Scanner::read_root()
    {
      // ensure that the markup stack is empty
      if(!_markups.empty())
        throw InternalError(__func__, __FILE__, __LINE__, "markup stack is not empty");

      // make sure that we did not yet read a line
      if(_have_read_root)
        throw InternalError(__func__, __FILE__, __LINE__, "root markup already read");

      // try to read first line
      if(!read_next_line())
        throw_syntax("Root markup is missing");

      // ensure that it is a markup
      if(!scan_markup())
        throw_syntax("First line is not a markup");

      // ensure that it is not a bogus markup
      if(_cur_is_closed)
        throw_syntax("Invalid closed root markup");
      if(_cur_is_termin)
        throw_syntax("Invalid root terminator markup");

      // we have read a root markup
      _have_read_root = true;
    }

    void Scanner::set_root_parser(std::shared_ptr<MarkupParser> parser)
    {
      // avoid bogus
      if(parser == nullptr)
        throw InternalError(__func__, __FILE__, __LINE__, "root parser must not be nullptr");

      // make sure that we do not have any markups yet
      if(!_markups.empty())
        throw InternalError(__func__, __FILE__, __LINE__, "root parser already exists");

      // make sure that we did read a root parser
      if(!_have_read_root)
        throw InternalError(__func__, __FILE__, __LINE__, "root markup not read yet");

      // push our parser to the markup stack
      _markups.push_back(MarkupInfo(_cur_iline, _cur_name, parser));

      // try to create the root node
      //_markups.back().parser()->create(_cur_iline, _cur_sline, _cur_attribs, false);
      create_top_parser();
    }

    void Scanner::scan()
    {
      // ensure that we have a root parser
      if(_markups.empty())
        throw InternalError(__func__, __FILE__, __LINE__, "root parser does not yet exist");

      // scan through all lines
      while(read_next_line())
      {
        // is it a comment?
        if(_cur_sline.starts_with("<!--"))
        {
          if (!_cur_sline.ends_with("-->"))
            throw_syntax("Missing '-->' at end of comment");
          else
            continue;
        }

        // is it an XML markup line?
        if(scan_markup())
          process_markup();
        else
          process_content();

        // all markups processed?
        if(_markups.empty())
          return;
      }

      // ensure that all _markups have been closed
      if(!_markups.empty())
      {
        String msg = String("Expected '</") + _markups.back().name() + ">' but found end-of-file";
        throw SyntaxError(_cur_iline, "", msg);
      }
    }

    void Scanner::scan(std::shared_ptr<MarkupParser> root_parser)
    {
      read_root();
      set_root_parser(root_parser);
      scan();
    }

    bool Scanner::read_next_line()
    {
      while(_stream.good())
      {
        // read another line
        std::getline(_stream, _cur_sline);
        _cur_sline.trim_me();
        ++_cur_iline;

        // non-empty line?
        if(!_cur_sline.empty())
          return true;
      }

      // no more lines
      return false;
    }

    bool Scanner::scan_markup()
    {
      // basic assumptions
      _cur_is_markup = _cur_is_closed = _cur_is_termin = false;
      _cur_name.clear();
      _cur_attribs.clear();

      // XML markers?
      bool xhead = _cur_sline.starts_with('<');
      bool xtail = _cur_sline.ends_with('>');

      // Not a markup?
      if((!xhead) && (!xtail))
        return false;

      // Invalid combination?
      if(xhead && !xtail)
        throw_syntax("Expected '>' at end of line");
      if(!xhead && xtail)
        throw_syntax("Expected '<' at beginning of line");

      // both xhead and xtail are true: we have a markup here
      _cur_is_markup = true;

      // extract line data; that's everything between '<' and '>'
      String sdata = String(_cur_sline.substr(1, _cur_sline.length()-2)).trim();

      // make sure that it is not empty
      if(sdata.empty())
        throw_syntax("Empty markup");

      // check whether we have additional '<' or '>' inside the line; we do not support that
      if(sdata.find('<') != sdata.npos)
        throw_syntax("Invalid '<' in line");
      if (sdata.find('>') != sdata.npos)
        throw_syntax("Invalid '>' in line");

      // closed or a terminator?
      _cur_is_termin = (sdata.front() == '/');
      _cur_is_closed = (sdata.back() == '/');

      // make sure we have no bogus here
      if(_cur_is_termin && _cur_is_closed)
        throw_syntax("Invalid markup");

      // remove '/'
      if(_cur_is_termin)
        sdata.pop_front();
      if(_cur_is_closed)
        sdata.pop_back();

      // empty?
      if((sdata = sdata.trim()).empty())
        throw_syntax("Empty markup");

      // extract markup name
      {
        size_t n0 = sdata.find_first_of(sdata.whitespaces());
        if(n0 != sdata.npos)
        {
          // markup name followed by whitespace
          _cur_name = sdata.substr(0, n0);
          sdata = String(sdata.substr(n0)).trim();
        }
        else
        {
          _cur_name = sdata;
          sdata.clear();
        }

        // ensure that the name is not empty
        if(_cur_name.empty())
          throw_syntax("Empty markup name");

        // ensure that it contains only alphanumeric chars
        if(std::isalpha(_cur_name.front()) == 0)
          throw_syntax("Invalid first character in markup name");
        for(std::size_t i(1); i < _cur_name.size(); ++i)
        {
          if(std::isalnum(_cur_name.at(i)) == 0)
            throw_syntax(String("Invalid character '") + _cur_name.at(i) + "' in markup name");
        }
      }

      // is this a terminator?
      if(_cur_is_termin)
      {
        // make sure that we have no trailing data
        if(!sdata.empty())
          throw_syntax("Invalid terminator");

        // okay, it's a terminator
        return true;
      }

      // extract the attributes
      while(!sdata.empty())
      {
        // try to find a '='
        size_t p = sdata.find_first_of('=');
        if(p == sdata.npos)
          throw_syntax("Invalid attribute list");

        // extract attribute name
        String akey = String(sdata.substr(0, p)).trim();
        sdata = String(sdata.substr(p+1)).trim();

        // no key name?
        if(akey.empty())
          throw_syntax("Attribute key name missing");

        // ensure that it contains only alphanumeric chars
        if(std::isalpha(akey.front()) == 0)
          throw_syntax("Invalid first character in attribute name '" + akey + "'");
        for(std::size_t i(1); i < akey.size(); ++i)
        {
          if(std::isalnum(akey.at(i)) == 0)
            throw_syntax(String("Invalid character '") + akey.at(i) + "' in attribute name '" + akey + "'");
        }

        // parse attribute
        if(sdata.empty() || sdata.front() != '"')
          throw_syntax("Expected '\"'");

        // find next '"'
        size_t q = sdata.find_first_of('"', 1);
        if(q == sdata.npos)
          throw_syntax("Missing '\"'");

        // extract attribute value
        String aval = String(sdata.substr(1, q-1)).trim();

        // push attribute
        _cur_attribs.emplace(akey, aval);

        // remove attribute
        sdata = String(sdata.substr(q+1)).trim();
      }

      // okay
      return true;
    }

    void Scanner::process_markup()
    {
      // not a markup?
      if(!_cur_is_markup)
        throw InternalError(__func__, __FILE__, __LINE__, "invalid markup process call");

      // is it a terminator?
      if(_cur_is_termin)
      {
        // ensure that we have something to terminate
        if(_markups.empty())
          throw_syntax("Unexpected terminator");

        // make sure the terminator matches the top marker
        if(_cur_name != _markups.back().name())
          throw_syntax(String("Expected '</") + _markups.back().name() + String(">'"));

        // destroy parser
        _markups.back().parser()->close(_cur_iline, _cur_sline);

        // pop the last marker
        _markups.pop_back();

        // we're done
        return;
      }

      // try to create a parser
      std::shared_ptr<MarkupParser> parser = _markups.back().parser()->markup(_cur_iline, _cur_sline, _cur_name);

      // make sure this markup is supported
      if(parser == nullptr)
        throw_grammar("Unsupported markup");

      // push new markup
      _markups.push_back(MarkupInfo(_cur_iline, _cur_name, parser));

      // create top markup parser
      create_top_parser();

      // pop markup if closed
      if(_cur_is_closed)
      {
        _markups.back().parser()->close(_cur_iline, _cur_sline);
        _markups.pop_back();
      }
    }

    void Scanner::create_top_parser()
    {
      if(_markups.empty())
        throw InternalError(__func__, __FILE__, __LINE__, "Empty markup stack");
      if(!_cur_is_markup)
        throw InternalError(__func__, __FILE__, __LINE__, "no markup to create");

      // get the parser attributes
      std::map<String, bool> parser_attribs;
      if(_markups.back().parser()->attribs(parser_attribs))
      {
        // first of all, check whether all arguments are valid
        for (auto it = _cur_attribs.begin(); it != _cur_attribs.end(); ++it)
        {
          if(parser_attribs.find(it->first) == parser_attribs.end())
          {
            throw_grammar(String("Unexpected attribute '") + it->first + "'");
          }
        }

        // now let's check if all mandatory attributes are given
        for(auto jt = parser_attribs.begin(); jt != parser_attribs.end(); ++jt)
        {
          if(jt->second)
          {
            if(_cur_attribs.find(jt->first) == _cur_attribs.end())
            {
              throw_grammar(String("Mandatory attribute '") + jt->first + "' missing");
            }
          }
        }
      }

      // create the markup
      _markups.back().parser()->create(_cur_iline, _cur_sline, _cur_name, _cur_attribs, _cur_is_closed);
    }

    void Scanner::process_content()
    {
      // process content
      if(!_markups.back().parser()->content(_cur_iline, _cur_sline))
        throw_grammar("Invalid content line");
    }
  } // namespace Xml
} // namespace FEAT
