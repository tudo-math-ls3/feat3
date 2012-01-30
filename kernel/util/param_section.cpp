#include <kernel/util/param_section.hpp>

#include <fstream>
#include <stack>

namespace FEAST
{
  ParamSection::ParamSection()
  {
  }

  ParamSection::~ParamSection()
  {
    // delete all sub-sections
    SectionMap::iterator it = _sections.begin();
    SectionMap::iterator jt = _sections.end();
    for(; it != jt ; ++it)
    {
      delete (*it).second;
    }
  }

  bool ParamSection::add_entry(String key, String value, bool replace)
  {
    // try to insert the key-value-pair
    std::pair<ValueMap::iterator, bool> rtn = _values.insert(std::pair<String, String>(key,value));
    if(!rtn.second)
    {
      // insertion failed, i.e. there already exists a pair with that key - replace it?
      if(!replace)
      {
        // Do not replace, so return false.
        return false;
      }
      // okay, replace the existing value string
      (*rtn.first).second = value;
    }
    // insertion / replacement successful
    return true;
  }

  ParamSection* ParamSection::add_section(String name)
  {
    // try to find the entry
    SectionMap::iterator it = _sections.find(name);

    // if it has not been found
    if(it != _sections.end())
    {
      return (*it).second;
    }

    // if it was found, create a new section
    ParamSection* sub_section = new ParamSection();
    _sections.insert(std::pair<String,ParamSection*>(name, sub_section));
    return sub_section;
  }

  void ParamSection::erase_section(String name)
  {
    SectionMap::iterator it = _sections.find(name);
    if(it != _sections.end())
    {
      _sections.erase (it);
    }
  }


  void ParamSection::erase_entry(String key)
  {
    ValueMap::iterator it = _values.find(key);
    if(it != _values.end())
    {
      _values.erase (it);
    }
  }

  void ParamSection::parse(String filename)
  {
    // try to open the file
    std::ifstream ifs(filename.c_str(), std::ios::in);

    // if something went wrong
    if(!ifs.is_open())
    {
      throw FileNotFound(filename);
    }

    //parsing
    try
    {
      parse(ifs);
      ifs.close();
    }
    catch(ParamSection::SyntaxError& exc)
    {
      // If the exception does not contain a filename, we'll recycle the exception
      // and include our filename now.
      if(exc.get_filename().empty())
      {
        throw(SyntaxError(exc.message(), filename));
      }
      else
      {
        throw exc;
      }
    }
  }


  void ParamSection::parse(std::istream& ifs)
  {
    // Initialize stack etc.
    String s, t;
    String::size_type found;
    bool sec = false; // marks if a section was declared
    int line = 0; //counts 'lines'

    ParamSection* current = this; // current section
    std::stack<ParamSection*> stack; // stores the 'opened' sections

    stack.push(current);

    while(!ifs.eof() && ifs.good())  // while the end is not reached and nothing went wrong
    {
      // get a line
      getline(ifs, s);
      ++line;

      // erasing comments
      found = s.find("#");
      if(found!= std::string::npos)
      {
        s = s.substr(0, found);
      }

      // erasing whitespaces
      s.trim();

      // if its empty
      if(s.size() == 0)
      {
        continue;
      }

      // if there is a "&", get the next line, too
      while (s.size()>1 && s.at(s.size()-1) == '&')
      {
        s.erase(s.length()-1);
        getline(ifs, t);
        ++line;
        found = t.find("#");
        if(found != std::string::npos)
        {
          t = t.substr(0, found);
        }
        t.trim();
        s = s + t;
      }

      // if its a section
      if((s.at(0) == '[') && (s.at(s.size()-1) == ']'))
      {
        // removing brackets and empty spaces
        s.erase(0, 1);
        s.erase(s.length()-1);
        s.trim();

        // if there is nothing left
        if(s.size() == 0)
        {
          throw SyntaxError("Missing section name in line " +stringify(line));
        }
        // adding section
        sec = true;
        current = current->add_section(s);
      }

      // if its a value
      else if(s.find("=") != std::string::npos)
      {
        sec = false;
        // parting the line
        found = s.find("=");
        t = s.substr(found + 1);
        s.erase(found);
        s.trim();
        t.trim();

        // if there is nothing left
        if(s.size() == 0)
        {
          throw SyntaxError("Missing key in line " +stringify(line));
        }

        current->add_entry (s, t);
      }

      // if its an "include@"
      else if(s.size() > 8 && s.substr(0,8) == "include@")
      {
        sec = false;
        // erasing "include@"
        t = s.substr(8);
        current->parse(t);
      }

      else if(s  == "{")
      {
        // if the new section was declared
        if(sec)
        {
          stack.push(current);
        }
        // if it was not or there is something in the line above, that is not valid
        else
        {
          throw SyntaxError("Wrong '{' or '}' or section initialization in line " +stringify(line) +" or above!");
        }
      }

      else if(s  == "}")
      {
        sec = false;
        stack.pop();
        // if the stack is not empty, erase the top entry
        if(stack.size() > 0)
        {
          current=stack.top();
        }
        else
        {
          throw SyntaxError("Missing '{' or '}' in line " +stringify(line) +" or above!");
        }
      }

      // if its something else
      else
      {
        throw SyntaxError("Unknown format in line " +stringify(line) +" : " + s);
      }
    } // end while

    // in the end, the stack must contain nothing but the root section
    if(stack.size() != 1)
    {
      throw SyntaxError("Missing '{' or '}' in line " +stringify(line) +" or above!");
    }
  }

  void ParamSection::merge(const ParamSection& section, bool replace)
  {
    ValueMap::const_iterator valiter(section._values.begin());
    ValueMap::const_iterator valend(section._values.end());

    // merging _values of the two sections
    for(; valiter != valend; ++valiter)
    {
      add_entry(valiter->first, valiter->second, replace);
    }

    SectionMap::const_iterator seciter(section._sections.begin());
    SectionMap::const_iterator secend(section._sections.end());

    // merging _sections of the two ParamSections
    for(; seciter != secend; ++seciter)
    {
      // get ParamSection pointer to the section corresponding to iter and merge again
      add_section(seciter->first)->merge(*seciter->second, replace);
    }
  }

  ////////////////////////////////////////////////////////
  /////////// Output functions (into file etc.) //////////
  ////////////////////////////////////////////////////////

  std::pair<String, bool> ParamSection::get_entry(String key) const
  {
    ValueMap::const_iterator iter;
    // try to find the entry
    iter = _values.find(key);
    // if it was not found
    if(iter == _values.end())
    {
      return std::make_pair("", false);
    }
    return std::make_pair(iter->second, true);
  }


  const ParamSection* ParamSection::get_section(String name) const
  {
    SectionMap::const_iterator iter = _sections.find(name);
    if(iter == _sections.end())
    {
      return nullptr;
    }
    return iter->second;
  }


  ParamSection* ParamSection::get_section(String name)
  {
    SectionMap::iterator iter = _sections.find(name);
    if(iter == _sections.end())
    {
      return nullptr;
    }
    return iter->second;
  }

  void ParamSection::dump(String filename) const
  {
    // open stream
    std::ofstream ofs(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!ofs.is_open())
    {
      throw FileError("Failed to create '" + filename +"' for dumping!");
    }

    // dump into stream
    dump(ofs);

    // close stream
    ofs.close();
  }

  void ParamSection::dump(std::ostream& os, String::size_type indent) const
  {
    // prefix string: 2*indent spaces
    String prefix(2*indent, ' ');

    // dump values
    ValueMap::const_iterator vit = _values.begin();
    ValueMap::const_iterator vend = _values.end();
    for( ; vit != vend ; ++vit)
    {
      os << prefix << (*vit).first << " = " << (*vit).second << std::endl;
    }

    // dump subsections
    SectionMap::const_iterator sit = _sections.begin();
    SectionMap::const_iterator send = _sections.end();
    for( ; sit != send ; ++sit)
    {
      // section name and opening brace
      os << prefix << "[" << (*sit).first << "]" << std::endl;
      os << prefix << "{" << std::endl;

      // dump subsection with increased indent
      (*sit).second->dump(os, indent + 1);

      // closing brace, including section name as comment
      os << prefix << "} # end of [" << (*sit).first << "]" << std::endl;
    }
  }


} //namespace FEAST
