#include <kernel/util/property_map.hpp>
#include <kernel/util/dist_file_io.hpp>

#include <fstream>
#include <stack>

namespace FEAT
{
  PropertyMap::PropertyMap(PropertyMap* parent) :
    _parent(parent)
  {
  }

  PropertyMap::~PropertyMap()
  {
    // delete all sub-sections
    SectionMap::iterator it(_sections.begin());
    SectionMap::iterator jt(_sections.end());
    for(; it != jt ; ++it)
    {
      if((*it).second != nullptr)
      {
        delete (*it).second;
      }
    }
  }

  bool PropertyMap::add_entry(String key, String value, bool replace)
  {
    // try to insert the key-value-pair
    std::pair<EntryMap::iterator, bool> rtn(_values.insert(std::make_pair(key,value)));
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

  PropertyMap* PropertyMap::add_section(String name)
  {
    // try to find the entry
    SectionMap::iterator it(_sections.find(name));

    // if it has been found
    if(it != _sections.end())
    {
      return (*it).second;
    }

    // if it was not found, create a new section
    PropertyMap* sub_section = new PropertyMap(this);
    _sections.insert(std::make_pair(name, sub_section));
    return sub_section;
  }

  bool PropertyMap::erase_entry(String key)
  {
    EntryMap::iterator it(_values.find(key));
    if(it != _values.end())
    {
      _values.erase(it);
      return true;
    }
    return false;
  }

  bool PropertyMap::erase_section(String name)
  {
    SectionMap::iterator it(_sections.find(name));
    if(it != _sections.end())
    {
      _sections.erase(it);
      return true;
    }
    return false;
  }

  PropertyMap* PropertyMap::query_section(String path)
  {
    // We're too lazy to copy-&-paste the code below, so we'll do some casting...
    return const_cast<PropertyMap*>(static_cast<const PropertyMap*>(this)->query_section(path));
  }

  const PropertyMap* PropertyMap::query_section(String path) const
  {
    // first of all, trim the input of leading and trailing slashes
    path.trim_me("/");

    // we need to split up the input path into:
    // the leading section name
    String sec_name(path);
    // the remainder of the path
    String rem_path("");

    // try to find the first slash inside our path
    std::size_t p = path.find('/');
    if(p != path.npos)
    {
      // split our path at the separator
      sec_name = path.substr(0, p).trim();
      rem_path = path.substr(p+1).trim();
    }

    // try to find the leading section
    const PropertyMap* lead = nullptr;
    if(sec_name == "!")
      lead = this->get_root();
    else if(sec_name == "~")
      lead = this->get_parent();
    else
      lead = this->get_sub_section(sec_name);

    // do we have a remaining path?
    if(rem_path.empty())
      return lead;
    else if(lead != nullptr)
      return lead->query_section(rem_path);
    else
      return nullptr;
  }

  std::pair<String, bool> PropertyMap::query(String key_path) const
  {
    // first of all, trim the input of leadin and trailing slashes
    key_path.trim_me("/");

    // try to find the last slash inside our path
    std::size_t p = key_path.rfind('/');

    // no path separator found?
    if(p == key_path.npos)
    {
      // interpret the path as a simple name then
      return get_entry(key_path);
    }

    // okay, split the input path at the slash; this gives us
    // the leading section path
    String sec_path = key_path.substr(0, p).trim();
    // the key name
    String key_name = key_path.substr(p+1).trim();

    // try to query the section
    const PropertyMap* sec = query_section(sec_path);

    // no luck finding the section?
    if(sec == nullptr)
      return std::make_pair(String(), false);
    else
      return sec->get_entry(key_name);
  }

  String PropertyMap::query(String key_path, String default_value) const
  {
    // use the other query function to search for the value
    std::pair<String, bool> sb = this->query(key_path);

    // if the key was found, return its value, otherwise return default value
    return sb.second ? sb.first : default_value;
  }

  std::pair<String, bool> PropertyMap::get_entry(String key) const
  {
    EntryMap::const_iterator iter(_values.find(key));
    if(iter == _values.end())
    {
      return std::make_pair("", false);
    }
    return std::make_pair(iter->second, true);
  }

  const PropertyMap* PropertyMap::get_sub_section(String name) const
  {
    SectionMap::const_iterator iter(_sections.find(name));
    if(iter == _sections.end())
    {
      return nullptr;
    }
    return iter->second;
  }

  PropertyMap* PropertyMap::get_sub_section(String name)
  {
    SectionMap::iterator iter(_sections.find(name));
    if(iter == _sections.end())
    {
      return nullptr;
    }
    return iter->second;
  }

  void PropertyMap::merge(const PropertyMap& section, bool replace)
  {
    EntryMap::const_iterator valiter(section._values.begin());
    EntryMap::const_iterator valend(section._values.end());

    // merging _values of the two sections
    for(; valiter != valend; ++valiter)
    {
      add_entry(valiter->first, valiter->second, replace);
    }

    SectionMap::const_iterator seciter(section._sections.begin());
    SectionMap::const_iterator secend(section._sections.end());

    // merging _sections of the two PropertyMaps
    for(; seciter != secend; ++seciter)
    {
      // get PropertyMap pointer to the section corresponding to iter and merge again
      add_section(seciter->first)->merge(*seciter->second, replace);
    }
  }

  void PropertyMap::read(const Dist::Comm& comm, String filename, bool replace)
  {
    std::stringstream is;
    DistFileIO::read_common(is, filename, comm);
    read(is, replace);
  }

  void PropertyMap::read(std::istream& ifs, bool replace)
  {
    // a stack to keep track of all currently open sections
    std::stack<PropertyMap*> stack;

    // a pointer to the currently used section; pushed onto the bottom of the stack
    PropertyMap* current = this;
    stack.push(current);

    // the string containing the current line
    String line;
    line.reserve(256);

    // other auxiliary variables
    String key, value;
    String::size_type found;
    int cur_line = 0;

    // an enumeration for keeping track what was read
    enum LastRead
    {
      None,       // nothing read until now
      Entry,      // last entity was a key-value pair
      Section,    // last entity was a section marker
      BraceOpen,  // last entity was an open curly brace
      BraceClose, // last entity was a closed curly brace
      Include     // last entity was an include
    };

    // keep track of what was the last thing we read
    LastRead last_read = None;

    // loop over all lines until we reach the end of the file
    while(!ifs.eof() && ifs.good())
    {
      // get a line
      getline(ifs, line);
      ++cur_line;

      // trim whitespaces; continue with next line if the current one is empty
      if(line.trim_me().empty())
      {
        continue;
      }

      // check for special keywords; these begin with an @ sign
      /*if(line.front() == '@')
      {
        // if its an "@include "
        if(line.compare(0, 9, "@include ") == 0)
        {
          // an include may appear anywhere in a file
          last_read = Include;

          // erase "@include "
          line.erase(0, 9);

          // parse the included file
          current->parse(line.trim(), replace);

          // okay
          continue;
        }
        else
        {
          throw SyntaxError("Unknown keyword in line " + stringify(cur_line) + " : " + line);
        }
      }
      else*/
      if((found = line.find('#')) != std::string::npos)
      {
        // erase comment
        line.erase(found);

        // trim whitespaces; continue with next line if the current one is empty
        if(line.trim_me().empty())
        {
          continue;
        }
      }

      // if its a section
      if((line.front() == '[') && (line.back() == ']'))
      {
        // a section may appear anywhere
        last_read = Section;

        // removing brackets and trimming
        line.pop_front();
        line.pop_back();

        // the section name must not be empty
        if(line.trim_me().empty())
        {
          throw SyntaxError("Missing section name in line " + stringify(cur_line));
        }

        // adding section
        current = stack.top()->add_section(line);
      }

      // if its a key-value pair
      else if((found = line.find('=')) != std::string::npos)
      {
        // an entry must not appear after a closing brace
        if(last_read == BraceClose)
        {
          throw SyntaxError("Unexpected key-value pair in line " + stringify(cur_line) + " after '}'");
        }
        last_read = Entry;

        // extract key
        key = line.substr(0, found);

        // the key must be not empty after trimming
        if(key.trim_me().empty())
        {
          throw SyntaxError("Missing key in line " + stringify(cur_line));
        }

        // extract value string; add it if it's empty
        value = String(line.substr(found + 1)).trim();
        if(value.empty())
        {
          current->add_entry(key, value, replace);
          continue;
        }

        // save current line number
        int old_line = cur_line;

        // check for line continuation
        while(value.back() == '&')
        {
          // erase the '&' character
          value.pop_back();

          // read next non empty line
          line.clear();
          while(line.empty() && !ifs.eof() && ifs.good())
          {
            // read next line
            getline(ifs, line);
            ++cur_line;

            // erase comments
            if((found = line.find('#')) != std::string::npos)
            {
              line.erase(found);
            }

            // erase whitespaces
            line.trim_me();
          }

          // do we have a non-empty string or an eof?
          if(line.empty())
          {
            // file ended before we could find a line continuation
            throw SyntaxError("Only empty lines from line " + stringify(old_line) + " on for line continuation");
          }

          // append line to value string
          value += line;

          // we'll continue the while-loop as the line might be further continued
        }

        // add the key-value pair
        current->add_entry(key, value, replace);
      }

      else if(line == "{")
      {
        // an open curly brace is allowed only after a section marker
        if(last_read == Section)
        {
          // push section onto the stack
          stack.push(current);
          last_read = BraceOpen;
        }
        // if it was not or there is something in the line above, that is not valid
        else
        {
          throw SyntaxError("Unexpected '{' in line " + stringify(cur_line));
        }
      }

      else if(line == "}")
      {
        // a closed curly brace is allowed anywhere except at the beginning of the file; in addition
        // to that the stack must contain at least 2 entries as the root section always is the bottom
        // entry of the stack and cannot be removed.
        if((last_read != None) && (stack.size() > 1))
        {
          // remove top section from stack
          stack.pop();
          // and get the section underneath it
          current = stack.top();
          last_read = BraceClose;
        }
        else
        {
          throw SyntaxError("Unexpected '}' in line " + stringify(cur_line));
        }
      }

      // if its something else
      else
      {
        throw SyntaxError("Unknown format in line " +stringify(cur_line) + " : " + line);
      }
    } // end while

    // in the end, the stack must contain nothing but the root section
    if(stack.size() > 1)
    {
      throw SyntaxError("Missing '}' in line " + stringify(cur_line) + " or above!");
    }
  }

  void PropertyMap::write(const Dist::Comm& comm, String filename) const
  {
    // only rank 0 writes
    if(comm.rank() > 0)
      return;

    // open stream
    std::ofstream ofs(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!ofs.is_open())
    {
      throw FileError("Failed to create '" + filename +"' for dumping!");
    }

    // dump into stream
    write(ofs);

    // close stream
    ofs.close();
  }

  void PropertyMap::write(std::ostream& os, String::size_type indent) const
  {
    // prefix string: 2*indent spaces
    String prefix(2*indent, ' ');

    // dump values
    EntryMap::const_iterator vit(_values.begin());
    EntryMap::const_iterator vend(_values.end());
    for( ; vit != vend ; ++vit)
    {
      os << prefix << (*vit).first << " = " << (*vit).second << std::endl;
    }

    // dump subsections
    SectionMap::const_iterator sit(_sections.begin());
    SectionMap::const_iterator send(_sections.end());
    for( ; sit != send ; ++sit)
    {
      // section name and opening brace
      os << prefix << "[" << (*sit).first << "]" << std::endl;
      os << prefix << "{" << std::endl;

      // dump subsection with increased indent
      (*sit).second->write(os, indent + 1);

      // closing brace, including section name as comment
      os << prefix << "} # end of [" << (*sit).first << "]" << std::endl;
    }
  }
} //namespace FEAT
