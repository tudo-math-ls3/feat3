// includes, FEAST
#include <kernel/util/mesh_streamer.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <fstream>

namespace FEAST
{
  /// \cond internal
  // read next line from stream
  String read_next_line(std::istream& ifs, Index& cur_line)
  {
    String line;
    if(ifs.eof() && ifs.good())
    {
      throw SyntaxError("Unexpected end of file at line " + stringify(cur_line));
    }
    getline(ifs, line);
    ++cur_line;
    return line.trim();
  }
  /// \endcond

  // default constructor
  MeshStreamer::MeshStreamer() :
    _num_charts(0),
    _num_meshparts(0),
    _info(""),
    _root_mesh_node(nullptr),
    _num_parsed_charts(0),
    _num_parsed_meshparts(0)
  {
    CONTEXT("MeshStreamer::MeshStreamer()");
  }

  MeshStreamer::MeshStreamer(const String& filename) :
    _num_charts(0),
    _num_meshparts(0),
    _info(""),
    _root_mesh_node(nullptr),
    _num_parsed_charts(0),
    _num_parsed_meshparts(0)
  {
    CONTEXT("MeshStreamer::MeshStreamer()");
    parse_mesh_file(filename);
  }

  MeshStreamer::~MeshStreamer()
  {
    CONTEXT("MeshStreamer::~MeshStreamer()");
    if(_root_mesh_node != nullptr)
    {
      delete _root_mesh_node;
    }
  }

  // returns the sub mesh parent specified by "parent_name"
  MeshStreamer::MeshNode* MeshStreamer::_find_meshpart_parent(String parent_name)
  {
    CONTEXT("MeshStreamer::_find_meshpart_parent");
    if(parent_name == "root")
      return _root_mesh_node;
    if(_root_mesh_node != nullptr)
      return _root_mesh_node->find_meshpart(parent_name);
    return nullptr;
  }

  void MeshStreamer::parse_multiple_files(std::deque<String>& filenames)
  {
    CONTEXT("MeshStreamer::parse_multiple_files(std::deque<String>&)");

    if(filenames.size() == 0)
      throw InternalError("No files specified for parsing!");

    // loop over all filenames
    for(auto it = filenames.begin(); it != filenames.end(); ++it)
    {
      parse_mesh_file(*it);
    }

    // Sanity checks
    if(_root_mesh_node == nullptr)
      throw InternalError("No root mesh found in files!");
    if(_num_parsed_charts != _num_charts)
      throw InternalError("Expected "+stringify(_num_charts)+" charts, but parsed "+stringify(_num_parsed_charts));
    if(_num_parsed_meshparts!= _num_meshparts)
      throw InternalError("Expected "+stringify(_num_meshparts)+" meshparts, but parsed "+stringify(_num_parsed_meshparts));
  }

  // parses the FEAST- mesh file given by filename
  void MeshStreamer::parse_mesh_file(String filename)
  {
    CONTEXT("MeshStreamer::parse_mesh_file(String)");

    // try to open the file
    std::ifstream ifs(filename.c_str(), std::ios::in);

    // if something went wrong
    if(!ifs.is_open())
    {
      throw FileNotFound(filename);
    }

    // parsing
    try
    {
      parse_mesh_file(ifs);
      ifs.close();
    }
    catch(SyntaxError& exc)
    {
      // If the exception does not contain a filename, we'll recycle the exception and include our filename now.
      if(exc.get_filename().empty())
      {
        throw(SyntaxError(exc.message(), filename));
      }
      else
      {
        throw;
      }
    }
  } // MeshStreamer::parse_mesh_file(String filename)


  // parses the mesh data stream given by ifs
  void MeshStreamer::parse_mesh_file(std::istream& ifs)
  {
    CONTEXT("MeshStreamer::parse_mesh_file(std::ifstream&)");

    // a string containing the current line
    String line = "42";
    line.reserve(256);

    // auxiliary variable that counts the lines
    Index cur_line = 0;

    line = read_next_line(ifs,cur_line);

    // first line must be "<feat_domain_file>"
    if(line != "<feat_domain_file>" || ifs.eof() || !ifs.good())
    {
      throw SyntaxError("Unknown file format. Expected <feat_domain_file>.");
    }
    line = read_next_line(ifs, cur_line);

    // header chunk
    if(!(line == "<header>"))
    {
      throw SyntaxError("Expected <header> in line " + stringify(cur_line));
    }
    cur_line = _parse_header_section(cur_line, ifs);

    line = read_next_line(ifs, cur_line);

    if(line == "<info>")
    {
      _info = _parse_info_section(cur_line, ifs);
      line = read_next_line(ifs, cur_line);
    }

    // loop over all lines until we reach the end of the mesh file
    while(!ifs.eof() && ifs.good())
    {
      // mesh chunk
      if(line == "<mesh>")
      {
        if(_root_mesh_node != nullptr)
          throw InternalError("Encountered second <mesh> section in line "+stringify(cur_line)+", only one is allowed!");
        cur_line = _parse_mesh_section(cur_line, false, ifs);

        // insert basic root mesh information
        _root_mesh_node->mesh_data.name = "root";
        _root_mesh_node->mesh_data.parent = "none";
        _root_mesh_node->mesh_data.chart = "none";

        // get a line
        line = read_next_line(ifs,cur_line);

        // trim whitespaces; throw error if the current one is empty
        if(line.trim_me().empty())
          throw SyntaxError("No empty lines allowed in line " + stringify(cur_line));

      }
      // if it is a meshpart chunk
      else if(line == "<meshpart>")
      {
        _num_parsed_meshparts++;
        cur_line = _parse_mesh_section(cur_line, true, ifs);
        line = read_next_line(ifs,cur_line);
      }

      // If it is a chart chunk
      else if(line == "<chart>")
      {
        _num_parsed_charts++;

        // Create new ChartContainer and fill it
        MeshStreamer::ChartContainer chart_container;
        cur_line = chart_container._parse_chart_section(cur_line, ifs);
        charts.push_back(chart_container);

        // Get first line after the chart chunk
        line = read_next_line(ifs,cur_line);
      }

      // if it is the end
      else if(line == "</feat_domain_file>")
        return;
      else
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    } // end while

    // if the file did not end properly
    if(line != "</feat_domain_file>" && !(!ifs.eof() && ifs.good()))
      throw SyntaxError("Reached end of file but expected </feat_domain_file> at line " + stringify(cur_line));

  } // MeshStreamer::parse_mesh_file(std::istream& ifs)


  // parses header-section streams
  Index MeshStreamer:: _parse_header_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::_parse_header_section");

    // string the current line is saved in
    String line("");
    std::vector<String> line_vec;
    // For checking the version number
    Index version(0);
    // This will hold all the lines in the chart header section
    std::map<String, String> my_data;

    // Get first line of the chart's header section
    line = read_next_line(ifs,cur_line);

    // Parse the whole header section into the map
    while(!ifs.eof() && ifs.good() && !(line == "</header>"))
    {
      line.split_by_charset(line_vec);

      if(line_vec.size() != 2)
        throw SyntaxError("Line "+stringify(cur_line)+" does not contain exactly 2 tokens!");

      // Check if the token is already present. If not, add it.
      const auto& tmp(my_data.find(line_vec.at(0)));
      if(tmp == my_data.end())
        my_data.insert(std::pair<String, String>(line_vec.at(0), line_vec.at(1)));
      else
        throw SyntaxError("Duplicate identifier "+line_vec.at(0)+" found in line "+stringify(cur_line));

      line = read_next_line(ifs,cur_line);

      // Note that this can be done by using emplace, but GCC of version < 4.8.0 does not implement this
      //const auto& tmp = my_data.emplace(line_vec.at(0), line_vec.at(1));
      //if(!tmp.second)
      //  throw SyntaxError("Duplicate identifier "+line_vec.at(0)+" found in line "+stringify(cur_line));
      //line = read_next_line(ifs,cur_line);
    }

    // Check if the map contains the required version token
    const auto& version_it(my_data.find("version"));
    if(version_it == my_data.end())
      throw SyntaxError("No version found in chart header section!");

    version_it->second.parse(version);

    if(version != 1)
      throw InternalError("Only version 1 is supported at the moment, but " + version_it->second + " found!");
    my_data.erase(version_it);

    // Check if the header contains a number of meshparts in the file
    const auto& num_is_meshpartes_it(my_data.find("meshparts"));
    if(num_is_meshpartes_it != my_data.end())
    {
      Index tmp(0);
      num_is_meshpartes_it->second.parse(tmp);
      _num_meshparts += tmp;
      my_data.erase(num_is_meshpartes_it);
    }

    // Check if the header contains a number of meshparts in the file
    const auto& num_charts_it(my_data.find("charts"));
    if(num_charts_it != my_data.end())
    {
      Index tmp(0);
      num_charts_it->second.parse(tmp);
      _num_charts += tmp;
      my_data.erase(num_charts_it);
    }

    // Check if the map contains anything else; this should not happen
    if(my_data.size() > 0)
    {
      String msg("");
      for(auto& it:my_data)
        msg += (it.first + " " + it.second + " ");

      throw InternalError("File header section contained unrecognised entries: "+msg);
    }

    return cur_line;

  } // Index MeshStreamer::_parse_header_section(Index cur_line, std::istream& ifs)

  // parses/ignores the given info-section stream
  String MeshStreamer::_parse_info_section(Index& cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::_parse_info_section");

    // string the current line is saved in
    String line;
    String infoline("");

    while(1)
    {
      // get a line
      line = read_next_line(ifs,cur_line);

      if(line != "</info>" && !ifs.eof() && ifs.good())
      {
        infoline+=line+"\n";
      }
      else
      {
        break;
      }
    }
    infoline.trim_me();
    return infoline;
  } // Index MeshStreamer::_parse_info_section(Index cur_line, std::istream& ifs)

  // parses the given mesh-section stream
  Index MeshStreamer::_parse_mesh_section(Index cur_line, bool is_meshpart, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::_parse_mesh_section");

    // create a new mesh node
    MeshNode* mesh_node = new MeshNode();
    if(!is_meshpart)
    {
      ASSERT_(_root_mesh_node == nullptr);
      _root_mesh_node = mesh_node;
    }

    cur_line = (mesh_node->mesh_data)._parse_mesh_section(cur_line, is_meshpart, ifs);

    if(is_meshpart)
    {
      MeshNode* parent = _find_meshpart_parent((mesh_node->mesh_data).parent);
      ASSERT_(parent != nullptr);
      parent->meshpart_map.insert(std::make_pair((mesh_node->mesh_data).name, mesh_node));
    }

    // return number of lines read so far
    return cur_line;

  } // MeshStreamer::parse_mesh_section(Index cur_line, std::istream& ifs)

  // inserts is_meshpart into the tree structure
  void MeshStreamer::_insert_meshpart(MeshNode* mesh_node)
  {
    CONTEXT("MeshStreamer::_insert_is_meshpart");

    if((mesh_node->mesh_data).parent == "root")
    {
      ASSERT_(_root_mesh_node != nullptr);
      _root_mesh_node->meshpart_map.insert(std::make_pair((mesh_node->mesh_data).name, mesh_node));
    }
    else
    {
      MeshNode* parent = _root_mesh_node->find_meshpart((mesh_node->mesh_data).parent);
      ASSERT_(parent != nullptr);
      parent->meshpart_map.insert(std::make_pair((mesh_node->mesh_data).name, mesh_node));
    }
    _num_meshparts += mesh_node->get_num_meshparts_below() + 1;
  } // MeshStreamer::_insert_is_meshpart(MeshNode* mesh_node)

  // deletes is_meshpart from the tree structure
  void MeshStreamer::_delete_meshpart(MeshNode* mesh_node)
  {
    CONTEXT("MeshStreamer::_delete_meshpart");
    _delete_meshpart(mesh_node->mesh_data.name);
  } // MeshStreamer::_delete_meshpart(MeshNode* mesh_node)

  // deletes is_meshpart from the tree structure
  void MeshStreamer::_delete_meshpart(String name)
  {
    CONTEXT("MeshStreamer::_delete_meshpart");
    if(name == "root")
    {
      ASSERT_(_root_mesh_node != nullptr);
      _num_meshparts = 0;
      delete _root_mesh_node;
    }
    else
    {
      MeshNode* mesh_node = _root_mesh_node->find_meshpart(name);
      ASSERT_(mesh_node != nullptr);
      _num_meshparts -= (mesh_node->get_num_meshparts_below() + 1);

      MeshNode* parent;
      if((mesh_node->mesh_data).parent == "root")
      {
        parent = _root_mesh_node;
      }
      else
      {
        parent = _root_mesh_node->find_meshpart((mesh_node->mesh_data).parent);
      }

      MeshNode::MeshpartMap::iterator it(parent->meshpart_map.begin()), jt(parent->meshpart_map.end());
      for(; it != jt; ++it)
      {
        if(it->first.compare_no_case(name) == 0)
        {
          parent->meshpart_map.erase(it);
          break;
        }
      }
      delete mesh_node;
    }
  } // MeshStreamer::_delete_meshpart(String name)

  // returns the number of is_meshpartes
  Index MeshStreamer::get_num_meshparts() const
  {
    CONTEXT("MeshStreamer::get_num_meshparts()");
    return _num_meshparts;
  }

  // returns the global mesh information
  String MeshStreamer::get_info() const
  {
    CONTEXT("MeshStreamer::get_info()");
    return _info;
  }


  // returns a pointer to the MeshDataContainer specified by "name"
  MeshStreamer::MeshDataContainer* MeshStreamer::get_mesh(String name)
  {
    CONTEXT("MeshStreamer::get_mesh()");
    if(_root_mesh_node == nullptr)
      return nullptr;
    if(name.compare_no_case("root") == 0)
      return &_root_mesh_node->mesh_data;
    MeshNode* node = _root_mesh_node->find_meshpart(name);
    return (node != nullptr) ? &node->mesh_data : nullptr;
  }

  // writes the data into a file specified by "filename"
  void MeshStreamer::write_mesh_file(String filename) const
  {
    CONTEXT("MeshStreamer::write_mesh_file()");
    std::ofstream ofs (filename.c_str());
    // if something went wrong
    if(!ofs.is_open())
    {
      throw FileNotFound(filename);
    }

    write_mesh_file(ofs);

    ofs.close();
  }

  // writes the data into the ostream
  void MeshStreamer::write_mesh_file(std::ostream& ofs) const
  {
    CONTEXT("MeshStreamer::write_mesh_file()");

    // FILE section
    ofs << "<feat_domain_file>" << std::endl;

    // HEADER section
    ofs << "<header>" << std::endl;
    ofs << " version " << "1" << std::endl;
    ofs << " meshparts " << _num_meshparts << std::endl;
    ofs << " charts " << _num_charts << std::endl;
    // END OF HEADER section
    ofs << "</header>" << std::endl;

    if(!_info.empty())
    {
      // INFO section
      ofs << "<info>" << std::endl;
      ofs << " " <<_info << std::endl;
      // END OF INFO section
      ofs << "</info>" << std::endl;
    }

    // drop mesh data
    _root_mesh_node->write(ofs,false);

    // Write chart data
    for(auto& it : charts)
      it.write(ofs);


    // END OF FILE section
    ofs << "</feat_domain_file>";
  }

  String MeshStreamer::BaseContainer::_parse_info_section(Index& cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::BaseContainer::_parse_info_section");

    // string the current line is saved in
    String line;
    String infoline("");

    while(1)
    {
      // get a line
      line = read_next_line(ifs,cur_line);

      if(line != "</info>" && !ifs.eof() && ifs.good())
      {
        infoline+=line+"\n";
      }
      else
      {
        break;
      }
    }
    infoline.trim_me();
    return infoline;
  } // Index MeshStreamer::BaseContainer::_parse_info_section(Index cur_line, std::istream& ifs)

  // parses a counts subchunk
  Index MeshStreamer::BaseContainer::_parse_counts_chunk(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::BaseContainer::_parse_counts_section");

    // string the current line is saved in
    String line;
    std::vector<String> line_vec;

    // default setting: zero
    for(int i(0); i < 4; ++i)
      num_entities[i] = Index(0);

    if(!(slices.size() == 0))
    {
      Index temp_holder;
      line = read_next_line(ifs,cur_line);
      line.split_by_charset(line_vec);

      if(!(line_vec.at(0) == "slices"))
      {
        throw SyntaxError("Unknown file format. Expected slices in line " + stringify(cur_line));
      }
      (line_vec.at(1)).parse(temp_holder);
      slices.at(0) = temp_holder;
      num_entities[0] = temp_holder+1;
      for(Index Iter(2);Iter < line_vec.size();++Iter)
      {
        (line_vec.at(Iter)).parse(temp_holder);
        slices.push_back(temp_holder);
        num_entities[0] = num_entities[0]*(temp_holder + 1);
      }
      if(!(slices.size() == coord_per_vertex))
      {
        throw SyntaxError("Number of slices does not match the dimension of the mesh in line " + stringify(cur_line));
      }

      line = read_next_line(ifs,cur_line);
      line.split_by_charset(line_vec);

      if(!(line == "</counts>"))
      {
        throw SyntaxError("Unknown file format. Expected </counts> in line " + stringify(cur_line));
      }
    }
    else
    {
      while(!ifs.eof() && ifs.good())
      {
        // get a line
        line = read_next_line(ifs,cur_line);
        line.split_by_charset(line_vec);

        // vertex number
        if(line_vec[0] == "verts")
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing vertex number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(num_entities[0]);
        }

        // edge number
        else if(line_vec[0] == "edges")
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing edge number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(num_entities[1]);
        }

        // trias / quads number
        else if((line_vec[0] == "trias") || (line_vec[0] == "quads") )
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing triangle or quad number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(num_entities[2]);
        }

        // tetra /hexa number
        else if((line_vec[0] == "tetras") || (line_vec[0] == "hexas") )
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing tetra / hexa number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(num_entities[3]);
        }

        // if it is the end of the counts sub chunk
        else if(line == "</counts>")
        {
          break;
        }
        else
        {
          throw SyntaxError("Unknown file format in line " + stringify(cur_line));
        }
      }
    }
    // return number of read lines
    return cur_line;
  } // MeshStreamer::BaseContainer::_parse_counts_chunk


  // parses a parent-indices chunk
  Index MeshStreamer::BaseContainer::_parse_parents_chunk(Index cur_line, std::istream& ifs, String line)
  {
    CONTEXT("MeshStreamer::BaseContainer::_parse_parents_chunk");

    line.trim_me();
    std::vector<Index> idx;
    Index dim = 0, shape = 0;

    // get the dimension of the entity
    if(line == "<vert_idx>")
    {
      dim = 0;
    }
    else if (line == "<edge_idx>")
    {
      dim = 1;
    }
    else if (line == "<tria_idx>")
    {
      dim = 2;
      shape = 1;
    }
    else if (line == "<quad_idx>")
    {
      dim = 2;
      shape = 2;
    }
    else if (line == "<tetra_idx>")
    {
      dim = 3;
      shape = 1;
    }
    else if (line == "<hexa_idx>")
    {
      dim = 3;
      shape = 2;
    }
    else
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    while(!ifs.eof() && ifs.good())
    {

      // get a line
      line = read_next_line(ifs,cur_line);

      Index current_value;

      // if it is the end of the parent index chunk
      if(line.find("idx") != String::npos)
      {
        switch (dim)
        {
          case 0:
            if (! (line == "</vert_idx>"))
            {
              throw SyntaxError("Unknown file format in line " + stringify(cur_line));
            }
            break;
          case 1:
            if (! (line == "</edge_idx>"))
            {
              throw SyntaxError("Unknown file format in line " + stringify(cur_line));
            }
            break;
          case 2:
            if ( shape == 1)
            {
              if (! (line == "</tria_idx>"))
              {
                throw SyntaxError("Unknown file format in line " + stringify(cur_line));
              }
            }
            if ( shape == 2)
            {
              if (! (line == "</quad_idx>"))
              {
                throw SyntaxError("Unknown file format in line " + stringify(cur_line));
              }
            }
            break;

          case 3:
            if ( shape == 1)
            {
              if (! (line == "</tetra_idx>"))
              {
                throw SyntaxError("Unknown file format in line " + stringify(cur_line));
              }
            }
            if ( shape == 2)
            {
              if (! (line == "</hexa_idx>"))
              {
                throw SyntaxError("Unknown file format in line " + stringify(cur_line));
              }
            }
            break;
        }
        parent_indices[dim] = idx;
        break;
      }
      else
      {
        if(!line.parse(current_value))
        {
          throw SyntaxError("Unknown file format in line " + stringify(cur_line));
        }
        idx.push_back(current_value);
      }
    }

    // Sanity check
    if(idx.size() != num_entities[dim])
      throw InternalError("Parsed number of parents of dimension "+stringify(dim)+" does not match num_entities!");

    // return number of lines read so far
    return cur_line;
  } // MeshStreamer::BaseContainer::_parse_parents_chunk

  Index MeshStreamer::MeshDataContainer::_parse_mesh_header_section(Index cur_line, bool is_meshpart, std::istream& ifs)
  {
    // string the current line is saved in
    String line("");
    std::vector<String> line_vec;
    // This will hold all the lines in the mesh header section
    std::map<String, String> my_data;

    // Get first line of the mesh's header section
    line = read_next_line(ifs,cur_line);

    // Parse the whole header section into the map
    while(!ifs.eof() && ifs.good() && !(line == "</header>"))
    {
      line.split_by_charset(line_vec);

      if(line_vec.size() != 2)
        throw SyntaxError("Line "+stringify(cur_line)+" does not contain exactly 2 tokens!");

      // Check if the token is already present. If not, add it.
      const auto& tmp(my_data.find(line_vec.at(0)));
      if(tmp == my_data.end())
        my_data.insert(std::pair<String, String>(line_vec.at(0), line_vec.at(1)));
      else
        throw SyntaxError("Duplicate identifier "+line_vec.at(0)+" found in line "+stringify(cur_line));

      line = read_next_line(ifs,cur_line);

      // Note that this can be done by using emplace, but GCC of version < 4.8.0 does not implement this
      //const auto& tmp = my_data.emplace(line_vec.at(0), line_vec.at(1));
      //if(!tmp.second)
      //  throw SyntaxError("Duplicate identifier "+line_vec.at(0)+" found in line "+stringify(cur_line));
      //line = read_next_line(ifs,cur_line);
    }

    // Check if the map contains the required type token
    const auto& type_it(my_data.find("type"));
    if(type_it == my_data.end())
        throw SyntaxError("No type information found in mesh header section!");
    else
    {
      mesh_type = convert_mesh_type(type_it->second);

      if(mesh_type == mt_structured)
        slices.push_back(0);

      my_data.erase(type_it);
    }

    // Check if the map contains the required shape token
    const auto& shape_it(my_data.find("shape"));
    if(shape_it == my_data.end())
      throw SyntaxError("No shape information found in mesh header section!");
    else
    {
      shape_type = convert_shape_type(shape_it->second);
      my_data.erase(shape_it);
    }

    // Check if the map contains the attributes token
    auto attributes_it(my_data.find("attribute_sets"));
    // If the information is not there, default to 0
    if(attributes_it == my_data.end())
      attribute_count = 0;
    else
    {
      attributes_it->second.parse(attribute_count);
      my_data.erase(attributes_it);
    }


    // If this is not a meshpart, we need the number of coordinates per vertex
    if(!is_meshpart)
    {
      // Check if the map contains the coords token
      auto coords_it(my_data.find("coords"));
      if(coords_it == my_data.end())
        throw SyntaxError("No coords information found in mesh header section!");
      else
      {
        coords_it->second.parse(coord_per_vertex);
        my_data.erase(coords_it);
      }
    }
    // If this is a meshpart, the number of coordinates per vertex is optional, but we need
    // name and parent information, and chart information might be present
    else
    {
      // Check if the map contains the coords token
      auto coords_it(my_data.find("coords"));
      // If the information is not there, default to 0
      if(coords_it == my_data.end())
        coord_per_vertex = 0;
      else
      {
        coords_it->second.parse(coord_per_vertex);
        my_data.erase(coords_it);
      }

      // Check if the map contains the chart token
      auto chart_it(my_data.find("chart"));
      // If this is a meshpart, we check for the chart token
      if(chart_it != my_data.end())
      {
        chart = chart_it->second;
        my_data.erase(chart_it);
      }

      // Check if the map contains the name token
      auto name_it(my_data.find("name"));
      if(name_it == my_data.end())
        throw SyntaxError("No name information found in meshpart header section!");
      else
      {
        name = name_it->second;
        my_data.erase(name_it);
      }

      // Check if the map contains the parent token
      auto parent_it(my_data.find("parent"));
      // If this is a meshpart, we check for the parent token
      if(parent_it == my_data.end())
        throw SyntaxError("No parent information found in meshpart header section!");
      else
      {
        parent = parent_it->second;
        my_data.erase(parent_it);
      }

    } // if(!is_meshpart)

    // Only Hypercube shapes are allowed for structured meshes
    if(mesh_type == mt_structured &&
        (shape_type == MeshDataContainer::st_tria || shape_type == MeshDataContainer::st_tetra))
        {
          throw SyntaxError("Unsupported shape_type for structured mesh, in line " + stringify(cur_line));
        }

    // Check if the map contains anything else; this should not happen
    if(my_data.size() > 0)
    {
      String msg("");
      for(auto& it:my_data)
        msg += (it.first + " " + it.second + " ");

      throw InternalError("File header section contained unrecognised entries: "+msg);
    }

    return cur_line;

  } //MeshStreamer::MeshDataContainer::_parse_mesh_header_section(Index, std::istream&)

  // parses the given mesh-section stream
  Index MeshStreamer::MeshDataContainer::_parse_mesh_section(Index cur_line, bool is_meshpart, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::MeshDataContainer::_parse_mesh_section");

    bool coordfile = false, adjfile = false;

    // string the current line is saved in
    String line;
    std::vector<String> line_vec;

    // if it is the root mesh
    String break_line = is_meshpart ? "</meshpart>" : "</mesh>";

    coord_per_vertex = 0;

    // get a line
    line = read_next_line(ifs,cur_line);
    if(!(line == "<header>"))
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));

    _parse_mesh_header_section(cur_line, is_meshpart, ifs);

    // Get the first line after the end of the header section
    line = read_next_line(ifs,cur_line);

    // if it is an info sub chunk
    if(line == "<info>")
    {
      info = _parse_info_section(cur_line, ifs);
      line = read_next_line(ifs,cur_line);
    } // info sub-chunk

    // if it is the counts sub chunk
    if(line == "<counts>")
    {
      cur_line = _parse_counts_chunk(cur_line, ifs);
      line = read_next_line(ifs,cur_line);
    }
    else if (!adjfile)
    {
      throw SyntaxError("Unknown file format. Expected <counts> in line " + stringify(cur_line));
    }// counts sub-chunk

    // if it is a coords sub chunk
    if(line == "<coords>")
    {
      cur_line = _parse_coords_chunk(cur_line, ifs);
      line = read_next_line(ifs,cur_line);
    }
    else if(!coordfile && coord_per_vertex!=0)
    {
      throw SyntaxError("Unknown file format. Expected <coords> in line " + stringify(cur_line));
    }// coords sub-chunk

    while(!ifs.eof() && ifs.good())
    {
      // if it is an adjacency sub chunk
      if(line.find('@') != std::string::npos)
      {
        if(mesh_type == mt_structured)
        {
          throw SyntaxError("No adjacency chunk allowed for structured meshes in line " + stringify(cur_line));
        }
        cur_line = _parse_adjacency_chunk(cur_line, ifs, line);
      }// adjacencies sub-chunk

      // if it is a parent index chunk (of a meshpart)
      else if((line == "<vert_idx>"||
               line == "<edge_idx>"||
               line == "<tria_idx>"||
               line == "<quad_idx>"||
               line == "<tetra_idx>"||
               line == "<hexa_idx>") && is_meshpart)
      {
        cur_line = _parse_parents_chunk(cur_line, ifs, line);
      } // parent index sub chunk
      else if(line=="<attribute>")
      {
        cur_line = _parse_attribute_chunk(cur_line, ifs);
      }

      // if it is the end of the mesh section
      else if(line == break_line)
      {
        break;
      }
      else
      {
        throw SyntaxError("Unknown file format. Expected </mesh> or </meshpart> in line " + stringify(cur_line));
      }

      line = read_next_line(ifs,cur_line);
    } // while


    if(line != break_line && !(!ifs.eof() && ifs.good()))
    {
      throw SyntaxError("Reached end of file but expected " + break_line + " at line " + stringify(cur_line));
    }

    // Sanity check: Does the number of parsed attributes match the specified number?
    Index real_num_attributes(0);
    for(Index i(0); i < 4; ++i)
      real_num_attributes += Index((this->attributes[i]).size());

    // Sanity check
    if(real_num_attributes != this->attribute_count)
      throw InternalError("Parsed " + stringify(real_num_attributes) + " attribute sets but " + stringify(this->attribute_count) + " were specified!");

    // return number of lines read so far
    return cur_line;

  } // MeshStreamer::MeshDataContainer::parse_mesh_section(Index cur_line, std::istream& ifs)

  // parses the coord file given by filepath
  void MeshStreamer::MeshDataContainer::parse_coord_file(String filename)
  {
    CONTEXT("MeshStreamer::MeshDataContainer::parse_coord_file(String filename)");

    // try to open the file
    std::ifstream ifs(filename.c_str(), std::ios::in);

    // if something went wrong
    if(!ifs.is_open())
    {
      throw FileNotFound(filename);
    }

    // parsing
    try
    {
      parse_coord_file(ifs);
      ifs.close();
    }
    catch(SyntaxError& exc)
    {
      // If the exception does not contain a filename, we'll recycle the exception and include our filename now.
      if(exc.get_filename().empty())
      {
        throw(SyntaxError(exc.message(), filename));
      }
      else
      {
        throw;
      }
    }
  } // MeshStreamer::MeshDataContainer::parse_coord_file(String filename, MeshStreamer::MeshDataContainer *mesh)

  // parses the coord-file stream given by ifs
  void MeshStreamer::MeshDataContainer::parse_coord_file(std::istream& ifs)
  {
    CONTEXT("MeshStreamer::MeshDataContainer::parse_coord_file(std::ifstream&)");

    // a string containing the current line
    String line = "42";
    std::vector<String> line_vec;
    line.reserve(256);

    // auxiliary variables
    String coord;
    Index cur_line = 0;

    line = read_next_line(ifs,cur_line);

    if(line != "<feast_coord_file>" && !ifs.eof() && ifs.good())
    {
      throw SyntaxError("Unknown coordfile format. Expected <feast_coord_file> in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);

    // if it is the header chunk
    if(!(line == "<header>"))
    {
      throw SyntaxError("Unknown coordfile format. Expected <header> in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "version")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing version number in coordfile in line " + stringify(cur_line));
      }
      if(!(line_vec[1] == "1"))
      {
        throw SyntaxError("Wrong coord_version in coordfile. Coord_version should be 1 in line " + stringify(cur_line));
      }
    }
    else
    {
      throw SyntaxError("Unknown coordfile format. Expected version in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "verts")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing number of vertices in coordfile in line " + stringify(cur_line));
      }
      line_vec[1].parse(num_entities[0]);
    }
    else
    {
      throw SyntaxError("Unknown coordfile format. Expected verts in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "coords")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing coordinate number in coordfile in line " + stringify(cur_line));
      }
      line_vec[1].parse(coord_per_vertex);
    }
    else
    {
      throw SyntaxError("Unknown coordfile format. Expected coords in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);

    if(!(line == "</header>"))
    {
      throw SyntaxError("Unknown coordfile format. Expected </header> in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);

    if(line == "<info>")
    {
      info = _parse_info_section(cur_line, ifs);
      line = read_next_line(ifs,cur_line);
    } // info sub chunk

    if(line == "<coords>")
    {
      cur_line = _parse_coords_chunk(cur_line, ifs);
    }
    else
    {
      throw SyntaxError("Unknown coordfile format. Expected <coords> in line " + stringify(cur_line));
    }// coords sub-chunk

    line = read_next_line(ifs,cur_line);

    if(!(line == "</feast_coord_file>"))
    {
      throw SyntaxError("Unknown coordfile format. Expected </feast_coord_file> in line " + stringify(cur_line));
    }
  } // MeshStreamer::MeshDataContainer::parse_coord_file(std::istream& ifs, MeshStreamer::MeshDataContainer *mesh)

  // parses the adjacency file given by filename
  void MeshStreamer::MeshDataContainer::parse_adjacency_file(String filename)
  {
    CONTEXT("MeshStreamer::MeshDataContainer::parse_adjacency_file(String filename)");

    // try to open the file
    std::ifstream ifs(filename.c_str(), std::ios::in);

    // if something went wrong
    if(!ifs.is_open())
    {
      throw FileNotFound(filename);
    }

    // parsing
    try
    {
      parse_adjacency_file(ifs);
      ifs.close();
    }
    catch(SyntaxError& exc)
    {
      // If the exception does not contain a filename, we'll recycle the exception and include our filename now.
      if(exc.get_filename().empty())
      {
        throw(SyntaxError(exc.message(), filename));
      }
      else
      {
        throw;
      }
    }
  } // MeshStreamer::MeshDataContainer::parse_adjacency_file(String filename)


  // parses the adjacency-file-stream ifs
  void MeshStreamer::MeshDataContainer::parse_adjacency_file(std::istream& ifs)
  {
    CONTEXT("MeshStreamer::MeshDataContainer::parse_adjacency_file(std::ifstream&)");

    // a string containing the current line
    String line = "42";
    line.reserve(256);
    std::vector<String> line_vec;

    // (auxiliary) variables
    String adj, face, shape, break_line;
    Index cur_line = 0;

    line = read_next_line(ifs,cur_line);
    // ignore everything until the adjacency file begins
    if(line != "<feast_adjacency_file>" && !ifs.eof() && ifs.good())
    {
      throw SyntaxError("Beginning of adjacency file not found");
    }

    line = read_next_line(ifs,cur_line);

    if(!(line == "<header>"))
    {
      throw SyntaxError("Unknown adjacencyfile format. Expected <header> in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "version")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing version number in adjacencyfile in line " + stringify(cur_line));
      }

      if(!(line_vec[1] == "1"))
      {
        throw SyntaxError("Wrong adjacency_version. adjacency_version should be 1 in adjacencyfile in line " + stringify(cur_line));
      }
    }
    else
    {
      throw SyntaxError("Unknown adjacencyfile format. Expected version in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "type")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing type in adjacencyfile in line " + stringify(cur_line));
      }
      mesh_type = convert_mesh_type(line_vec[1]);
    }
    else
    {
      throw SyntaxError("Unknown adjacencyfile format. Expected type in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "shape")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing shape type in adjacencyfile in line " + stringify(cur_line));
      }
      shape_type = convert_shape_type(line_vec[1]);
    }
    else
    {
      throw SyntaxError("Unknown adjancencyfile format. Expected shape in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);

    if(!(line == "</header>"))
    {
      throw SyntaxError("Unknown adjacencyfile format. Expected </header> in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);

    // if it is an info chunk
    if(line == "<info>")
    {
      _parse_info_section(cur_line, ifs);

      line = read_next_line(ifs,cur_line);
    }

      // if it is the counts sub chunk
    if(line == "<counts>")
    {
      cur_line = _parse_counts_chunk(cur_line, ifs);
    }
    else
    {
      throw SyntaxError("Unknown adjacencyfile format. Expected <counts> in line " + stringify(cur_line));
    }// counts sub-chunk

    while(!ifs.eof() && ifs.good())
    {
      line = read_next_line(ifs,cur_line);

      // if it is an adjacency sub chunk
      if(line.find('@') != std::string::npos)
      {
        cur_line =  _parse_adjacency_chunk(cur_line, ifs, line);
      } // adjacencies sub-chunk

      // if it is the end of the file
      else if(line == "</feast_adjacency_file>")
      {
        return;
      }
      else
      {
        throw SyntaxError("Unknown adjacencyfile format. Expected </feast_adjacency_file> in line " + stringify(cur_line));
      }
    }

    if(!(line == "</feast_adjacency_file>"))
    {
      throw SyntaxError("Unknown adjacencyfile format. Expected </feast_adjacency_file> in line " + stringify(cur_line));
    }

  } // MeshStreamer::MeshDataContainer::parse_adjacency_file(std::istream& ifs, MeshStreamer::MeshDataContainer *mesh)

  /// Parses an attribute chunk
  Index MeshStreamer::MeshDataContainer::_parse_attribute_chunk(Index cur_line, std::istream& ifs)
  {
    Index dimension(~Index(0));
    String my_identifier("");
    int value_dim(0);
    Index value_count(0);

    String line("");
    std::vector<String> line_vec;

    line = read_next_line(ifs,cur_line);
    if(line!="<header>")
      throw SyntaxError("Attribute section does not start with header section in line " + stringify(cur_line));

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);
    if(line_vec[0] == "dimension")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing dimension number in line " + stringify(cur_line));
      }

      (line_vec.at(1)).parse(dimension);
    }
    else
      throw SyntaxError("Expected dimension in line " + stringify(cur_line));

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);
    if(line_vec[0] == "name")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing name indentifier in line " + stringify(cur_line));
      }

      (line_vec.at(1)).parse(my_identifier);
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);
    if(line_vec[0] == "value_dim")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing value_dim count in line " + stringify(cur_line));
      }

      (line_vec.at(1)).parse(value_dim);
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);
    if(line_vec[0] == "value_count")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing value_count count in line " + stringify(cur_line));
      }

      (line_vec.at(1)).parse(value_count);
    }

    line = read_next_line(ifs,cur_line);
    if(line!="</header>")
      throw SyntaxError("Expected end of attribute header in line " + stringify(cur_line));

    // Init AttributeContainer with the parsed sizes
    AttributeContainer my_attribute_container(my_identifier, value_dim, value_count);

    line = read_next_line(ifs,cur_line);
    if(line!="<values>")
      throw SyntaxError("Expected <values> in line " + stringify(cur_line));

    while(!ifs.eof() && ifs.good())
    {
      // get a line
      line = read_next_line(ifs,cur_line);

      // if it is the end of the values sub chunk
      if(line == "</values>")
      {
        break;
      }
      else
      {
        // auxiliary variables
        AttributeContainer::ValueType current_value;
        AttributeContainer::ValueVec v;

        // separate by " "
        line.split_by_charset(line_vec);

        // parse the substrings
        for(Index i(0); i < line_vec.size(); ++i)
        {
          if ( !(line_vec[i]).parse(current_value) )
          {
            throw SyntaxError("Wrong value format in line " + stringify(cur_line));
          }
          v.push_back(current_value);
        }

        // if the number of entries does not match the coord_per_vertex variable
        if (v.size() != Index(value_dim))
        {
          throw SyntaxError("Entry does not match value_dim count in line " + stringify(cur_line));
        }

        // add to the coordinate stack
        my_attribute_container.values.push_back(v);
      }
    }
    line = read_next_line(ifs,cur_line);
    if(line!="</attribute>")
      throw SyntaxError("Expected </attribute> in line " + stringify(cur_line));

    // Sanity check: Does the number of parsed values match the specified number?
    if(my_attribute_container.values.size() != my_attribute_container.value_count)
      throw InternalError("Parsed " + stringify(my_attribute_container.values.size()) + " values but " + stringify(my_attribute_container.value_count) + " were specified");

    // Add the temporay AttributeContainer to the MeshDataContainer's attributes
    (this->attributes[dimension]).push_back(my_attribute_container);

    return cur_line;

  } // MeshStreamer::MeshDataContainer::_parse_attribute_chunk

  // parses a coords subchunk
  Index MeshStreamer::MeshDataContainer::_parse_coords_chunk(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::MeshDataContainer::_parse_coords_section");

    // string the current line is saved in
    String line;

    while(!ifs.eof() && ifs.good())
    {
      // get a line
      line = read_next_line(ifs,cur_line);

      // if it is the end of the coord sub chunk
      if(line == "</coords>")
      {
        break;
      }
      else
      {
        // auxiliary variables
        std::vector<String> line_vec;
        std::vector<double> v;
        double current_value;

        // separate by " "
        line.split_by_charset(line_vec);

        // parse the substrings
        for(Index i(0); i < line_vec.size(); ++i)
        {
          if ( !(line_vec[i]).parse(current_value) )
          {
            throw SyntaxError("Wrong coordinate format in line " + stringify(cur_line));
          }
          v.push_back(current_value);
        }

        // if the number of entries does not match the coord_per_vertex variable
        if (v.size() != coord_per_vertex)
        {
          throw SyntaxError("Wrong coordinate format in line " + stringify(cur_line));
        }

        // add to the coordinate stack
        (coords).push_back(v);
      }
    }
    // Sanity check
    if(coords.size() != num_entities[0])
      throw InternalError("Number of coordinates parsed does not match with num_entites[0]");
    // current number of lines
    return cur_line;
  } //MeshStreamer::MeshDataContainer::_parse_coords_chunk


  // parses an adjacency subchunk
  Index MeshStreamer::MeshDataContainer::_parse_adjacency_chunk(Index cur_line, std::istream& ifs, String line)
  {
    CONTEXT("MeshStreamer::MeshDataContainer::_parse_adjacency_section");

    // various auxiliary variables
    Index shape_dim = 0;
    String shape;
    String::size_type found;
    String face;

    // get position of the '@'
    found = line.find('@');

    // separate face and shape
    line.trim_me();
    face = line.substr(0, found);
    shape = line.substr(found + 1);

    face.pop_front();
    shape.pop_back();

    face.trim_me();
    shape.trim_me();

    // if it is a "vert@shape" data
    if(face == "vert")
    {
      // get the shape dimension
      if(shape == "edge")
      {
        shape_dim = 1;
      }
      else if(shape == "quad" || shape == "tria")
      {
        shape_dim = 2;
      }
      else if(shape == "hexa" || shape == "tetra")
      {
        shape_dim = 3;
      }
      else
      {
        throw SyntaxError("Unknown format. Invalid shape type in line " + stringify(cur_line));
      }

      // adjacency vector
      std::vector<std::vector<Index> > a_stack;

      while(!ifs.eof() && ifs.good())
      {
        // get a line
        line = read_next_line(ifs,cur_line);

        // if it is the end of the sub chunk
        if(line.find('@') != std::string::npos)
        {
          adjacencies[0][shape_dim] = a_stack;
          break;
        }
        else
        {
          // auxiliary variables
          std::vector<String> line_vec;
          std::vector<Index> a;
          Index current_value;

          // separate by " "
          line.split_by_charset(line_vec);

          // parse the substrings
          for(Index i(0); i < line_vec.size(); ++i)
          {
            line_vec[i].trim_me();
            if( !(line_vec[i]).empty())
            {
              if( !(line_vec[i]).parse(current_value) )
              {
                throw SyntaxError("Wrong adjacency format in line " + stringify(cur_line));
              }
              a.push_back(current_value);
            }
          }

          // add the vector to the adjacency stack
          a_stack.push_back(a);
        }
      }
    }
    // if it is something else
    else
    {
      while(!ifs.eof() && ifs.good())
      {
        // get a line
        line = read_next_line(ifs,cur_line);

        // if it is the end of the sub chunk
        if(line.find('@') != std::string::npos)
        {
          break;
        }
      }
    }
    if(adjacencies[0][shape_dim].size() != num_entities[shape_dim])
      throw InternalError("Number of adjacencies for dimension "+stringify(shape_dim)+" does not match num_entities!");
    // return number of lines read so far
    return cur_line;
  } // MeshStreamer::MeshDataContainer::_parse_adjacency_chunk

  // converts meshtype (MeshType to String)
  String MeshStreamer::MeshDataContainer::convert_mesh_type(const MeshStreamer::MeshDataContainer::MeshType mesh_type_in)
  {
    if(mesh_type_in == mt_conformal)
    {
      return "conformal";
    }
    else if (mesh_type_in == mt_structured)
    {
      return "structured";
    }
    else
    {
      // todo error anpassen
      throw InternalError("Unknown meshtype");
    }
  }// String MeshStreamer::MeshDataContainer::convert_mesh_type(const MeshType mesh_type) const

  // converts meshtype (String to MeshType)
  MeshStreamer::MeshDataContainer::MeshType MeshStreamer::MeshDataContainer::convert_mesh_type(const String mesh_type_in)
  {
    if(mesh_type_in == "conformal")
    {
      return mt_conformal;
    }
    else if(mesh_type_in == "structured")
    {
      return mt_structured;
    }
    else
    {
      throw InternalError("Unknown meshtype");
    }
  } // MeshType MeshStreamer::MeshDataContainer::convert_mesh_type(const String mesh_type) const

  // converts shapetype (ShapeType to String)
  String MeshStreamer::MeshDataContainer::convert_shape_type(const MeshStreamer::MeshDataContainer::ShapeType shape_type_in)
  {
    if(shape_type_in == MeshDataContainer::st_vert)
    {
      return "vertex";
    }
    if(shape_type_in == MeshDataContainer::st_edge)
    {
      return "edge";
    }
    else if(shape_type_in == MeshDataContainer::st_tria)
    {
      return "tria";
    }
    else if(shape_type_in == MeshDataContainer::st_quad)
    {
      return "quad";
    }
    else if(shape_type_in == MeshDataContainer::st_hexa)
    {
      return "hexa";
    }
    else if(shape_type_in == MeshDataContainer::st_tetra)
    {
      return "tetra";
    }
    else
    {
      throw InternalError("Unknown shapetype");
    }
  } // String MeshStreamer::MeshDataContainer::convert_shape_type(const MeshStreamer::MeshDataContainer::ShapeType  shape_type) const

  // converts shapetype (String to ShapeType)
  MeshStreamer::MeshDataContainer::ShapeType MeshStreamer::MeshDataContainer::convert_shape_type(const String shape_type_in)
  {
    if(shape_type_in == "vertex")
    {
      return st_vert;
    }
    if(shape_type_in == "edge")
    {
      return st_edge;
    }
    else if(shape_type_in == "tria")
    {
      return st_tria;
    }
    else if(shape_type_in == "quad")
    {
      return st_quad;
    }
    else if(shape_type_in == "hexa")
    {
      return st_hexa;
    }
    else if(shape_type_in == "tetra")
    {
      return st_tetra;
    }
    else
    {
      throw InternalError("Unknown shapetype: "+ shape_type_in);
    }
  } // ShapeType MeshStreamer::MeshDataContainer::convert_shape_type(const String shape_type) const

  // Writes the mesh data of this mesh and all is_meshpartes related to this mesh
  // into the output stream.
  void MeshStreamer::MeshNode::write(std::ostream& ofs, bool is_meshpart) const
  {
    mesh_data.write(ofs, is_meshpart);

    // loop through all is_meshpartes related to this mesh and drop the data as well
    MeshpartMap::const_iterator it(meshpart_map.begin()), jt(meshpart_map.end());
    for(; it != jt; ++it)
      it->second->write(ofs, true);
  }

  // Writes the mesh data of this mesh and all is_meshpartes related to this mesh
  // into the output stream.
  void MeshStreamer::MeshDataContainer::write(std::ostream& ofs, bool is_meshpart) const
  {
    // choose the right tag
    if (is_meshpart)
      ofs << "<meshpart>" << std::endl;
    else
      ofs << "<mesh>" << std::endl;

    // header section
    ofs << " <header>" << std::endl;

    if (is_meshpart)
    {
      ofs << "  name " << name << std::endl;
      ofs << "  parent " << parent << std::endl;
      if ( !chart.empty() )
        ofs << "  chart " << chart << std::endl;
    }

    ofs << "  type " << convert_mesh_type(mesh_type) << std::endl;
    ofs << "  shape " << convert_shape_type(shape_type) << std::endl;
    if(coord_per_vertex != 0)
      ofs << "  coords " << coord_per_vertex << std::endl;
    if( attribute_count > 0)
      ofs << "  attribute_sets " << attribute_count << std::endl;
    ofs << " </header>" << std::endl;

    // info section
    if( !info.empty() )
    {
      ofs << " <info>" << std::endl;
      ofs << "  " << info << std::endl;
      ofs << " </info>" << std::endl;
    }

    // count section
    ofs << " <counts>" << std::endl;
    if(mesh_type == MeshDataContainer::mt_conformal)
    {
      if( num_entities[0] != 0)
        ofs << "  verts " << num_entities[0] << std::endl;
      if( num_entities[1] != 0)
        ofs << "  edges " << num_entities[1] << std::endl;

      if( num_entities[2] != 0)
      {
        if(shape_type == MeshDataContainer::st_quad || shape_type == MeshDataContainer::st_hexa)
          ofs << "  quads " << num_entities[2] << std::endl;
        else if(shape_type == MeshDataContainer::st_tria || shape_type == MeshDataContainer::st_tetra)
          ofs << "  trias " << num_entities[2] << std::endl;
        else
          throw InternalError("Mesh of type "+stringify(shape_type)+" has num_entities of dim 2!");
      }

      if( num_entities[3] != 0)
      {
        if(shape_type == MeshDataContainer::st_tetra)
          ofs << "  tetras " << num_entities[3] << std::endl;
        else if( shape_type == MeshDataContainer::st_hexa)
          ofs << "  hexas " << num_entities[3] << std::endl;
        else
          throw InternalError("Mesh of type "+stringify(shape_type)+" has num_entities of dim 3!");
      }
    }
    else if( mesh_type == MeshDataContainer::mt_structured)
    {
      ofs << "  slices";
      for (Index i(0); i < slices.size() ; ++i)
      {
        ofs << " " << slices.at(i);
      }
      ofs << std::endl;
    }
    ofs << " </counts>" << std::endl;

    // coord section (if any)
    if(coord_per_vertex > 0)
    {
      ofs << " <coords>" << std::endl;
      for (Index i(0); i < num_entities[0]; ++i)
      {
        for (Index j(0); j < coord_per_vertex ; ++j)
        {
          ofs  << "  " << scientify((coords[i])[j]);
        }
        ofs << std::endl;
      }
      ofs << " </coords>" << std::endl;
    }

    // adjacency section
    if(mesh_type == MeshDataContainer::mt_conformal)
    {
      if (!( (adjacencies[0][1]).size() == 0 ))
      {
        ofs << " <vert@edge>" << std::endl;
        for (Index i(0); i < (adjacencies[0][1]).size() ; ++i)
        {
          ofs << "  " << ((adjacencies[0][1])[i])[0];
          ofs << " " << ((adjacencies[0][1])[i])[1] << std::endl;
        }
        ofs << " </vert@edge>" << std::endl;
      }

      if (!( (adjacencies[0][2]).size() == 0 ))
      {
        if ( ((adjacencies[0][2])[0]).size() == 3 )
        {
          ofs << " <vert@tria>" << std::endl;
          for (Index i(0); i < (adjacencies[0][2]).size() ; ++i)
          {
            ofs << "  " << ((adjacencies[0][2])[i])[0];
            ofs << " " << ((adjacencies[0][2])[i])[1];
            ofs << " " << ((adjacencies[0][2])[i])[2] << std::endl;
          }
          ofs << " </vert@tria>" << std::endl;
        }
        else if ( ((adjacencies[0][2])[0]).size() == 4 )
        {
          ofs << " <vert@quad>" << std::endl;
          for (Index i(0); i < (adjacencies[0][2]).size() ; ++i)
          {
            ofs << "  " << ((adjacencies[0][2])[i])[0];
            ofs << " " << ((adjacencies[0][2])[i])[1];
            ofs << " " << ((adjacencies[0][2])[i])[2];
            ofs << " " << ((adjacencies[0][2])[i])[3] << std::endl;
          }
          ofs << " </vert@quad>" << std::endl;
        }
      }

      if (!( (adjacencies[0][3]).size() == 0 ))
      {
        if ( ((adjacencies[0][3])[0]).size() == 4 )
        {
          ofs << " <vert@tetra>" << std::endl;
          for (Index i(0); i < (adjacencies[0][3]).size() ; ++i)
          {
            ofs << "  " << ((adjacencies[0][3])[i])[0];
            ofs << " " << ((adjacencies[0][3])[i])[1];
            ofs << " " << ((adjacencies[0][3])[i])[2];
            ofs << " " << ((adjacencies[0][3])[i])[3] << std::endl;
          }
          ofs << " </vert@tetra>" << std::endl;
        }
        else if ( ((adjacencies[0][3])[0]).size() == 8 )
        {
          ofs << " <vert@hexa>" << std::endl;
          for (Index i(0); i < (adjacencies[0][3]).size() ; ++i)
          {
            ofs << "  " << ((adjacencies[0][3])[i])[0];
            ofs << " " << ((adjacencies[0][3])[i])[1];
            ofs << " " << ((adjacencies[0][3])[i])[2];
            ofs << " " << ((adjacencies[0][3])[i])[3];
            ofs << " " << ((adjacencies[0][3])[i])[4];
            ofs << " " << ((adjacencies[0][3])[i])[5];
            ofs << " " << ((adjacencies[0][3])[i])[6];
            ofs << " " << ((adjacencies[0][3])[i])[7] << std::endl;
          }
          ofs << " </vert@hexa>" << std::endl;
        }
      }
    }

    // if it is a meshpart, add the parent indices
    if (is_meshpart)
    {
      ofs << " <vert_idx>" << std::endl;
      for (Index i(0); i < num_entities[0]; ++i)
      {
        ofs << "  " << (parent_indices[0])[i] << std::endl;
      }
      ofs << " </vert_idx>" << std::endl;

      if ( !(num_entities[1]== 0) )
      {
        ofs << " <edge_idx>" << std::endl;
        for (Index i(0); i < num_entities[1] ; ++i)
        {
          ofs << "  " << (parent_indices[1])[i] << std::endl;
        }
        ofs << " </edge_idx>" << std::endl;
      }

      if ( !(num_entities[2] == 0) )
      {
        if(shape_type == MeshDataContainer::st_tria|| shape_type == MeshDataContainer::st_tetra)
        {
          ofs << " <tria_idx>" << std::endl;
          for (Index i(0); i < num_entities[2]; ++i)
          {
            ofs << "  " << (parent_indices[2])[i] << std::endl;
          }
          ofs << " </tria_idx>" << std::endl;
        }
        else if(shape_type == MeshDataContainer::st_quad || shape_type == MeshDataContainer::st_hexa)
        {
          ofs << " <quad_idx>" << std::endl;
          for (Index i(0); i < num_entities[2]; ++i)
          {
            ofs << "  " << (parent_indices[2])[i] << std::endl;
          }
          ofs << " </quad_idx>" << std::endl;
        }
        else
          throw InternalError("Mesh of type "+stringify(shape_type)+" has num_entities of dim 2!");
      }

      if ( !(num_entities[3]== 0) )
      {
        if(shape_type == MeshDataContainer::st_tetra)
        {
          ofs << " <tetra_idx>" << std::endl;
          for (Index i(0); i < num_entities[3]; ++i)
          {
            ofs << "  " << (parent_indices[3])[i] << std::endl;
          }
          ofs << " </tetra_idx>" << std::endl;
        }
        else if(shape_type == MeshDataContainer::st_hexa)
        {
          ofs << " <hexa_idx>" << std::endl;
          for (Index i(0); i < num_entities[3]; ++i)
          {
            ofs << "  " << (parent_indices[3])[i] << std::endl;
          }
          ofs << " </hexa_idx>" << std::endl;
        }
        else
          throw InternalError("Mesh of type "+stringify(shape_type)+" has num_entities of dim 3!");
      }
    }

    if(attribute_count > 0)
    {
      Index attributes_written(0);
      for(Index attribute_dim(0); attribute_dim < 4; ++attribute_dim)
      {
        for(auto& it:attributes[attribute_dim])
        {
          attributes_written++;
          // Write one attribute
          ofs << " <attribute>" << std::endl;
          // Write header
          ofs << "  <header>" << std::endl;
          ofs << "   dimension " << attribute_dim << std::endl;
          ofs << "   name " << it.identifier<< std::endl;
          ofs << "   value_dim " << it.value_dim << std::endl;
          ofs << "   value_count " << it.values.size() << std::endl;
          ofs << "  </header>" << std::endl;
          // Write values
          ofs << "  <values>" << std::endl;
          // Iterate over all value vectors
          for(auto& val_vec_it:it.values)
          {
            ASSERT(val_vec_it.size() == Index(it.value_dim), "Expected value of dimension " + stringify(it.value_dim) + " but got " + stringify(val_vec_it.size()));
            ofs << "   ";
            // Write all member of the value vector
            for(auto& val_it:val_vec_it)
              ofs << " " << val_it << std::endl;
          }
          ofs << "  </values>" << std::endl;
          ofs << " </attribute>" << std::endl;
        }
      }
      // Sanity check
      if(attributes_written != attribute_count)
        throw InternalError("MeshDataContainer is supposed to contain " + stringify(attribute_count) + " attributes but " + stringify(attributes_written) + " were written");
    }

    // choose the right tag to end the mesh section
    if (is_meshpart)
      ofs<<"</meshpart>"<<std::endl;
    else
      ofs<<"</mesh>"<<std::endl;

  }// write_mesh_data

  /// Parse header subsection of the chart section
  Index MeshStreamer::ChartContainer::_parse_chart_header_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::ChartContainer::_parse_chart_header_section");

    // string the current line is saved in
    String line("");
    std::vector<String> line_vec;

    // This will hold all the lines in the chart header section
    std::map<String, String> my_data;

    // Get first line of the chart's header section
    line = read_next_line(ifs,cur_line);

    // Parse the whole header section into the map
    while(!ifs.eof() && ifs.good() && !(line == "</header>"))
    {
      line.split_by_charset(line_vec);

      if(line_vec.size() != 2)
        throw SyntaxError("Line "+stringify(cur_line)+" does not contain exactly 2 tokens!");

      // Check if the token is already present. If not, add it.
      const auto& tmp(my_data.find(line_vec.at(0)));
      if(tmp == my_data.end())
        my_data.insert(std::pair<String, String>(line_vec.at(0), line_vec.at(1)));
      else
        throw SyntaxError("Duplicate identifier "+line_vec.at(0)+" found in line "+stringify(cur_line));

      line = read_next_line(ifs,cur_line);

      // Note that this can be done by using emplace, but GCC of version < 4.8.0 does not implement this
      //const auto& tmp = my_data.emplace(line_vec.at(0), line_vec.at(1));
      //if(!tmp.second)
      //  throw SyntaxError("Duplicate identifier "+line_vec.at(0)+" found in line "+stringify(cur_line));
      //line = read_next_line(ifs,cur_line);
    }

    // Check if the map contains the required name token
    const auto& name_it(my_data.find("name"));
    if(name_it == my_data.end())
      throw SyntaxError("No name found in chart header section!");
    name = name_it->second;
    my_data.erase(name_it);

    // Check if the map contains the required type token
    const auto& type_it(my_data.find("type"));
    if(type_it == my_data.end())
      throw SyntaxError("No type found in chart header section!");
    type = type_it->second;
    my_data.erase(type_it);

    // Check if the map contains anything else; this should not happen
    if(my_data.size() > 0)
    {
      String msg("");
      for(auto& it:my_data)
        msg += (it.first + " " + it.second + " ");

      throw InternalError("Chart header section contained unrecognised entries: "+msg);
    }

    return cur_line;

  } // MeshStreamer::ChartContainer::_parse_chart_header_section

  /// Parses chart section
  Index MeshStreamer::ChartContainer::_parse_chart_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::ChartContainer::_parse_chart_section");

    // string the current line is saved in
    String line("");
    std::vector<String> line_vec;

    line = read_next_line(ifs,cur_line);

    // The first line has to be the header
    if( !(line == "<header>") )
      throw SyntaxError("Expected <header> in line "+stringify(cur_line));
    else
      _parse_chart_header_section(cur_line, ifs);

    line = read_next_line(ifs,cur_line);

    // Optional: Info section
    if( (line == "<info>") )
    {
      line = read_next_line(ifs,cur_line);
      info = line.trim();
      line = read_next_line(ifs,cur_line);
      if(line != "</info>")
        throw SyntaxError("Expected end of chart info in line "+stringify(cur_line));
      line = read_next_line(ifs,cur_line);
    }

    start_of_data = cur_line;

    if(type == "discrete")
    {
      cur_line = mesh_data._parse_mesh_section(cur_line, false, ifs);
      line = read_next_line(ifs,cur_line);
      if(line != "</chart>")
        throw SyntaxError("Expected end of chart section in line "+stringify(cur_line));
    }
    else
    {
      // Parse rest of the chart section into the data deque
      while(!ifs.eof() && ifs.good() && !(line == "</chart>"))
      {
        data.push_back(line.trim());
        line = read_next_line(ifs,cur_line);
      }
    }

    return cur_line;
  } // MeshStreamer::ChartContainer::_parse_chart_section

  void MeshStreamer::ChartContainer::write(std::ostream& ofs) const
  {
    CONTEXT("MeshStreamer::ChartContainer::write()");

    // FILE section
    ofs << "<chart>" << std::endl;

    // HEADER section
    ofs << "<header>" << std::endl;
    ofs << " name " << name << std::endl;
    ofs << " type " << type << std::endl;
    // END OF HEADER section
    ofs << "</header>" << std::endl;

    if(!info.empty())
    {
      // INFO section
      ofs << "<info>" << std::endl;
      ofs << " " << info << std::endl;
      // END OF INFO section
      ofs << "</info>" << std::endl;
    }

    if(type == "discrete")
      mesh_data.write(ofs, false);
    else
    {
      for(auto& it : data)
        ofs << it << std::endl;
    }
    ofs << "</chart>" << std::endl;


  }

} //namespace FEAST
