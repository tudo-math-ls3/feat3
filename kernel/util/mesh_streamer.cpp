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
    _num_submeshes(0),
    _num_cellsets(0),
    _chart_path(""),
    _root_mesh_node(nullptr)
  {
    CONTEXT("MeshStreamer::MeshStreamer()");
  }

  // default destructor
  MeshStreamer::~MeshStreamer()
  {
    CONTEXT("MeshStreamer::~MeshStreamer()");
    if(_root_mesh_node != nullptr)
    {
      delete _root_mesh_node;
    }
  }

  // returns the parent mesh specified by "parent_name"
  MeshStreamer::CellSetParent* MeshStreamer::_find_cell_set_parent(String parent_name)
  {
    CONTEXT("MeshStreamer::_find_cell_set_parent");
    if(parent_name == "root")
      return _root_mesh_node;
    if(_root_mesh_node != nullptr)
      return _root_mesh_node->find_cell_set_parent(parent_name);
    return nullptr;
  }

  // returns the sub mesh parent specified by "parent_name"
  MeshStreamer::MeshNode* MeshStreamer::_find_sub_mesh_parent(String parent_name)
  {
    CONTEXT("MeshStreamer::_find_sub_mesh_parent");
    if(parent_name == "root")
      return _root_mesh_node;
    if(_root_mesh_node != nullptr)
      return _root_mesh_node->find_sub_mesh(parent_name);
    return nullptr;
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

    // first line must be "<feast_mesh_file>"
    if(line != "<feast_mesh_file>" || ifs.eof() || !ifs.good())
    {
      throw SyntaxError("Unknown file format. Expected <feast_mesh_file>.");
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

    // mesh chunk
    if(!(line == "<mesh>"))
    {
      throw SyntaxError("Expected <mesh> in line " + stringify(cur_line));
    }

    cur_line = _parse_mesh_section(cur_line, false, ifs);

    // insert basic root mesh information
    _root_mesh_node->mesh_data.name = "root";
    _root_mesh_node->mesh_data.parent = "none";
    _root_mesh_node->mesh_data.chart = "none";

    // loop over all lines until we reach the end of the mesh file
    while(!ifs.eof() && ifs.good())
    {

      // get a line
      line = read_next_line(ifs,cur_line);

      // trim whitespaces; throw error if the current one is empty
      if(line.trim_me().empty())
      {
        throw SyntaxError("No empty lines allowed in line " + stringify(cur_line));
      }

      // if it is a submesh chunk
      else if(line == "<submesh>")
      {
        cur_line = _parse_mesh_section(cur_line, true, ifs);
      }

      // if it is a cellset chunk
      else if(line == "<cellset>")
      {
        cur_line = _parse_cellset_section(cur_line, ifs);
      }

      // if it is the end
      else if(line == "</feast_mesh_file>")
      {
        return;
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    } // end while

    // if the file did not end properly

    if(line != "</feast_mesh_file>" && !(!ifs.eof() && ifs.good()))
    {
      throw SyntaxError("Reached end of file but expected </feast_mesh_file> at line " + stringify(cur_line));
    }

  } // MeshStreamer::parse_mesh_file(std::istream& ifs)


  // parses header-section streams
  Index MeshStreamer:: _parse_header_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::_parse_header_section");

    // string the current line is saved in
    String line;
    std::vector<String> line_vec;

    // get a line
    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);
    // version
    if(line_vec[0] == "version")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing version number in line " + stringify(cur_line));
      }
      if(!(line_vec.at(1) == "1") )
      {
        throw SyntaxError("Wrong version. Version should be 1, line " + stringify(cur_line));
      }
    }
    else
    {
      throw SyntaxError("Expected version, line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);
    // chart file
    if(line_vec[0] == "chart_file")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing chart file path in line " + stringify(cur_line));
      }
      _chart_path = line_vec.at(1);

      line = read_next_line(ifs,cur_line);
      line.split_by_charset(line_vec);
    }

    // submeshes
    if(line_vec[0] == "submeshes")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing number of submeshes in line " + stringify(cur_line));
      }
      (line_vec.at(1)).parse(_num_submeshes);
    }
    else
    {
      throw SyntaxError("Expected number of submeshes, line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    // cellsets
    if(line_vec[0] == "cellsets")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing number of cellsets in line " + stringify(cur_line));
      }
      (line_vec.at(1)).parse(_num_cellsets);
    }
    else
    {
      throw SyntaxError("Expected number of cellsets, line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);

    // if it is the end of the header section
    if(!(line == "</header>"))
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    // return the number of lines
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

  // parses the cellset-section stream ifs
  Index MeshStreamer::_parse_cellset_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::_parse_cellset_section");

    // (auxiliary) variables
    String line;

    // create a new cell-set node
    CellSetNode* cell_set_node = new CellSetNode();

    // cellset data container of the root mesh
    CellSetContainer& current_set(cell_set_node->cell_set);

    cur_line = current_set._parse_cellset_section(cur_line, ifs);

    CellSetParent* parent = _find_cell_set_parent(current_set.parent);
    ASSERT_(parent != nullptr);
    parent->cell_set_map.insert(std::make_pair(current_set.name, cell_set_node));

    // return number of lines read so far
    return cur_line;
  } // Index MeshStreamer::_parse_cellset_section(Index cur_line, std::istream& ifs)


  // parses the given mesh-section stream
  Index MeshStreamer::_parse_mesh_section(Index cur_line, bool submesh, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::_parse_mesh_section");

    // create a new mesh node
    MeshNode* mesh_node = new MeshNode();
    if(!submesh)
    {
      ASSERT_(_root_mesh_node == nullptr);
      _root_mesh_node = mesh_node;
    }

    cur_line = (mesh_node->mesh_data)._parse_mesh_section(cur_line, submesh, ifs);

    if(submesh)
    {
      MeshNode* parent = _find_sub_mesh_parent((mesh_node->mesh_data).parent);
      ASSERT_(parent != nullptr);
      parent->sub_mesh_map.insert(std::make_pair((mesh_node->mesh_data).name, mesh_node));
    }

    // return number of lines read so far
    return cur_line;

  } // MeshStreamer::parse_mesh_section(Index cur_line, std::istream& ifs)

  // inserts submesh into the tree structure
  void MeshStreamer::_insert_sub_mesh(MeshNode* mesh_node)
  {
    CONTEXT("MeshStreamer::_insert_submesh");

    if((mesh_node->mesh_data).parent == "root")
    {
      ASSERT_(_root_mesh_node != nullptr);
      _root_mesh_node->sub_mesh_map.insert(std::make_pair((mesh_node->mesh_data).name, mesh_node));
    }
    else
    {
      MeshNode* parent = _root_mesh_node->find_sub_mesh((mesh_node->mesh_data).parent);
      ASSERT_(parent != nullptr);
      parent->sub_mesh_map.insert(std::make_pair((mesh_node->mesh_data).name, mesh_node));
    }
    _num_submeshes += mesh_node->get_num_sub_meshes_below() + 1;
  } // MeshStreamer::_insert_submesh(MeshNode* mesh_node)

  // deletes submesh from the tree structure
  void MeshStreamer::_delete_sub_mesh(MeshNode* mesh_node)
  {
    CONTEXT("MeshStreamer::_delete_sub_mesh");
    _delete_sub_mesh(mesh_node->mesh_data.name);
  } // MeshStreamer::_delete_sub_mesh(MeshNode* mesh_node)

  // deletes submesh from the tree structure
  void MeshStreamer::_delete_sub_mesh(String name)
  {
    CONTEXT("MeshStreamer::_delete_sub_mesh");
    if(name == "root")
    {
      ASSERT_(_root_mesh_node != nullptr);
      _num_submeshes = 0;
      delete _root_mesh_node;
    }
    else
    {
      MeshNode* mesh_node = _root_mesh_node->find_sub_mesh(name);
      ASSERT_(mesh_node != nullptr);
      _num_submeshes -= (mesh_node->get_num_sub_meshes_below() + 1);

      MeshNode* parent;
      if((mesh_node->mesh_data).parent == "root")
      {
        parent = _root_mesh_node;
      }
      else
      {
        parent = _root_mesh_node->find_sub_mesh((mesh_node->mesh_data).parent);
      }

      MeshNode::SubMeshMap::iterator it(parent->sub_mesh_map.begin()), jt(parent->sub_mesh_map.end());
      for(; it != jt; ++it)
      {
        if(it->first.compare_no_case(name) == 0)
        {
          parent->sub_mesh_map.erase(it);
          break;
        }
      }
      delete mesh_node;
    }
  } // MeshStreamer::_delete_sub_mesh(String name)

  // returns the chart path
  String MeshStreamer::get_chart_path() const
  {
    CONTEXT("MeshStreamer::get_chart_path()");
    return _chart_path;
  }


  // returns the number of submeshes
  Index MeshStreamer::get_num_submeshes() const
  {
    CONTEXT("MeshStreamer::get_num_submeshes()");
    return _num_submeshes;
  }


  // returns the number of cellsets
  Index MeshStreamer::get_num_cellsets() const
  {
    CONTEXT("MeshStreamer::get_num_cellsets()");
    return _num_cellsets;
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
    MeshNode* node = _root_mesh_node->find_sub_mesh(name);
    return (node != nullptr) ? &node->mesh_data : nullptr;
  }

  // returns a pointer to the CellSetContainer specified by "name"
  MeshStreamer::CellSetContainer* MeshStreamer::get_cell_set(String name)
  {
    CONTEXT("MeshStreamer::get_cell_set()");
    if(_root_mesh_node == nullptr)
      return nullptr;
    CellSetNode* node = _root_mesh_node->find_cell_set(name);
    return (node != nullptr) ? &node->cell_set : nullptr;
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
    ofs << "<feast_mesh_file>" << std::endl;

    // HEADER section
    ofs << "<header>" << std::endl;
    ofs << " version " << "1" << std::endl;
    if ( !_chart_path.empty() )
      ofs << " chart_file " << _chart_path << std::endl;
    ofs << " submeshes " << _num_submeshes << std::endl;
    ofs << " cellsets " << _num_cellsets << std::endl;
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

    // END OF FILE section
    ofs << "</feast_mesh_file>";
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
    vertex_count = 0;
    edge_count = 0;
    quad_count = 0;
    tria_count = 0;
    tetra_count = 0;
    hexa_count = 0;

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
      vertex_count = temp_holder+1;
      for(Index Iter(2);Iter < line_vec.size();++Iter)
      {
        (line_vec.at(Iter)).parse(temp_holder);
        slices.push_back(temp_holder);
        vertex_count = vertex_count*(temp_holder + 1);
      }
      if(!(slices.size() == coord_per_vertex))
      {
        throw SyntaxError("Number of slices missmatches the dimension of the mesh in line " + stringify(cur_line));
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
          (line_vec.at(1)).parse(vertex_count);
        }

        // edge number
        else if(line_vec[0] == "edges")
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing edge number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(edge_count);
        }

        // trias number
        else if(line_vec[0] == "trias")
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing triangle number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(tria_count);
        }

        // quad number
        else if(line_vec[0] == "quads")
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing quad number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(quad_count);
        }

        // tetra number
        else if(line_vec[0] == "tetras")
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing tetra number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(tetra_count);
        }

        // hexa number
        else if(line_vec[0] == "hexas")
        {
          if(!(line_vec.size() == 2))
          {
            throw SyntaxError("Missing hexa number in line " + stringify(cur_line));
          }
          (line_vec.at(1)).parse(hexa_count);
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

    // return number of lines read so far
    return cur_line;
  } // MeshStreamer::BaseContainer::_parse_parents_chunk


  // parses the cellset-section stream ifs
  Index MeshStreamer::CellSetContainer::_parse_cellset_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::CellSetContainer::_parse_cellset_section");

    // (auxiliary) variables
    String line;
    std::vector<String> line_vec;

    // get a line
    line = read_next_line(ifs,cur_line);
    if(!(line=="<header>"))
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "name")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing name in line " + stringify(cur_line));
      }
      name = line_vec.at(1);
    }
    else
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    // if it is the parent
    if(line_vec[0] == "parent")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing parent in line " + stringify(cur_line));
      }
      parent = line_vec.at(1);
    }
    else
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);

    // end of the header section
    if(!(line == "</header>"))
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);

    if(line == "<info>")
    {
      info = _parse_info_section(cur_line, ifs);

      line = read_next_line(ifs,cur_line);
    } // info sub-chunk

    // if it is the counts sub chunk
    if(line == "<counts>")
    {
      cur_line = _parse_counts_chunk(cur_line, ifs);
    } // counts sub-chunk
    else
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    while(!ifs.eof() && ifs.good())
    {
      line = read_next_line(ifs,cur_line);
      // if it is a parent index chunk
      if(line == "<vert_idx>" ||
         line == "<edge_idx>" ||
         line == "<tria_idx>" ||
         line == "<quad_idx>" ||
         line == "<tetra_idx>" ||
         line == "<hexa_idx>")
      {
        cur_line = _parse_parents_chunk(cur_line, ifs, line);
      } // parent index sub chunk

      // if it is the end of the cellset section
      else if(line == "</cellset>")
      {
        break;
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    } // while

    if(line != "</cellset>" && !(!ifs.eof() && ifs.good()))
    {
      throw SyntaxError("Reached end of file but expected </cellset> at line " + stringify(cur_line));
    }
    // return number of lines read so far
    return cur_line;
  } // Index MeshStreamer::CellSetContainer::_parse_cellset_section(Index cur_line, std::istream& ifs)


  // parses the given mesh-section stream
  Index MeshStreamer::MeshDataContainer::_parse_mesh_section(Index cur_line, bool submesh, std::istream& ifs)
  {
    CONTEXT("MeshStreamer::MeshDataContainer::_parse_mesh_section");

    bool coordfile = false, adjfile = false;

    // string the current line is saved in
    String line;
    std::vector<String> line_vec;

    // if it is the root mesh
    String break_line = submesh ? "</submesh>" : "</mesh>";

    coord_per_vertex = 0;

    // get a line
    line = read_next_line(ifs,cur_line);
    if(!(line == "<header>"))
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(submesh)
    {
      if(line_vec[0] == "name")
      {
        if(!(line_vec.size() == 2))
        {
          throw SyntaxError("Missing name in line " + stringify(cur_line));
        }
        name = line_vec.at(1);
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }

      line = read_next_line(ifs,cur_line);
      line.split_by_charset(line_vec);

      if(line_vec[0] == "parent")
      {
        if(!(line_vec.size() == 2))
        {
          throw SyntaxError("Missing parent in line " + stringify(cur_line));
        }
        parent = line_vec.at(1);
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }

      line = read_next_line(ifs,cur_line);
      line.split_by_charset(line_vec);

      if(line_vec[0] == "chart")
      {
        if(!(line_vec.size() == 2))
        {
          throw SyntaxError("Missing chart name in line " + stringify(cur_line));
        }
        chart = line_vec.at(1);

        line = read_next_line(ifs,cur_line);
        line.split_by_charset(line_vec);
      }
    }

    // if it is the type
    if(line_vec[0] == "type")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing type in line " + stringify(cur_line));
      }
      mesh_type = convert_mesh_type(line_vec.at(1));
      if(mesh_type == mt_structured)
      {
        slices.push_back(0);
      }
    }
    else
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "shape")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing shape type in line " + stringify(cur_line));
      }
      shape_type = convert_shape_type(line_vec.at(1));
      if(mesh_type == mt_structured &&
         (shape_type == st_tria || shape_type == st_tetra ||
          shape_type == st_tria_quad || shape_type == st_tetra_hexa))
      {
        throw SyntaxError("Unsupported shape_type for structured mesh, in line " + stringify(cur_line));
      }
    }
    else
    {
      throw SyntaxError("Unknown file format in line " + stringify(cur_line));
    }

    line = read_next_line(ifs,cur_line);
    line.split_by_charset(line_vec);

    if(line_vec[0] == "coord_file")
    {
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing coordinate file path in line " + stringify(cur_line));
      }
      // parse the coord file
      parse_coord_file(line_vec.at(1));
      coordfile = true;
      line = read_next_line(ifs,cur_line);
      line.split_by_charset(line_vec);
    }

    if(!coordfile)
    {
      if(line_vec[0] == "coords")
      {
        if(!(line_vec.size() == 2))
        {
          throw SyntaxError("Missing coordinate number in line " + stringify(cur_line));
        }
        (line_vec.at(1)).parse(coord_per_vertex);
        line = read_next_line(ifs,cur_line);
        line.split_by_charset(line_vec);
      }
      else
      {
        throw SyntaxError("Missing coordinate number in line " + stringify(cur_line));
      }
    }
    else
    {
      if(line_vec[0] == "coords")
      {
        if(!(line_vec.size() == 2))
        {
          throw SyntaxError("Missing coordinate number in line " + stringify(cur_line));
        }
        Index coord_per_vertex_cmp;
        (line_vec.at(1)).parse(coord_per_vertex_cmp);
        if(coord_per_vertex_cmp != coord_per_vertex)
        {
          throw SyntaxError("Coordinate number missmatch in line " + stringify(cur_line));
        }
        line = read_next_line(ifs,cur_line);
        line.split_by_charset(line_vec);
      }
    }

    // if it is the adjacency file path
    if(line_vec[0] == "adj_file")
    {
      if(mesh_type == mt_structured)
      {
        throw SyntaxError("No adjacency file supported for structured mesh in line " + stringify(cur_line));
      }
      if(!(line_vec.size() == 2))
      {
        throw SyntaxError("Missing adjacency file path in line " + stringify(cur_line));
      }
      parse_adjacency_file(line_vec.at(1));

      line = read_next_line(ifs,cur_line);
      line.split_by_charset(line_vec);
      adjfile = true;
    }

    if(line_vec[0] != "</header>")
    {
      throw SyntaxError("Unknown file format. Expected </header> in line " + stringify(cur_line));
    }

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
    else if(!coordfile)
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

      // if it is a parent index chunk (of a submesh)
      else if((line == "<vert_idx>"||
               line == "<edge_idx>"||
               line == "<tria_idx>"||
               line == "<quad_idx>"||
               line == "<tetra_idx>"||
               line == "<hexa_idx>") && submesh)
      {
        cur_line = _parse_parents_chunk(cur_line, ifs, line);
      } // parent index sub chunk

      // if it is the end of the mesh section
      else if(line == break_line)
      {
        break;
      }
      else
      {
        throw SyntaxError("Unknown file format. Expected </mesh> or </submesh> in line " + stringify(cur_line));
      }

      line = read_next_line(ifs,cur_line);
    } // while

    if(line != break_line && !(!ifs.eof() && ifs.good()))
    {
      throw SyntaxError("Reached end of file but expected " + break_line + " at line " + stringify(cur_line));
    }

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
      line_vec[1].parse(vertex_count);
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
    if(shape_type_in == st_edge)
    {
      return "edge";
    }
    else if(shape_type_in == st_tria)
    {
      return "tria";
    }
    else if(shape_type_in == st_quad)
    {
      return "quad";
    }
    else if(shape_type_in == st_hexa)
    {
      return "hexa";
    }
    else if(shape_type_in == st_tetra)
    {
      return "tetra";
    }
    else if(shape_type_in == st_tria_quad)
    {
      return "tria-quad";
    }
    else if(shape_type_in == st_tetra_hexa)
    {
      return "tetra-hexa";
    }
    else
    {
      throw InternalError("Unknown shapetype");
    }
  } // String MeshStreamer::MeshDataContainer::convert_shape_type(const MeshStreamer::MeshDataContainer::ShapeType  shape_type) const

  // converts shapetype (String to ShapeType)
  MeshStreamer::MeshDataContainer::ShapeType MeshStreamer::MeshDataContainer::convert_shape_type(const String shape_type_in)
  {
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
    else if(shape_type_in == "tria-quad")
    {
      return st_tria_quad;
    }
    else if(shape_type_in == "tetra-hexa")
    {
      return st_tetra_hexa;
    }
    else
    {
      // todo error anpassen
      throw InternalError("Unknown shapetype");
    }
  } // ShapeType MeshStreamer::MeshDataContainer::convert_shape_type(const String shape_type) const


  // writes the stored cell set data into the output stream.
  void MeshStreamer::CellSetNode::write(std::ostream& ofs) const
  {
    ofs << "<cellset>" << std::endl;

    // header section
    ofs << " <header>" << std::endl;
    ofs << "  name " << cell_set.name << std::endl;
    ofs << "  parent " << cell_set.parent <<std::endl;
    ofs << " </header>" << std::endl;

    // info section
    if( !cell_set.info.empty() )
    {
      ofs << " <info>" << std::endl;
      ofs << "  " << cell_set.info << std::endl;
      ofs << " </info>" << std::endl;
    }
    // count section
    ofs << " <counts>" << std::endl;
    if(cell_set.slices.size() == 0)
    {
    if( cell_set.vertex_count != 0)
      ofs << "  verts " << cell_set.vertex_count << std::endl;
    if( cell_set.edge_count != 0)
      ofs << "  edges " << cell_set.edge_count << std::endl;
    if( cell_set.quad_count != 0)
      ofs << "  quads " << cell_set.quad_count << std::endl;
    if( cell_set.tria_count != 0)
      ofs << "  trias " << cell_set.tria_count << std::endl;
    if( cell_set.tetra_count != 0)
      ofs << "  tetras " << cell_set.tetra_count << std::endl;
    if( cell_set.hexa_count != 0)
      ofs << "  hexas " << cell_set.hexa_count << std::endl;
    }
    else
    {
      ofs << "  slices";
      for (Index i(0); i < cell_set.slices.size() ; ++i)
      {
        ofs  << "   " << cell_set.slices.at(i);
      }
      ofs << std::endl;
    }
    ofs << " </counts>" << std::endl;

    // parent indices
    ofs << " <vert_idx>" << std::endl;
    for (Index i(0); i < cell_set.vertex_count ; ++i)
    {
      ofs << "  " << (cell_set.parent_indices[0])[i] << std::endl;
    }
    ofs << " </vert_idx>" << std::endl;

    if ( !(cell_set.edge_count == 0) )
    {
      ofs << " <edge_idx>" << std::endl;
      for (Index i(0); i < cell_set.edge_count ; ++i)
      {
        ofs << "  " << (cell_set.parent_indices[1])[i] << std::endl;
      }
      ofs << " </edge_idx>" << std::endl;
      }

    if ( !(cell_set.tria_count == 0) )
    {
      ofs << " <tria_idx>" << std::endl;
      for (Index i(0); i < cell_set.tria_count ; ++i)
      {
        ofs << "  " << (cell_set.parent_indices[2])[i] << std::endl;
      }
      ofs << " </tria_idx>" << std::endl;
    }
    if ( !(cell_set.quad_count == 0) )
    {
      ofs << " <quad_idx>" << std::endl;
      for (Index i(0); i < cell_set.quad_count ; ++i)
      {
        ofs << "  " << (cell_set.parent_indices[2])[i] << std::endl;
      }
      ofs << " </quad_idx>" << std::endl;
    }

    if ( !(cell_set.tetra_count == 0) )
    {
      ofs << " <tetra_idx>" << std::endl;
      for (Index i(0); i < cell_set.tetra_count ; ++i)
      {
        ofs << "  " << (cell_set.parent_indices[3])[i] << std::endl;
      }
      ofs << "  </tetra_idx>" << std::endl;
    }

    if ( !(cell_set.hexa_count == 0) )
    {
      ofs << " <hexa_idx>" << std::endl;
      for (Index i(0); i < cell_set.hexa_count ; ++i)
      {
        ofs << "  " << (cell_set.parent_indices[3])[i] << std::endl;
      }
      ofs << " </hexa_idx>" << std::endl;
    }
    ofs << "</cellset>" << std::endl;
  }//write_cell_set_data

  // Writes the mesh data of this mesh and all submeshes related to this mesh
  // into the output stream.
  void MeshStreamer::MeshNode::write(std::ostream& ofs, bool submesh) const
  {
    // choose the right tag
    if (submesh)
      ofs << "<submesh>" << std::endl;
    else
      ofs << "<mesh>" << std::endl;

    // header section
    ofs << " <header>" << std::endl;
    if (submesh)
    {
      ofs << "  name " << mesh_data.name << std::endl;
      ofs << "  parent " << mesh_data.parent << std::endl;
      if ( !mesh_data.chart.empty() )
        ofs << "  chart " << mesh_data.chart << std::endl;
    }
    ofs << "  type " << mesh_data.convert_mesh_type(mesh_data.mesh_type) << std::endl;
    ofs << "  shape " << mesh_data.convert_shape_type(mesh_data.shape_type) << std::endl;
    if(mesh_data.coord_per_vertex!=0)
      ofs << "  coords " << mesh_data.coord_per_vertex << std::endl;
    ofs << " </header>" << std::endl;

    // info section
    if( !mesh_data.info.empty() )
    {
      ofs << " <info>" << std::endl;
      ofs << mesh_data.info << std::endl;
      ofs << " </info>" << std::endl;
    }

    // count section
    ofs << " <counts>" << std::endl;
    if(mesh_data.mesh_type == MeshDataContainer::mt_conformal)
    {
    if( mesh_data.vertex_count != 0)
      ofs << "  verts " << mesh_data.vertex_count << std::endl;
    if( mesh_data.edge_count != 0)
      ofs << "  edges " << mesh_data.edge_count << std::endl;
    if( mesh_data.quad_count != 0)
      ofs << "  quads " << mesh_data.quad_count << std::endl;
    if( mesh_data.tria_count != 0)
      ofs << "  trias " << mesh_data.tria_count << std::endl;
    if( mesh_data.tetra_count != 0)
      ofs << "  tetras " << mesh_data.tetra_count << std::endl;
    if( mesh_data.hexa_count != 0)
      ofs << "  hexas " << mesh_data.hexa_count << std::endl;
    }
    else if( mesh_data.mesh_type == MeshDataContainer::mt_structured)
    {
      ofs << "  slices";
      for (Index i(0); i < mesh_data.slices.size() ; ++i)
      {
        ofs << " " << mesh_data.slices.at(i);
      }
      ofs << std::endl;
    }
    ofs << " </counts>" << std::endl;

    // coord section
    ofs << " <coords>" << std::endl;
    for (Index i(0); i < mesh_data.vertex_count ; ++i)
    {
      for (Index j(0); j < mesh_data.coord_per_vertex ; ++j)
      {
        ofs  << "  " << scientify((mesh_data.coords[i])[j]);
      }
      ofs << std::endl;
    }
    ofs << " </coords>" << std::endl;

    // adjacency section
    if(mesh_data.mesh_type == MeshDataContainer::mt_conformal)
    {
      ofs << " <vert@edge>" << std::endl;
      for (Index i(0); i < (mesh_data.adjacencies[0][1]).size() ; ++i)
      {
        ofs << "  " << ((mesh_data.adjacencies[0][1])[i])[0];
        ofs << " " << ((mesh_data.adjacencies[0][1])[i])[1] << std::endl;
      }
      ofs << " </vert@edge>" << std::endl;

      if (!( (mesh_data.adjacencies[0][2]).size() == 0 ))
      {
        if ( ((mesh_data.adjacencies[0][2])[0]).size() == 3 )
        {
          ofs << " <vert@tria>" << std::endl;
          for (Index i(0); i < (mesh_data.adjacencies[0][2]).size() ; ++i)
          {
            ofs << "  " << ((mesh_data.adjacencies[0][2])[i])[0];
            ofs << " " << ((mesh_data.adjacencies[0][2])[i])[1];
            ofs << " " << ((mesh_data.adjacencies[0][2])[i])[2] << std::endl;
          }
          ofs << " </vert@tria>" << std::endl;
        }
        else if ( ((mesh_data.adjacencies[0][2])[0]).size() == 4 )
        {
          ofs << " <vert@quad>" << std::endl;
          for (Index i(0); i < (mesh_data.adjacencies[0][2]).size() ; ++i)
          {
            ofs << "  " << ((mesh_data.adjacencies[0][2])[i])[0];
            ofs << " " << ((mesh_data.adjacencies[0][2])[i])[1];
            ofs << " " << ((mesh_data.adjacencies[0][2])[i])[2];
            ofs << " " << ((mesh_data.adjacencies[0][2])[i])[3] << std::endl;
          }
          ofs << " </vert@quad>" << std::endl;
        }
      }

      if (!( (mesh_data.adjacencies[0][3]).size() == 0 ))
      {
        if ( ((mesh_data.adjacencies[0][3])[0]).size() == 4 )
        {
          ofs << " <vert@tetra>" << std::endl;
          for (Index i(0); i < (mesh_data.adjacencies[0][3]).size() ; ++i)
          {
            ofs << "  " << ((mesh_data.adjacencies[0][3])[i])[0];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[1];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[2];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[3] << std::endl;
          }
          ofs << " </vert@tetra>" << std::endl;
        }
        else if ( ((mesh_data.adjacencies[0][3])[0]).size() == 8 )
        {
          ofs << " <vert@hexa>" << std::endl;
          for (Index i(0); i < (mesh_data.adjacencies[0][3]).size() ; ++i)
          {
            ofs << "  " << ((mesh_data.adjacencies[0][3])[i])[0];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[1];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[2];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[3];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[4];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[5];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[6];
            ofs << " " << ((mesh_data.adjacencies[0][3])[i])[7] << std::endl;
          }
          ofs << " </vert@hexa>" << std::endl;
        }
      }
    }

    // if it is a submesh, add the parent indices
    if (submesh)
    {
      ofs << " <vert_idx>" << std::endl;
      for (Index i(0); i < mesh_data.vertex_count ; ++i)
      {
        ofs << "  " << (mesh_data.parent_indices[0])[i] << std::endl;
      }
      ofs << " </vert_idx>" << std::endl;

      if ( !(mesh_data.edge_count == 0) )
      {
        ofs << " <edge_idx>" << std::endl;
        for (Index i(0); i < mesh_data.edge_count ; ++i)
        {
          ofs << "  " << (mesh_data.parent_indices[1])[i] << std::endl;
        }
        ofs << " </edge_idx>" << std::endl;
      }

      if ( !(mesh_data.tria_count == 0) )
      {
        ofs << " <tria_idx>" << std::endl;
        for (Index i(0); i < mesh_data.tria_count ; ++i)
        {
          ofs << "  " << (mesh_data.parent_indices[2])[i] << std::endl;
        }
        ofs << " </tria_idx>" << std::endl;
      }

      if ( !(mesh_data.quad_count == 0) )
      {
        ofs << " <quad_idx>" << std::endl;
        for (Index i(0); i < mesh_data.quad_count ; ++i)
        {
          ofs << "  " << (mesh_data.parent_indices[2])[i] << std::endl;
        }
        ofs << " </quad_idx>" << std::endl;
      }

      if ( !(mesh_data.tetra_count == 0) )
      {
        ofs << " <tetra_idx>" << std::endl;
        for (Index i(0); i < mesh_data.tetra_count ; ++i)
        {
          ofs << "  " << (mesh_data.parent_indices[3])[i] << std::endl;
        }
        ofs << " </tetra_idx>" << std::endl;
      }

      if ( !(mesh_data.hexa_count == 0) )
      {
        ofs << " <hexa_idx>" << std::endl;
        for (Index i(0); i < mesh_data.hexa_count ; ++i)
        {
          ofs << "  " << (mesh_data.parent_indices[3])[i] << std::endl;
        }
        ofs << " </hexa_idx>" << std::endl;
      }
    }

    // choose the right tag to end the mesh section
    if (submesh)
      ofs<<"</submesh>"<<std::endl;
    else
      ofs<<"</mesh>"<<std::endl;

    // loop through all submeshes related to this mesh and drop the data as well
     SubMeshMap::const_iterator it(sub_mesh_map.begin()), jt(sub_mesh_map.end());
     for(; it != jt; ++it)
     {
       it->second->write(ofs, true);
     }
      // loop through all cell sets related to this mesh and drop the data as well
     CellSetMap::const_iterator it_cell(cell_set_map.begin()), jt_cell(cell_set_map.end());
     for(; it_cell != jt_cell; ++it_cell)
     {
       it_cell->second->write(ofs);
     }
   }// write_mesh_data

} //namespace FEAST
