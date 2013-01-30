// includes, FEAST
#include <kernel/util/mesh_reader.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <fstream>

namespace FEAST
{

  // default constructor
  MeshReader::MeshReader() :
    _version(0),
    _num_submeshes(0),
    _num_cellsets(0),
    _chart_path(""),
    _root_mesh_node(nullptr)
  {
    CONTEXT("MeshReader::MeshReader()");
  }

  // default destructor
  MeshReader::~MeshReader()
  {
    CONTEXT("MeshReader::~MeshReader()");
    if(_root_mesh_node != nullptr)
    {
      delete _root_mesh_node;
    }
  }

  // returns the parent mesh specified by "parent_name"
  MeshReader::CellSetParent* MeshReader::_find_cell_set_parent(String parent_name)
  {
    CONTEXT("MeshReader::_find_cell_set_parent");
    if(parent_name == "root")
      return _root_mesh_node;
    if(_root_mesh_node != nullptr)
      return _root_mesh_node->find_cell_set_parent(parent_name);
    return nullptr;
  }

  // returns the sub mesh parent specified by "parent_name"
  MeshReader::MeshNode* MeshReader::_find_sub_mesh_parent(String parent_name)
  {
    CONTEXT("MeshReader::_find_sub_mesh_parent");
    if(parent_name == "root")
      return _root_mesh_node;
    if(_root_mesh_node != nullptr)
      return _root_mesh_node->find_sub_mesh(parent_name);
    return nullptr;
  }

  // parses the FEAST- mesh file given by filename
  void MeshReader::parse_mesh_file(String filename)
  {
    CONTEXT("MeshReader::parse_mesh_file(String)");

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
        throw exc;
      }
    }
  } // MeshReader::parse_mesh_file(String filename)


  // parses the mesh data stream given by ifs
  void MeshReader::parse_mesh_file(std::istream& ifs)
  {
    CONTEXT("MeshReader::parse_mesh_file(std::ifstream&)");

    // a string containing the current line
    String line = "42";
    line.reserve(256);

    // auxiliary variable that counts the lines
    Index cur_line = 0;

    // ignore everything until the mesh file begins
    while(line != "<feast_mesh_file>" && !ifs.eof() && ifs.good())
    {
      getline(ifs, line);
      ++cur_line;
      line.trim_me();
    }

    // if the beginning of the file was not found or something went wrong
    if(ifs.eof() || !ifs.good())
    {
      throw SyntaxError("Beginning of mesh file not found");
    }

    // loop over all lines until we reach the end of the mesh file
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

      // if it is the header chunk
      if(line == "<header>")
      {
        cur_line = _parse_header_section(cur_line, ifs);
      }

      // if it is an info chunk
      else if(line == "<info>")
      {
        cur_line = _parse_info_section(cur_line, ifs);
      }

      // if it is a mesh chunk
      else if(line == "<mesh>")
      {
        cur_line = _parse_mesh_section(cur_line, false, ifs);

        // insert basic root mesh information
        _root_mesh_node->mesh_data.name = "root";
        _root_mesh_node->mesh_data.parent = "none";
        _root_mesh_node->mesh_data.chart = "none";
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
        break;
      }
    } // end while

    // if the file did not end properly
    if(line != "</feast_mesh_file>")
    {
      throw SyntaxError("Unexpected ending of file in line " + stringify(cur_line));
    }
  } // MeshReader::parse_mesh_file(std::istream& ifs)


  // parses header-section streams
  Index MeshReader:: _parse_header_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshReader::parse_header_section");

    // string the current line is saved in
    String line;

    // while everything is fine
    while(!ifs.eof() && ifs.good())
    {
      // get a line
      getline(ifs, line);
      ++cur_line;
      line.trim_me();

      // if it is the version
      if(line.compare(0, 8, "version ") == 0)
      {
        // erase "version "
        line.erase(0, 8);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing version number in line " + stringify(cur_line));
        }
        line.parse(_version);
      }

      // if it is the chart-file-path
      else if(line.compare(0, 11, "chart_file ") == 0)
      {
        // erase "chart_file "
        line.erase(0, 11);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing chart file path in line " + stringify(cur_line));
        }
        _chart_path = line;
      }

      // if it is the number of submeshes
      else if(line.compare(0, 10, "submeshes ") == 0)
      {
        // erase "submeshes "
        line.erase(0, 10);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing number of submeshes in line " + stringify(cur_line));
        }
        line.parse(_num_submeshes);
      }

      // if it is the number of cellsets
      else if(line.compare(0, 9, "cellsets ") == 0)
      {
        // erase "cellsets "
        line.erase(0, 9);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing number of cellsets in line " + stringify(cur_line));
        }
        line.parse(_num_cellsets);
      }

      // if it is the end of the header section
      else if(line == "</header>")
      {
        break;
      }

      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    }

    // return the number of lines
    return cur_line;

  } // Index MeshReader::parse_header_section(Index cur_line, std::istream& ifs)


  // parses/ignores the given info-section stream
  Index MeshReader::_parse_info_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshReader::parse_info_section");

    // string the current line is saved in
    String line;

    // skip the info block
    while(line != "</info>" && !ifs.eof() && ifs.good())
    {
      // get a line
      getline(ifs, line);
      ++cur_line;
      line.trim_me();
    }
    return cur_line;
  } // Index MeshReader::parse_info_section(Index cur_line, std::istream& ifs)


  // parses the given mesh-section stream
  Index MeshReader::_parse_mesh_section(Index cur_line, bool submesh, std::istream& ifs)
  {
    CONTEXT("MeshReader::parse_mesh_section");

    // string the current line is saved in
    String line;

    // if it is the root mesh
    String break_line = submesh ? "</submesh>" : "</mesh>";

    // create a new mesh node
    MeshNode* mesh_node = new MeshNode();
    if(!submesh)
    {
      ASSERT_(_root_mesh_node == nullptr);
      _root_mesh_node = mesh_node;
    }

    // mesh data container of the current mesh
    MeshDataContainer& current_mesh(mesh_node->mesh_data);

    // while everything is fine
    while(!ifs.eof() && ifs.good())
    {
      // get a line
      getline(ifs, line);
      ++cur_line;
      line.trim_me();

      // if it is the header sub chunk
      if(line.compare(0, 8, "<header>") == 0)
      {
        while(!ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line;
          line.trim_me();

          // if it is the type
          if(line.compare(0, 5, "type ") == 0)
          {
            // erase "type "
            line.erase(0, 5);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing type in line " + stringify(cur_line));
            }
            current_mesh.mesh_type = line;
          }

          // if it is the shape
          else if(line.compare(0, 6, "shape ") == 0)
          {
            // erase "shape "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing shape type in line " + stringify(cur_line));
            }
            current_mesh.shape_type = line;
          }

          // if it is the coord-file path
          else if(line.compare(0, 11, "coord_file ") == 0)
          {
            // erase "coord_file "
            line.erase(0, 11);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing coordinate file path in line " + stringify(cur_line));
            }
            current_mesh.coord_path = line;
            // parse the coord file
            parse_coord_file(line, &current_mesh);
          }

          // if it is the number of coordinates per vertex
          else if(line.compare(0, 7, "coords ") == 0)
          {
            // erase "coords "
            line.erase(0, 7);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing coordinate number in line " + stringify(cur_line));
            }
            line.parse(current_mesh.coord_per_vertex);
          }

          // if it is the adjacency file path
          else if(line.compare(0, 9, "adj_file ") == 0)
          {
            // erase "adj_file "
            line.erase(0, 9);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing adjacency file path in line " + stringify(cur_line));
            }
            current_mesh.adj_path = line;
            parse_adjacency_file(line, &current_mesh);
          }

          // if it is the name and the mesh is a submesh
          else if(line.compare(0, 5, "name ") == 0 && submesh)
          {
            // erase "name "
            line.erase(0, 5);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing name in line " + stringify(cur_line));
            }
            current_mesh.name = line;
          }

          // if it is the parent (of a submesh)
          else if(line.compare(0, 7, "parent ") == 0 && submesh)
          {
            // erase "parent "
            line.erase(0, 7);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing parent in line " + stringify(cur_line));
            }
            current_mesh.parent = line;
          }

          // if it is the chart
          else if(line.compare(0, 6, "chart ") == 0 && submesh)
          {
            // erase "chart "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing chart name in line " + stringify(cur_line));
            }
            current_mesh.chart = line;
          }

          // if it is the end of the header file
          else if(line == "</header>")
          {
            break;
          }

          else
          {
            throw SyntaxError("Unknown file format in line " + stringify(cur_line));
          }
        }
      } // header sub-chunk

      // if it is an info sub chunk
      else if(line.compare(0, 6, "<info>") == 0)
      {
        cur_line = _parse_info_section(cur_line, ifs);
      } // info sub-chunk

      // if it is the counts sub chunk
      else if(line.compare(0, 8, "<counts>") == 0)
      {
        cur_line = _parse_counts_chunk(cur_line, ifs, &current_mesh);
      } // counts sub-chunk

      // if it is a coords sub chunk
      else if(line.compare(0, 8, "<coords>") == 0)
      {
        cur_line = _parse_coords_chunk(cur_line, ifs, &current_mesh);
      } // coords sub-chunk

      // if it is an adjacency sub chunk
      else if(line.find('@') != std::string::npos)
      {
        cur_line = _parse_adjacency_chunk(cur_line, ifs, &current_mesh, line);
      } // adjacencies sub-chunk

      // if it is a parent index chunk (of a submesh)
      else if((line.compare(0, 10, "<vert_idx>") == 0 ||
              line.compare(0, 10, "<edge_idx>") == 0 ||
              line.compare(0, 10, "<tria_idx>") == 0 ||
              line.compare(0, 10, "<quad_idx>") == 0 ||
              line.compare(0, 11, "<tetra_idx>") == 0 ||
              line.compare(0, 10, "<hexa_idx>") == 0) && submesh)
      {
        cur_line = _parse_parents_chunk(cur_line, ifs, &current_mesh, line);
      } // parent index sub chunk

      // if it is the end of the mesh section
      else if(line == break_line)
      {
        if(submesh)
        {
          MeshNode* parent = _find_sub_mesh_parent(current_mesh.parent);
          ASSERT_(parent != nullptr);
          parent->sub_mesh_map.insert(std::make_pair(current_mesh.name, mesh_node));
        }
        break;
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    } // while

    // return number of lines read so far
    return cur_line;

  } // MeshReader::parse_mesh_section(Index cur_line, std::istream& ifs)


  // parses the cellset-section stream ifs
  Index MeshReader::_parse_cellset_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshReader::parse_cellset_section");

    // (auxiliary) variables
    String line;

    // create a new cell-set node
    CellSetNode* cell_set_node = new CellSetNode();

    // cellset data container of the root mesh
    CellSetContainer& current_set(cell_set_node->cell_set);

    while(!ifs.eof() && ifs.good())
    {
      // get a line
      getline(ifs, line);
      ++cur_line;
      line.trim_me();

      // if it is the header sub chunk
      if(line.compare(0, 8, "<header>") == 0)
      {
        while(!ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line;
          line.trim_me();

          // if it is the name
          if(line.compare(0, 5, "name ") == 0)
          {
            // erase "name "
            line.erase(0, 5);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing name in line " + stringify(cur_line));
            }
            current_set.name = line;
          }

          // if it is the parent
          else if(line.compare(0, 7, "parent ") == 0)
          {
            // erase "parent "
            line.erase(0, 7);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing parent in line " + stringify(cur_line));
            }
            current_set.parent = line;
          }
          // if it is the end of the header section
          else if(line == "</header>")
          {
            break;
          }

          else
          {
            throw SyntaxError("Unknown file format in line " + stringify(cur_line));
          }
        }
      } // header sub-chunk

      // if it is an info sub chunk
      else if(line.compare(0, 6, "<info>") == 0)
      {
        cur_line = _parse_info_section(cur_line, ifs);
      } // info sub-chunk

      // if it is the counts sub chunk
      else if(line.compare(0, 8, "<counts>") == 0)
      {
        cur_line = _parse_counts_chunk(cur_line, ifs, &current_set);
      } // counts sub-chunk

      // if it is a parent index chunk
      else if(line.compare(0, 10, "<vert_idx>") == 0 ||
              line.compare(0, 10, "<edge_idx>") == 0 ||
              line.compare(0, 10, "<tria_idx>") == 0 ||
              line.compare(0, 10, "<quad_idx>") == 0 ||
              line.compare(0, 11, "<tetra_idx>") == 0 ||
              line.compare(0, 10, "<hexa_idx>") == 0)
      {
        cur_line = _parse_parents_chunk(cur_line, ifs, &current_set, line);
      } // parent index sub chunk

      // if it is the end of the cellset section
      else if(line == "</cellset>")
      {
        CellSetParent* parent = _find_cell_set_parent(current_set.parent);
        ASSERT_(parent != nullptr);
        parent->cell_set_map.insert(std::make_pair(current_set.name, cell_set_node));
        break;
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    } // while

    // return number of lines read so far
    return cur_line;
  } // Index MeshReader::parse_cellset_section(Index cur_line, std::istream& ifs)


  // parses the coord file given by filepath
  void MeshReader::parse_coord_file(String filename, MeshReader::MeshDataContainer *mesh)
  {
    CONTEXT("MeshReader::parse_coord_file(String filename)");

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
      parse_coord_file(ifs, mesh);
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
        throw exc;
      }
    }
  } // MeshReader::parse_coord_file(String filename, MeshReader::MeshDataContainer *mesh)


  // parses the coord-file stream given by ifs
  void MeshReader::parse_coord_file(std::istream& ifs, MeshReader::MeshDataContainer *mesh)
  {
    CONTEXT("MeshReader::parse_coord_file(std::ifstream&)");

    // a string containing the current line
    String line = "42";
    line.reserve(256);

    // auxiliary variables
    String coord;
    Index cur_line = 0;

    // ignore everything until the coord file begins
    while(line != "<feast_coord_file>" && !ifs.eof() && ifs.good())
    {
      getline(ifs, line);
      ++cur_line;
      line.trim_me();
    }

    // if the beginning of the file was not found or something went wrong
    if(ifs.eof() || !ifs.good())
    {
      throw SyntaxError("Beginning of coord file not found");
    }

    // loop over all lines until we reach the end of the coord file
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

      // if it is the header chunk
      if(line == "<header>")
      {
        while(!ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line;
          line.trim_me();

          if(line.compare(0, 8, "version ") == 0)
          {
            // erase "version "
            line.erase(0, 8);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing version number in line " + stringify(cur_line));
            }
            mesh->coord_version = line;
          }
          else if(line.compare(0, 6, "verts ") == 0)
          {
            // erase "verts "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing number of vertices in line " + stringify(cur_line));
            }
            line.parse(mesh->vertex_number);
          }
          else if(line.compare(0, 7, "coords ") == 0)
          {
            // erase "coords "
            line.erase(0, 7);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing coordinate number in line " + stringify(cur_line));
            }
            line.parse(mesh->coord_per_vertex);
          }
          else if(line == "</header>")
          {
            break;
          }
          else
          {
            throw SyntaxError("Unknown file format in line " + stringify(cur_line));
          }
        } // while
      } // header sub chunk

      // if it is an info chunk
      else if(line == "<info>")
      {
        cur_line = _parse_info_section(cur_line, ifs);
      } // info sub chunk

      // if it is a coord chunk
      else if(line.compare(0, 8, "<coords>") == 0)
      {
        cur_line = _parse_coords_chunk(cur_line, ifs, mesh);
      } // coords sub-chunk

      // if it is the end
      else if(line == "</feast_coord_file>")
      {
        break;
      }
      else
      {
        throw SyntaxError("Unknown format in line " + stringify(cur_line));
      }
    } // end while

    // if the file did not end properly
    if(line != "</feast_coord_file>")
    {
      throw SyntaxError("Unexpected file ending in line " + stringify(cur_line));
    }
  } // MeshReader::parse_coord_file(std::istream& ifs, MeshReader::MeshDataContainer *mesh)


  // parses the adjacency file given by filename
  void MeshReader::parse_adjacency_file(String filename, MeshReader::MeshDataContainer *mesh)
  {
    CONTEXT("MeshReader::parse_adjacency_file(String filename)");

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
      parse_adjacency_file(ifs, mesh);
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
        throw exc;
      }
    }
  } // MeshReader::parse_adjacency_file(String filename)


  // parses the adjacency-file-stream ifs
  void MeshReader::parse_adjacency_file(std::istream& ifs, MeshReader::MeshDataContainer *mesh)
  {
    CONTEXT("MeshReader::parse_adjacency_file(std::ifstream&)");

    // a string containing the current line
    String line = "42";
    line.reserve(256);

    // (auxiliary) variables
    String adj, face, shape, break_line;
    Index cur_line = 0;

    // ignore everything until the adjacency file begins
    while(line != "<feast_adjacency_file>" && !ifs.eof() && ifs.good())
    {
      getline(ifs, line);
      ++cur_line;
      line.trim_me();
    }

    // if the beginning of the file was not found or something went wrong
    if(ifs.eof() || !ifs.good())
    {
      throw SyntaxError("Beginning of adjacency file not found");
    }

    // loop over all lines until we reach the end of the adjacency file
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

      // if it is the header chunk
      if(line == "<header>")
      {
        while(!ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line;
          line.trim_me();

          if(line.compare(0, 8, "version ") == 0)
          {
            // erase "version "
            line.erase(0, 8);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing version number in line " + stringify(cur_line));
            }
            mesh->adjacency_version = line;
          }
          else if(line.compare(0, 5, "type ") == 0)
          {
            // erase "type "
            line.erase(0, 5);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing type in line " + stringify(cur_line));
            }
            mesh->mesh_type = line;
          }
          else if(line.compare(0, 6, "shape ") == 0)
          {
            // erase "shape "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing shape type in line " + stringify(cur_line));
            }
            mesh->shape_type = line;
          }
          // if it is the end of the header section
          else if(line == "</header>")
          {
            break;
          }
          else
          {
            throw SyntaxError("Unknown file format in line " + stringify(cur_line));
          }
        } // while
      } // header sub chunk

      // if it is an info chunk
      else if(line == "<info>")
      {
        cur_line = _parse_info_section(cur_line, ifs);
      }

      // if it is the counts sub chunk
      else if(line.compare(0, 8, "<counts>") == 0)
      {
        cur_line = _parse_counts_chunk(cur_line, ifs, mesh);
      } // counts sub-chunk

      // if it is an adjacency sub chunk
      else if(line.find('@') != std::string::npos)
      {
        cur_line =  _parse_adjacency_chunk(cur_line, ifs, mesh, line);
      } // adjacencies sub-chunk

      // if it is the end of the file
      else if(line == "</feast_adjacency_file>")
      {
        break;
      }
      else
      {
        throw SyntaxError("Unknown format in line " + stringify(cur_line));
      }
    } // end while

    // if the file did not end properly
    if(line != "</feast_adjacency_file>")
    {
      throw SyntaxError("Unexpected file ending in line " + stringify(cur_line));
    }
  } // MeshReader::parse_adjacency_file(std::istream& ifs, MeshReader::MeshDataContainer *mesh)


  // parses a counts subchunk
  Index MeshReader::_parse_counts_chunk(Index cur_line, std::istream& ifs, BaseContainer* container)
  {
    CONTEXT("MeshReader::_parse_counts_section");

    // string the current line is saved in
    String line;

    // default setting: zero
    container->vertex_number = 0;
    container->edge_number = 0;
    container->quad_number = 0;
    container->tria_number = 0;
    container->tetra_number = 0;
    container->hexa_number = 0;

    while(!ifs.eof() && ifs.good())
    {
      // get a line
      getline(ifs, line);
      ++cur_line;
      line.trim_me();

      // vertex number
      if(line.compare(0, 6, "verts ") == 0)
      {
        // erase "verts "
        line.erase(0, 6);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing vertex number in line " + stringify(cur_line));
        }
        line.parse(container->vertex_number);
      }

      // edge number
      else if(line.compare(0, 6, "edges ") == 0)
      {
        // erase "edges "
        line.erase(0, 6);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing edge number in line " + stringify(cur_line));
        }
        line.parse(container->edge_number);
      }

      // trias number
      else if(line.compare(0, 6, "trias ") == 0)
      {
        // erase "trias "
        line.erase(0, 6);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing triangle number in line " + stringify(cur_line));
        }
        line.parse(container->tria_number);
      }

      // quad number
      else if(line.compare(0, 6, "quads ") == 0)
      {
        // erase "quads "
        line.erase(0, 6);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing quad number in line " + stringify(cur_line));
        }
        line.parse(container->quad_number);
      }

      // tetra number
      else if(line.compare(0, 7, "tetras ") == 0)
      {
        // erase "tetras "
        line.erase(0, 7);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing tetra number in line " + stringify(cur_line));
        }
        line.parse(container->tetra_number);
      }

      // hexa number
      else if(line.compare(0, 6, "hexas ") == 0)
      {
        // erase "hexas "
        line.erase(0, 6);
        line.trim_me();
        if(line.empty())
        {
          throw SyntaxError("Missing hexa number in line " + stringify(cur_line));
        }
        line.parse(container->hexa_number);
      }
      // if it is the end of the counts sub chunk
      else if(line == "</counts>")
      {
        break;
      }
      // if it is an empty line
      else if(line.empty())
      {
        continue;
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    }
    // return number of read lines
    return cur_line;
  } // MeshReader::_parse_counts_chunk


  // parses a coords subchunk
  Index MeshReader::_parse_coords_chunk(Index cur_line, std::istream& ifs, MeshDataContainer* container)
  {
    CONTEXT("MeshReader::_parse_coords_section");

    // string the current line is saved in
    String line;

    while(!ifs.eof() && ifs.good())
    {
      // get a line
      getline(ifs, line);
      ++cur_line;
      line.trim_me();

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
        line.fork_by_string(line_vec, " ");

        // parse the substrings
        for(Index i(0); i < line_vec.size(); ++i)
        {
          line_vec[i].trim_me();
          if( !(line_vec[i]).empty())
          {
            line_vec[i].trim_me();
            if ( !(line_vec[i]).parse(current_value) )
            {
              throw SyntaxError("Wrong coordinate format in line " + stringify(cur_line));
            }
            v.push_back(current_value);
          }
        }

        // if the number of entries does not match the coord_per_vertex variable
        if (v.size() != container->coord_per_vertex)
        {
          throw SyntaxError("Wrong coordinate format in line " + stringify(cur_line));
        }

        // add to the coordinate stack
        (container->coords).push_back(v);
      }
    }
    // current number of lines
    return cur_line;
  } //MeshReader::_parse_coords_chunk


  // parses an adjacency subchunk
  Index MeshReader::_parse_adjacency_chunk(Index cur_line, std::istream& ifs, MeshDataContainer* container, String line)
  {
    CONTEXT("MeshReader::_parse_adjacency_section");

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
        throw SyntaxError("Unknown format in line " + stringify(cur_line));
      }

      // adjacency vector
      std::vector<std::vector<Index> > a_stack;

      while(!ifs.eof() && ifs.good())
      {
        // get a line
        getline(ifs, line);
        ++cur_line;
        line.trim_me();

        // if it is the end of the sub chunk
        if(line.find('@') != std::string::npos)
        {
          container->adjacencies[0][shape_dim] = a_stack;
          break;
        }
        else
        {
          // auxiliary variables
          std::vector<String> line_vec;
          std::vector<Index> a;
          Index current_value;

          // separate by " "
          line.fork_by_string(line_vec, " ");

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
        getline(ifs, line);
        ++cur_line;

        // if it is the end of the sub chunk
        if(line.find('@') != std::string::npos)
        {
          break;
        }
      }
    }
    // return number of lines read so far
    return cur_line;
  } // MeshReader::_parse_adjacency_chunk


  // parses a parent-indices chunk
  Index MeshReader::_parse_parents_chunk(Index cur_line, std::istream& ifs, BaseContainer* container, String line)
  {
    CONTEXT("MeshReader::_parse_parents_chunk");

    std::vector<Index> idx;
    Index dim = 0;

    // get the dimension of the entity
    if(line.compare(0, 10, "<vert_idx>") == 0)
    {
     dim = 0;
    }
    else if (line.compare(0, 10, "<edge_idx>") == 0)
    {
      dim = 1;
    }
    else if (line.compare(0, 10, "<tria_idx>") == 0 || line.compare(0, 10, "<quad_idx>") == 0)
    {
      dim = 2;
    }
    else if (line.compare(0, 11, "<tetra_idx>") == 0 || line.compare(0, 10, "<hexa_idx>") == 0)
    {
      dim = 3;
    }

    while(!ifs.eof() && ifs.good())
    {

      // get a line
      getline(ifs, line);
      ++cur_line;
      line.trim_me();

      Index current_value;

      // if it is the end of the parent index chunk
      if(line.find("idx") != String::npos)
      {
        container->parent_indices[dim] = idx;
        break;
      }
      else
      {
        line.parse(current_value);
        idx.push_back(current_value);
      }
    }

    // return number of lines read so far
    return cur_line;
  } // MeshReader::_parse_parents_chunk


  // returns the version
  Index MeshReader::get_version() const
  {
    CONTEXT("MeshReader::get_version()");
    return _version;
  }

  // returns the chart path
  String MeshReader::get_chart_path() const
  {
    CONTEXT("MeshReader::get_chart_path()");
    return _chart_path;
  }

  // returns the number of submeshes
  Index MeshReader::get_num_submeshes() const
  {
    CONTEXT("MeshReader::get_number_of_submeshes()");
    return _num_submeshes;
  }

  // returns the number of cellsets
  Index MeshReader::get_num_cellsets() const
  {
    CONTEXT("MeshReader::get_number_of_cellsets()");
    return _num_cellsets;
  }

  // returns a pointer to the MeshDataContainer specified by "name"
  MeshReader::MeshDataContainer* MeshReader::get_mesh(String name)
  {
    CONTEXT("MeshReader::get_mesh()");
    if(_root_mesh_node == nullptr)
      return nullptr;
    if(name.compare_no_case("root") == 0)
      return &_root_mesh_node->mesh_data;
    MeshNode* node = _root_mesh_node->find_sub_mesh(name);
    return (node != nullptr) ? &node->mesh_data : nullptr;
  }

  // returns a pointer to the CellSetContainer specified by "name"
  MeshReader::CellSetContainer* MeshReader::get_cell_set(String name)
  {
    CONTEXT("MeshReader::get_mesh()");
    if(_root_mesh_node == nullptr)
      return nullptr;
    CellSetNode* node = _root_mesh_node->find_cell_set(name);
    return (node != nullptr) ? &node->cell_set : nullptr;
  }
} //namespace FEAST
