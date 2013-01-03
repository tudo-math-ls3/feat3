// includes, FEAST
#include <kernel/util/mesh_reader.hpp>

// includes, system
#include <fstream>

namespace FEAST
{

  // default constructor
  MeshReader::MeshReader()
  {
    CONTEXT("MeshReader::MeshReader()");
  }


  // default destructor
  MeshReader::~MeshReader()
  {
    CONTEXT("MeshReader::~MeshReader()");
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
        cur_line = parse_header_section(cur_line, ifs);
      }

      // if it is an info chunk
      else if(line == "<info>")
      {
        cur_line = parse_info_section(cur_line, ifs);
      }

      // if it is a mesh chunk
      else if(line == "<mesh>")
      {
        cur_line = parse_mesh_section(cur_line, "root", ifs);

        // insert basic root mesh information
        (meshes.back()).name = "root";
        (meshes.back()).parent = "none";
        (meshes.back()).chart = "none";
      }

      // if it is a submesh chunk
      else if(line == "<submesh>")
      {
        cur_line = parse_mesh_section(cur_line, "sub", ifs);
      }

      // if it is a cellset chunk
      else if(line == "<cellset>")
      {
        cur_line = parse_cellset_section(cur_line, ifs);
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
  Index MeshReader::parse_header_section(Index cur_line, std::istream& ifs)
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
        version = line;
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
        chart_path = line;
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
        number_of_submeshes = atoi(line.c_str());
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
        number_of_cellsets =  atoi(line.c_str());
      }
      // if it is the end of the header section
      else if(line == "</header>")
      {
        return cur_line;
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    }

    // pointless default case to avoid warnings
    return 42;

  } // Index MeshReader::parse_header_section(Index cur_line, std::istream& ifs)


  // parses/ignores the given info-section stream
  Index MeshReader::parse_info_section(Index cur_line, std::istream& ifs)
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
  Index MeshReader::parse_mesh_section(Index cur_line, String flag, std::istream& ifs)
  {
    CONTEXT("MeshReader::parse_mesh_section");

    // string the current line is saved in
    String line;

    // coordinate counter
    Index coord_counter;

    // various auxiliary variables
    Index face_dim = 0;
    Index shape_dim = 0;
    String::size_type found;
    String coord, adj, face, shape, break_line;

    // if it is the root mesh
    if(flag == "root")
    {
      break_line = "</mesh>";
    }
    // if it is a submesh
    else if(flag == "sub")
    {
      break_line = "</submesh>";
    }

    // mesh data container of the current mesh
    MeshDataContainer current_mesh;

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
            current_mesh.coord_per_vertex = atoi(line.c_str());
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
          else if(line.compare(0, 5, "name ") == 0 && flag == "sub")
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
          else if(line.compare(0, 7, "parent ") == 0 && flag == "sub")
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
          else if(line.compare(0, 6, "chart ") == 0 && flag == "sub")
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
        // skip the info block
        while(line != "</info>" && !ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          line.trim_me();
          ++cur_line;
        }
      } // info sub-chunk

      // if it is the counts sub chunk
      else if(line.compare(0, 8, "<counts>") == 0)
      {

        // default setting: zero
        current_mesh.vertex_number = 0;
        current_mesh.edge_number = 0;
        current_mesh.quad_number = 0;
        current_mesh.tria_number = 0;
        current_mesh.tetra_number = 0;
        current_mesh.hexa_number = 0;

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
            current_mesh.vertex_number = atoi(line.c_str());
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
            current_mesh.edge_number = atoi(line.c_str());
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
            current_mesh.tria_number = atoi(line.c_str());
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
            current_mesh.quad_number = atoi(line.c_str());
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
            current_mesh.tetra_number = atoi(line.c_str());
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
            current_mesh.hexa_number = atoi(line.c_str());
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
      } // counts sub-chunk

      // if it is a coords sub chunk
      else if(line.compare(0, 8, "<coords>") == 0)
      {
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
            // coordinate counter and current coord-vector
            coord_counter = 0;
            std::vector<double> v;

            // while everything is fine
            while(!line.empty())
            {
              // separating coordinate
              found = line.find(' ');
              coord = line.substr(0, found);
              if(found != std::string::npos)
              {
                line = line.substr(found + 1);
              }
              else
              {
                line = "";
              }
              line.trim_me();
              coord.trim_me();
              // add to the coordinate vector
              v.push_back(atof(coord.c_str()));
              ++coord_counter;
            }

            // if the number of coordinates does not match the coord_per_vertex variable
            if (coord_counter != current_mesh.coord_per_vertex)
            {
              throw SyntaxError("Wrong coordinate format in line " + stringify(cur_line));
            }

            // add to the coordinate stack
            (current_mesh.coords).push_back(v);

          }
        }
      } // coords sub-chunk

      // if it is an adjacency sub chunk
      else if((found = line.find('@')) != std::string::npos)
      {
        // separate face and shape
        line.trim_me();
        face = line.substr(0, found);
        shape = line.substr(found + 1);

        face.pop_front();
        shape.pop_back();

        face.trim_me();
        shape.trim_me();

        // get the face-dimension
        if(face == "vert")
        {
          face_dim = 0;
        }
        else if(face == "edge")
        {
          face_dim = 1;
        }
        else if(face == "quad" || face == "tria")
        {
          face_dim = 2;
        }
        else if(face == "hexa" || face == "tetra")
        {
          face_dim = 3;
        }
        else
        {
          throw SyntaxError("Unknown format in line " + stringify(cur_line));
        }

        // get the shape dimension
        if(shape == "vert")
        {
          shape_dim = 0;
        }
        else if(shape == "edge")
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
            current_mesh.adjacencies[face_dim][shape_dim] = a_stack;
            break;
          }
          else
          {
            // current adjacency vector
            std::vector<Index> a;
            while(!line.empty())
            {
              // separating coordinate
              found = line.find(' ');
              adj = line.substr(0, found);
              if(found != std::string::npos)
              {
                line = line.substr(found + 1);
              }
              else
              {
                line = "";
              }
              line.trim_me();
              adj.trim_me();
              // add to the adjacency vector
              a.push_back(atoi(adj.c_str()));
            }
            // add the vector to the adjacency stack
            a_stack.push_back(a);
          }
        }
      } // adjacencies sub-chunk

      // if it is a parent index chunk (of a submesh)
      else if((line.compare(0, 10, "<vert_idx>") == 0 ||
              line.compare(0, 10, "<edge_idx>") == 0 ||
              line.compare(0, 10, "<tria_idx>") == 0 ||
              line.compare(0, 10, "<quad_idx>") == 0 ||
              line.compare(0, 11, "<tetra_idx>") == 0 ||
              line.compare(0, 10, "<hexa_idx>") == 0) && flag == "sub")
      {

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

          // if it is the end of the parent index chunk
          if(line.find("idx") != String::npos)
          {
            current_mesh.parent_indices[dim] = idx;
            break;
          }
          else
          {
            idx.push_back(atoi(line.c_str()));
          }
        }
      } // parent index sub chunk

      // if it is the end of the mesh section
      else if(line == break_line)
      {
        meshes.push_back(current_mesh);
        return cur_line;
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    } // while

    // pointless default case to avoid warnings
    return 42;

  } // MeshReader::parse_mesh_section(Index cur_line, std::istream& ifs)


  // parses the cellset-section stream ifs
  Index MeshReader::parse_cellset_section(Index cur_line, std::istream& ifs)
  {
    CONTEXT("MeshReader::parse_cellset_section");

    // (auxiliary) variables
    String line;

    // cellset data container of the root mesh
    CellSetContainer current_set;

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
        // skip the info block
        while(line != "</info>" && !ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          line.trim_me();
          ++cur_line;
        }
      } // info sub-chunk

      // if it is the counts sub chunk
      else if(line.compare(0, 8, "<counts>") == 0)
      {
        // default case: zero
        current_set.vertex_number = 0;
        current_set.edge_number = 0;
        current_set.quad_number = 0;
        current_set.tria_number = 0;
        current_set.tetra_number = 0;
        current_set.hexa_number = 0;

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
            current_set.vertex_number = atoi(line.c_str());
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
            current_set.edge_number = atoi(line.c_str());
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
            current_set.tria_number = atoi(line.c_str());
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
            current_set.quad_number = atoi(line.c_str());
          }

          // tetra number
          else if(line.compare(0, 7, "tetras ") == 0)
          {
            // erase "hexas "
            line.erase(0, 7);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing tetra number in line " + stringify(cur_line));
            }
            current_set.tetra_number = atoi(line.c_str());
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
            current_set.hexa_number = atoi(line.c_str());
          }
          // if it is the end of the sub chunk
          else if(line == "</counts>")
          {
            break;
          }
          // if the line is empty
          else if(line.empty())
          {
            continue;
          }
          else
          {
            throw SyntaxError("Unknown file format in line " + stringify(cur_line));
          }
        }
      } // counts sub-chunk

      // if it is a parent index chunk
      else if(line.compare(0, 10, "<vert_idx>") == 0 ||
              line.compare(0, 10, "<edge_idx>") == 0 ||
              line.compare(0, 10, "<tria_idx>") == 0 ||
              line.compare(0, 10, "<quad_idx>") == 0 ||
              line.compare(0, 11, "<tetra_idx>") == 0 ||
              line.compare(0, 10, "<hexa_idx>") == 0)
      {

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

          // if it is the end of the parent index chunk
          if(line.find("idx") != std::string::npos)
          {
            current_set.parent_indices[dim] = idx;
            break;
          }
          else
          {
            idx.push_back(atoi(line.c_str()));
          }
        }
      } // parent index sub chunk

      // if it is the end of the cellset section
      else if(line == "</cellset>")
      {
        cell_sets.push_back(current_set);
        return cur_line;
      }
      else
      {
        throw SyntaxError("Unknown file format in line " + stringify(cur_line));
      }
    } // while

    // pointless default case to avoid warnings
    return 42;

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
    Index coord_counter;
    String::size_type found;
    String coord;
    Index cur_line2 = 0;

    // ignore everything until the coord file begins
    while(line != "<feast_coord_file>" && !ifs.eof() && ifs.good())
    {
      getline(ifs, line);
      ++cur_line2;
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
      ++cur_line2;

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
          ++cur_line2;
          line.trim_me();

          if(line.compare(0, 8, "version ") == 0)
          {

            // erase "version "
            line.erase(0, 8);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing version number in line " + stringify(cur_line2));
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
              throw SyntaxError("Missing number of vertices in line " + stringify(cur_line2));
            }
            mesh->vertex_number = atoi(line.c_str());
          }
          else if(line.compare(0, 7, "coords ") == 0)
          {
            // erase "coords "
            line.erase(0, 7);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing coordinate number in line " + stringify(cur_line2));
            }
            mesh->coord_per_vertex = atoi(line.c_str());
          }
          else if(line == "</header>")
          {
            break;
          }
          else
          {
            throw SyntaxError("Unknown file format in line " + stringify(cur_line2));
          }
        } // while
      } // header sub chunk

      // if it is an info chunk
      else if(line == "<info>")
      {
        while(line != "</info>" && !ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line2;
          line.trim_me();
        }
      }

      // if it is a coord chunk
      else if(line.compare(0, 8, "<coords>") == 0)
      {
        while(!ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line2;
          line.trim_me();

          // if it is the end of the sub chunk
          if(line == "</coords>")
          {
            break;
          }
          else
          {
            coord_counter = 0;
            std::vector<double> v;
            while(!line.empty())
            {
              // separating coordinate
              found = line.find(' ');
              coord = line.substr(0, found);
              if(found != std::string::npos)
              {
                line = line.substr(found + 1);
              }
              else
              {
                line = "";
              }
              line.trim_me();
              coord.trim_me();
              // add to the coordinate vector
              v.push_back(atof(coord.c_str()));
              ++coord_counter;
            }

            if (coord_counter != mesh->coord_per_vertex)
            {
              throw SyntaxError("Wrong coordinate format in line " + stringify(cur_line2));
            }

            // add to the coordinate stack
            (mesh->coords).push_back(v);
          }
        }
      } // coords sub-chunk

      // if it is the end
      else if(line == "</feast_coord_file>")
      {
        break;
      }
      else
      {
        throw SyntaxError("Unknown format in line " + stringify(cur_line2));
      }
    } // end while

    // if the file did not end properly
    if(line != "</feast_coord_file>")
    {
      throw SyntaxError("Unexpected file ending in line " + stringify(cur_line2));
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
    Index face_dim = 0;
    Index shape_dim = 0;
    String::size_type found;
    String adj, face, shape, break_line;
    Index cur_line2 = 0;

    // ignore everything until the adjacency file begins
    while(line != "<feast_adjacency_file>" && !ifs.eof() && ifs.good())
    {
      getline(ifs, line);
      ++cur_line2;
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
      ++cur_line2;

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
          ++cur_line2;
          line.trim_me();

          if(line.compare(0, 8, "version ") == 0)
          {
            // erase "version "
            line.erase(0, 8);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing version number in line " + stringify(cur_line2));
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
              throw SyntaxError("Missing type in line " + stringify(cur_line2));
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
              throw SyntaxError("Missing shape type in line " + stringify(cur_line2));
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
            throw SyntaxError("Unknown file format in line " + stringify(cur_line2));
          }
        } // while
      } // header sub chunk

      // if it is an info chunk
      else if(line == "<info>")
      {
        while(line != "</info>" && !ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line2;
          line.trim_me();
        }
      }

      // if it is the counts sub chunk
      else if(line.compare(0, 8, "<counts>") == 0)
      {

        // default setting: zero
        mesh->vertex_number = 0;
        mesh->edge_number = 0;
        mesh->quad_number = 0;
        mesh->tria_number = 0;
        mesh->tetra_number = 0;
        mesh->hexa_number = 0;

        while(!ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line2;
          line.trim_me();

          // vertex number
          if(line.compare(0, 6, "verts ") == 0)
          {
            // erase "verts "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing vertex number in line " + stringify(cur_line2));
            }
            mesh->vertex_number = atoi(line.c_str());
          }

          // edge number
          else if(line.compare(0, 6, "edges ") == 0)
          {
            // erase "edges "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing edge number in line " + stringify(cur_line2));
            }
            mesh->edge_number = atoi(line.c_str());
          }

          // trias number
          else if(line.compare(0, 6, "trias ") == 0)
          {
            // erase "trias "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing triangle number in line " + stringify(cur_line2));
            }
            mesh->tria_number = atoi(line.c_str());
          }

          // quad number
          else if(line.compare(0, 6, "quads ") == 0)
          {
            // erase "quads "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing quad number in line " + stringify(cur_line2));
            }
            mesh->quad_number = atoi(line.c_str());
          }

          // tetra number
          else if(line.compare(0, 7, "tetras ") == 0)
          {
            // erase "tetras "
            line.erase(0, 7);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing tetra number in line " + stringify(cur_line2));
            }
            mesh->tetra_number = atoi(line.c_str());
          }

          // hexa number
          else if(line.compare(0, 6, "hexas ") == 0)
          {
            // erase "hexas "
            line.erase(0, 6);
            line.trim_me();
            if(line.empty())
            {
              throw SyntaxError("Missing hexa number in line " + stringify(cur_line2));
            }
            mesh->hexa_number = atoi(line.c_str());
          }
          // if it is the end of the sub chunk
          else if(line == "</counts>")
          {
            break;
          }
          else if(line.empty())
          {
            continue;
          }
          else
          {
            throw SyntaxError("Unknown file format in line " + stringify(cur_line2));
          }
        }
      } // counts sub-chunk

      // if it is an adjacency sub chunk
      else if((found = line.find('@')) != std::string::npos)
      {
        // separate face and shape
        line.trim_me();
        face = line.substr(0, found);
        shape = line.substr(found + 1);

        face.pop_front();
        shape.pop_back();

        face.trim_me();
        shape.trim_me();

        // get the face dimension
        if(face == "vert")
        {
          face_dim = 0;
        }
        else if(face == "edge")
        {
          face_dim = 1;
        }
        else if(face == "quad" || face == "tria")
        {
          face_dim = 2;
        }
        else if(face == "hexa" || face == "tetra")
        {
          face_dim = 3;
        }
        else
        {
          throw SyntaxError("Unknown format in line " + stringify(cur_line2));
        }

        // get the shape dimension
        if(shape == "vert")
        {
          shape_dim = 0;
        }
        else if(shape == "edge")
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
          throw SyntaxError("Unknown format in line " + stringify(cur_line2));
        }

        std::vector<std::vector<Index> > a_stack;

        while(!ifs.eof() && ifs.good())
        {
          // get a line
          getline(ifs, line);
          ++cur_line2;
          line.trim_me();

          // if it is the end of the sub chunk
          if(line.find('@') != std::string::npos)
          {
            mesh->adjacencies[face_dim][shape_dim] = a_stack;
            break;
          }
          else
          {
            std::vector<Index> a;
            while(!line.empty())
            {
              // separating coordinates
              found = line.find(' ');
              adj = line.substr(0, found);
              if(found != std::string::npos)
              {
                line = line.substr(found + 1);
              }
              else
              {
                line = "";
              }
              line.trim_me();
              adj.trim_me();
              // add to the adjacency vector
              a.push_back(atoi(adj.c_str()));
            }
            // add to the adjacency list
            a_stack.push_back(a);
          }
        }
      } // adjacencies sub-chunk

      // if it is the end of the file
      else if(line == "</feast_adjacency_file>")
      {
        break;
      }
      else
      {
        throw SyntaxError("Unknown format in line " + stringify(cur_line2));
      }
    } // end while

    // if the file did not end properly
    if(line != "</feast_adjacency_file>")
    {
      throw SyntaxError("Unexpected file ending in line " + stringify(cur_line2));
    }
  } // MeshReader::parse_adjacency_file(std::istream& ifs, MeshReader::MeshDataContainer *mesh)

  // returns the version
  String MeshReader::get_version()
  {
    CONTEXT("MeshReader::get_version()");
    return version;
  }

  // returns the chart path
  String MeshReader::get_chart_path()
  {
    CONTEXT("MeshReader::get_chart_path()");
    return chart_path;
  }

  // returns the number of submeshes
  Index MeshReader::get_number_of_submeshes()
  {
    CONTEXT("MeshReader::get_number_of_submeshes()");
    return number_of_submeshes;
  }

  // returns the number of cellsets
  Index MeshReader::get_number_of_cellsets()
  {
    CONTEXT("MeshReader::get_number_of_cellsets()");
    return number_of_cellsets;
  }

  // returns the mesh with the name "mesh_name"
  std::pair<MeshReader::MeshDataContainer, bool> MeshReader::get_mesh(String mesh_name)
  {
    CONTEXT("MeshReader::get_mesh()");

    // trim the input name
    mesh_name.trim_me();

    std::vector<MeshDataContainer>::const_iterator iter;
    for(iter = meshes.begin(); iter != meshes.end() ; ++iter)
    {
      if(iter->name == mesh_name)
      {
        return std::make_pair(*iter, true);
      }
    }
    return std::make_pair(*meshes.begin(), false);
  } // MeshReader::get_mesh(String mesh_name)

  // returns the cellset with the name "cellset_name"
  std::pair<MeshReader::CellSetContainer, bool> MeshReader::get_cellset(String cellset_name)
  {
    CONTEXT("MeshReader::get_cellset()");

    // trim the input name
    cellset_name.trim_me();

    std::vector<CellSetContainer>::const_iterator iter;
    for(iter = cell_sets.begin(); iter !=cell_sets.end() ; ++iter)
    {
      if(iter->name == cellset_name)
      {
        return std::make_pair(*iter, true);
      }
    }
    return std::make_pair(*cell_sets.begin(), false);
  } // MeshReader::get_cellset(String cellset_name)

  bool MeshReader::no_meshes()
  {
    CONTEXT("MeshReader::no_meshes()");
    return meshes.empty();
  }

  bool MeshReader::no_cellsets()
  {
    CONTEXT("MeshReader::no_cellsets()");
    return cell_sets.empty();
  }

} //namespace FEAST
