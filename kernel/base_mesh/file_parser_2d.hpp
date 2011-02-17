#pragma once
#ifndef KERNEL_BM_FILE_PARSER_2D_HPP
#define KERNEL_BM_FILE_PARSER_2D_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
//#include <vector>   // for std::vector
#include <math.h>   // for sin, cos

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/file_reader_ascii.hpp>
// COMMENT_HILMAR: wenn Fusion mit MPI vollendet, dann nutze das hier:
//#include <kernel/error_handler.hpp>
#include <kernel/base_mesh/base_mesh_2d.hpp>
//#include <kernel/base_mesh/vertex.hpp>
//#include <kernel/base_mesh/cell.hpp>
//#include <kernel/base_mesh/cell_1d_edge.hpp>
//#include <kernel/base_mesh/cell_2d_quad.hpp>
//#include <kernel/base_mesh/cell_2d_tri.hpp>

namespace FEAST
{
  namespace BaseMesh
  {
    /**
    * \brief base mesh file parser for 2D files
    *
    * This class parses 2D mesh files in FEAST1 format. It is friend of BaseMesh2D to have full access to BaseMesh2D's
    * private members.
    *
    * \note Currently, we only have 2D files, so I began with a 2D file reader only. I didn't care yet how to
    * distinguish/unite 1D, 2D and 3D. However we do it... we should try to avoid too much code duplication. Most
    * probably, there will be only *one* file parser for 1D, 2D and 3D files...
    *
    * \note This is only a very rudimentary file parser which will be rewritten when we have defined the new file
    * format. It only supports some features for basic test purposes.
    *
    * \note FEAST1 file format does not distinguish world dimension and space dimension. When it does in future,
    * this class needs the corresponding template parameter world_dim_ (similar to class BaseMesh2D).
    *
    * \author Dominik Goeddeke
    * \author Hilmar Wobker
    */
    class FileParser2D
    {
      /// shortcuts various cell types to save typing of template parameters
      typedef Vertex<2> Vertex_;
      typedef Edge<2, 2> Edge_;
//      typedef Tri<2, 2> Tri_;
      typedef Quad<2, 2> Quad_;
//      typedef Cell<2, 2, 2> Cell_;


    private:

      /// pointer to the base mesh
      BaseMesh2D<2>* _bm;
      /// number boundaries in the mesh file
      unsigned int _num_boundaries;
      /// array for storing number of segments per boundary
      unsigned int* _num_segments;
      /// array for storing the segment types (line or circle)
      unsigned int** _segment_type;
      /// cartesian coordinates of the segments' start vertices, in case of circles: centre
      double*** _start_vertex;
      /// the segments' vectors from start vertex to the next vertex, in case of circles: radius + dummy
      double*** _vector_to_end_vertex;
      /// the segments' start and end parameter values
      double*** _param_int;
      /// the segments' start and end angle in case of circle sections
      double*** _circle_section;

      /// number of vertices
      unsigned int _num_vertices;
      /// number of edges
      unsigned int _num_edges;
      /// number of cells
      unsigned int _num_cells;

      /// reads FEAST_GEOM to extract vertex coordinates
      void _parse_geom(FileReaderASCII* mesh_file)
      {
        mesh_file->read("FEAST_GEOM");
        mesh_file->read("2D");
        int version_major;
        mesh_file->read(version_major);
        int version_minor;
        mesh_file->read(version_minor);

        // number of boundary components
        mesh_file->read(_num_boundaries);

        // allocate arrays
        _num_segments = new unsigned int[_num_boundaries];
        _segment_type = new unsigned int*[_num_boundaries];
        _param_int = new double**[_num_boundaries];
        _start_vertex = new double**[_num_boundaries];
        _vector_to_end_vertex = new double**[_num_boundaries];
        _circle_section = new double**[_num_boundaries];

        // loop over boundary components
        for(unsigned int ibc(0) ; ibc < _num_boundaries ; ++ibc)
        {
          // read number of segments
          mesh_file->read(_num_segments[ibc]);
          // segment type
          _segment_type[ibc] = new unsigned int[_num_segments[ibc]];
          // cartesian coordinates of the segments' start vertices, in case of circles: centre
          _start_vertex[ibc] = new double*[_num_segments[ibc]];
          // the segments' vectors from start vertex to the next vertex, in case of circles: radius + dummy
          _vector_to_end_vertex[ibc] = new double*[_num_segments[ibc]];
          // the segments' start and end parameter values
          _param_int[ibc] = new double*[_num_segments[ibc]];
          // start and end angle of circle section
          _circle_section[ibc] = new double*[_num_segments[ibc]];

          for(unsigned int iseg(0) ; iseg < _num_segments[ibc] ; ++iseg)
          {
            _start_vertex[ibc][iseg] = new double[2];
            _vector_to_end_vertex[ibc][iseg] = new double[2];
            _param_int[ibc][iseg] = new double[2];
            // read segment type
            mesh_file->read(_segment_type[ibc][iseg]);
//std::cout << "segment " << iseg << ", type " << _segment_type[ibc][iseg] << std::endl;
            switch(_segment_type[ibc][iseg])
            {
              // line
              case 0:
              {
                // read start vertex
                mesh_file->read(_start_vertex[ibc][iseg]);
                // read vector from start to end vertex
                mesh_file->read(_vector_to_end_vertex[ibc][iseg]);
                // read parameter values
                mesh_file->read(_param_int[ibc][iseg]);
//std::cout << "start vert: " << _start_vertex[ibc][iseg][0] << " " << _start_vertex[ibc][iseg][1] << std::endl;
//std::cout << "vector    : " << _vector_to_end_vertex[ibc][iseg][0] << " " << _vector_to_end_vertex[ibc][iseg][1] << std::endl;
//std::cout << "param     : " << _param_int[ibc][iseg][0] << " " << _param_int[ibc][iseg][1] << std::endl;
                break;
              }
              // circle
              case 1:
              {
//                throw InternalError("Circle segments are not supported in this file parser.");
                _circle_section[ibc][iseg] = new double[2];
                // read centre
                mesh_file->read(_start_vertex[ibc][iseg]);
                // read radius and dummy
                mesh_file->read(_vector_to_end_vertex[ibc][iseg]);
                // read start and end angle
                mesh_file->read(_circle_section[ibc][iseg]);
                // read parameter values
                mesh_file->read(_param_int[ibc][iseg]);
//std::cout << "centre : " << _start_vertex[ibc][iseg][0] << " " << _start_vertex[ibc][iseg][1] << std::endl;
//std::cout << "radius : " << _vector_to_end_vertex[ibc][iseg][0] << std::endl;
//std::cout << "section: " << _circle_section[ibc][iseg][0] << " " << _circle_section[ibc][iseg][1] << std::endl;
//std::cout << "param  : " << _param_int[ibc][iseg][0] << " " << _param_int[ibc][iseg][1] << std::endl;
                break;
              }
              default:
              {
                throw InternalError("Unknown/unhandled segment type found.");
              }
            }
          }
        }
      } // _parse_geom()

      /// reads and skips FEAST_FGEOM parts
      void _parse_fgeom(FileReaderASCII* mesh_file)
      {
        mesh_file->read("FEAST_FGEOM");
        mesh_file->read("2D");
        int version_major;
        mesh_file->read(version_major);
        int version_minor;
        mesh_file->read(version_minor);

        // number of boundary components
        unsigned int num_boundaries;
        mesh_file->read(num_boundaries);

        // loop over boundary components
        for(unsigned int ibc(0) ; ibc < num_boundaries ; ++ibc)
        {
          unsigned int num_segments;
          mesh_file->read(num_segments);
          for(unsigned int iseg(0) ; iseg < num_segments ; ++iseg)
          {
            int seg_type;
            mesh_file->read(seg_type);

            switch(seg_type)
            {
              // line
              case 0:
              {
                //FGEOM: skip three lines
                mesh_file->skip_line();
                mesh_file->skip_line();
                mesh_file->skip_line();
                break;
              }
              // circle
              case 1:
              {
                // FGEOM: skip four lines
                mesh_file->skip_line();
                mesh_file->skip_line();
                mesh_file->skip_line();
                mesh_file->skip_line();
                break;
              }
              default:
              {
                throw InternalError("Unknown/unhandled boundary type found.");
              }
            }
          }
        }
      } // _parse_fgeom()


      /// reads MESH portion of FEAST files into this BaseMesh
      void _parse_mesh(FileReaderASCII * mesh_file)
      {
        // header
        mesh_file->read("FEAST_MESH");
        mesh_file->read("2D");
        int version_major;
        mesh_file->read(version_major);
        int version_minor;
        mesh_file->read(version_minor);
        int type;
        mesh_file->read(type);
        // number of entities
        mesh_file->read(_num_vertices);
        mesh_file->read(_num_cells);
        mesh_file->read(_num_edges);

        // read and set vertices
        for(unsigned int ivertex = 0 ; ivertex < _num_vertices ; ++ivertex)
        {
          double coords[2];
          double x, y;
          unsigned int itype;
          mesh_file->read(x, y, itype);
          // distinguish inner vertices and boundary vertices
          if(itype == 0)
          {
            // inner vertex
//std::cout << "found inner vertex " << x << " " << y << std::endl;
            coords[0] = x;
            coords[1] = y;
          }
          else
          {
            unsigned int ibc(itype - 1);
            assert(ibc >= 0);
            assert(ibc < _num_boundaries);
            // search parameter interval
            for(unsigned int iseg = 0 ; iseg < _num_segments[ibc] ; ++iseg)
            {
              double pstart = _param_int[ibc][iseg][0];
              double pend = _param_int[ibc][iseg][1];
              if(x >= pstart && x <= pend)
              {
                // interval found, distinguish segment type
                switch(_segment_type[ibc][iseg])
                {
                  // line segment
                  case 0:
                  {
                    double denom = pend - pstart;
                    assert (denom > 0);
                    double factor = (x - pstart)/denom;
                    // compute cartesian coordinates
                    coords[0] = _start_vertex[ibc][iseg][0] + factor * _vector_to_end_vertex[ibc][iseg][0];
                    coords[1] = _start_vertex[ibc][iseg][1] + factor * _vector_to_end_vertex[ibc][iseg][1];
//std::cout << "found boundary vertex on line seg: " << coords[0] << " " << coords[1] << std::endl;
                    break;
                  }
                  // circle segment
                  case 1:
                  {
                    double radius = _vector_to_end_vertex[ibc][iseg][0];
                    double denom = pend - pstart;
                    assert (denom > 0);
                    double angle = (  (pend - x)   * _circle_section[ibc][iseg][0]
                                    + (x - pstart) * _circle_section[ibc][iseg][1])/denom;
                    // compute cartesian coordinates
                    coords[0] = _start_vertex[ibc][iseg][0] + radius * cos(angle);
                    coords[1] = _start_vertex[ibc][iseg][1] + radius * sin(angle);
//std::cout << "found boundary vertex on circle seg: " << coords[0] << " " << coords[1] << std::endl;
                  }
                }
                // interval found, exit for loop
                break;
              }
            }
          }
          // add vertex to base mesh
          _bm->_add(new Vertex_(coords));
        }

        // read  and set edges
        for(unsigned int iedge = 0 ; iedge < _num_edges ; ++iedge)
        {
          unsigned int v0, v1;
          mesh_file->read(v0, v1);
          // substract 1 from 1-based to 0-based indexing
          _bm->_add(new Edge_(_bm->vertex(v0-1), _bm->vertex(v1-1), 0));
        }

        // vector for buffering neighbourhood information that can not be immediately set
        std::vector<int*> buffered_neighbourhood;

        // read and set cells
        for(unsigned int icell = 0 ; icell < _num_cells ; ++icell)
        {
          // type of this cell
          unsigned int type;
          mesh_file->read(type);
          switch(type)
          {
            case(0):
            {
              // quad

              // set number of vertices and edges
              unsigned int nverts = 4;
              unsigned int nedges = 4;

              // when reading vertex and edge indices, map them from standard ccw-numbering to new numbering scheme:
              //      v3---e2----v2                v2---e1----v3
              //      |           |                |           |
              //      |           |                |           |
              //      e3         e1     ------>    e2         e3
              //      |           |                |           |
              //      |           |                |           |
              //      v0---e0----v1                v0---e0----v1
              //
              // replace verts: 0 -> 0, 1 -> 1, 2 -> 3, 3 -> 2
              // replace edges: 0 -> 0, 1 -> 3, 2 -> 1, 3 -> 2

              // read vertex indices
              unsigned int* v = new unsigned int[nverts];
              mesh_file->read(v[0], v[1], v[3], v[2]);

              // read edge indices
              unsigned int* e = new unsigned int[nedges];
              mesh_file->read(e[0], e[3], e[1], e[2]);

              // substract 1 from 1-based to 0-based indexing
              for(unsigned int ivert = 0 ; ivert < nverts ; ++ivert)
              {
                v[ivert]--;
              }
              for(unsigned int iedge = 0 ; iedge < nedges ; ++iedge)
              {
                e[iedge]--;
              }
//std::cout << "Vertices: " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;
//std::cout << "Edges   : " << e[0] << " "<< e[1] << " " << e[2] << " " << e[3] << std::endl;

              _bm->_add(new Quad_(_bm->vertex(v[0]), _bm->vertex(v[1]), _bm->vertex(v[2]), _bm->vertex(v[3]),
                                  _bm->edge(e[0]), _bm->edge(e[1]), _bm->edge(e[2]), _bm->edge(e[3]), 0));
              delete [] v;
              delete [] e;
              // read edge neighbours
              int eneigh[4];
              mesh_file->read(eneigh[0], eneigh[3], eneigh[1], eneigh[2]);
              // substract 1 from 1-based to 0-based indexing
              for(unsigned int iedge = 0 ; iedge < nedges ; ++iedge)
              {
                eneigh[iedge]--;
              }
//std::cout << "edge neighs: " << eneigh[0] << " "<< eneigh[1] << " " << eneigh[2] << " " << eneigh[3] << std::endl;

              // set edge neighbourhood information
              for(unsigned int iedge = 0 ; iedge < nedges ; ++iedge)
              {
                if(eneigh[iedge] >= 0)
                {
                  if((unsigned int)eneigh[iedge] < icell)
                  {
                    // when the index of the neighbour cell is less than that of the current cell, the neighbour cell
                    // already exists and can be set as neighbour
                    _bm->cell(icell)->add_neighbour(SDIM_EDGE, iedge, _bm->cell(eneigh[iedge]));
                  }
                  else
                  {
                    //otherwise, the information has to be buffered and neighbourhood has to be set later on
                    int* buf = new int[4];
                    buf[0] = icell; buf[1] = SDIM_EDGE; buf[2] = iedge; buf[3] = eneigh[iedge];
                    buffered_neighbourhood.push_back(buf);
                  }
                }
              }

              // read vertex neighbours
              // Not only map from standard ccw-numbering to new numbering scheme here, but at the same time shift
              // the vertex neighbours one to the left. Why? Since for some strange reason
              //   0 0 4 0
              // does not mean "cell 4 is vertex neighbour at the third vertex", but "cell 4 is vertex neighbour at the
              // vertex which is the *end* vertex of the third edge" (i.e. the fourth vertex). To get what one would
              // expect, the line would have to be changed to
              //   0 0 0 4
              // (which is equivalent to shifting the vertex neighbours one to the left).
              // Original:            v0, v1, v2, v3
              // change of numbering: v0, v1, v3, v2
              // shilt:               v1, v3, v2, v0
              int vneigh[4];
              mesh_file->read(vneigh[1], vneigh[3], vneigh[2], vneigh[0]);
//std::cout << "vert neighs: ";
              for(unsigned int ivert = 0 ; ivert < nverts ; ++ivert)
              {
                if(vneigh[ivert] > 0)
                {
//std::cout << vneigh[ivert]-1 << " ";
                  if((unsigned int)vneigh[ivert]-1 < icell)
                  {
                    // when the index of the neighbour cell is less than that of the current cell, the neighbour cell
                    // already exists and can be set as neighbour
                    _bm->cell(icell)->add_neighbour(SDIM_VERTEX, ivert, _bm->cell(vneigh[ivert]-1));
                  }
                  else
                  {
                    //otherwise, the information has to be buffered and neighbourhood has to be set later on
                    int* buf = new int[4];
                    buf[0] = icell; buf[1] = SDIM_VERTEX; buf[2] = ivert; buf[3] = vneigh[ivert]-1;
                    buffered_neighbourhood.push_back(buf);
                  }
                }
                else if(vneigh[ivert] < 0)
                {
//std::cout << "[";
                  // special case: negative number means several vertex neighbours
                  unsigned int num_extra = -vneigh[ivert];
                  int* vdneigh = new int[num_extra];
                  // lacking a proper "line tokeniser", perform this intermediate hack
                  switch(num_extra)
                  {
                    case 1:
                    {
                      mesh_file->read(vdneigh[0]);
                      break;
                    }
                    case 2:
                    {
                      mesh_file->read(vdneigh[0], vdneigh[1]);
                      break;
                    }
                    case 3:
                    {
                      mesh_file->read(vdneigh[0], vdneigh[1], vdneigh[2]);
                      break;
                    }
                    case 4:
                    {
                      mesh_file->read(vdneigh[0], vdneigh[1], vdneigh[2], vdneigh[3]);
                      break;
                    }
                    default:
                    {
                      throw InternalError("More than 4 vertex neighbours cannot be parsed.");
                    }
                  }
                  for(unsigned int ineigh = 0 ; ineigh < num_extra ; ++ineigh)
                  {
//std::cout << vdneigh[ineigh]-1 << " ";
                    if((unsigned int)vdneigh[ineigh]-1 < icell)
                    {
                      // when the index of the neighbour cell is less than that of the current cell, the neighbour cell
                      // already exists and can be set as neighbour
                      _bm->cell(icell)->add_neighbour(SDIM_VERTEX, ivert, _bm->cell(vdneigh[ineigh]-1));
                    }
                    else
                    {
                      //otherwise, the information has to be buffered and neighbourhood has to be set later on
                      int* buf = new int[4];
                      buf[0] = icell; buf[1] = SDIM_VERTEX; buf[2] = ivert; buf[3] = vdneigh[ineigh]-1;
                      buffered_neighbourhood.push_back(buf);
                    }
                  }
//std::cout << "] ";
                  delete [] vdneigh;
                }
//                else
//                {
//std::cout << "x ";
//                }
              } // loop over point neighbours
//std::cout << std::endl;

              // finally, skip the line with refinement level etc.
              int ref_level, ref_mode;
              double ref_factor1, ref_factor2, ref_factor3;
              mesh_file->read(ref_level, ref_mode, ref_factor1, ref_factor2, ref_factor3);
              // and continue with next cell
              break;
            } // case(0) (quad)

            case(1):
            {
              // triangle
              throw InternalError("Parsing of triangles not implemented yet.");
            }

            default:
            {
              throw InternalError("Unknown cell type <" + StringUtils::stringify(type) + "> found.");
            }
          } // switch(type)
        } // for(unsigned int icell = 0 ; icell < _num_cells ; ++icell)

        // now set buffered neighbourhood information
        for(unsigned int ipos = 0 ; ipos < buffered_neighbourhood.size() ; ++ipos)
        {
          int* buf = buffered_neighbourhood[ipos];
          _bm->cell(buf[0])->add_neighbour((subdim)buf[1], buf[2], _bm->cell(buf[3]));
          delete [] buf;
        }
        buffered_neighbourhood.clear();
      } // _parse_mesh

/*
COMMENT_HILMAR: will be adapted later
      /// parses partition part of FEAST files
      void _parse_partition(FileReaderASCII * mesh_file)
      {
        // header
        mesh_file->read("FEAST_PART");
        mesh_file->read("2D");
        int version_major;
        mesh_file->read(version_major);
        int version_minor;
        mesh_file->read(version_minor);
        // number of entities
        unsigned int num_cells;
        mesh_file->read(num_cells);
        if(num_cells != _num_cells)
          throw InternalError("Number of BaseMeshCells from mesh file does not match that of partition file.");
        mesh_file->read(_num_partitions);
        // TODO: THIS IS A HACK UNTIL HILMAR AND I DECIDE WHAT DTO DO WITH THE BASEMESH
        if(_num_cells != _num_partitions)
          throw InternalError("Sorry. Hilmar and Dominik currently only support one basemeshcell per process.");

        // cells
        for(unsigned int icell = 0 ; icell < _num_cells; ++icell)
        {
          unsigned int pb, mb;
          mesh_file->read(pb, mb);
          _cells.at(icell)->set_parallel_block(pb-1);
          _cells.at(icell)->set_matrix_block(mb-1);
        }
      } // _parse_partition
*/

    public:

      /* ***************************
      * constructors & destructors *
      *****************************/
//COMMENT_HILMAR: DTOR and clean up missing!

      /* *****************
      * member functions *
      *******************/
      /**
      * \brief function for reading a 2D mesh in FEAST1 format
      *
      * \author Dominik Goeddeke
      * \author Hilmar Wobker
      */
      inline void parse(std::string const& file_name, BaseMesh2D<2>* bm)
      {
        FileReaderASCII* mesh_file = new FileReaderASCII(file_name, '#', true);

        // set base mesh pointer
        _bm = bm;

        // file header
        mesh_file->read("FEAST");
        mesh_file->read("2D");
        int version_major;
        mesh_file->read(version_major);
        int version_minor;
        mesh_file->read(version_minor);
        if(version_major != 3 && version_minor != 0)
        {
          throw InternalError("Only file format 3.0 (inline) is currently supported");
        }
        // description of this base mesh
        std::string saux;
        mesh_file->read(saux);

        // remaining parts of this file
        _parse_geom(mesh_file);
        _parse_fgeom(mesh_file);
        _parse_mesh(mesh_file);
//        _parse_partition(mesh_file);
        // HACK: do not read in FEAST_BC and FEAST_PROP, instead, just close the file and get on with life

        // clean up file reader
        delete mesh_file;
      } // parse
    }; // class FileParser2D
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BM_FILE_PARSER_2D_HPP
