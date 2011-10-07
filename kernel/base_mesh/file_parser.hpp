/* GENERAL_REMARK_BY_HILMAR:
 * Parsing of part and prop files is missing (see issue00012).
 * Generally, Peter wanted to take a deeper look at the base mesh implementation.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_BASE_MESH_FILE_PARSER_HPP
#define KERNEL_BASE_MESH_FILE_PARSER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/file_reader_ascii.hpp>
#include <kernel/logger.hpp>
#include <kernel/base_mesh/cell_3d_hexa.hpp>
#include <kernel/base_mesh/bm.hpp>

// includes, system
#include <iostream> // for std::ostream
#include <math.h>   // for sin, cos
#include <assert.h>

namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
  namespace BaseMesh
  {
    /**
    * \brief template class for base mesh file parser
    *
    * To be specialised w.r.t. space_dim_
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Dominik Goeddeke
    * \author Hilmar Wobker
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class FileParser
    {
    };


    /**
    * \brief template class for base mesh file parser, specialisation for 1D files
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \note As long as no real file parser is implemented, we only manually create a 1D base mesh here.
    *
    * \author Dominik Goeddeke
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class FileParser<1, world_dim_>
    {
      /// shortcut to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;

      /// shortcut to save typing of template parameters
      typedef Edge<1, world_dim_> Edge_;

    private:

      /// pointer to the base mesh
      BM<1, world_dim_>* _bm;

    public:

      /* *****************
      * member functions *
      *******************/
      /**
      * \brief function for reading a 1D mesh in FEAST1 format
      *
      * \param[in] file_name
      * name of the mesh file
      *
      * \param[in,out] bm
      * base mesh to be built from the mesh file
      *
      * \note As long as no real file parser is implemented, we only manually create a 1D base mesh here.
      *
      * \author Dominik Goeddeke
      * \author Hilmar Wobker
      */
      inline void parse(String const& file_name, BM<1, world_dim_>* bm)
      {
        CONTEXT("BaseMesh::FileParser::parse()");
        ASSERT(world_dim_ == 1, "Currently, world dim " + stringify(world_dim_) + " must be 1.");

        String s = "No file parser implemented yet!\n";
        s += "A manual 1D base mesh is created (v0 is at (0,0) and all edges are unit length):\n";
        s += "v0--------e0--------v1--------e1--------v2--------e2--------v3\n";
        Logger::log(s, Logger::master);

        // TODO: remove the following line as soon as a real file parser is implemented
        Logger::log("Write file name " + file_name + " to prevent compiler warning about unused variable.");

        // set base mesh pointer
        _bm = bm;

        // Base mesh example consisting four vertices and three edges
        //   v0---e0---v1---e1---v2---e2---v3

        // create the four vertices
        // v0
        double coords[1];
        coords[0] = 0.0;

        // shortcut for _subcells object
        Subcells<1, 1, world_dim_>* sc = &_bm->_subcells;

        // add vertex to base mesh
        sc->_add(new Vertex_(coords));

        // v1
        coords[0] = 1.0;
        // add vertex to base mesh
        sc->_add(new Vertex_(coords));

        // v2
        coords[0] = 2.0;
        // add vertex to base mesh
        sc->_add(new Vertex_(coords));

        // v3
        coords[0] = 3.0;
        // add vertex to base mesh
        sc->_add(new Vertex_(coords));

        // create the three edges (=cells)
        // e0
        _bm->_add(new Edge_(sc->vertex(0), sc->vertex(1), 0));
        // e1
        _bm->_add(new Edge_(sc->vertex(1), sc->vertex(2), 0));
        // e2
        _bm->_add(new Edge_(sc->vertex(2), sc->vertex(3), 0));

        // set neighbourhood information
        _bm->cell(0)->add_neighbour(SDIM_VERTEX, 1, _bm->cell(1));
        _bm->cell(1)->add_neighbour(SDIM_VERTEX, 0, _bm->cell(0));
        _bm->cell(1)->add_neighbour(SDIM_VERTEX, 1, _bm->cell(2));
        _bm->cell(2)->add_neighbour(SDIM_VERTEX, 0, _bm->cell(1));
      } // parse
    }; // FileParser<1, world_dim_>



    /**
    * \brief template class for base mesh file parser, specialisation for 2D files
    *
    * This class parses 2D mesh files in FEAST1 format. It is friend of class BM to have full access to base mesh's
    * private members.
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \note Currently, we only have 2D files, so we implemented this 2D file reader only. Before doing the same for
    * 1D and 3D, one should think about how to avoid too much code duplication.
    *
    * \note This is only a very rudimentary file parser which will be rewritten when we have defined the new file
    * format (adapted the old one, resp.). It only supports some features for basic test purposes.
    *
    * \note FEAST1 file format does not distinguish world dimension and space dimension. The world_dim_ template
    * parameter actually has to be equal to 2 here.
    *
    * \author Dominik Goeddeke
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class FileParser<2, world_dim_>
    {
      /// shortcut to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;
      /// shortcut to save typing of template parameters
      typedef Edge<2, world_dim_> Edge_;
      /// shortcut to save typing of template parameters
      typedef Quad<2, world_dim_> Quad_;


    private:

      /* *****************
      * member variables *
      *******************/
      /// pointer to the base mesh
      BM<2, world_dim_>* _bm;
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

      /* *****************
      * member functions *
      *******************/
      /// reads FEAST_GEOM to extract vertex coordinates
      void _parse_geom(FileReaderASCII* mesh_file)
      {
        CONTEXT("BaseMesh::FileParser::_parse_geom()");
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
                _circle_section[ibc][iseg] = nullptr;
                // read start vertex
                mesh_file->read2(_start_vertex[ibc][iseg]);
                // read vector from start to end vertex
                mesh_file->read2(_vector_to_end_vertex[ibc][iseg]);
                // read parameter values
                mesh_file->read2(_param_int[ibc][iseg]);
//std::cout << "start vert: " << _start_vertex[ibc][iseg][0] << " " << _start_vertex[ibc][iseg][1] << std::endl;
//std::cout << "vector    : " << _vector_to_end_vertex[ibc][iseg][0] << " " << _vector_to_end_vertex[ibc][iseg][1] << std::endl;
//std::cout << "param     : " << _param_int[ibc][iseg][0] << " " << _param_int[ibc][iseg][1] << std::endl;
                break;
              }
              // circle
              case 1:
              {
                _circle_section[ibc][iseg] = new double[2];
                // read centre
                mesh_file->read2(_start_vertex[ibc][iseg]);
                // read radius and dummy
                mesh_file->read2(_vector_to_end_vertex[ibc][iseg]);
                // read start and end angle
                mesh_file->read2(_circle_section[ibc][iseg]);
                // read parameter values
                mesh_file->read2(_param_int[ibc][iseg]);
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
        CONTEXT("BaseMesh::FileParser::_parse_fgeom()");
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


      /// reads MESH portion of FEAST files into this base mesh
      void _parse_mesh(FileReaderASCII * mesh_file)
      {
        CONTEXT("BaseMesh::FileParser::_parse_mesh()");
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

        // shortcut for _subcells object
        Subcells<2, 2, world_dim_>* sc = &_bm->_subcells;

        // read and set vertices
        for(unsigned int ivertex = 0 ; ivertex < _num_vertices ; ++ivertex)
        {
          double coords[2] = {0};
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
            ASSERT(itype >= 1, "itype " + stringify(itype) + " must be >= 1.");
            unsigned int ibc(itype - 1);
            ASSERT(ibc < _num_boundaries, "Boundary index " + stringify(ibc) + " must not exceed number of boundaries"
                   + stringify(_num_boundaries) + ".");
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
                    break;
                  }
                  default:
                  {
                    throw InternalError("Unknown/unhandled segment type found.");
                  }
                }
                // interval found, exit for loop
                break;
              }
            }
          }
          // add vertex to base mesh
          sc->_add(new Vertex_(coords));
        }

        // read and set edges
        for(unsigned int iedge = 0 ; iedge < _num_edges ; ++iedge)
        {
          unsigned int v0, v1;
          mesh_file->read(v0, v1);
          ASSERT(v0 >= 1, "Vertex index " + stringify(v0) + " must be >=1.");
          ASSERT(v1 >= 1, "Vertex index " + stringify(v1) + " must be >=1.");
          // substract 1 from 1-based to 0-based indexing
          sc->_add(new Edge_(sc->vertex(v0-1), sc->vertex(v1-1), 0));
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

              _bm->_add(new Quad_(sc->vertex(v[0]), sc->vertex(v[1]), sc->vertex(v[2]), sc->vertex(v[3]),
                                  sc->edge(e[0]), sc->edge(e[1]), sc->edge(e[2]), sc->edge(e[3]), 0));
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
                      ASSERT(vdneigh[ineigh] >= 1, "Index " + stringify(vdneigh[ineigh]) + " of the "
                             + stringify(ineigh) + "-th neighbour must be >=1.");
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
              throw InternalError("Unknown cell type <" + stringify(type) + "> found.");
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

//COMMENT_HILMAR: see issue00012
//      /// parses partition part of FEAST files
//      void _parse_partition(FileReaderASCII * mesh_file)
//      {
//        CONTEXT("BaseMesh::FileParser::_parse_partition()");
//        // header
//        mesh_file->read("FEAST_PART");
//        mesh_file->read("2D");
//        int version_major;
//        mesh_file->read(version_major);
//        int version_minor;
//        mesh_file->read(version_minor);
//        // number of entities
//        unsigned int num_cells;
//        mesh_file->read(num_cells);
//        if(num_cells != _num_cells)
//          throw InternalError("Number of base mesh cells from mesh file does not match that of partition file.");
//        mesh_file->read(_num_partitions);
//        // TODO: THIS IS A HACK UNTIL HILMAR AND I DECIDE WHAT TO DO WITH THE BASE MESH
//        if(_num_cells != _num_partitions)
//          throw InternalError("Sorry. Hilmar and Dominik currently only support one base mesh cell per process.");
//
//        // cells
//        for(unsigned int icell = 0 ; icell < _num_cells; ++icell)
//        {
//          unsigned int pb, mb;
//          mesh_file->read(pb, mb);
//          _cells.at(icell)->set_parallel_block(pb-1);
//          _cells.at(icell)->set_matrix_block(mb-1);
//        }
//      } // _parse_partition

    public:

      /* ***************************
      * constructors & destructors *
      *****************************/
      /// DTOR
      ~FileParser<2, world_dim_>()
      {
        CONTEXT("BaseMesh::FileParser::~FileParser()");
        // loop over boundary components
        for(unsigned int ibc(0) ; ibc < _num_boundaries ; ++ibc)
        {
          delete [] _segment_type[ibc];

          for(unsigned int iseg(0) ; iseg < _num_segments[ibc] ; ++iseg)
          {
            delete [] _start_vertex[ibc][iseg];
            delete [] _vector_to_end_vertex[ibc][iseg];
            delete [] _param_int[ibc][iseg];
            if(_circle_section[ibc][iseg] != nullptr)
            {
              delete [] _circle_section[ibc][iseg];
            }
          }
          delete [] _param_int[ibc];
          delete [] _start_vertex[ibc];
          delete [] _vector_to_end_vertex[ibc];
          delete [] _circle_section[ibc];
        }
        delete [] _segment_type;
        delete [] _param_int;
        delete [] _start_vertex;
        delete [] _vector_to_end_vertex;
        delete [] _circle_section;
        delete [] _num_segments;
      }


      /* *****************
      * member functions *
      *******************/
      /**
      * \brief function for reading a 2D mesh in FEAST1 format
      *
      * \param[in] file_name
      * name of the mesh file
      *
      * \param[in,out] bm
      * base mesh to be built from the mesh file
      *
      * \author Dominik Goeddeke
      * \author Hilmar Wobker
      */
      inline void parse(String const& file_name, BM<2, world_dim_>* bm)
      {
        CONTEXT("BaseMesh::FileParser::parse()");
        // FEAST1 file format does not distinguish world dimension and space dimension. So, for now the world_dim_
        // template parameter actually has to be equal to 2 here.
        ASSERT(world_dim_ == 2, "Currently, world dim " + stringify(world_dim_) + " must be 2.");

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
        String saux;
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





    /**
    * \brief template class for base mesh file parser, specialisation for 3D files
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \note As long as no real file parser is implemented, we only manually create a 3D base mesh here.
    *
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class FileParser<3, world_dim_>
    {
      /// shortcut to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;

      /// shortcut to save typing of template parameters
      typedef Edge<3, world_dim_> Edge_;

      /// shortcut to save typing of template parameters
      typedef Quad<3, world_dim_> Quad_;

      /// shortcut to save typing of template parameters
      typedef Hexa<3, world_dim_> Hexa_;

    private:

      /// pointer to the base mesh
      BM<3, world_dim_>* _bm;

    public:

      /* *****************
      * member functions *
      *******************/
      /**
      * \brief function for reading a 3D mesh in FEAST1 format
      *
      * \param[in] file_name
      * name of the mesh file
      *
      * \param[in,out] bm
      * base mesh to be built from the mesh file
      *
      * \note As long as no real file parser is implemented, we only manually create a 3D base mesh here.
      *
      * \author Dominik Goeddeke
      * \author Hilmar Wobker
      */
      inline void parse(String const& file_name, BM<3, world_dim_>* bm)
      {
        CONTEXT("BaseMesh::FileParser::parse()");
        ASSERT(world_dim_ == 3, "Currently, world dim " + stringify(world_dim_) + " must be 3.");

        Logger::log("No file parser implemented yet! A manual 3D base mesh is created.", Logger::master);

        // TODO: remove the following line as soon as a real file parser is implemented
        Logger::log("Write file name " + file_name + " to prevent compiler warning about unused variable.");

        // set base mesh pointer
        _bm = bm;

        // shortcut for _subcells object
        Subcells<3, 3, world_dim_>* sc = &_bm->_subcells;

        /*
        * Test mesh consisting of 4 hexas:
        * View direction (0,0,-1) (from the front), (x,y) \in [0,2]x[0,2]:
        *   plane z=0 (back face):        plane z=1 (center):           plane z=2 (front):
        *
        * (0,2,0)             (2,2,0)   (0,2,1)             (2,2,1)   (0,2,2)             (2,2,2)
        *                                 13---26---14---27---15        05---08---06---09---07
        *                                  |         |         |         |         |         |
        *                                 23   13   24   14   25        05   01   06   02   07      y
        *                                  |         |         |         |         |         |      /\
        *             18---35---19        10---21---11---22---12        02---03---03---04---04      |
        *              |         |                   |         |                   |         |      |
        *             33   17   34                  19   12   20                  01   00   02      |
        *              |         |                   |         |                   |         |       --------> x
        *             16---32---17                  08---18---09                  00---00---01
        * (0,0,0)             (2,0,0)   (0,0,1)             (2,0,1)   (0,0,2)             (2,0,2)
        *
        * View direction (-1,0,0) (from the right), (y,z) \in [0,2]x[2,0]:
        *   plane x=0 (right):          plane x=1 (center):         plane x=2 (left):
        *
        * (0,2,2)             (0,2,0)   (1,2,2)             (1,2,0)   (2,2,2)             (2,2,0)
        *   05---15---13                  06---16---14                  07---17---15
        *    |         |                   |         |                   |         |
        *   05   07   23                  06   08   24                  07   09   25                          y
        *    |         |                   |         |                   |         |                          /\
        *   02---12---10                  03---13---11---30---18        04---14---12---31---19                |
        *                                  |         |         |         |         |         |                |
        *                                 01   03   19   16   33        02   04   20   17   34                |
        *                                  |         |         |         |         |         |     z <--------
        *                                 00---19---08---28---16        01---11---09---29---17
        * (0,0,2)             (0,0,0)   (1,0,2)             (1,0,0)   (2,0,2)             (2,0,0)
        *
        * View direction (0,-1,0) (from the top), (x,z) \in [0,2]x[2,0]:
        *   plane y=0 (bottom):         plane y=1 (center):         plane y=2 (top):
        *
        * (0,0,0)             (2,0,0)   (0,1,0)             (2,1,0)   (0,2,0)             (2,2,0)
        *             16---32---17                  18---35---19
        *              |         |                   |         |
        *             28   15   29                  30   18   31                                     --------> x
        *              |         |                   |         |                                    |
        *             08---18---09        10---21---11---22---12        13---26---14---27---15      |
        *              |         |         |         |         |         |         |         |      |
        *             10   02   11        12   05   13   06   14        15   10   16   11   17      \/
        *              |         |         |         |         |         |         |         |      z
        *             00---00---01        02---03---03---04---04        05---08---06---09---07
        * (0,0,2)             (2,0,2)   (0,1,2)             (2,1,2)   (0,2,2)             (2,2,2)
        *
        */
        // first, a bunch of vertices
        // v0 = (1,0,2)
        double coords[3] = {1,0,2};
        sc->_add(new Vertex_(coords));

        // v1 = (2,0,2)
        coords[0] = 2;
        sc->_add(new Vertex_(coords));

        // v2 = (0,1,2)
        coords[0] = 0;
        coords[1] = 1;
        sc->_add(new Vertex_(coords));

        // v3 = (1,1,2)
        coords[0] = 1;
        sc->_add(new Vertex_(coords));

        // v4 = (2,1,2)
        coords[0] = 2;
        sc->_add(new Vertex_(coords));

        // v5 = (0,2,2)
        coords[0] = 0;
        coords[1] = 2;
        sc->_add(new Vertex_(coords));

        // v6 = (1,2,2)
        coords[0] = 1;
        sc->_add(new Vertex_(coords));

        // v7 = (2,2,2)
        coords[0] = 2;
        sc->_add(new Vertex_(coords));

        // v8 = (1,0,1)
        coords[0] = 1;
        coords[1] = 0;
        coords[2] = 1;
        sc->_add(new Vertex_(coords));

        // v9 = (2,0,1)
        coords[0] = 2;
        sc->_add(new Vertex_(coords));

        // v10 = (0,1,1)
        coords[0] = 0;
        coords[1] = 1;
        sc->_add(new Vertex_(coords));

        // v11 = (1,1,1)
        coords[0] = 1;
        sc->_add(new Vertex_(coords));

        // v12 = (2,1,1)
        coords[0] = 2;
        sc->_add(new Vertex_(coords));

        // v13 = (0,2,1)
        coords[0] = 0;
        coords[1] = 2;
        sc->_add(new Vertex_(coords));

        // v14 = (1,2,1)
        coords[0] = 1;
        sc->_add(new Vertex_(coords));

        // v15 = (2,2,1)
        coords[0] = 2;
        sc->_add(new Vertex_(coords));

        // v16 = (1,0,0)
        coords[0] = 1;
        coords[1] = 0;
        coords[2] = 0;
        sc->_add(new Vertex_(coords));

        // v17 = (2,0,0)
        coords[0] = 2;
        sc->_add(new Vertex_(coords));

        // v18 = (1,1,0)
        coords[0] = 1;
        coords[1] = 1;
        sc->_add(new Vertex_(coords));

        // v19 = (2,1,0)
        coords[0] = 2;
        sc->_add(new Vertex_(coords));

        // then, a bunch of edges
        // e0
        sc->_add(new Edge_(sc->vertex(0), sc->vertex(1), 0));
        // e1
        sc->_add(new Edge_(sc->vertex(0), sc->vertex(3), 0));
        // e2
        sc->_add(new Edge_(sc->vertex(1), sc->vertex(4), 0));
        // e3
        sc->_add(new Edge_(sc->vertex(2), sc->vertex(3), 0));
        // e4
        sc->_add(new Edge_(sc->vertex(3), sc->vertex(4), 0));
        // e5
        sc->_add(new Edge_(sc->vertex(2), sc->vertex(5), 0));
        // e6
        sc->_add(new Edge_(sc->vertex(3), sc->vertex(6), 0));
        // e7
        sc->_add(new Edge_(sc->vertex(4), sc->vertex(7), 0));
        // e8
        sc->_add(new Edge_(sc->vertex(5), sc->vertex(6), 0));
        // e9
        sc->_add(new Edge_(sc->vertex(6), sc->vertex(7), 0));
        // e10
        sc->_add(new Edge_(sc->vertex(0), sc->vertex(8), 0));
        // e11
        sc->_add(new Edge_(sc->vertex(1), sc->vertex(9), 0));
        // e12
        sc->_add(new Edge_(sc->vertex(2), sc->vertex(10), 0));
        // e13
        sc->_add(new Edge_(sc->vertex(3), sc->vertex(11), 0));
        // e14
        sc->_add(new Edge_(sc->vertex(4), sc->vertex(12), 0));
        // e15
        sc->_add(new Edge_(sc->vertex(5), sc->vertex(13), 0));
        // e16
        sc->_add(new Edge_(sc->vertex(6), sc->vertex(14), 0));
        // e17
        sc->_add(new Edge_(sc->vertex(7), sc->vertex(15), 0));
        // e18
        sc->_add(new Edge_(sc->vertex(8), sc->vertex(9), 0));
        // e19
        sc->_add(new Edge_(sc->vertex(8), sc->vertex(11), 0));
        // e20
        sc->_add(new Edge_(sc->vertex(9), sc->vertex(12), 0));
        // e21
        sc->_add(new Edge_(sc->vertex(10), sc->vertex(11), 0));
        // e22
        sc->_add(new Edge_(sc->vertex(11), sc->vertex(12), 0));
        // e23
        sc->_add(new Edge_(sc->vertex(10), sc->vertex(13), 0));
        // e24
        sc->_add(new Edge_(sc->vertex(11), sc->vertex(14), 0));
        // e25
        sc->_add(new Edge_(sc->vertex(12), sc->vertex(15), 0));
        // e26
        sc->_add(new Edge_(sc->vertex(13), sc->vertex(14), 0));
        // e27
        sc->_add(new Edge_(sc->vertex(14), sc->vertex(15), 0));
        // e28
        sc->_add(new Edge_(sc->vertex(8), sc->vertex(16), 0));
        // e29
        sc->_add(new Edge_(sc->vertex(9), sc->vertex(17), 0));
        // e30
        sc->_add(new Edge_(sc->vertex(11), sc->vertex(18), 0));
        // e31
        sc->_add(new Edge_(sc->vertex(12), sc->vertex(19), 0));
        // e32
        sc->_add(new Edge_(sc->vertex(16), sc->vertex(17), 0));
        // e33
        sc->_add(new Edge_(sc->vertex(16), sc->vertex(18), 0));
        // e34
        sc->_add(new Edge_(sc->vertex(17), sc->vertex(19), 0));
        // e35
        sc->_add(new Edge_(sc->vertex(18), sc->vertex(19), 0));

        // now faces
        // f0
        //sc->_add(new Quad_(sc->vertex(0), sc->vertex(1), sc->vertex(3), sc->vertex(4),
        //                   sc->edge(0), sc->edge(4), sc->edge(1), sc->edge(2), 0));
        // deliberately, use different local numbering
        sc->_add(new Quad_(sc->vertex(1), sc->vertex(4), sc->vertex(0), sc->vertex(3),
                           sc->edge(2), sc->edge(1), sc->edge(0), sc->edge(4), 0));
        // f1
        sc->_add(new Quad_(sc->vertex(2), sc->vertex(3), sc->vertex(5), sc->vertex(6),
                           sc->edge(3), sc->edge(8), sc->edge(5), sc->edge(6), 0));
        // f2
        sc->_add(new Quad_(sc->vertex(3), sc->vertex(4), sc->vertex(6), sc->vertex(7),
                           sc->edge(4), sc->edge(9), sc->edge(6), sc->edge(7), 0));
        // f3
        //sc->_add(new Quad_(sc->vertex(0), sc->vertex(1), sc->vertex(8), sc->vertex(9),
        //                   sc->edge(0), sc->edge(18), sc->edge(10), sc->edge(11), 0));
        // deliberately, use different local numbering
        sc->_add(new Quad_(sc->vertex(1), sc->vertex(9), sc->vertex(0), sc->vertex(8),
                           sc->edge(11), sc->edge(10), sc->edge(0), sc->edge(18), 0));
        // f4
        sc->_add(new Quad_(sc->vertex(0), sc->vertex(8), sc->vertex(3), sc->vertex(11),
                           sc->edge(10), sc->edge(13), sc->edge(1), sc->edge(19), 0));
        // f5
        sc->_add(new Quad_(sc->vertex(1), sc->vertex(9), sc->vertex(4), sc->vertex(12),
                           sc->edge(11), sc->edge(14), sc->edge(2), sc->edge(20), 0));
        // f6
        sc->_add(new Quad_(sc->vertex(2), sc->vertex(3), sc->vertex(10), sc->vertex(11),
                           sc->edge(3), sc->edge(21), sc->edge(12), sc->edge(13), 0));
        // f7
        //sc->_add(new Quad_(sc->vertex(3), sc->vertex(4), sc->vertex(11), sc->vertex(12),
        //                   sc->edge(4), sc->edge(22), sc->edge(13), sc->edge(14), 0));
        // deliberately, use different local numbering
        sc->_add(new Quad_(sc->vertex(3), sc->vertex(11), sc->vertex(4), sc->vertex(12),
                           sc->edge(13), sc->edge(14), sc->edge(4), sc->edge(22), 0));
        // f8
        sc->_add(new Quad_(sc->vertex(2), sc->vertex(10), sc->vertex(5), sc->vertex(13),
                           sc->edge(12), sc->edge(15), sc->edge(5), sc->edge(23), 0));
        // f9
        sc->_add(new Quad_(sc->vertex(3), sc->vertex(11), sc->vertex(6), sc->vertex(14),
                           sc->edge(13), sc->edge(16), sc->edge(6), sc->edge(24), 0));
        // f10
        sc->_add(new Quad_(sc->vertex(4), sc->vertex(12), sc->vertex(7), sc->vertex(15),
                           sc->edge(14), sc->edge(17), sc->edge(7), sc->edge(25), 0));
        // f11
        sc->_add(new Quad_(sc->vertex(5), sc->vertex(6), sc->vertex(13), sc->vertex(14),
                           sc->edge(8), sc->edge(26), sc->edge(15), sc->edge(16), 0));
        // f12
        sc->_add(new Quad_(sc->vertex(6), sc->vertex(7), sc->vertex(14), sc->vertex(15),
                           sc->edge(9), sc->edge(27), sc->edge(16), sc->edge(17), 0));
        // f13
        sc->_add(new Quad_(sc->vertex(8), sc->vertex(9), sc->vertex(11), sc->vertex(12),
                           sc->edge(18), sc->edge(22), sc->edge(19), sc->edge(20), 0));
        // f14
        sc->_add(new Quad_(sc->vertex(10), sc->vertex(11), sc->vertex(13), sc->vertex(14),
                           sc->edge(21), sc->edge(26), sc->edge(23), sc->edge(24), 0));
        // f15
        sc->_add(new Quad_(sc->vertex(11), sc->vertex(12), sc->vertex(14), sc->vertex(15),
                           sc->edge(22), sc->edge(27), sc->edge(24), sc->edge(25), 0));
        // f16
        sc->_add(new Quad_(sc->vertex(8), sc->vertex(9), sc->vertex(16), sc->vertex(17),
                           sc->edge(18), sc->edge(32), sc->edge(28), sc->edge(29), 0));
        // f17
        sc->_add(new Quad_(sc->vertex(8), sc->vertex(16), sc->vertex(11), sc->vertex(18),
                           sc->edge(28), sc->edge(30), sc->edge(19), sc->edge(33), 0));
        // f18
        sc->_add(new Quad_(sc->vertex(9), sc->vertex(17), sc->vertex(12), sc->vertex(19),
                           sc->edge(29), sc->edge(31), sc->edge(20), sc->edge(34), 0));
        // f19
        sc->_add(new Quad_(sc->vertex(11), sc->vertex(12), sc->vertex(18), sc->vertex(19),
                           sc->edge(22), sc->edge(35), sc->edge(30), sc->edge(31), 0));
        // f20
        sc->_add(new Quad_(sc->vertex(16), sc->vertex(17), sc->vertex(18), sc->vertex(19),
                           sc->edge(32), sc->edge(35), sc->edge(33), sc->edge(34), 0));

        // finally, cells
        // c0
        _bm->_add(new Hexa_(sc->vertex(0), sc->vertex(1), sc->vertex(8), sc->vertex(9),
                            sc->vertex(3), sc->vertex(4), sc->vertex(11), sc->vertex(12),
                            sc->edge(0), sc->edge(18), sc->edge(4), sc->edge(22), sc->edge(10), sc->edge(11),
                            sc->edge(13), sc->edge(14), sc->edge(1), sc->edge(2), sc->edge(19), sc->edge(20),
                            sc->face(3), sc->face(7), sc->face(0), sc->face(13), sc->face(4), sc->face(5), 0));
        // c1
        //_bm->_add(new Hexa_(sc->vertex(2), sc->vertex(3), sc->vertex(10), sc->vertex(11),
        //                sc->vertex(5), sc->vertex(6), sc->vertex(13), sc->vertex(14),
        //                sc->edge(3), sc->edge(21), sc->edge(8), sc->edge(26), sc->edge(12), sc->edge(13),
        //                sc->edge(15), sc->edge(16), sc->edge(5), sc->edge(6), sc->edge(23), sc->edge(24),
        //                sc->face(6), sc->face(11), sc->face(1), sc->face(14), sc->face(8), sc->face(9), 0));
        // deliberately, use different local numbering
        _bm->_add(new Hexa_(sc->vertex(3), sc->vertex(6), sc->vertex(11), sc->vertex(14),
                            sc->vertex(2), sc->vertex(5), sc->vertex(10), sc->vertex(13),
                            sc->edge(6), sc->edge(24), sc->edge(5), sc->edge(23), sc->edge(13), sc->edge(16),
                            sc->edge(12), sc->edge(15), sc->edge(3), sc->edge(8), sc->edge(21), sc->edge(26),
                            sc->face(9), sc->face(8), sc->face(1), sc->face(14), sc->face(6), sc->face(11), 0));
        // c2
        _bm->_add(new Hexa_(sc->vertex(3), sc->vertex(4), sc->vertex(11), sc->vertex(12),
                            sc->vertex(6), sc->vertex(7), sc->vertex(14), sc->vertex(15),
                            sc->edge(4), sc->edge(22), sc->edge(9), sc->edge(27), sc->edge(13), sc->edge(14),
                            sc->edge(16), sc->edge(17), sc->edge(6), sc->edge(7), sc->edge(24), sc->edge(25),
                            sc->face(7), sc->face(12), sc->face(2), sc->face(15), sc->face(9), sc->face(10), 0));
        // c3
        //_bm->_add(new Hexa_(sc->vertex(8), sc->vertex(9), sc->vertex(16), sc->vertex(17),
        //                    sc->vertex(11), sc->vertex(12), sc->vertex(18), sc->vertex(19),
        //                    sc->edge(18), sc->edge(32), sc->edge(22), sc->edge(35), sc->edge(28), sc->edge(29),
        //                    sc->edge(30), sc->edge(31), sc->edge(19), sc->edge(20), sc->edge(33), sc->edge(34),
        //                    sc->face(16), sc->face(19), sc->face(13), sc->face(20), sc->face(17), sc->face(18), 0));
        // deliberately, use different local numbering
        _bm->_add(new Hexa_(sc->vertex(12), sc->vertex(11), sc->vertex(19), sc->vertex(18),
                            sc->vertex(9), sc->vertex(8), sc->vertex(17), sc->vertex(16),
                            sc->edge(22), sc->edge(35), sc->edge(18), sc->edge(32), sc->edge(31), sc->edge(30),
                            sc->edge(29), sc->edge(28), sc->edge(20), sc->edge(19), sc->edge(34), sc->edge(33),
                            sc->face(19), sc->face(16), sc->face(13), sc->face(20), sc->face(18), sc->face(17), 0));

        // neighbourhood

        // face neighbours
        _bm->cell(0)->add_neighbour(SDIM_FACE, 1, _bm->cell(2));
        _bm->cell(0)->add_neighbour(SDIM_FACE, 3, _bm->cell(3));
        _bm->cell(1)->add_neighbour(SDIM_FACE, 0, _bm->cell(2));
        _bm->cell(2)->add_neighbour(SDIM_FACE, 0, _bm->cell(0));
        _bm->cell(2)->add_neighbour(SDIM_FACE, 4, _bm->cell(1));
        _bm->cell(3)->add_neighbour(SDIM_FACE, 2, _bm->cell(0));

        // edge neighbours
        _bm->cell(0)->add_neighbour(SDIM_EDGE, 6, _bm->cell(1));
        _bm->cell(1)->add_neighbour(SDIM_EDGE, 4, _bm->cell(0));
        _bm->cell(2)->add_neighbour(SDIM_EDGE, 1, _bm->cell(3));
        _bm->cell(3)->add_neighbour(SDIM_EDGE, 0, _bm->cell(2));

        // vertex neighbours
        _bm->cell(1)->add_neighbour(SDIM_VERTEX, 2, _bm->cell(3));
        _bm->cell(3)->add_neighbour(SDIM_VERTEX, 1, _bm->cell(1));
      } // parse
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // KERNEL_BASE_MESH_FILE_PARSER_HPP

/*
Manually created 2D base mesh, maybe still useful. But before it can be used, it must be adapted.
      /// test CTOR, generates a hard-wired test mesh
      BaseMesh2D()
      {
        // Base mesh example consisting of three quads and two tris:
        //    v0---e0---v1---e1---v2 \.
        //    |          |         |    \.
        //   e2    c0   e3   c1   e4  c2  \ e5
        //    |          |         |         \.
        //    v3---e6---v4---e7---v5----e8---v6
        //                      /  |         |
        //                  e9/ c3 e10  c4  e11
        //                /        |         |
        //              v7---e12--v8---e13---v9

        // create the ten vertices
        // v0
        Vertex_* v = new Vertex_();
        v->set_coord(0, 0.0);
        v->set_coord(1, 1.0);
        _add(v);

        // v1
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, 1.0);
        _add(v);

        // v2
        v = new Vertex_();
        v->set_coord(0, 2.0);
        v->set_coord(1, 1.0);
        _add(v);

        // v3
        v = new Vertex_();
        v->set_coord(0, 0.0);
        v->set_coord(1, 0.0);
        _add(v);

        // v4
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, 0.0);
        _add(v);

        // v5
        v = new Vertex_();
        v->set_coord(0, 2.0);
        v->set_coord(1, 0.0);
        _add(v);

        // v6
        v = new Vertex_();
        v->set_coord(0, 3.0);
        v->set_coord(1, 0.0);
        _add(v);

        // v7
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, -1.0);
        _add(v);

        // v8
        v = new Vertex_();
        v->set_coord(0, 2.0);
        v->set_coord(1, -1.0);
        _add(v);

        // v9
        v = new Vertex_();
        v->set_coord(0, 3.0);
        v->set_coord(1, -1.0);
        _add(v);


        // create the 14 edges
        // (just to ease manual sanity checks, always use the vertex of smaller global index as start vertex)
        // e0
        Edge_* e = new Edge_(_vertices[0], _vertices[1], 0);
        _add(e);
        // e1
        e = new Edge_(_vertices[1], _vertices[2], 0);
        _add(e);
        // e2
        e = new Edge_(_vertices[0], _vertices[3], 0);
        _add(e);
        // e3
        e = new Edge_(_vertices[1], _vertices[4], 0);
        _add(e);
        // e4
        e = new Edge_(_vertices[2], _vertices[5], 0);
        _add(e);
        // e5
        e = new Edge_(_vertices[2], _vertices[6], 0);
        _add(e);
        // e6
        e = new Edge_(_vertices[3], _vertices[4], 0);
        _add(e);
        // e7
        e = new Edge_(_vertices[4], _vertices[5], 0);
        _add(e);
        // e8
        e = new Edge_(_vertices[5], _vertices[6], 0);
        _add(e);
        // e9
        e = new Edge_(_vertices[5], _vertices[7], 0);
        _add(e);
        // e10
        e = new Edge_(_vertices[5], _vertices[8], 0);
        _add(e);
        // e11
        e = new Edge_(_vertices[6], _vertices[9], 0);
        _add(e);
        // e12
        e = new Edge_(_vertices[7], _vertices[8], 0);
        _add(e);
        // e13
        e = new Edge_(_vertices[8], _vertices[9], 0);
        _add(e);

        // create quad cell c0
        Quad_* quad =
          new Quad_(_vertices[3], _vertices[4], _vertices[0], _vertices[1],
                    _edges[6], _edges[0], _edges[2], _edges[3], 0);
        _add(quad);

        // create quad cell c1
        quad = new Quad_(_vertices[4], _vertices[5], _vertices[1], _vertices[2],
                         _edges[7], _edges[1], _edges[3], _edges[4], 0);
        _add(quad);

        // create tri cell c2
        Tri_* tri = new Tri_(_vertices[5], _vertices[6], _vertices[2],
                             _edges[8], _edges[5], _edges[4], 0);
        _add(tri);
        // create tri cell c3
        tri = new Tri_(_vertices[7], _vertices[8], _vertices[5],
                       _edges[12], _edges[10], _edges[9], 0);
        _add(tri);

        // create quad cell c4
        quad = new Quad_(_vertices[8], _vertices[9], _vertices[5], _vertices[6],
                         _edges[13], _edges[8], _edges[10], _edges[11], 0);
        _add(quad);

        // set neighbourhood information (emulated file parser part 2)
        _cells[0]->add_neighbour(SDIM_EDGE, 3, _cells[1]);

        _cells[1]->add_neighbour(SDIM_EDGE, 2, _cells[0]);
        _cells[1]->add_neighbour(SDIM_EDGE, 3, _cells[2]);
        _cells[1]->add_neighbour(SDIM_VERTEX, 1, _cells[4]);
        _cells[1]->add_neighbour(SDIM_VERTEX, 1, _cells[3]);

        _cells[2]->add_neighbour(SDIM_EDGE, 2, _cells[1]);
        _cells[2]->add_neighbour(SDIM_EDGE, 0, _cells[4]);
        _cells[2]->add_neighbour(SDIM_VERTEX, 0, _cells[3]);

        _cells[3]->add_neighbour(SDIM_EDGE, 1, _cells[4]);
        _cells[3]->add_neighbour(SDIM_VERTEX, 2, _cells[2]);
        _cells[3]->add_neighbour(SDIM_VERTEX, 2, _cells[1]);

        _cells[4]->add_neighbour(SDIM_EDGE, 1, _cells[2]);
        _cells[4]->add_neighbour(SDIM_EDGE, 2, _cells[3]);
        _cells[4]->add_neighbour(SDIM_VERTEX, 2, _cells[1]);
      }
*/
