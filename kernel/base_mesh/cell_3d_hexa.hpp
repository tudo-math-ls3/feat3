#pragma once
#ifndef KERNEL_BASE_MESH_CELL_3D_HEXA_HPP
#define KERNEL_BASE_MESH_CELL_3D_HEXA_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <typeinfo>  // for typeid()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/base_mesh/vertex.hpp>
#include <kernel/base_mesh/cell.hpp>
#include <kernel/base_mesh/cell_data_validation.hpp>
#include <kernel/base_mesh/cell_1d_edge.hpp>
#include <kernel/base_mesh/cell_2d_quad.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    /**
    * \brief 3D base mesh cell of type hexa
    *
    * numbering scheme:
    *                                                                 back
    *                             6----03----7                        6----03----7                  6              7
    *                    top     /         /                          |          |       left    /  |    right  /  |
    * coord-system            06    1   07                            |          |            06    |        07    |
    *  z                   /         /                 front          10    3   11          /      10      /      11
    * /\      y          4----02----5                  4----02----5   |          |        4         |    5         |
    *  |    /                                          |          |   |          |        |    4    |    |    5    |
    *  | /                        2----01----3         |          |   2----01----3        |         2    |         3
    *  o------> x      bottom    /         /           08    2   09                       08     /       09     /
    *                         04    0   05             |          |                       |   04         |   05
    *                      /         /                 |          |                       | /            | /
    *                    0----00----1                  0----00----1                       0              1
    *
    * vertex coordinates of standard hexa [0,1]x[0,1]x[0,1]:
    *   v0: (0, 0, 0)
    *   v1: (1, 0, 0)
    *   v2: (0, 1, 0)
    *   v3: (1, 1, 0)
    *   v4: (0, 0, 1)
    *   v5: (1, 0, 1)
    *   v6: (0, 1, 1)
    *   v7: (1, 1, 1)
    *
    * edges:
    *   in x-direction    in y-direction    in z-direction
    *   e0: (v0, v1)      e4: (v0, v2)       e8: (v0, v4)
    *   e1: (v2, v3)      e5: (v1, v3)       e9: (v1, v5)
    *   e2: (v4, v5)      e6: (v4, v6)      e10: (v2, v6)
    *   e3: (v6, v7)      e7: (v5, v7)      e11: (v3, v7)
    *
    * faces:
    *   f0: (v0, v1, v2, v3),  (e0, e1, e4, e5)     (bottom, z=0)
    *   f1: (v4, v5, v6, v7),  (e2, e3, e6, e7)     (top, z=1)
    *   f2: (v0, v1, v4, v5),  (e0, e2,  e8,  e9)   (front, y=0)
    *   f3: (v2, v3, v6, v7),  (e1, e3, e10, e11)   (back, y=1)
    *   f4: (v0, v2, v4, v6),  (e4, e6, e8, e10)    (left, x=0)
    *   f5: (v1, v3, v5, v7),  (e5, e7, e9, e11)    (right, x=1)
    *
    * "Orientation in the hexa" means that edge and face vertices are traversed by increasing local vertex index.
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    */
// COMMENT_HILMAR: Um Code-Redundanz zu vermeiden, koennte wir ueberlegen, eine weitere Klasse Cell3D einzufuehren,
// die von Cell<3, space_dim_, world_dim_> erbt, und von der dann wieder um Tet und Hexa erben.
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Hexa
      : public Cell<3, space_dim_, world_dim_>
    {
      /// shortcuts to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;
      typedef Edge<space_dim_, world_dim_> Edge_;
      typedef Cell<1, space_dim_, world_dim_> Cell_1D_;
      typedef Quad<space_dim_, world_dim_> Quad_;
      typedef Cell<2, space_dim_, world_dim_> Cell_2D_;

    private:
      /// vertices of the hexa
      Vertex_* _vertices[8];

      /// edges of the hexa
      Cell_1D_* _edges[12];

      /// faces of the hexa
      Cell_2D_* _faces[6];

      /// relation between local quad numbering and face numbering within the hexa for each face (see struct Numbering)
      unsigned char _face_numbering[6];


      /// returns index (w.r.t. to hexa numbering) of the start vertex (iv=0) or the end vertex (iv=1) of edge iedge
      inline unsigned char _edge_vertex(unsigned char iedge, unsigned char iv)
      {
        CONTEXT("BaseMesh::Hexa::_edge_vertex()");
        assert(iedge < num_edges());
        assert(iv < 2);
        // the index is inquired from the fixed numbering scheme stored in Numbering::hexa_edge_vertices
        return Numbering::hexa_edge_vertices[iedge][iv];
      }


      /// returns true when the orientation of the edge coincides with its orientation within the hexa
      inline bool _edge_has_correct_orientation(unsigned char iedge)
      {
        CONTEXT("BaseMesh::Hexa::_edge_has_correct_orientation()");
        assert(iedge < num_edges());
        // return true when the edge's start vertex within the hexa is local vertex 0 within the edge structure
        return (vertex(_edge_vertex(iedge,0)) == edge(iedge)->vertex(0));
      }


      /**
      * \brief determines the relation between local face numbering and face numbering within the hexa
      *
      * This function determines the relation between local face numbering and face numbering within the hexa. It is
      * called from the constructor of the hexa. The relation between the numberings is stored for each face in the
      * array _face_numbering. For details see the description inside the struct Numbering.
      */
      inline void _determine_face_numbering()
      {
        CONTEXT("BaseMesh::Hexa::_determine_face_numbering()");
        try
        {
          for(unsigned char iface(0) ; iface < num_faces() ; ++iface)
          {
            // init the array entry
            _face_numbering[iface] = 42;
            for(unsigned char ivert(0) ; ivert < face(iface)->num_vertices() ; ++ivert)
            {
              // inquire whether the vertex with index ivert within the face structure is the first vertex of this face
              // within the hex structure
              if(face(iface)->vertex(ivert) == vertex(Numbering::hexa_face_vertices[iface][0]))
              {
                // If so then we either have relation ivert or ivert+4. Determine which of the two by comparing the
                // orientations. Within the face structure, we have to use the corresponding function
                // next_vertex_ccw(..), while in the hex structure we know that the vertex with index 1 is the ccw-next
                // of vertex with index 0.
                if(face(iface)->next_vertex_ccw(ivert) == vertex(Numbering::hexa_face_vertices[iface][1]))
                {
                  // same orientation
                  _face_numbering[iface] = ivert;
                }
                else if(face(iface)->previous_vertex_ccw(ivert) == vertex(Numbering::hexa_face_vertices[iface][1]))
                {
                  // opposite orientation
                  _face_numbering[iface] = ivert+4;
                }
                else
                {
                  throw InternalError("Something is wrong with the numbering of the "
                                      + StringUtils::stringify((int)iface) + "-th face with index "
                                      + face(iface)->print_index() + " in hexa " + this->print_index() + ".");
                }
              }
            } // for(unsigned char ivert(0) ; ivert < face(iface)->num_vertices() ; ++ivert)

            if (_face_numbering[iface] == 42)
            {
              throw InternalError("Vertex " + vertex(Numbering::hexa_face_vertices[iface][0])->print_index()
                                  + " not found in " + StringUtils::stringify((int)iface) + "-th face with index "
                                  + face(iface)->print_index() + " in hexa " + this->print_index() + ".");
            }
          } // for(unsigned char iface(0) ; iface < num_faces() ; ++iface)
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      }


    public:
      /// CTOR
      Hexa(
        Vertex_* v0, Vertex_* v1, Vertex_* v2, Vertex_* v3, Vertex_* v4, Vertex_* v5, Vertex_* v6, Vertex_* v7,
        Cell_1D_* e0, Cell_1D_* e1, Cell_1D_* e2, Cell_1D_* e3, Cell_1D_* e4, Cell_1D_* e5,
        Cell_1D_* e6, Cell_1D_* e7,Cell_1D_* e8, Cell_1D_* e9, Cell_1D_* e10, Cell_1D_* e11,
        Cell_2D_* f0, Cell_2D_* f1,Cell_2D_* f2, Cell_2D_* f3, Cell_2D_* f4, Cell_2D_* f5,
        unsigned char ref_level)
        : Cell<3, space_dim_, world_dim_>(ref_level)
      {
        CONTEXT("BaseMesh::Hexa::Hexa()");
        _vertices[0] = v0;
        _vertices[1] = v1;
        _vertices[2] = v2;
        _vertices[3] = v3;
        _vertices[4] = v4;
        _vertices[5] = v5;
        _vertices[6] = v6;
        _vertices[7] = v7;
        _edges[0] = e0;
        _edges[1] = e1;
        _edges[2] = e2;
        _edges[3] = e3;
        _edges[4] = e4;
        _edges[5] = e5;
        _edges[6] = e6;
        _edges[7] = e7;
        _edges[8] = e8;
        _edges[9] = e9;
        _edges[10] = e10;
        _edges[11] = e11;
        // assure that the edges are in fact of type Edge<space_dim_, world_dim_>, and not "only"
        // of type Cell<1, space_dim_, world_dim_>
        for(int i(0) ; i < 12 ; ++i)
        {
          assert(typeid(*_edges[i]) == typeid(Edge<space_dim_, world_dim_>));
        }
        _faces[0] = f0;
        _faces[1] = f1;
        _faces[2] = f2;
        _faces[3] = f3;
        _faces[4] = f4;
        _faces[5] = f5;
        // assure that the faces are in fact of type Quad_, and not "only" of type Cell_2D_
        for(int i(0) ; i < 6 ; ++i)
        {
          assert(typeid(*_faces[i]) == typeid(Quad_));
        }

        unsigned char num_subitems_per_subdim[3] = {8, 12, 6};
        this->_set_num_subitems_per_subdim(3, num_subitems_per_subdim);
        this->_init_neighbours();
// COMMENT_HILMAR: Eigentlich haette ich das lieber in die Konstruktoren-Liste gepackt, also sowas in der Art:
//    : CellData<3, space_dim_, world_dim_>({8,12,6})
// (was nicht kompiliert). Wie kann man denn on-the-fly ein Array anlegen und durchreichen?

        // determine and store the relation between local face numbering and face numbering within the hexa
        _determine_face_numbering();
      }


      /// returns number of vertices
      inline unsigned char num_vertices() const
      {
        CONTEXT("BaseMesh::Hexa::num_vertices()");
        return 8;
      }


      /// returns vertex at given index
      inline Vertex_* vertex(unsigned char const index) const
      {
        CONTEXT("BaseMesh::Hexa::vertex()");
        assert(index < num_vertices());
        return _vertices[index];
      }


      /// returns number of edges
      inline unsigned char num_edges() const
      {
        CONTEXT("BaseMesh::Hexa::num_edges()");
        return 12;
      }


      /// returns edge at given index
      inline Cell_1D_* edge(unsigned char const index) const
      {
        CONTEXT("BaseMesh::Hexa::edge()");
        assert(index < num_edges());
        return _edges[index];
      }


      /// returns number of faces
      inline unsigned char num_faces() const
      {
        CONTEXT("BaseMesh::Hexa::num_faces()");
        return 6;
      }


      /// returns face at given index
      inline Cell_2D_* face(unsigned char const index) const
      {
        CONTEXT("BaseMesh::Hexa::face()");
        assert(index < num_faces());
        return _faces[index];
      }


      /**
      * \brief subdivision routine splitting a hex and storing parent/child information
      *
      * \todo New centre vertex of the hexa has to be computed properly!
      */
      inline void subdivide(SubdivisionData<3, space_dim_, world_dim_>* subdiv_data)
      {
        CONTEXT("BaseMesh::Hexa::subdivide()");
        try
        {
          // assure that this cell has not been divided yet
          if(!this->active())
          {
            throw InternalError("Hexa " + this->print_index() + " is already subdivided!");
          }

          this->set_subdiv_data(subdiv_data);

          // clear all vectors of created entities in the SubdivisionData object
          this->subdiv_data()->clear_created();

          if(this->subdiv_data()->type == NONCONFORM_SAME_TYPE)
          {
            // Perform a nonconform subdivision without changing the cell type (1 hexa --> 8 hexas).

            /**
            * \brief vertices that this action creates and/or reuses (12 + 6 + 1)
            *
            * numbering scheme: new vertices
            *   - on existing edges: new_vertex_index = old_edge_index (positions 0-11 in array new_vertices[])
            *   - in face centres: new_vertex_index = 12 + old_face_index (positions 12-17)
            *   - in centre of the hexa: new_vertex_index = greatest new index (position 18)
            */
            Vertex_* new_vertices[19];

            /**
            * \brief edges that this action creates and/or reuses (12*2 + 6*4 + 6)
            *
            * numbering scheme: new edges
            *   - children of existing edges: new_edge_index = 2*(old_edge_index) + {0,1}
            *     (corresponding to edge vertex index within the hex numbering) (positions 0 - 23 in array new_edges[])
            *   - on faces: new_edge_index = 24 + 4*old_face_index + {0,1,2,3}
            *     (corresponding to face edge index within the hex numbering) (positions 24 - 47)
            *   - towards centre of the hexa: new_edge_index = 48 + old_face_index (positions 48-53)
            */
            Cell_1D_* new_edges[54];

            /**
            * \brief faces that this action creates and/or reuses (4*6 + 12)
            *
            * numbering scheme: new faces
            *   - children of existing faces: new_face_index = 4*(old_face_index) + {0,1,2,3}
            *     (corresponding to face vertex index within the hex numbering) (positions 0 - 23 in array new_faces[])
            *   - in the interior: new_face_index = 24 + old_edge_index (positions 24 - 35)
            *     (each new interior face bisects one old edge)
            */
            Cell_2D_* new_faces[36];

            // numbering of 8 new hexas: new_hexa_index = old_vertex_index
            // (each new hexa can be associated with one vertex of the old hexa)

            // loop over all faces and split them eventually, creating new vertices, edges and faces on the way
            for(unsigned char iface(0) ; iface < num_faces() ; ++iface)
            {
              // if edge has no children, create them
              if (face(iface)->active())
              {

                // subdivide face
                SubdivisionData<2, space_dim_, world_dim_>* subdiv_data
                  = new SubdivisionData<2, space_dim_, world_dim_>(NONCONFORM_SAME_TYPE);
                face(iface)->subdivide(subdiv_data);

                // add the created vertices, edges and faces to the vector of created vertices/edges/faces
                for(unsigned int i(0) ; i < face(iface)->subdiv_data()->created_vertices.size() ; ++i)
                {
                  this->subdiv_data()->created_vertices.push_back(face(iface)->subdiv_data()->created_vertices[i]);
                }

                for(unsigned int i(0) ; i < face(iface)->subdiv_data()->created_edges.size() ; ++i)
                {
                  this->subdiv_data()->created_edges.push_back(face(iface)->subdiv_data()->created_edges[i]);
                }

                for(unsigned int i(0) ; i < face(iface)->subdiv_data()->created_cells.size() ; ++i)
                {
                  this->subdiv_data()->created_faces.push_back(face(iface)->subdiv_data()->created_cells[i]);
                }
                // debug output
                Logger::log("Face " + StringUtils::stringify((int)iface) + ":"
                            + StringUtils::stringify(face(iface)->subdiv_data()->created_vertices.size())
                            + " vertices, " + StringUtils::stringify(face(iface)->subdiv_data()->created_edges.size())
                            + " edges, " + StringUtils::stringify(face(iface)->subdiv_data()->created_cells.size())
                            + " faces created.\n");
              }
            } // for(unsigned char iface(0) ; iface < num_faces() ; ++iface)


            // add vertices lying in the centres of the old edges to the array of new vertices (indices 0-11)
            // (they have already been pushed to the subdivision data structure)
            for (unsigned char iedge(0) ; iedge < num_edges() ; ++iedge)
            {
              // exploit that the vertex shared by the edge children is stored as second vertex within the structure of
              // both edge children
              new_vertices[iedge] = edge(iedge)->child(0)->vertex(1);
            }

            // add vertices lying in the centres of the faces to the array of new vertices (indices 12-17)
            // (they have already been pushed to the subdivision data structure)
            for (unsigned char iface(0) ; iface < num_faces() ; ++iface)
            {
              // exploit that the centre vertex of the face has index i in child i (no matter if it has been
              // created in this subdivision step or already by the neighbour hex) (see quad subdivision routine)
              new_vertices[12 + iface] = face(iface)->child(0)->vertex(0);
            }

            // create the centre vertex of the hexa
  // COMMENT_HILMAR: For the time being simply compute the centre vertex of the hexa as average of the eight corner
  // vertices until we find out, what is the best way of computing this point correctly.
            double p[world_dim_];
            for(unsigned char i(0) ; i < world_dim_ ; ++i)
            {
              p[i] = 0;
              for(int j(0) ; j < num_vertices()-1 ; ++j)
              {
                p[i] += vertex(j)->coord(i);
              }
              p[i] /= num_vertices();
            }
            new_vertices[18] = new Vertex<world_dim_>(p);
            this->subdiv_data()->created_vertices.push_back(new_vertices[18]);

            // add edges being children of the old edges to the array of new edges (indices 0-23)
            // (they have already been pushed to the subdivision data structure)
            for (unsigned char iedge(0) ; iedge < num_edges() ; ++iedge)
            {
              // inquire whether the internal edge orientation equals its orientation within the hexa
              if (_edge_has_correct_orientation(iedge))
              {
                new_edges[2*iedge] = edge(iedge)->child(0);
                new_edges[2*iedge+1] = edge(iedge)->child(1);
              }
              else
              {
                new_edges[2*iedge] = edge(iedge)->child(1);
                new_edges[2*iedge+1] = edge(iedge)->child(0);
              }
            }

            // add edges lying in the interior of the faces to the array of new edges (indices 24-47)
            // (they have already been pushed to the subdivision data structure)
            for (unsigned char iface(0) ; iface < num_faces() ; ++iface)
            {
              assert(face(iface)->num_children() == 4);
              for (unsigned char iedge(0) ; iedge < 4 ; ++iedge)
              {
                // Exploit that the edge from edge iedge of the parent face towards the centre vertex of the subdivided
                // face has index iedge in the two children adjacent to this edge (see quad subdivision routine).
                //
                // Example (see figure): The interior edge marked by "xx", which is connected to parent edge pe3, can be
                // accessed as local edge e3 of child c1 (v1 being the start vertex of pe3) or of child e3 (v3 being the
                // end vertex of edge pe3).
                // v2           pe1          v3
                //     ---------------------
                //     |   e3    |   e2    |
                //     |         |         |
                //     |e0 c2  e1|e1 c3  e0|
                //     |         |         |
                //     |   e2    |   e3    |
                // pe2 ----------|-----xx--- pe3
                //     |   e2    |   e3    |
                //     |         |         |
                //     |e1 c0  e0|e0 c1  e1|
                //     |         |         |
                //     |   e3    |   e2    |
                //     ---------------------
                // v0           pe0          v1

                // The face numbering has to be mapped to the local face structure via the array
                // Numbering::quad_to_quad_mappings_edges[][]. So, how to get the correct edge?
                // First, get the index of the current edge within the face structure
                unsigned char iedge_in_face = Numbering::quad_to_quad_mappings_edges[_face_numbering[iface]][iedge];
                // Then, get the index of the start vertex of this edge within the face structure (one could also take
                // the end vertex). This index is also the index of the child which we use to get the desired edge.
                unsigned char ivert_start = Numbering::quad_edge_vertices[iedge_in_face][0];
                // Now, the desired edge is the iedge_in_face-th edge of the ivert_start-th child of the quad.
                new_edges[24 + 4*iface + iedge] = face(iface)->child(ivert_start)->edge(iedge_in_face);
              }
            }

            // create edges from face centres to centre vertex of the hexa, add them to the array of new edges
            // (indices 48-53) and to the subdivision data structure
            for (unsigned char iface(0) ; iface < num_faces() ; ++iface)
            {
              new_edges[48 + iface] = new Edge_(new_vertices[12 + iface], new_vertices[18], 0);
              this->subdiv_data()->created_edges.push_back(new_edges[48 + iface]);
            }

            // add children of faces of the parent hexa to the array of new faces (indices 0-23)
            // (they have already been pushed to the subdivision data structure)
            for (unsigned char iface(0) ; iface < num_faces() ; ++iface)
            {
              assert(face(iface)->num_children() == 4);
              for (unsigned char ivert(0) ; ivert < 4 ; ++ivert)
              {
                new_faces[4*iface + ivert]
                  = face(iface)->child(Numbering::quad_to_quad_mappings_vertices[_face_numbering[iface]][ivert]);

              }
            }

            // create interior faces, add them to the array of new faces (indices 24-35) and to the subdivision data
            // structure (each new interior face can be associated with one edge of the parent hexa)
            // building rules:
            //   a) interior face in new_faces[24 + iedge] is associated with edge iedge of the parent hexa,
            //   b) interior numbering of the faces analogue to that of quad subdivision (as if the four quads building
            //      one "slice" were children of a non-existing parent quad)

            // four quads building the "slice" between the left and right face of the parent hexa
            new_faces[24] = new Quad_(new_vertices[18], new_vertices[12], new_vertices[14], new_vertices[0],
                                      new_edges[48], new_edges[32], new_edges[50], new_edges[24], 0);
            new_faces[25] = new Quad_(new_vertices[12], new_vertices[18], new_vertices[1], new_vertices[15],
                                      new_edges[48], new_edges[36], new_edges[25], new_edges[51], 0);
            new_faces[26] = new Quad_(new_vertices[14], new_vertices[2], new_vertices[18], new_vertices[13],
                                      new_edges[33], new_edges[49], new_edges[50], new_edges[28], 0);
            new_faces[27] = new Quad_(new_vertices[3], new_vertices[15], new_vertices[13], new_vertices[18],
                                      new_edges[37], new_edges[49], new_edges[29], new_edges[51], 0);
            // four quads building the "slice" between the front and back face of the parent hexa
            new_faces[28] = new Quad_(new_vertices[18], new_vertices[12], new_vertices[16], new_vertices[4],
                                      new_edges[48], new_edges[40], new_edges[52], new_edges[26], 0);
            new_faces[29] = new Quad_(new_vertices[12], new_vertices[18], new_vertices[5], new_vertices[17],
                                      new_edges[48], new_edges[44], new_edges[27], new_edges[53], 0);
            new_faces[30] = new Quad_(new_vertices[16], new_vertices[6], new_vertices[18], new_vertices[13],
                                      new_edges[41], new_edges[49], new_edges[52], new_edges[30], 0);
            new_faces[31] = new Quad_(new_vertices[7], new_vertices[17], new_vertices[13], new_vertices[18],
                                      new_edges[45], new_edges[49], new_edges[31], new_edges[53], 0);
            // four quads building the "slice" between the bottom and top face of the parent hexa
            new_faces[32] = new Quad_(new_vertices[18], new_vertices[14], new_vertices[16], new_vertices[8],
                                      new_edges[50], new_edges[42], new_edges[52], new_edges[34], 0);
            new_faces[33] = new Quad_(new_vertices[14], new_vertices[18], new_vertices[9], new_vertices[17],
                                      new_edges[50], new_edges[46], new_edges[35], new_edges[53], 0);
            new_faces[34] = new Quad_(new_vertices[16], new_vertices[10], new_vertices[18], new_vertices[15],
                                      new_edges[43], new_edges[51], new_edges[52], new_edges[38], 0);
            new_faces[35] = new Quad_(new_vertices[11], new_vertices[17], new_vertices[15], new_vertices[18],
                                      new_edges[47], new_edges[51], new_edges[39], new_edges[53], 0);
            for (unsigned char i(24) ; i < 24 + num_edges() ; ++i)
            {
              this->subdiv_data()->created_faces.push_back(new_faces[i]);
            }

            // set number of children to 8
            this->_set_num_children(8);

            // finally, create and add new hexas
            // building rules:
            //   - vertex i of the parent hexa is vertex 0 of the i-th child hexa
            //   - new_faces[i] is face 0 in i-th child hexa.
            //   - edge 0 of each child lies either on face 2 (front) or face 3 (back) of the parent hexa
            //  ==> face 0 of each child lies either on face 0 (bottom) or face 1 (top) of the parent hexa
            //  ==> faces 0/2/4 of each child lie on the boundary of the parent hexa (i.e. are children of parent hexa
            //      faces), faces 1/3/5 are interior faces
            //  ==> face 1 of child i is face 1 of child i+4, i=0,1,2,3 (between top and bottom children)
            //  ==> face 3 of child i is face 3 of child i+2, i=0,1,4,5 (between front and back children)
            //  ==> face 5 of child i is face 5 of child i+1, i=0,2,4,6 (between left and right children)

            // child 0 (front, bottom, left)
            _set_child(0, new Hexa(vertex(0), new_vertices[0], new_vertices[4], new_vertices[12],
                                   new_vertices[8], new_vertices[14], new_vertices[16], new_vertices[18],
                                   new_edges[0], new_edges[26], new_edges[34], new_edges[52],
                                   new_edges[8], new_edges[24], new_edges[42], new_edges[50],
                                   new_edges[16], new_edges[32], new_edges[40], new_edges[48],
                                   new_faces[0], new_faces[32], new_faces[8],
                                   new_faces[28], new_faces[16], new_faces[24],
                                   this->refinement_level()+1));
            // child 1 (front, bottom, right)
            _set_child(1, new Hexa(vertex(1), new_vertices[0], new_vertices[5], new_vertices[12],
                                   new_vertices[9], new_vertices[14], new_vertices[17], new_vertices[18],
                                   new_edges[1], new_edges[27], new_edges[35], new_edges[53],
                                   new_edges[10], new_edges[24], new_edges[46], new_edges[50],
                                   new_edges[18], new_edges[32], new_edges[44], new_edges[48],
                                   new_faces[1], new_faces[33], new_faces[9],
                                   new_faces[29], new_faces[20], new_faces[24],
                                   this->refinement_level()+1));
            // child 2 (back, bottom, left)
            _set_child(2, new Hexa(vertex(2), new_vertices[1], new_vertices[4], new_vertices[12],
                                   new_vertices[10], new_vertices[15], new_vertices[16], new_vertices[18],
                                   new_edges[2], new_edges[26], new_edges[38], new_edges[52],
                                   new_edges[9], new_edges[25], new_edges[43], new_edges[51],
                                   new_edges[20], new_edges[36], new_edges[40], new_edges[48],
                                   new_faces[2], new_faces[34], new_faces[12],
                                   new_faces[28], new_faces[17], new_faces[25],
                                   this->refinement_level()+1));
            // child 3 (back, bottom, right)
            _set_child(3, new Hexa(vertex(3), new_vertices[1], new_vertices[5], new_vertices[12],
                                   new_vertices[11], new_vertices[15], new_vertices[17], new_vertices[18],
                                   new_edges[3], new_edges[27], new_edges[39], new_edges[53],
                                   new_edges[11], new_edges[25], new_edges[47], new_edges[51],
                                   new_edges[22], new_edges[36], new_edges[44], new_edges[48],
                                   new_faces[3], new_faces[35], new_faces[13],
                                   new_faces[29], new_faces[21], new_faces[25],
                                   this->refinement_level()+1));
            // child 4 (front, top, left)
            _set_child(4, new Hexa(vertex(4), new_vertices[2], new_vertices[6], new_vertices[13],
                                   new_vertices[8], new_vertices[14], new_vertices[16], new_vertices[18],
                                   new_edges[4], new_edges[30], new_edges[34], new_edges[52],
                                   new_edges[12], new_edges[28], new_edges[42], new_edges[50],
                                   new_edges[17], new_edges[33], new_edges[41], new_edges[49],
                                   new_faces[4], new_faces[32], new_faces[10],
                                   new_faces[30], new_faces[18], new_faces[26],
                                   this->refinement_level()+1));
            // child 5 (front, top, right)
            _set_child(5, new Hexa(vertex(5), new_vertices[2], new_vertices[7], new_vertices[13],
                                   new_vertices[9], new_vertices[14], new_vertices[17], new_vertices[18],
                                   new_edges[5], new_edges[31], new_edges[35], new_edges[53],
                                   new_edges[14], new_edges[28], new_edges[46], new_edges[50],
                                   new_edges[19], new_edges[33], new_edges[45], new_edges[49],
                                   new_faces[5], new_faces[33], new_faces[11],
                                   new_faces[31], new_faces[22], new_faces[26],
                                   this->refinement_level()+1));
            // child 6 (back, top, left)
            _set_child(6, new Hexa(vertex(6), new_vertices[3], new_vertices[6], new_vertices[13],
                                   new_vertices[10], new_vertices[15], new_vertices[16], new_vertices[18],
                                   new_edges[6], new_edges[30], new_edges[38], new_edges[52],
                                   new_edges[13], new_edges[29], new_edges[43], new_edges[51],
                                   new_edges[21], new_edges[37], new_edges[41], new_edges[49],
                                   new_faces[6], new_faces[34], new_faces[14],
                                   new_faces[30], new_faces[19], new_faces[27],
                                   this->refinement_level()+1));
            // child 7 (back, top, right)
            _set_child(7, new Hexa(vertex(7), new_vertices[3], new_vertices[7], new_vertices[13],
                                   new_vertices[11], new_vertices[15], new_vertices[17], new_vertices[18],
                                   new_edges[7], new_edges[31], new_edges[39], new_edges[53],
                                   new_edges[15], new_edges[29], new_edges[47], new_edges[51],
                                   new_edges[23], new_edges[37], new_edges[45], new_edges[49],
                                   new_faces[7], new_faces[35], new_faces[15],
                                   new_faces[31], new_faces[23], new_faces[27],
                                   this->refinement_level()+1));

            // add the hexas to the vector of new created cells
            for (unsigned char i(0) ; i < this->num_children() ; ++i)
            {
              this->child(i)->set_parent(this);
              this->subdiv_data()->created_cells.push_back(this->child(i));
            }

            // set internal neighbourhood (external neighbourhood is set outside this function)
            // (in case space_dim_ > 3, an empty dummy function is called; see CellData)
            // face neighbours
            this->child(0)->add_neighbour(SDIM_FACE, 1, this->child(4));
            this->child(0)->add_neighbour(SDIM_FACE, 3, this->child(2));
            this->child(0)->add_neighbour(SDIM_FACE, 5, this->child(1));
            this->child(1)->add_neighbour(SDIM_FACE, 1, this->child(5));
            this->child(1)->add_neighbour(SDIM_FACE, 3, this->child(3));
            this->child(1)->add_neighbour(SDIM_FACE, 5, this->child(0));
            this->child(2)->add_neighbour(SDIM_FACE, 1, this->child(6));
            this->child(2)->add_neighbour(SDIM_FACE, 3, this->child(0));
            this->child(2)->add_neighbour(SDIM_FACE, 5, this->child(3));
            this->child(3)->add_neighbour(SDIM_FACE, 1, this->child(7));
            this->child(3)->add_neighbour(SDIM_FACE, 3, this->child(1));
            this->child(3)->add_neighbour(SDIM_FACE, 5, this->child(2));
            this->child(4)->add_neighbour(SDIM_FACE, 1, this->child(0));
            this->child(4)->add_neighbour(SDIM_FACE, 3, this->child(6));
            this->child(4)->add_neighbour(SDIM_FACE, 5, this->child(5));
            this->child(5)->add_neighbour(SDIM_FACE, 1, this->child(1));
            this->child(5)->add_neighbour(SDIM_FACE, 3, this->child(7));
            this->child(5)->add_neighbour(SDIM_FACE, 5, this->child(4));
            this->child(6)->add_neighbour(SDIM_FACE, 1, this->child(2));
            this->child(6)->add_neighbour(SDIM_FACE, 3, this->child(4));
            this->child(6)->add_neighbour(SDIM_FACE, 5, this->child(7));
            this->child(7)->add_neighbour(SDIM_FACE, 1, this->child(3));
            this->child(7)->add_neighbour(SDIM_FACE, 3, this->child(5));
            this->child(7)->add_neighbour(SDIM_FACE, 5, this->child(6));
            // edge neighbours
            this->child(0)->add_neighbour(SDIM_EDGE, 3, this->child(6));
            this->child(0)->add_neighbour(SDIM_EDGE, 7, this->child(5));
            this->child(0)->add_neighbour(SDIM_EDGE, 11, this->child(3));
            this->child(1)->add_neighbour(SDIM_EDGE, 3, this->child(7));
            this->child(1)->add_neighbour(SDIM_EDGE, 7, this->child(4));
            this->child(1)->add_neighbour(SDIM_EDGE, 11, this->child(2));
            this->child(2)->add_neighbour(SDIM_EDGE, 3, this->child(4));
            this->child(2)->add_neighbour(SDIM_EDGE, 7, this->child(7));
            this->child(2)->add_neighbour(SDIM_EDGE, 11, this->child(1));
            this->child(3)->add_neighbour(SDIM_EDGE, 3, this->child(5));
            this->child(3)->add_neighbour(SDIM_EDGE, 7, this->child(6));
            this->child(3)->add_neighbour(SDIM_EDGE, 11, this->child(0));
            this->child(4)->add_neighbour(SDIM_EDGE, 3, this->child(2));
            this->child(4)->add_neighbour(SDIM_EDGE, 7, this->child(1));
            this->child(4)->add_neighbour(SDIM_EDGE, 11, this->child(7));
            this->child(5)->add_neighbour(SDIM_EDGE, 3, this->child(3));
            this->child(5)->add_neighbour(SDIM_EDGE, 7, this->child(0));
            this->child(5)->add_neighbour(SDIM_EDGE, 11, this->child(6));
            this->child(6)->add_neighbour(SDIM_EDGE, 3, this->child(0));
            this->child(6)->add_neighbour(SDIM_EDGE, 7, this->child(3));
            this->child(6)->add_neighbour(SDIM_EDGE, 11, this->child(5));
            this->child(7)->add_neighbour(SDIM_EDGE, 3, this->child(1));
            this->child(7)->add_neighbour(SDIM_EDGE, 7, this->child(2));
            this->child(7)->add_neighbour(SDIM_EDGE, 11, this->child(4));
            // vertex neighbours
            this->child(0)->add_neighbour(SDIM_VERTEX, 7, this->child(7));
            this->child(1)->add_neighbour(SDIM_VERTEX, 7, this->child(6));
            this->child(2)->add_neighbour(SDIM_VERTEX, 7, this->child(5));
            this->child(3)->add_neighbour(SDIM_VERTEX, 7, this->child(4));
            this->child(4)->add_neighbour(SDIM_VERTEX, 7, this->child(3));
            this->child(5)->add_neighbour(SDIM_VERTEX, 7, this->child(2));
            this->child(6)->add_neighbour(SDIM_VERTEX, 7, this->child(1));
            this->child(7)->add_neighbour(SDIM_VERTEX, 7, this->child(0));
          }
          else
          {
            throw InternalError("Wrong type of subdivision in hexa " + this->print_index()
                                + ". Currently, only subdivision NONCONFORM_SAME_TYPE is supported.");
          }
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      } // subdivide()


      /**
      * \brief validate this cell
      *
      * \param[in] stream
      * stream validation info is written into
      */
      inline void validate(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Hexa::validate()");
        try
        {
          if(space_dim_ == 3)
          {
            stream << "Validating hexaeder " << this->print_index() << "..." << std::endl;
          }

          std::string s = "Hexaeder " + this->print_index() + ": ";

          // validate that all vertices, edges and faces are set
          for(unsigned char ivert(0) ; ivert < num_vertices() ; ++ivert)
          {
            if (vertex(ivert) == nullptr)
            {
              s += "Vertex " + StringUtils::stringify((int)ivert) + " is null.\n";
              throw new InternalError(s);
            }
          }
          for(unsigned char iedge(0) ; iedge < num_edges() ; ++iedge)
          {
            if (edge(iedge) == nullptr)
            {
              s += "Edge " + StringUtils::stringify((int)iedge) + " is null.\n";
              throw new InternalError(s);
            }
          }
          for(unsigned char iface(0) ; iface < num_faces() ; ++iface)
          {
            if (face(iface) == nullptr)
            {
              s += "Face " + StringUtils::stringify((int)iface) + " is null.\n";
              throw new InternalError(s);
            }
          }

          // validate children numbering
          if(!this->active())
          {
            if(this->subdiv_data()->type == NONCONFORM_SAME_TYPE)
            {
              // check whether the centre vertex of the hexaeder is vertex 7 in each child
              for(unsigned char ichild(0) ; ichild < this->num_children()-1 ; ++ichild)
              {
                if (this->child(ichild)->vertex(7) != this->child(ichild+1)->vertex(7))
                {
                  s += "Vertex 7 in child " + StringUtils::stringify((int)ichild) + " or child "
                       + StringUtils::stringify(ichild+1) + " is not the centre vertex of the parent hexa!\n";
                  throw new InternalError(s);
                }
              }

              // check whether vertex i of the parent hexa is vertex 0 of child i
              for(unsigned char ichild(0) ; ichild < this->num_children() ; ++ichild)
              {
                if (this->child(ichild)->vertex(0) != this->vertex(ichild))
                {
                  s += "Vertex 0 in child " + StringUtils::stringify((int)ichild) + " is not vertex "
                       + StringUtils::stringify(ichild+1) + " in the parent hexa!\n";
                  throw new InternalError(s);
                }
              }

              // check whether
              //   - face 0 of children 0-3 (4-7) is a child of face 0 (face 1) of the parent hexa
              //   - face 2 of children 0,1,4,5 (2,3,6,7) is a child of face 2 (face 3) of the parent hexa
              //   - face 4 of children 0,2,4,6 (1,3,5,7) is a child of face 4 (face 5) of the parent hexa
              for(unsigned char ifacelocal(0) ; ifacelocal <= 4 ; ifacelocal +=2)
              {
                for(unsigned char ichild(0) ; ichild < this->num_children() ; ++ichild)
                {
                  bool face_found = false;
                  unsigned char face_index;
                  switch(ifacelocal)
                  {
                    case(0):
                    {
                      // map 0-3 to 0 and 4-7 to 1
                      face_index = ichild/4;
                      break;
                    }
                    case(2):
                    {
                      // map 0,1,4,5 to 2 and 2,3,6,7 to 3
                      face_index = 2 + (ichild%4)/2;
                      break;
                    }
                    case(4):
                    {
                      // map 0,2,4,6 to 4 and 1,3,5,7 to 5
                      face_index = 4 + ichild%2;
                      break;
                    }
                  }
                  for(unsigned char ifacechild(0) ; ifacechild < this->face(face_index)->num_children() ; ++ifacechild)
                  {
                    if(this->child(ichild)->face(ifacelocal) == this->face(face_index)->child(ifacechild))
                    {
                      face_found = true;
                      break;
                    }
                  }
                  if(!face_found)
                  {
                    s += "Face " + StringUtils::stringify((int)ifacelocal) + " of child "
                         + StringUtils::stringify((int)ichild) + " is not a child of face "
                         + StringUtils::stringify((int)face_index) + " of the parent hexa!\n";
                    throw new InternalError(s);
                  }
                }
              }

              // check whether
              //   - face 1 of child i equals face 1 of child i+4, i=0,1,2,3 (between top and bottom children)
              //   - face 3 of child i equals face 3 of child i+2, i=0,1,4,5 (between front and back children)
              //   - face 5 of child i equals face 5 of child i+1, i=0,2,4,6 (between left and right children)
              unsigned char children[4] = {0,1,2,3};
              for(unsigned char ifacelocal(1) ; ifacelocal <= 5 ; ifacelocal +=2)
              {
                unsigned char ioffset;
                switch(ifacelocal)
                {
                  case(1):
                  {
                    ioffset = 4;
                    break;
                  }
                  case(3):
                  {
                    ioffset = 2;
                    children[1] = 1; children[2] = 4; children[3] = 5;
                    break;
                  }
                  case(5):
                  {
                    ioffset = 1;
                    children[1] = 2; children[3] = 6;
                    break;
                  }
                }

                for(unsigned char ichild(0) ; ichild < 4 ; ++ichild)
                {
                  if(   this->child(children[ichild])->face(ifacelocal)
                     != this->child(children[ichild] + ioffset)->face(ifacelocal))
                  {
                    s += "Face " + StringUtils::stringify((int)ifacelocal) + " of child "
                         + StringUtils::stringify((int)children[ichild]) + " is not equal to face "
                         + StringUtils::stringify((int)ifacelocal) + " of child "
                         + StringUtils::stringify(children[ichild] + ioffset) + "!\n";
                    throw new InternalError(s);
                  }
                }
              }
            }
          }

          // validate subitems (here: faces and edges)
          for(unsigned char iface(0) ; iface < num_faces() ; ++iface)
          {
            face(iface)->validate(stream);
          }
          // validate subitems (here: faces and edges)
          for(unsigned char iedge(0) ; iedge < num_edges() ; ++iedge)
          {
            edge(iedge)->validate(stream);
          }

          // validate parent-child relations
          this->validate_history(stream);

          // validate neighbours
          if (this->active())
          {
            CellDataValidation<3, space_dim_, world_dim_>::validate_neighbourhood(this, stream);
          }
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      }


      /// print information about this quad
      inline void print(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Hexa::print()");
        stream << "Hexa";
        this->print_index(stream);

        stream << ": [V ";
        _vertices[0]->print_index(stream);
        for(int i(1) ; i < num_vertices() ; ++i)
        {
          stream << ", ";
          _vertices[i]->print_index(stream);
        }
        stream << "]"<< std::endl << "    [E ";
        _edges[0]->print_index(stream);
        if(_edge_has_correct_orientation(0))
        {
          stream << "(+)";
        }
        else
        {
          stream << "(-)";
        }
        for(int i(1) ; i < num_edges() ; ++i)
        {
          stream << ", ";
          _edges[i]->print_index(stream);
          if(_edge_has_correct_orientation(i))
          {
            stream << "(+)";
          }
          else
          {
            stream << "(-)";
          }
        }
        stream << "]"<< std::endl << "    [F ";
        _faces[0]->print_index(stream);
        stream << "(" << (int)_face_numbering[0] << ")";
        for(int i(1) ; i < num_faces() ; ++i)
        {
          stream << ", ";
          _faces[i]->print_index(stream);
          stream << "(" << (int)_face_numbering[i] << ")";
        }
        stream << "]" << std::endl << "    ";
        Cell<3, space_dim_, world_dim_>::print_history(stream);
        // print neighbourhood information (if there is any)
        CellData<3, space_dim_, world_dim_>::print(stream);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_3D_HEXA_HPP
