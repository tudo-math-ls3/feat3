#pragma once
#ifndef KERNEL_BASE_MESH_CELL_3D_HEXA_HPP
#define KERNEL_BASE_MESH_CELL_3D_HEXA_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <typeinfo>  // for typeid()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_vertex.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_cell_1d_edge.hpp>
#include <kernel/base_mesh_cell_2d_quad.hpp>

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
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    * \author Peter Zajac
    */
// COMMENT_HILMAR: Um Code-Redundanz zu vermeiden, koennte wir ueberlegen, eine weitere Klasse Cell3D einzufuehren,
// die von Cell<3, space_dim_, world_dim_> erbt, und von der dann wieder um Tet und Hexa erben.
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Hexa
      : public Cell<3, space_dim_, world_dim_>
    {
      /// shortcut for type Vertex<world_dim_>
      typedef Vertex<world_dim_> Vertex_;

      /// shortcut for type Cell<1, space_dim_, world_dim_>
      typedef Edge<space_dim_, world_dim_> Edge_;

      /// shortcut for type Cell<1, space_dim_, world_dim_>
      typedef Cell<1, space_dim_, world_dim_> Cell_1D_;

      /// shortcut for type Quad<space_dim_, world_dim_>
      typedef Quad<space_dim_, world_dim_> Quad_;

      /// shortcut for type Cell<2, space_dim_, world_dim_>
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
        assert(iedge < num_edges());
        assert(iv < 2);
        // the index is inquired from the fixed numbering scheme stored in Numbering::hexa_edge_vertices
        return Numbering::hexa_edge_vertices[iedge][iv];
      }

      /// returns true when the orientation of the edge coincides with its orientation within the hexa
      inline bool _edge_has_correct_orientation(unsigned char iedge)
      {
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
              // orientations. Within the face structure, we have to use the corresponding function next_vertex_ccw(..),
              // while in the hex structure we know that the vertex with index 1 is the ccw-next of vertex with index 0.
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
                std::cerr << "Something is wrong with the numbering of face " << (int)iface << " in hexa." << std::endl;
                exit(1);
              }
            }
          } // for(unsigned char ivert(0) ; ivert < face(iface)->num_vertices() ; ++ivert)

          if (_face_numbering[iface] == 42)
          {
            std::cerr << "Vertex not found in face " << (int)iface << " in hexa." << std::endl;
            exit(1);
          }
          //std::cout << "face " << (int)iface << "(" << (int)face(iface)->index() << "), numb "
          //          << (int)_face_numbering[iface] << std::endl;
        } // for(unsigned char iface(0) ; iface < num_faces() ; ++iface)
      }


    public:
      /// CTOR
      Hexa(Vertex_* v0, Vertex_* v1, Vertex_* v2, Vertex_* v3, Vertex_* v4, Vertex_* v5, Vertex_* v6, Vertex_* v7,
           Cell_1D_* e0, Cell_1D_* e1, Cell_1D_* e2, Cell_1D_* e3, Cell_1D_* e4, Cell_1D_* e5,
           Cell_1D_* e6, Cell_1D_* e7,Cell_1D_* e8, Cell_1D_* e9, Cell_1D_* e10, Cell_1D_* e11,
           Cell_2D_* f0, Cell_2D_* f1,Cell_2D_* f2, Cell_2D_* f3, Cell_2D_* f4, Cell_2D_* f5)
      {
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

        unsigned char num_subcells_per_subdimension[3] = {8, 12, 6};
        this->_init_neighbours(3, num_subcells_per_subdimension);
// COMMENT_HILMAR: Eigentlich haette ich das lieber in die Konstruktoren-Liste gepackt, also sowas in der Art:
//    : CellData<3, space_dim_, world_dim_>({8,12,6})
// (was nicht kompiliert). Wie kann man denn on-the-fly ein Array anlegen und durchreichen?

        // determine and store the relation between local face numbering and face numbering within the hexa
        _determine_face_numbering();
      }


      /// returns number of vertices
      inline unsigned char num_vertices() const
      {
        return 8;
      }


      /// returns vertex at given index
      inline Vertex_* vertex(unsigned char const index) const
      {
        assert(index < num_vertices());
        return _vertices[index];
      }


      /// returns number of edges
      inline unsigned char num_edges() const
      {
        return 12;
      }


      /// returns edge at given index
      inline Cell_1D_* edge(unsigned char const index) const
      {
        assert(index < num_edges());
        return _edges[index];
      }


      /// returns number of faces
      inline unsigned char num_faces() const
      {
        return 6;
      }


      /// returns face at given index
      inline Cell_2D_* face(unsigned char const index) const
      {
        assert(index < num_faces());
        return _faces[index];
      }


      /// subdivision routine splitting a hex and storing parent/child information
// COMMENT_HILMAR: this is currently hard-wired to splitting the hexa into eight hexas. Later, this is parameterised
// via the information in the SubdivisionData object.
      inline void subdivide(SubdivisionData<3, space_dim_, world_dim_>& subdiv_data)
      {
        // assure that this cell has not been divided yet
        if(!this->active())
        {
          std::cerr << "Hexa " << this->index() << " is already subdivided! Aborting program.";
          exit(1);
        }

        // clear all vectors of created entities in the SubdivisionData object
        subdiv_data.clear_created();

        /**
        * \brief vertices that this action creates and/or reuses (12 + 6 + 1)
        *
        * numbering scheme: new vertices
        *   - on existing edges: new_vertex_index = old_edge_index (indices 0-11)
        *   - in face centres: new_vertex_index = 12 + old_face_index (indices 12-17)
        *   - in centre of the hexa: new_vertex_index = greatest new index (index 18)
        */
        Vertex_* new_vertices[19];

        /**
        * \brief edges that this action creates and/or reuses (12*2 + 6*4 + 6)
        *
        * numbering scheme: new edges
        *   - as children of existing edges: new_edge_index = 2*(old_edge_index) + {0,1}
        *     (increasing with old_vertex_index) (indices 0 - 23)
        *   - on faces: new_edge_index = 24 + 4*old_face_index + {0,1,2,3} (indices 24 - 47)
        *   - towards centre of the hexa: new_edge_index = 48 + old_face_index (indices 48-53)
        */
        Cell_1D_* new_edges[54];

        /**
        * \brief faces that this action creates and/or reuses (4*6 + 12)
        *
        * numbering scheme: new faces
        *   - as children of existing faces: new_face_index = 4*(old_face_index) + {0,1,2,3}
        *     (increasing with old_vertex_index) (indices 0 - 23)
        *   - in the interior: new_face_index = 24 + old_edge_index (indices 24 - 35)
        *     (each new interior face bisects one old edge)
        */
        Cell_2D_* new_faces[36];

        // numbering of 8 new hexas: new_hexa_index = old_vertex_index
        // (each new hexa can be associated with one vertex of the old hexa)

//COMMENT_HILMAR: brauchen wir das ueberhaupt noch?
//        // store old active-mask of each edge, because it gets overwritten once edges are subdivided
//        bool old_edge_active_mask[12];
//        // store old active-mask of each face, because it gets overwritten once faces are subdivided
//        bool old_face_active_mask[6];

        SubdivisionData<2, space_dim_, world_dim_> subdiv_data_face;

        // loop over all faces and split them eventually, creating new vertices, edges and faces on the way
        for(unsigned char iface(0) ; iface < num_faces() ; ++iface)
        {
          // if edge has no children, create them
          if (face(iface)->active())
          {
//            // store old active mask
//            old_face_active_mask[iface] = true;

            // subdivide face
            face(iface)->subdivide(subdiv_data_face);

            // add the created vertices, edges and faces to the vector of created vertices/edges/faces
// COMMENT_HILMAR: temporarily ignoring orientation here!
            for(unsigned int i(0) ; i < subdiv_data_face.created_vertices.size() ; ++i)
            {
              subdiv_data.created_vertices.push_back(subdiv_data_face.created_vertices[i]);
            }

            for(unsigned int i(0) ; i < subdiv_data_face.created_edges.size() ; ++i)
            {
              subdiv_data.created_edges.push_back(subdiv_data_face.created_edges[i]);
            }

            for(unsigned int i(0) ; i < subdiv_data_face.created_cells.size() ; ++i)
            {
              subdiv_data.created_faces.push_back(subdiv_data_face.created_cells[i]);
            }
std::cout << "face " << (int)iface << ":" << std::endl;
std::cout << subdiv_data_face.created_vertices.size() << " vertices";
std::cout << ", " << subdiv_data_face.created_edges.size() << " edges";
std::cout << ", " << subdiv_data_face.created_cells.size() << " faces created." << std::endl;
          }
          else // edge has children, reuse them
          {
//            // store old active mask
//            old_face_active_mask[iface] = false;
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
          // exploit that the centre vertex of the face has index 3 in all its children (no matter if it has been
          // created in this subdivision step or already by the neighbour hex) (see quad subdivision routine)
          new_vertices[12 + iface] = face(iface)->child(0)->vertex(3);
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
        subdiv_data.created_vertices.push_back(new_vertices[18]);

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
          for (unsigned char iedge(0) ; iedge < 4 ; ++iedge)
          {
            // Exploit that the edge from edge iedge of the parent face towards the centre vertex of the subdivided
            // face has index 3 in child iedge (see quad subdivision routine). The face numbering has to be mapped
            // to the local face structure via the array Numbering::quad_to_quad_mappings_edges[][].
            new_edges[24 + 4*iface + iedge]
              = face(iface)->child(Numbering::quad_to_quad_mappings_edges[_face_numbering[iface]][iedge])->edge(3);
//std::cout << "index = " << 24 + 4*iface + iedge << ", edge hex = " << (int) iedge << ", edge local = "
//  << (int)Numbering::quad_to_quad_mappings_edges[_face_numbering[iface]][iedge] << std::endl;
          }
        }

        // create edges from face centres to centre vertex of the hexa, add them to the array of new edges
        // (indices 48-53) and to the subdivision data structure
        for (unsigned char iface(0) ; iface < num_faces() ; ++iface)
        {
          new_edges[48 + iface] = new Edge_(new_vertices[12 + iface], new_vertices[18]);
          subdiv_data.created_edges.push_back(new_edges[48 + iface]);
        }


        // now faces
        // ...


// COMMENT_HILMAR: code below not adapted yet!!
//        // set number of children to 8
//        this->_set_num_children(8);
//
//
//        // finally, create and add new quads
//        _set_child(0, new Quad_(vertex(0), new_vertices[0], new_vertices[4], new_vertices[3],
//                               new_edges[0], new_edges[8], new_edges[11], new_edges[7]));
//        _set_child(1, new Quad_(new_vertices[0], vertex(1), new_vertices[1], new_vertices[4],
//                               new_edges[1], new_edges[2], new_edges[9], new_edges[8]));
//        _set_child(2, new Quad_(new_vertices[4], new_vertices[1], vertex(2), new_vertices[2],
//                               new_edges[9], new_edges[3], new_edges[4], new_edges[10]));
//        _set_child(3, new Quad_(new_vertices[3], new_vertices[4], new_vertices[2], vertex(3),
//                               new_edges[11], new_edges[10], new_edges[5], new_edges[6]));
//
//        for (unsigned char i(0) ; i < 4 ; ++i)
//        {
//          this->child(i)->set_parent(this);
//          subdiv_data.created_cells.push_back(this->child(i));
//        }
      } // subdivide()


      /// print information about this quad
      inline void print(std::ostream& stream)
      {
        stream << "Hexa";
        Item::print(stream);
        stream << ": [";

        for(int i(0) ; i < num_faces() ; ++i)
        {
          stream << "F" << _faces[i]->index() << "(" << _face_numbering[i] << ")";
          if(i < num_faces()-1)
          {
            stream << ", ";
          }
          else
          {
            stream << "]";
          }
        }
        Cell<3, space_dim_, world_dim_>::print_history(stream);
        // print neighbourhood information (if there is any)
        CellData<3, space_dim_, world_dim_>::print(stream);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_3D_HEXA_HPP
