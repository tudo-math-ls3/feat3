#pragma once
#ifndef KERNEL_BASE_MESH_CELL_2D_TRI_HPP
#define KERNEL_BASE_MESH_CELL_2D_TRI_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_vertex.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_cell_1d_edge.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    /**
    * \brief 2D base mesh cell of type triangle
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    * \author Peter Zajac
    */
// COMMENT_HILMAR: Um Code-Redundanz zu vermeiden, koennte wir ueberlegen, eine weitere Klasse Cell2D einzufuehren,
// die von Cell<2, space_dim_, world_dim_> erbt, und von der dann wieder um Quad und Tri erben. Darin koennte
// man zum Beispiel die Funktion _edge_has_correct_orientation() implementieren.
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Tri
      : public Cell<2, space_dim_, world_dim_>
    {
      /// shortcuts for type Vertex<world_dim_>
      typedef Vertex<world_dim_> Vertex_;

      /// shortcuts for type Edge<space_dim_, world_dim_>
      typedef Edge<space_dim_, world_dim_> Edge_;

    private:
      /// vertices of the triangle
      Vertex_* _vertices[3];

      /// edges of the triangle
      Edge_* _edges[3];


      /// returns true when edge with local index iedge has the same orientation as the quad
      inline bool _edge_has_correct_orientation(unsigned char iedge)
      {
        // the orientation of the edge is correct (i.e. the same as that of the quad), when its start vertex within
        // the quad (i.e. the vertex with the same local index) is local vertex 0 within the edge structure
        return (vertex(iedge) == edge(iedge)->vertex(0));
      }


    public:
      /// CTOR
      Tri(Vertex_* v0, Vertex_* v1, Vertex_* v2, Edge_* e0, Edge_* e1, Edge_* e2)
      {
        _vertices[0] = v0;
        _vertices[1] = v1;
        _vertices[2] = v2;
        _edges[0] = e0;
        _edges[1] = e1;
        _edges[2] = e2;

        unsigned char temp[2] = {3,3};
        this->_init_neighbours(2, temp);
// TODO: Eigentlich haette ich das lieber in die Konstruktoren-Liste gepackt, also sowas in der Art:
//    : CellData<2, space_dim_, world_dim_>({4,4})
// (was nicht kompiliert). Wie kann man denn on-the-fly ein Array anlegen und durchreichen?
      }


      /// returns number of vertices
      inline unsigned char num_vertices() const
      {
        return 3;
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
        return 3;
      }


      /// returns edge at given index
      inline Edge_* edge(unsigned char const index) const
      {
        assert(index < num_edges());
        return _edges[index];
      }


      /// subdivision routine splitting a tri and storing parent/child information
      inline void subdivide(SubdivisionData<2, space_dim_, world_dim_>& subdiv_data)
      {
        // assure that this cell has not been divided yet
        if(!this->active())
        {
          std::cerr << "Tri " << this->index() << " is already subdivided! Aborting program.";
          exit(1);
        }

        // clear all vectors of created entities in the SubdivisionData object
        subdiv_data.clear_created();

        std::cerr << "Subdivision for triangles not implemented yet!!" << std::endl;
        // TODO: subdivion code
      }


      /// print information about this tri
      inline void print(std::ostream& stream)
      {
        stream << "Tri";
        Item::print(stream);
        stream << ": [";

        for(int i(0) ; i < 3 ; ++i)
        {
          stream << "E" << _edges[i]->index();
          if(_edge_has_correct_orientation(i))
          {
            stream << "(+)";
          }
          else
          {
            stream << "(-)";
          }
          if(i<3)
          {
            stream << ", ";
          }
          else
          {
            stream << "]";
          }
        }
        Cell<2, space_dim_, world_dim_>::print_history(stream);
        CellData<2, space_dim_, world_dim_>::print(stream);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_2D_TRI_HPP
