#pragma once
#ifndef KERNEL_BASE_MESH_CELL_1D_EDGE_HPP
#define KERNEL_BASE_MESH_CELL_1D_EDGE_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_vertex.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    /**
    * \brief 1D base mesh cell of type edge
    *
    * \note
    * Be aware of the difference between Cell<1, space_dim_, world_dim_> and Edge<space_dim_, world_dim_>. It is
    * analogous to the difference between Cell<2, space_dim_, world_dim_> and Quad<space_dim_, world_dim_>
    * (Tri<space_dim_, world_dim_>, resp.). The base mesh only stores the general type! Actually, the distinction
    * between Cell<1, space_dim_, world_dim_> and Edge<space_dim_, world_dim_> is not really necessary, since there is
    * only one edge type. However, when one wants to avoid this class Edge, then one would have to write a template
    * specialisation of the class Cell<cell_dim_, ...> for cell_dim_ = 1.
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    * \author Peter Zajac
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Edge
      : public Cell<1, space_dim_, world_dim_>
    {
      /// shortcut for type Vertex<world_dim_>
      typedef Vertex<world_dim_> Vertex_;


    private:
      /// vertices of the edge
      Vertex_* _vertices[2];

    public:
      /// CTOR
      Edge(Vertex_* v0, Vertex_* v1)
      {
        _vertices[0] = v0;
        _vertices[1] = v1;

        unsigned char num_subcells_per_subdimension[1] = {2};
        this->_init_neighbours(1, num_subcells_per_subdimension);
      }


      /// returns number of vertices
      inline unsigned char num_vertices() const
      {
        return 2;
      }


      /// returns vertex at given index
      inline Vertex_* vertex(unsigned char const index) const
      {
        assert(index < num_vertices());
        return _vertices[index];
      }


      /// subdivision routine splitting an edge into two and storing parent/child information
      inline void subdivide(SubdivisionData<1, space_dim_, world_dim_>& subdiv_data)
      {
        // assure that this cell has not been divided yet
        if(!this->active())
        {
          std::cerr << "Edge " << this->index() << " is already subdivided! Aborting program.";
          exit(1);
        }

        // clear all vectors of created entities in the SubdivisionData object
        subdiv_data.clear_created();

        // create new vertex as mid point of this edge
        double p[world_dim_];
        for(int i(0) ; i < world_dim_ ; ++i)
        {
          p[i] = vertex(0)->coords(i) + 0.5*(vertex(1)->coords(i) - vertex(0)->coords(i) );
        }
        // create new vertex
        subdiv_data.created_vertex = new Vertex_(p);

        // create new edges (note the numbering of the vertices) and set them as children of this edge
        this->_set_num_children(2);
        _set_child(0, new Edge(vertex(0), subdiv_data.created_vertex));
        _set_child(1, new Edge(vertex(1), subdiv_data.created_vertex));
        // update the parent relationship
        this->child(0)->set_parent(this);
        this->child(1)->set_parent(this);

        // add new edges to the vector of created cells
        subdiv_data.created_cells.push_back(this->child(0));
        subdiv_data.created_cells.push_back(this->child(1));
      }

      /// print information about this edge
      inline void print(std::ostream& stream)
      {
        stream << "E";
        Item::print(stream);
        stream << ": [";
        _vertices[0]->print(stream);
        stream << ", ";
        _vertices[1]->print(stream);
        stream << "]";
        Cell<1, space_dim_, world_dim_>::print_history(stream);
        // print neighbourhood information (if there is any)
        CellData<1, space_dim_, world_dim_>::print(stream);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_1D_EDGE_HPP
