#pragma once
#ifndef KERNEL_BASE_MESH_CELL_1D_EDGE_HPP
#define KERNEL_BASE_MESH_CELL_1D_EDGE_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_cell_data_checker.hpp>
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
      Edge(
        Vertex_* v0, Vertex_* v1,
        unsigned char ref_level)
        : Cell<1, space_dim_, world_dim_>(ref_level)
      {
        _vertices[0] = v0;
        _vertices[1] = v1;

        unsigned char num_subitems_per_subdim[1] = {2};
        this->_set_num_subitems_per_subdim(1, num_subitems_per_subdim);
        this->_init_neighbours();
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
      inline void subdivide()
      {
        // assure that this cell has not been divided yet
        if(!this->active())
        {
          std::cerr << "Edge " << this->print_index() << " is already subdivided! Aborting program." << std::endl;
          exit(1);
        }

        if(!this->subdiv_data_initialised())
        {
          std::cerr << "Edge " << this->print_index()
                    << " cannot be subdivided! Set subdivision data first! Aborting program." << std::endl;
          exit(1);
        }

        if(this->subdiv_data()->type == CONFORM_SAME_TYPE)
        {
          // split the edge into two edges

          // clear all vectors of created entities in the SubdivisionData object
          this->subdiv_data()->clear_created();

          // create new vertex as mid point of this edge
          double p[world_dim_];
          for(int i(0) ; i < world_dim_ ; ++i)
          {
            p[i] = vertex(0)->coord(i) + 0.5*(vertex(1)->coord(i) - vertex(0)->coord(i) );
          }
          // create new vertex
          this->subdiv_data()->created_vertex = new Vertex_(p);

          // create new edges and set them as children of this edge
          // Note the numbering of the vertices: the new created vertex is the second one within the structure of both
          // edge children. This is exploited at some places and must not be changed.
          this->_set_num_children(2);
          _set_child(0, new Edge(vertex(0), this->subdiv_data()->created_vertex, this->refinement_level()+1));
          _set_child(1, new Edge(vertex(1), this->subdiv_data()->created_vertex, this->refinement_level()+1));
          // update the parent relationship
          this->child(0)->set_parent(this);
          this->child(1)->set_parent(this);

          // add new edges to the vector of created cells
          this->subdiv_data()->created_cells.push_back(this->child(0));
          this->subdiv_data()->created_cells.push_back(this->child(1));

          // set internal neighbourhood (external neighbourhood is set outside this function)
          // (in case space_dim_ > 1, an empty dummy function is called; see CellData)
          // vertex neighbours
          this->child(0)->add_neighbour(SDIM_VERTEX, 1, this->child(1));
          this->child(1)->add_neighbour(SDIM_VERTEX, 1, this->child(0));
        }
        else
        {
          std::cerr << "Wrong type of subdivision in edge " << this->print_index()
                    << ". There is only one type of subdivision for edges: CONFORM_SAME_TYPE. Aborting program."
                    << std::endl;
          exit(1);
        }
      } // subdivide()


      /// validate this cell
// COMMENT_HILMAR: will be done via exceptions
      void validate() const
      {
        try
        {
          if(space_dim_ == 1)
          {
            std::cout << "Validating edge " + this->print_index() + "\n";
          }

          std::string s;

          // validate that all vertices are set
          for(unsigned char ivert(0) ; ivert < num_vertices() ; ++ivert)
          {
            if (vertex(ivert) == nullptr)
            {
              s = "Edge " + this->print_index() + ": Vertex " + StringUtils::stringify((int)ivert) + " is null.\n";
              throw new InternalError(s);
            }
          }

          // validate children numbering
          if(!this->active())
          {
            if(this->subdiv_data()->type == CONFORM_SAME_TYPE)
            {
              // check whether the common vertex of the two children is set correctly
              if(this->child(0)->vertex(1) != this->child(1)->vertex(1))
              {
                s = "Edge " + this->print_index() + ": Shared vertex of the two children is not set correctly!\n";
                throw new InternalError(s);
              }
              // check whether the other vertices are set correctly
              if(this->child(0)->vertex(0) != vertex(0) || this->child(1)->vertex(0) != vertex(1))
              {
                s = "Edge " + this->print_index() + ": One of the end vertices of the children is not set correctly!\n";
                throw new InternalError(s);
              }
            }
          }
          // validate parent-child relations
          this->validate_history();

          if (this->active())
          {
            CellDataChecker<1, space_dim_, world_dim_>::check_neighbourhood(this);
          }
        }
        catch(InternalError* e)
        {
          std::cerr << e->message() << std::endl;
          exit(1);
        }
      }


      /// print information about this edge
      inline void print(std::ostream& stream)
      {
        stream << "E" << this->print_index() << ": [";
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
