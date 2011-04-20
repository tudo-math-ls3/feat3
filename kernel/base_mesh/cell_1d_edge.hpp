#pragma once
#ifndef KERNEL_BASE_MESH_CELL_1D_EDGE_HPP
#define KERNEL_BASE_MESH_CELL_1D_EDGE_HPP 1

// includes, system
#include <iostream> // for std::ostream

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/base_mesh/cell.hpp>
#include <kernel/base_mesh/cell_data_validation.hpp>
#include <kernel/base_mesh/vertex.hpp>

/// FEAST namespace
namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
  namespace BaseMesh
  {

    /**
    * \brief 1D base mesh cell of type edge
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
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
      /// shortcut to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;


    private:

      /// vertices of the edge
      Vertex_* _vertices[2];

    public:

      /**
      * \brief CTOR
      *
      * \param[in] v0, v1
      * vertices the edge is built of
      *
      * \param[in] ref_level
      * refinement level the edge is created on
      */
      Edge(
        Vertex_* v0,
        Vertex_* v1,
        unsigned char ref_level)
        : Cell<1, space_dim_, world_dim_>(ref_level)
      {
        CONTEXT("BaseMesh::Edge::Edge()");
        _vertices[0] = v0;
        _vertices[1] = v1;

        unsigned char num_subitems_per_subdim[1] = {2};
        this->_set_num_subitems_per_subdim(1, num_subitems_per_subdim);
        this->_init_neighbours();
      }


      /* **********
      * functions *
      ************/
      /**
      * \brief returns number of vertices
      *
      * \return number of vertices
      */
      inline unsigned char num_vertices() const
      {
        CONTEXT("BaseMesh::Edge::num_vertices()");
        return 2;
      }


      /**
      * \brief returns pointer to the vertex at given index
      *
      * \param[in] index
      * index of the vertex to be returned
      *
      * \return pointer to the vertex at given index
      */
      inline Vertex_* vertex(unsigned char const index) const
      {
        CONTEXT("BaseMesh::Edge::vertex()");
        ASSERT(index < num_vertices(), "Index " + stringify(index) + " must not exceed number of vertices "
               + stringify(num_vertices()) + ".");
        return _vertices[index];
      }


      /**
      * \brief subdivision routine splitting an edge into two and storing parent/child information
      *
      * \param[in,out] subdiv_data
      * pointer to the subdivision data object
      */
      inline void subdivide(SubdivisionData<1, space_dim_, world_dim_>* subdiv_data)
      {
        CONTEXT("BaseMesh::Edge::subdivide()");
        try
        {
          // assure that this cell has not been divided yet
          if(!this->active())
          {
            throw InternalError("Edge " + this->print_index() + " is already subdivided!");
          }

          this->set_subdiv_data(subdiv_data);

          // clear all vectors of created entities in the SubdivisionData object
          this->subdiv_data()->clear_created();

          if(this->subdiv_data()->type == CONFORM_SAME_TYPE)
          {
            // split the edge into two edges

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
            this->_set_child(0, new Edge(vertex(0), this->subdiv_data()->created_vertex, this->refinement_level()+1));
            this->_set_child(1, new Edge(vertex(1), this->subdiv_data()->created_vertex, this->refinement_level()+1));
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
            throw InternalError("Wrong type of subdivision in edge " + this->print_index()
                                + ". Currently, only subdivision CONFORM_SAME_TYPE is supported.");
          }
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      } // subdivide()


      /**
      * \brief validates this cell
      *
      * \param[in,out] stream
      * stream validation info is written into
      */
      void validate(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Edge::validate()");
        try
        {
          if(space_dim_ == 1)
          {
            stream << "Validating edge " + this->print_index() << "..." << std::endl;
          }

          std::string s = "Edge " + this->print_index() + ": ";

          // validate that all vertices are set
          for(unsigned char ivert(0) ; ivert < num_vertices() ; ++ivert)
          {
            if (vertex(ivert) == nullptr)
            {
              s += "Vertex " + stringify((int)ivert) + " is null.\n";
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
                s += "Shared vertex of the two children is not set correctly!\n";
                throw new InternalError(s);
              }
              // check whether the other vertices are set correctly
              if(this->child(0)->vertex(0) != vertex(0) || this->child(1)->vertex(0) != vertex(1))
              {
                s += "One of the end vertices of the children is not set correctly!\n";
                throw new InternalError(s);
              }
            }
          }
          // validate parent-child relations
          this->validate_history(stream);

          if (this->active())
          {
            CellDataValidation<1, space_dim_, world_dim_>::validate_neighbourhood(this, stream);
          }
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      }


      /**
      * \brief prints information about this edge to the given stream
      *
      * \param[in,out] stream
      * stream to write into
      */
      inline void print(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Edge::print()");
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
