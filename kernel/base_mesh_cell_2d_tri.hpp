#pragma once
#ifndef KERNEL_BASE_MESH_CELL_2D_TRI_HPP
#define KERNEL_BASE_MESH_CELL_2D_TRI_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <typeinfo>  // for typeid()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/base_mesh_vertex.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_cell_data_checker.hpp>
#include <kernel/base_mesh_cell_1d_edge.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    /**
    * \brief 2D base mesh cell of type triangle
    *
    * \author Hilmar Wobker
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

      /// shortcut for type Cell<1, space_dim_, world_dim_>
      typedef Cell<1, space_dim_, world_dim_> Cell_1D_;

    private:
      /// vertices of the triangle
      Vertex_* _vertices[3];

      /// edges of the triangle
      Cell_1D_* _edges[3];


      /// returns true when edge with local index iedge has the same orientation as the quad
      inline bool _edge_has_correct_orientation(unsigned char iedge)
      {
        // the orientation of the edge is correct (i.e. the same as that of the quad), when its start vertex within
        // the quad (i.e. the vertex with the same local index) is local vertex 0 within the edge structure
        return (vertex(iedge) == edge(iedge)->vertex(0));
      }


    public:
      /// CTOR
      Tri(
        Vertex_* v0, Vertex_* v1, Vertex_* v2,
        Cell_1D_* e0, Cell_1D_* e1, Cell_1D_* e2,
        unsigned char ref_level)
        : Cell<2, space_dim_, world_dim_>(ref_level)
      {
        _vertices[0] = v0;
        _vertices[1] = v1;
        _vertices[2] = v2;
        _edges[0] = e0;
        _edges[1] = e1;
        _edges[2] = e2;
        // assure that the edges are in fact of type Edge<space_dim_, world_dim_>, and not "only"
        // of type Cell<1, space_dim_, world_dim_>
        for(int i(0) ; i < 3 ; ++i)
        {
          assert(typeid(*_edges[i]) == typeid(Edge<space_dim_, world_dim_>));
        }

        unsigned char num_subitems_per_subdim[2] = {3,3};
        this->_set_num_subitems_per_subdim(2, num_subitems_per_subdim);
        this->_init_neighbours();

// TODO: Eigentlich haette ich das lieber in die Konstruktoren-Liste gepackt, also sowas in der Art:
//    : CellData<2, space_dim_, world_dim_>({3,3})
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
      inline Cell_1D_* edge(unsigned char const index) const
      {
        assert(index < num_edges());
        return _edges[index];
      }


      /// returns next vertex of vertex with given index w.r.t. to ccw ordering
      inline Vertex_* next_vertex_ccw(unsigned char const index) const
      {
        assert(index < num_vertices());
// TODO: to be implemented correctly!
        return nullptr;
      }


      /// returns previous vertex of vertex with given index w.r.t. to ccw ordering
      inline Vertex_* previous_vertex_ccw(unsigned char const index) const
      {
        assert(index < num_vertices());
// TODO: to be implemented correctly!
        return nullptr;
      }


      /// returns next edge of edge with given index w.r.t. to ccw ordering
      inline Cell_1D_* next_edge_ccw(unsigned char const index) const
      {
        assert(index < num_edges());
// TODO: to be implemented correctly!
        return nullptr;
      }


      /// returns previous edge of edge with given index w.r.t. to ccw ordering
      inline Cell_1D_* previous_edge_ccw(unsigned char const index) const
      {
        assert(index < num_edges());
// TODO: to be implemented correctly!
        return nullptr;
      }


      /// subdivision routine splitting a tri and storing parent/child information
      inline void subdivide()
      {
        // assure that this cell has not been divided yet
        if(!this->active())
        {
          std::cerr << "Tri ";
          this->print_index(std::cerr);
          std::cerr << " is already subdivided! Aborting program." << std::endl;
          exit(1);
        }

        if(!this->subdiv_data_initialised())
        {
          std::cerr << "Tri ";
          this->print_index(std::cerr);
          std::cerr << " cannot be subdivided! Initialise subdivision data first! Aborting program." << std::endl;
          exit(1);
        }

        // clear all vectors of created entities in the SubdivisionData object
        this->subdiv_data()->clear_created();

        // TODO: perform subdivision
        // ...
        std::cerr << "Subdivision for triangles not implemented yet!!" << std::endl;
        exit(1);
      } // subdivide()


      /// validates the cell
      inline void validate() const
      {
        try
        {
          if(space_dim_ == 2)
          {
            std::cout << "Validating triangle " << this->print_index() << std::endl;
          }

          std::string s = "Triangle " + this->print_index() + ": ";

          // validate that all vertices and edges are set
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

          // validate subitems (here: egdes)
          for(unsigned char iedge(0) ; iedge < num_edges() ; ++iedge)
          {
            edge(iedge)->validate();
          }

          // validate children numbering
          // TODO

          // validate parent-child relations
          this->validate_history();

          // validate neighbours
          if (this->active())
          {
            CellDataChecker<2, space_dim_, world_dim_>::check_neighbourhood(this);
          }
        }
        catch(InternalError* e)
        {
          std::cerr << e->message() << std::endl;
          exit(1);
        }
      }

      /// print information about this tri
      inline void print(std::ostream& stream)
      {
        stream << "Tri";
        this->print_index(stream);

        stream << ": [V ";
        _vertices[0]->print_index(stream);
        for(int i(1) ; i < num_vertices() ; ++i)
        {
          stream << ", ";
          _vertices[i]->print_index(stream);
        }
        stream << "] [";

        stream << "E ";
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
        stream << "] ";
        Cell<2, space_dim_, world_dim_>::print_history(stream);
        // print neighbourhood information (if there is any)
        CellData<2, space_dim_, world_dim_>::print(stream);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_2D_TRI_HPP
