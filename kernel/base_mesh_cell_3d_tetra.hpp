#pragma once
#ifndef KERNEL_BASE_MESH_CELL_3D_TETRA_HPP
#define KERNEL_BASE_MESH_CELL_3D_TETRA_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <typeinfo>  // for typeid()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_vertex.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_cell_data_checker.hpp>
#include <kernel/base_mesh_cell_1d_edge.hpp>
#include <kernel/base_mesh_cell_2d_tri.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    /**
    * \brief 3D base mesh cell of type tetra
    *
    * \author Hilmar Wobker
    */
// COMMENT_HILMAR: Um Code-Redundanz zu vermeiden, koennte wir ueberlegen, eine weitere Klasse Cell3D einzufuehren,
// die von Cell<3, space_dim_, world_dim_> erbt, und von der dann wieder um Tetra und Hexa erben.
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Tetra
      : public Cell<3, space_dim_, world_dim_>
    {
      /// shortcut for type Vertex<world_dim_>
      typedef Vertex<world_dim_> Vertex_;

      /// shortcut for type Cell<1, space_dim_, world_dim_>
      typedef Cell<1, space_dim_, world_dim_> Cell_1D_;

      /// shortcut for type Tri<space_dim_, world_dim_>
      typedef Tri<space_dim_, world_dim_> Tri_;

      /// shortcut for type Cell<2, space_dim_, world_dim_>
      typedef Cell<2, space_dim_, world_dim_> Cell_2D_;

    private:
      /// vertices of the tetra
      Vertex_* _vertices[4];

      /// edges of the tetra
      Cell_1D_* _edges[6];

      /// edges of the tetra
      Cell_2D_* _faces[4];


    public:
      /// CTOR
      Tetra(
        Vertex_* v0, Vertex_* v1, Vertex_* v2, Vertex_* v3,
        Cell_1D_* e0, Cell_1D_* e1, Cell_1D_* e2, Cell_1D_* e3, Cell_1D_* e4, Cell_1D_* e5,
        Cell_2D_* f0, Cell_2D_* f1,Cell_2D_* f2, Cell_2D_* f3,
        unsigned char ref_level)
        : Cell<3, space_dim_, world_dim_>(ref_level)
      {
        _vertices[0] = v0;
        _vertices[1] = v1;
        _vertices[2] = v2;
        _vertices[3] = v3;
        _edges[0] = e0;
        _edges[1] = e1;
        _edges[2] = e2;
        _edges[3] = e3;
        _edges[4] = e4;
        _edges[5] = e5;
        // assure that the edges are in fact of type Edge<space_dim_, world_dim_>, and not "only"
        // of type Cell<1, space_dim_, world_dim_>
        for(int i(0) ; i < 6 ; ++i)
        {
          assert(typeid(*_edges[i]) == typeid(Edge<space_dim_, world_dim_>));
        }
        _faces[0] = f0;
        _faces[1] = f1;
        _faces[2] = f2;
        _faces[3] = f3;
        // assure that the faces are in fact of type Tri_, and not "only" of type Cell_2D_
        for(int i(0) ; i < 4 ; ++i)
        {
          assert(typeid(*_faces[i]) == typeid(Tri_));
        }

        unsigned char num_subitems_per_subdim[3] = {4, 6, 4};
        this->_set_num_subitems_per_subdim(3, num_subitems_per_subdim);
        this->_init_neighbours();
// COMMENT_HILMAR: Eigentlich haette ich das lieber in die Konstruktoren-Liste gepackt, also sowas in der Art:
//    : CellData<3, space_dim_, world_dim_>({4,6,4})
// (was nicht kompiliert). Wie kann man denn on-the-fly ein Array anlegen und durchreichen?
      }


      /// returns number of vertices
      inline unsigned char num_vertices() const
      {
        return 4;
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
        return 6;
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
        return 4;
      }


      /// returns face at given index
      inline Cell_2D_* face(unsigned char const index) const
      {
        assert(index < num_faces());
        return _faces[index];
      }


      /// subdivision routine splitting a tetra and storing parent/child information
// COMMENT_HILMAR: this is currently hard-wired to splitting the tetra into 8 tetras (called 'standard partition' or
// '3D-Freudenthal-Bey partition'. Later, this is parameterised via the information in the SubdivisionData object.
      inline void subdivide(SubdivisionData<3, space_dim_, world_dim_>& subdiv_data)
      {
        // assure that this cell has not been divided yet
        if(!this->active())
        {
          std::cerr << "Tetra ";
          this->print_index(std::cerr);
          std::cerr << " is already subdivided! Aborting program." << std::endl;
          exit(1);
        }

        // to be implemented...

      } // subdivide()


      /// validates the cell
// COMMENT_HILMAR: will be done via exceptions
      inline void validate() const
      {
        // to be implemented
      }


      /// print information about this quad
      inline void print(std::ostream& stream)
      {
        stream << "Tetra";
        this->print_index(stream);
        stream << ": [";

        for(int i(0) ; i < num_faces() ; ++i)
        {
          stream << "F";
          _faces[i]->print_index(stream);
          if(i < num_faces()-1)
          {
            stream << ", ";
          }
          else
          {
            stream << "]";
          }
        }
        stream << std::endl << "    ";
        Cell<3, space_dim_, world_dim_>::print_history(stream);
        // print neighbourhood information (if there is any)
        CellData<3, space_dim_, world_dim_>::print(stream);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_3D_TETRA_HPP
