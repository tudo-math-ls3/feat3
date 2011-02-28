#pragma once
#ifndef KERNEL_BASE_MESH_CELL_3D_TETRA_HPP
#define KERNEL_BASE_MESH_CELL_3D_TETRA_HPP 1

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
#include <kernel/base_mesh/cell_2d_tri.hpp>

/// FEAST namespace
namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
  namespace BaseMesh
  {

    /**
    * \brief 3D base mesh cell of type tetra
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
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
      /// shortcut to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;

      /// shortcut to save typing of template parameters
      typedef Cell<1, space_dim_, world_dim_> Cell_1D_;

      /// shortcut to save typing of template parameters
      typedef Tri<space_dim_, world_dim_> Tri_;

      /// shortcut to save typing of template parameters
      typedef Cell<2, space_dim_, world_dim_> Cell_2D_;


    private:

      /* *****************
      * member variables *
      *******************/
      /// vertices of the tetra
      Vertex_* _vertices[4];

      /// edges of the tetra
      Cell_1D_* _edges[6];

      /// edges of the tetra
      Cell_2D_* _faces[4];


    public:

      /**
      * \brief CTOR
      *
      * \param[in] v0, v1, v2, v3
      * vertices the triangle is built of
      *
      * \param[in] e0, e1, e2, e3, e4, e5
      * edges the triangle is built of
      *
      * \param[in] f0, f1, f2, f3
      * faces the triangle is built of
      *
      * \param[in] ref_level
      * refinement level the triangle is created on
      */
      Tetra(
        Vertex_* v0, Vertex_* v1, Vertex_* v2, Vertex_* v3,
        Cell_1D_* e0, Cell_1D_* e1, Cell_1D_* e2, Cell_1D_* e3, Cell_1D_* e4, Cell_1D_* e5,
        Cell_2D_* f0, Cell_2D_* f1,Cell_2D_* f2, Cell_2D_* f3,
        unsigned char ref_level)
        : Cell<3, space_dim_, world_dim_>(ref_level)
      {
        CONTEXT("BaseMesh::Tetra::Tetra()");
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
        CONTEXT("BaseMesh::Tetra::num_vertices()");
        return 4;
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
        CONTEXT("BaseMesh::Tetra::vertex()");
        assert(index < num_vertices());
        return _vertices[index];
      }


      /**
      * \brief returns number of edges
      *
      * \return number of edges
      */
      inline unsigned char num_edges() const
      {
        CONTEXT("BaseMesh::Tetra::num_edges()");
        return 6;
      }


      /**
      * \brief returns pointer to the edge at given index
      *
      * \param[in] index
      * index of the edge to be returned
      *
      * \return pointer to the edge at given index
      */
      inline Cell_1D_* edge(unsigned char const index) const
      {
        CONTEXT("BaseMesh::Tetra::edge()");
        assert(index < num_edges());
        return _edges[index];
      }


      /**
      * \brief returns number of faces
      *
      * \return number of faces
      */
      inline unsigned char num_faces() const
      {
        CONTEXT("BaseMesh::Tetra::num_faces()");
        return 4;
      }


      /**
      * \brief returns pointer to the face at given index
      *
      * \param[in] index
      * index of the face to be returned
      *
      * \return pointer to the face at given index
      */
      inline Cell_2D_* face(unsigned char const index) const
      {
        CONTEXT("BaseMesh::Tetra::face()");
        assert(index < num_faces());
        return _faces[index];
      }


      /**
      * \brief subdivision routine splitting a tetra and storing parent/child information
      *
      * \param[in,out] subdiv_data
      * pointer to the subdivision data object
      */
      inline void subdivide(SubdivisionData<3, space_dim_, world_dim_>* subdiv_data)
      {
        CONTEXT("BaseMesh::Tetra::subdivide()");
        try
        {
          // assure that this cell has not been divided yet
          if(!this->active())
          {
              throw InternalError("Tetraeder " + this->print_index() + " is already subdivided!");
          }

          this->set_subdiv_data(subdiv_data);

          // clear all vectors of created entities in the SubdivisionData object
          this->subdiv_data()->clear_created();


          if(this->subdiv_data()->type == NONCONFORM_SAME_TYPE)
          {
            // Perform a nonconform subdivision without changing the cell type (1 tetra --> 8 tetras; called
            // 'standard partition' or '3D-Freudenthal-Bey partition'.
            // TODO: perform subdivision
            // ...
            throw InternalError("Subdivision for tetraeders not implemented yet!!");
          }
          else
          {
            throw InternalError("Wrong type of subdivision in tetra " + this->print_index()
                                + ". Currently, only subdivision NONCONFORM_SAME_TYPE is supported.");
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
      inline void validate(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Tetra::validate()");
        try
        {
          if(space_dim_ == 3)
          {
            stream << "Validating triangle " << this->print_index() << "..." << std::endl;
          }

          std::string s = "Triangle " + this->print_index() + ": ";

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

          // to be implemented

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


      /**
      * \brief prints information about this tetraeder to the given stream
      *
      * \param[in,out] stream
      * stream to write into
      */
      inline void print(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Tetra::print()");
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
