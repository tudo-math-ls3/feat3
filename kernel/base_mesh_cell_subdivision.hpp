#pragma once
#ifndef KERNEL_BASE_MESH_CELL_SUBDIV_HPP
#define KERNEL_BASE_MESH_CELL_SUBDIV_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_vertex.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    enum type_of_subdivision
    {
      // subdivision into cells of the same type, eventually leading to a nonconform mesh
      // (quad --> four quads, tri --> four tris)
      // -------------       -------------
      // |           |       |     |     |
      // |           |       |     |     |
      // |           | --->  -------------
      // |           |       |     |     |
      // |           |       |     |     |
      // -------------       -------------
      NONCONFORM_SAME_TYPE,

      // subdivision into cells of the other type, eventually leading to a nonconform mesh
      // quad --> six tris (dividing edges, adding one line between two opposing new vertices)
      // -------------       -------------
      // |           |       |    /|\    |
      // |           |       |  /  |  \  |
      // |           |  -->  |/    |    \|
      // |           |       |\    |    /|
      // |           |       |  \  |  /  |
      // |           |       |    \|/    |
      // -------------       -------------
      // tri --> four quads (dividing edges and adding vertex in triangle center)
      NONCONFORM_CHANGE_TYPE,

      // subdivision into cells of the same type, leading to a nonconform mesh
      // quad --> five quads (adding small quad within the quad, connect two vertices, resp.)
      // tri --> three tris (adding vertex in triangle center)
      CONFORM_SAME_TYPE,

      // subdivision into cells of the other type leading to a conform mesh
      // quad --> two tris (adding diagonal)
      // -------------       -------------
      // |           |       |          /|
      // |           |       |        /  |
      // |           |  -->  |      /    |
      // |           |       |    /      |
      // |           |       |  /        |
      // |           |       |/          |
      // -------------       -------------
      // quad --> four tris (adding both diagonals)
      // -------------       -------------
      // |           |       |\         /|
      // |           |       |  \     /  |
      // |           |  -->  |    \ /    |
      // |           |       |    / \    |
      // |           |       |  /     \  |
      // |           |       |/         \|
      // -------------       -------------
      // tri --> quads: not possible!
      CONFORM_CHANGE_TYPE
    };


    /// forward declaration of class BaseMesh::Cell
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Cell;


    /**
    * \brief class containing subdivision specific data (empty definition to be specialised by cell_dim_)
    *
    * The main purpose of this template design is to enable the usage of a common interface for the function
    *   void subdivide(SubdivisionData<cell_dim_, space_dim_, world_dim_>& subdiv_data)
    * such that it can be declared in class Cell<cell_dim_, space_dim_, world_dim_>.
    * On the one hand, the class will contain "return" vectors of entities that have been created during the subdivision
    * process, on the other hand it holds parameters that steer the subdivision (type of subdivision, anisotropy,
    * factors, ...). It's hard to define such data independently of the cell dimension, hence the class is specialised
    * via the cell dimension. Another advantage: The interface of the function subdivide(...) will never have
    * to be changed again.
    *
    * \author Hilmar Wobker
    */
// COMMENT_HILMAR: Wie und wo sollen die SubdivisionData-Objekte gespeichert werden? Als member von Cell<...>?
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    struct SubdivisionData
    {
    };


    /**
    * \brief subdivision specific data for 1D cells
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    struct SubdivisionData<1, space_dim_, world_dim_>
    {
      /// type of subdivision
      type_of_subdivision type;

      /// new vertex created during subdivision
      Vertex<world_dim_>* created_vertex;

      /**
      * \brief new cells created during subdivision
      *
      * For sake of simplicity also pointers to the created cells are stored here. Thus, they can simply be added
      * to the base mesh via BaseMesh1D.add_created_items(SubdivisionData<...>& subdiv_data).
      * The alternative would be to to access these cells as children of the just subdivided cell.
      */
      std::vector<Cell<1, space_dim_, world_dim_>*> created_cells;

      /// CTOR
      SubdivisionData<1, space_dim_, world_dim_>(type_of_subdivision t)
      {
        type = t;
      }

      /// clears all vectors of created entities
      inline void clear_created()
      {
        created_vertex = nullptr;
        created_cells.clear();
      }
    };


    /**
    * \brief subdivision specific data for 2D cells
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    struct SubdivisionData<2, space_dim_, world_dim_>
    {
      /// type of subdivision
      type_of_subdivision type;

      /// new vertices created during subdivision
      std::vector<Vertex<world_dim_>*> created_vertices;

      /// new edges created during subdivision
      std::vector<Cell<1, space_dim_, world_dim_>*> created_edges;

      /**
      * \brief new cells created during subdivision
      *
      * For sake of simplicity also pointers to the created cells are stored here. Thus, they can simply
      * be added to the base mesh via BaseMesh2D.add_created_items(SubdivisionData<...>& subdiv_data).
      * The alternative would be to to access these cells as children of the just subdivided cell.
      */
      std::vector<Cell<2, space_dim_, world_dim_>*> created_cells;

      /// CTOR
      SubdivisionData<2, space_dim_, world_dim_>(type_of_subdivision t)
      {
        type = t;
      }

      /// clears all vectors of created entities
      inline void clear_created()
      {
        created_vertices.clear();
        created_edges.clear();
        created_cells.clear();
      }
    };


    /**
    * \brief subdivision specific data for 3D cells
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    struct SubdivisionData<3, space_dim_, world_dim_>
    {
      /// type of subdivision
      type_of_subdivision type;

      /// new vertices created during subdivision
      std::vector<Vertex<world_dim_>*> created_vertices;

      /// new edges created during subdivision
      std::vector<Cell<1, space_dim_, world_dim_>*> created_edges;

      /// new faces created during subdivision
      std::vector<Cell<2, space_dim_, world_dim_>*> created_faces;

      /**
      * \brief new cells created during subdivision
      *
      * For sake of simplicity also pointers to the created cells are stored here. Thus, they can simply
      * be added to the base mesh via BaseMesh3D.add_created_items(SubdivisionData<...>& subdiv_data).
      * The alternative would be to to access these cells as children of the just subdivided cell.
      */
      std::vector<Cell<3, space_dim_, world_dim_>*> created_cells;

      /// CTOR
      SubdivisionData<3, space_dim_, world_dim_>(type_of_subdivision t)
      {
        type = t;
      }

      /// clears all vectors of created entities
      inline void clear_created()
      {
        created_vertices.clear();
        created_edges.clear();
        created_faces.clear();
        created_cells.clear();
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_SUBDIV_HPP
