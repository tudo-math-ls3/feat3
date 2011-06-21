/* GENERAL_REMARK_BY_HILMAR:
 * Generally, Peter wanted to take a deeper look at the base mesh implementation.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_BASE_MESH_CELL_SUBDIV_HPP
#define KERNEL_BASE_MESH_CELL_SUBDIV_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh/vertex.hpp>

/// FEAST namespace
namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
  namespace BaseMesh
  {

    /**
    * \brief type of subdivisions
    *
    * There are (at least) two motivations for subdividing base mesh cells:
    * \li change of shape (e.g. quad to triangle for better approximation of geometry; some special functionality
    *     may be available for a special element shape only)
    * \li more MPI processes needed
    * Depending on the situation, it may be
    * \li necessary to perform only conform subdivisions (do not split edges/faces)
    * \li allowed (or even wanted) to perform nonconform subdivisions (split edges/faces)
    *
    * \note These classifications are not complete yet. For example, transition elements are missing which may be
    * necessary to preserve or establish conformity of the mesh (e.g. prism elements for conform transitions from
    * hexas to tetras and vice versa). These classifications are closely connected to some "expert system" which
    * eventually decides where and how to subdivide. Until such an "expert system" is developed, the
    * classifications cannot be completed. The subdivision types below are just a starting point which may be useful
    * later.
    *
    * \note Possible strategies in such an "expert system" may be: "First try a conform subdivision. If this is not
    * possible (maybe due to badly shaped elements), then perform nonconform subdivision."
    * Note that for preserving or establishing the conformity of the mesh, the "expert system" has to examine the
    * neighbourhood of a cell. That means, whether a special subdivision leads to a conform mesh or not, does not
    * only depend on the cell to be subdivided itself! Consequently, one has to think about renaming the types below or
    * adding types of the form TO_QUAD_WITH_EDGE_BISECTION, TO_QUAD_WITHOUT_EDGE_BISECTION etc.
    * But that all depends on how this "expert system" will work. Hence, it doesn't make sense *now* to put too much
    * effort in these details.
    */
    enum type_of_subdivision
    {
      /**
      * \brief subdivision into cells of the same type, nonconform (splitting edges/faces)
      *
      * 2D: tri --> four tris
      * \verbatim
      *            /\                               /\
      *          /    \                           /    \
      *        /        \         --->          /________\
      *      /            \                   / \        / \
      *    /                \               /     \    /     \
      *  /                    \           /         \/         \
      * ------------------------         ------------------------
      * \endverbatim
      * 2D: quad --> four quads
      * \verbatim
      * -------------       -------------
      * |           |       |     |     |
      * |           |       |     |     |
      * |           |       |     |     |
      * |           | --->  -------------
      * |           |       |     |     |
      * |           |       |     |     |
      * |           |       |     |     |
      * -------------       -------------
      * \endverbatim
      * 3D: tetra --> eight tetras
      * 3D: hexa --> eight hexas
      */
      NONCONFORM_SAME_TYPE,

      /**
      * \brief subdivision into cells of the other type, eventually leading to a nonconform mesh
      *
      * 2D: tri --> three quads
      * \verbatim
      *            / \                               / \
      *          /     \                           /     \
      *        /         \         --->          / \     / \
      *      /             \                   /     \ /     \
      *    /                 \               /        |        \
      *  /                     \           /          |          \
      * -------------------------         -------------------------
      * \endverbatim
      * 2D: quad --> eight tris
      * \verbatim
      * -------------       -------------
      * |           |       |    /|    /|
      * |           |       |  /  |  /  |
      * |           |       |/    |/    |
      * |           | --->  -------------
      * |           |       |    /|    /|
      * |           |       |  /  |  /  |
      * |           |       |/    |/    |
      * -------------       -------------
      * \endverbatim
      * 3D: tetra --> hexas ???
      * 3D: hexa --> tetras ???
      */
      NONCONFORM_CHANGE_TYPE,

      /**
      * \brief subdivision into cells of the same type, leading to a nonconform mesh
      *
      * 2D: tri --> three tris (adding vertex in tri centre)
      * \verbatim
      *            / \                               /|\
      *          /     \                           /  |  \
      *        /         \         --->          /    |    \
      *      /             \                   /    __|__    \      (ugly ASCII art...)
      *    /                 \               /  ___/     \___  \
      *  /                     \           /  /               \  \
      * -------------------------         -------------------------
      * \endverbatim
      * 2D: quad --> five quads
      * \verbatim
      * ---------------      ---------------
      * |             |      | \         / |
      * |             |      |   \     /   |
      * |             |      |    -----    |
      * |             | ---> |    |   |    |
      * |             |      |    -----    |
      * |             |      |   /     \   |
      * |             |      | /         \ |
      * ---------------      ---------------
      * \endverbatim
      * 3D: tetra --> four tetras (adding vertex in tetra center; analogue to tri)
      * 3D: hexa --> seven hexas (small hexa in the centre, surrounded by six hexas; analogue to quad)
      */
      CONFORM_SAME_TYPE,

      /**
      * \brief subdivision into cells of the other type leading to a conform mesh
      *
      * 2D: tri --> quads: NOT POSSIBLE!
      * 2D: quad --> two tris (adding diagonal)
      * \verbatim
      * ---------------       ---------------
      * |             |       |            /|
      * |             |       |          /  |
      * |             |       |        /    |
      * |             | --->  |      /      |
      * |             |       |    /        |
      * |             |       |  /          |
      * |             |       |/            |
      * ---------------       ---------------
      * \endverbatim
      * 3D: tetra --> hexas: NOT POSSIBLE!
      * 3D: hexa --> tetras: NOT POSSIBLE!
      */
      CONFORM_CHANGE_TYPE,

      /**
      * \brief nonconform subdivision into triangles (only 2D)
      *
      * tri -->  four tris (see ASCII art of NONCONFORM_SAME_TYPE)
      * quad --> eight tris (see ASCII art of NONCONFORM_CHANGE_TYPE)
      */
      NONCONFORM_TO_TRI,

      /**
      * \brief nonconform subdivision into quads (only 2D)
      *
      * tri --> three quad (see ASCII art of NONCONFORM_CHANGE_TYPE)
      * quad --> four quads (see ASCII art of NONCONFORM_SAME_TYPE)
      */
      NONCONFORM_TO_QUAD,

      /**
      * \brief conform subdivision into trianges (only 2D)
      *
      * tri --> three tris (see ASCII art of CONFORM_SAME_TYPE)
      * quad --> two tris (see ASCII art of CONFORM_CHANGE_TYPE)
      */
      CONFORM_TO_TRI,

      /**
      * \brief conform subdivision into quads (only 2D)
      *
      * tri --> quads: NOT POSSIBLE!
      * quad --> five quads (see ASCII art of CONFORM_SAME_TYPE)
      */
      CONFORM_TO_QUAD,

      /**
      * \brief nonconform subdivision into tetraeders (only 3D)
      *
      * tetra --> eight tetras
      * hexa --> ??? tetras (not thought about yet...)
      */
      NONCONFORM_TO_TETRA,

      /**
      * \brief nonconform subdivision into hexaeders (only 3D)
      *
      * tetra --> ??? hexas (not thought about yet...)
      * hexa --> eight hexas
      */
      NONCONFORM_TO_HEXA,

      /**
      * \brief conform subdivision into tetraeders (only 3D)
      *
      * tetra --> four tetras (adding vertex in tetra center)
      * hexa --> tetras: NOT POSSIBLE!
      */
      CONFORM_TO_TETRA,

      /**
      * \brief conform subdivision into hexaeders (only 3D)
      *
      * tetra --> hexas: NOT POSSIBLE!
      * hexa --> seven hexas (small hexa in the centre, surrounded by six hexas; analogue to quad)
      */
      CONFORM_TO_HEXA
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
      SubdivisionData(type_of_subdivision t)
      {
        CONTEXT("BaseMesh::SubdivisionData::SubdivisionData()");
        type = t;
      }

      /// clears all vectors of created entities
      inline void clear_created()
      {
        CONTEXT("BaseMesh::SubdivisionData::clear_created()");
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
      SubdivisionData(type_of_subdivision t)
      {
        CONTEXT("BaseMesh::SubdivisionData::SubdivisionData()");
        type = t;
      }

      /// clears all vectors of created entities
      inline void clear_created()
      {
        CONTEXT("BaseMesh::SubdivisionData::clear_created()");
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
      SubdivisionData(type_of_subdivision t)
      {
        CONTEXT("BaseMesh::SubdivisionData::SubdivisionData()");
        type = t;
      }

      /// clears all vectors of created entities
      inline void clear_created()
      {
        CONTEXT("BaseMesh::SubdivisionData::clear_created()");
        created_vertices.clear();
        created_edges.clear();
        created_faces.clear();
        created_cells.clear();
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_SUBDIV_HPP
