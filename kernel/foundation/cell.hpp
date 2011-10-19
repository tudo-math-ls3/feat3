/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR below.
 * Generally, Peter wanted to take a deeper look at the base mesh implementation.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_BASE_MESH_CELL_HPP
#define KERNEL_BASE_MESH_CELL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/foundation/error_handler.hpp>

#include <kernel/foundation/cell_data.hpp>
#include <kernel/foundation/cell_subdivision.hpp>
#include <kernel/foundation/vertex.hpp>

// includes, system
#include <iostream> // for std::ostream
#include <vector>   // for std::vector

namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
  namespace BaseMesh
  {
    /**
    * \brief stores the fixed numbering schemes
    *
    * \author Hilmar Wobker
    */
    struct Numbering
    {
      /// indices of start and end vertex of the four edges in a quad
      static unsigned char quad_edge_vertices[][2];

      /**
      * \brief index of the next vertex in a quad w.r.t. ccw ordering ([1,3,0,2])
      *
      * Due to this numbering scheme of the quad
      *   2---1---3
      *   |       |
      *   2       3
      *   |       |
      *   0---0---1
      * the respectively next vertex w.r.t. ccw ordering is given by the mapping [1,3,0,2]
      */
      static unsigned char quad_next_vertex_ccw[];

      /**
      * \brief  index of the previous vertex in a quad w.r.t. ccw ordering ([2,0,3,1])
      *
      * See the detailed description of #quad_next_vertex_ccw.
      */
      static unsigned char quad_previous_vertex_ccw[];

      /**
      * \brief  index of the next edge in a quad w.r.t. ccw ordering ([3,2,0,1])
      *
      * See the detailed description of #quad_next_vertex_ccw.
      */
      static unsigned char quad_next_edge_ccw[];

      /**
      * \brief  index of the previous edge in a quad w.r.t. ccw ordering ([2,3,1,0])
      *
      * See the detailed description of #quad_next_vertex_ccw.
      */
      static unsigned char quad_previous_edge_ccw[];

      /// indices of start and end vertex of the twelve edges in a hexa
      static unsigned char hexa_edge_vertices[][2];

      /// indices of the four vertices of the six faces in a hexa
      static unsigned char hexa_face_vertices[][4];

      /// indices of the four edges of the six faces in a hexa
      static unsigned char hexa_face_edges[][4];

      /**
      * \brief quad-to-quad mappings of vertices
      *
      * On the one hand a quad face in a 3D cell has a certain numbering w.r.t. the numbering of the 3D cell, i.e.
      * the numbering of the 3D cell determines which are first, second, ... vertex/edge of the face. On the
      * other hand the quad is stored as a 2D cell with a certain numbering. These two numberings usually do not
      * coincide. There are eight possibilities how the two numberings can be related. The first four have the same
      * orientation as the reference numeration, i.e. they are only rotated, the last four have opposite orientation.
      * The eight possibilites lead to eight different mappings of vertices and edges, resp., which we denote by, e.g.,
      *   V1:2031
      *   (relation 2: vertex 0 (of the reference numbering) equals vert. 2 in the different numbering,
      *                vert. 1 equals vert. 0, ...)
      *   E7:3210
      *   (relation 7: edge 0 (of the reference numbering) equals edge 3 in the different numbering, ...).
      *
      * same orientation as reference                          opposite orientation
      *
      * relation 0   relation 1   relation 2   relation 3      relation 4   relation 5   relation 6   relation 7
      *
      * 2---1---3    0---2---2    3---3---1    1---0---0       1---3---3    3---1---2    0---0---1    2---2---0
      * |       |    |       |    |       |    |       |       |       |    |       |    |       |    |       |
      * 2       3    0       1    1       0    3       2       0       1    3       2    2       3    1       0
      * |       |    |       |    |       |    |       |       |       |    |       |    |       |    |       |
      * 0---0---1    1---3---3    2---2---0    3---1---2       0---2---2    1---0---0    2---1---3    3---3---1
      *  V0:0123      V1:1302      V2:2031      V3:3210         V4:0213      V5:1032      V6:2301      V7:3120
      *  E0:0123      E1:3201      E2:2310      E3:1032         E4:2301      E5:0132      E6:1023      E7:3210
      * (reference)
      *
      * When vertex i of the given numbering equals vertex 0 in the reference numbering, then we have either
      * relation i or relation i+4 (depending on the orientation).
      */
// COMMENT_HILMAR: Sieht hier jemand eine Moeglichkeit, dieses Mapping ohne explizites Abspeichern zu erhalten? Also
// durch Berechnung in Abhaengigkeit vom Knoten- bzw. Kantenindex? Ich seh's nicht...
// COMMENT_HILMAR: Mit dem Standard-ccw-Mapping waere das alles einfacher... da koennte man mit sowas aehnlichem wie
// (ivertex + ishift)%4 arbeiten.
      static unsigned char quad_to_quad_mappings_vertices[][4];

      /**
      * \brief quad-to-quad mappings of edges
      *
      * See the description of #quad_to_quad_mappings_vertices.
      */
      static unsigned char quad_to_quad_mappings_edges[][4];
    };





    /**
    * \brief dimension specific function interface for cell class (empty definition to be specialised by cell dimension)
    *
    * While vertices are common to all cell dimensions, edges and faces do not exist in all dimensions. So,
    * getter and setter routines for these entities cannot be provided in the general Cell class.
    * The alternative would be to specialise the complete Cell class by cell dimension.
    * Why do we need this interface already in the Cell class? Since we must be able to do something like
    *   face(iface)->child(0)->edge(1)
    * where child(0) is of type Cell<2, ...> (and not, e.g., of type Quad<...>!).
    *
    * \tparam cell_dim_
    * cell dimension (must be <= space dimension)
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class CellInterface
    {
    };




    /**
    * \brief function interface for 1D cells
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class CellInterface<1, space_dim_, world_dim_>
    {

    public:

      /// pure virtual function for returning number of vertices
      virtual unsigned char num_vertices() const = 0;


      /// pure virtual function for returning pointer to the vertex at given index
      virtual Vertex<world_dim_>* vertex(unsigned char const) const = 0;
    };




    /**
    * \brief function interface for 2D cells
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class CellInterface<2, space_dim_, world_dim_>
    {

    public:

      /// pure virtual function for returning number of vertices
      virtual unsigned char num_vertices() const = 0;


      /// pure virtual function for returning pointer to the vertex at given index
      virtual Vertex<world_dim_>* vertex(unsigned char const) const = 0;


      /// pure virtual function for returning number of edges
      virtual unsigned char num_edges() const = 0;


      /// pure virtual function for returning pointer to the edge at given index
      virtual Cell<1, space_dim_, world_dim_>* edge(unsigned char const) const = 0;


      /// pure virtual function for returning next vertex of vertex with given index w.r.t. to ccw ordering
      virtual Vertex<world_dim_>* next_vertex_ccw(unsigned char const) const = 0;


      /// pure virtual function for returning previous vertex of vertex with given index w.r.t. to ccw ordering
      virtual Vertex<world_dim_>* previous_vertex_ccw(unsigned char const) const = 0;


      /// pure virtual function for returning next edge of edge with given index w.r.t. to ccw ordering
      virtual Cell<1, space_dim_, world_dim_>* next_edge_ccw(unsigned char const) const = 0;


      /// pure virtual function for returning previous edge of edge with given index w.r.t. to ccw ordering
      virtual Cell<1, space_dim_, world_dim_>* previous_edge_ccw(unsigned char const) const = 0;
    };




    /**
    * \brief function interface for 3D cells
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class CellInterface<3, space_dim_, world_dim_>
    {

    public:

      /// pure virtual function for returning number of vertices
      virtual unsigned char num_vertices() const = 0;


      /// pure virtual function for returning pointer to the vertex at given index
      virtual Vertex<world_dim_>* vertex(unsigned char const) const = 0;


      /// pure virtual function for returning number of edges
      virtual unsigned char num_edges() const = 0;


      /// pure virtual function for returning pointer to the edge at given index
      virtual Cell<1, space_dim_, world_dim_>* edge(unsigned char const) const = 0;


      /// pure virtual function for returning number of faces
      virtual unsigned char num_faces() const = 0;


      /// pure virtual function for returning pointer to the face at given index
      virtual Cell<2, space_dim_, world_dim_>* face(unsigned char const) const = 0;
    };




    /**
    * \brief general base mesh cell class containing parent/child information
    *
    * Used for cells of maximum dimension (e.g., quads in a 2D world), but also for those of lower dimension
    * (e.g., edges in a 2D world). For the latter, however, the CellData class is empty such that no unnecessary
    * neighbourhood information etc. included.
    *
    * \tparam cell_dim_
    * cell dimension (must be <= space dimension)
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Cell
      : public CellData<cell_dim_, space_dim_, world_dim_>,
        public CellInterface<cell_dim_, space_dim_, world_dim_>
    {

    private:

      /// parent of this cell
      Cell* _parent;

      /// number of children (when zero, then the cell is called active)
      unsigned char _num_children;

      /// array of children of this cell
      Cell** _children;

      /// subdivision data
      SubdivisionData<cell_dim_, space_dim_, world_dim_>* _subdiv_data;

      /// refinement level of this cell
      unsigned char _refinement_level;


    protected:

      /**
      * \brief sets number of children and allocates the _children array
      *
      * This function may only be called, when new children are to be created. "Old" children have to be cleared before
      * via unset_children(...).
      * Rationale:
      * - It won't happen often that a cell creates and removes children over and over again.
      * - A situation like "I already have 4 children, and now I want to have 2 more" will not occur.
      * - Allocating/deallocating an array of less than 10 pointers does not really harm.
      *
      * \param[in] num
      * new number of children
      */
      inline void _set_num_children(unsigned char const num)
      {
        CONTEXT("BaseMesh::Cell::_set_num_children()");
        // this function must not be called when there *are* already children
        ASSERT(_children == nullptr, "Array of children must be the nullptr.");
        ASSERT(_num_children == 0, "Number of children " + stringify(_num_children) + " must be zero.");
        // and it must not be called to unset children (use unset_children() for that)
        ASSERT(num > 0,
               "This function must not be used to set number of children to zero. Use unset_children() instead.");
        _num_children = num;
        _children = new Cell*[_num_children];
        // nullify the new pointer array
        for (unsigned char i(0) ; i < _num_children ; ++i)
        {
          _children[i] = nullptr;
        }
      }


      /**
      * \brief sets child at given index
      *
      * The cell-type does not have to be the same (quads can have tris as children), so a general BaseMeshItem2D is
      * passed to the function.
      *
      * \param[in] index
      * position of the child in array #_children
      *
      * \param[in] e
      * cell pointer to be set as child
      */
      inline void _set_child(
        unsigned char const index,
        Cell* e)
      {
        CONTEXT("BaseMesh::Cell::_set_child()");
        ASSERT(index < num_children(), "Index " + stringify(index) + " must not exceed number of children "
               + stringify(num_children()) + " .");
        // ensure that this function is not used to unset children (use unset_children() for that)
        ASSERT(e != nullptr,
               "This function must not be used to set children to null. Use unset_children() instead.");
        _children[index] = e;
      }


      /**
      * \brief unsets all children
      *
      * Deletes the array _children and sets _num_children to zero. The actual destruction of the children is
      * triggered by the base mesh.
      *
      * \note Unsetting children should not be done via _set_child(i, nullptr).
      */
      inline void _unset_children()
      {
        CONTEXT("BaseMesh::Cell::_unset_children()");
        if(_children != nullptr)
        {
          delete [] _children;
          _children = nullptr;
          _num_children = 0;
        }
      }


      /**
      * \brief deletes the subdivision data object
      *
      * Deletes the subdivision data pointer #_subdiv_data when it is not the nullptr.
      */
      inline void _delete_subdiv_data()
      {
        CONTEXT("BaseMesh::Cell::_delete_subdiv_data()");
        if(_subdiv_data != nullptr)
        {
          _subdiv_data->clear_created();
          delete _subdiv_data;
          _subdiv_data = nullptr;
        }
      }


    public:

      /**
      * \brief CTOR
      *
      * \param[in] ref_level
      * level of refinement on which this cell is created
      */
      Cell(unsigned char ref_level)
        : _parent(nullptr),
          _num_children(0),
          _children(nullptr),
          _subdiv_data(nullptr),
          _refinement_level(ref_level)
      {
        CONTEXT("BaseMesh::Cell::Cell()");
        ASSERT(world_dim_ >= space_dim_, "Space dimension " + stringify(space_dim_)
               + " must not exceed world dimension " + stringify(world_dim_) + ".");
        ASSERT(space_dim_ >= cell_dim_, "Cell dimension " + stringify(space_dim_)
               + " must not exceed space dimension " + stringify(cell_dim_) + ".");
      }


      /// DTOR (automatically virtual since DTOR of base class Item is virtual)
      ~Cell()
      {
        CONTEXT("BaseMesh::Cell::~Cell()");
        _delete_subdiv_data();
        _unset_children();
        _parent = nullptr;
      }


      /**
      * \brief returns pointer to parent cell
      *
      * \return pointer to parent cell
      */
      inline Cell* parent() const
      {
        CONTEXT("BaseMesh::Cell::parent()");
        return _parent;
      }


      /**
      * \brief sets pointer to parent cell
      *
      * \param[in] par
      * pointer to parent cell
      */
      inline void set_parent(Cell* const par)
      {
        CONTEXT("BaseMesh::Cell::set_parent()");
        _parent = par;
      }


      /**
      * \brief returns number of children
      *
      * \return number of children
      */
      inline unsigned char num_children() const
      {
        CONTEXT("BaseMesh::Cell::num_children()");
        return _num_children;
      }


      /**
      * \brief returns pointer to child at given index
      *
      * \param[in] index of the child to be returned
      *
      * \return pointer to child at given index
      */
      inline Cell* child(unsigned char index) const
      {
        CONTEXT("BaseMesh::Cell::child()");
        ASSERT(index < num_children(), "Index " + stringify(index) + " must not exceed number of children "
               + stringify(num_children()) + ".");
        return _children[index];
      }


      /**
      * \brief returns refinement level this cell has been created at
      *
      * \return refinement level this cell has been created at
      */
      inline unsigned char refinement_level() const
      {
        CONTEXT("BaseMesh::Cell::refinement_level()");
        return _refinement_level;
      }


      /**
      * \brief returns pointer to subdivision data object
      *
      * \return pointer to subdivision data object
      */
      inline SubdivisionData<cell_dim_, space_dim_, world_dim_>* subdiv_data() const
      {
        CONTEXT("BaseMesh::Cell::subdiv_data()");
        return _subdiv_data;
      }


      /**
      * \brief sets pointer subdivision data object
      *
      * Assumes that the subdivision data pointer _subdiv_data is nullptr.
      *
      * \param[in] subdiv_data
      * pointer subdivision data object
      */
      inline void set_subdiv_data(SubdivisionData<cell_dim_, space_dim_, world_dim_>* subdiv_data)
      {
        CONTEXT("BaseMesh::Cell::set_subdiv_data()");
        try
        {
          if(_subdiv_data != nullptr)
          {
            throw InternalError("SubdivisionData object of cell " + this->print_index() + " is not null!");
          }
          _subdiv_data = subdiv_data;
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      }


      /**
      * \brief returns whether this cell is active (i.e. whether it does not have children)
      *
      * \return flag whether this cell is active
      *
      * \todo Falsch fuer Zellen niedrigerer Dimension! Bsp.: zwei benachbarte Quads
      * verschiedenen refinement levels, naemlich 0 und 1 --> gemeinsame Kante (Level 0) ist aktiv, obwohl sie Kinder
      * (Level 1) hat! Folge: Wir muessen den Aktivitaetsgrad entweder tatsaechlich als Variable abspeichern, oder
      * verzichten fuer Zellen niedrigerer Dimension ganz darauf. Das haette zur Folge, dass beim Verschicken des
      * aktuellen Gitters an die Rechenprozesse nicht einfach alle Zellen niedrigerer Dimension durchlaufen werden
      * koennen und alle aktiven als "zu verschicken" markiert werden koennen. Stattdessen muesste man ueber alle
      * aktiven cells hoechster Dimension die darin enthaltenen subcells als "zu verschicken" markieren.
      */
      inline bool active() const
      {
        CONTEXT("BaseMesh::Cell::active()");
        return (num_children() == 0);
      }


      /// pure virtual function for subdividing the cell
      virtual void subdivide(SubdivisionData<cell_dim_, space_dim_, world_dim_>*) = 0;


      /// pure virtual function for validating the cell
      virtual void validate(std::ostream&) const = 0;


      /// pure virtual function for printing the cell
      virtual void print(std::ostream&) const = 0;


      /**
      * \brief validates parent-child relations of this cell
      *
      * \param[in,out] stream
      * stream validation info is written into
      */
      void validate_history(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Cell::validate_history()");
        try
        {
          stream << "Validating history of cell " << this->print_index() << "..." << std::endl;

          String s = "Error in cell " + this->print_index() + ": ";

          // if this cell has children, do some respective validations
          if(!active())
          {
            for(unsigned char ichild(0) ; ichild < this->num_children() ; ++ichild)
            {
              // check whether parent is set correctly
              if(child(ichild)->parent() != this)
              {
                s += "Parent-child relation to the " + stringify((int)ichild) + "-th child not correctly set!\n";
                throw new InternalError(s);
              }
              // check whether refinement levels are correct
              if(this->child(ichild)->refinement_level() != this->refinement_level()+1)
              {
                s += "Refinement level of this cell or its " + stringify((int)ichild)
                     + "-th child is not correctly set!\n";
                throw new InternalError(s);
              }
            }
          }
          // if this cell has a parent, do some further validations
          if(parent() != nullptr)
          {
            if(parent()->active())
            {
              s += "Its parent is active, i.e. has no children, which is not possible!";
              throw new InternalError(s);
            }
            // validate that this cell can actually be found among the children of its parent
            bool child_found = false;
            for(unsigned char ichild(0) ; ichild < parent()->num_children() ; ++ichild)
            {
              if(parent()->child(ichild) == this)
              {
                child_found = true;
                break;
              }
            }
            if(!child_found)
            {
              s += "It can not be found among the children of its parent!";
              throw new InternalError(s);
            }
          }
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      }


      /**
      * \brief prints parent-child relations of this cell to the given stream
      *
      * \param[in,out] stream
      * stream to write into
      */
      inline void print_history(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Cell::print_history()");
        stream << "[parent: ";
        if (parent() != nullptr)
        {
          parent()->print_index(stream);
        }
        else
        {
          stream << "-";
        }
        if (num_children() > 0)
        {
          stream << ", children: ";
          child(0)->print_index(stream);
          for(int i(1) ; i < num_children() ; ++i)
          {
            stream << ", ";
            child(i)->print_index(stream);
          }
        }
        else
        {
          stream << ", children: -";
        }
        stream << "]";
      }
    }; // class Cell
  } // namespace BaseMesh
} // namespace FEAST

#endif // KERNEL_BASE_MESH_CELL_HPP
