#pragma once
#ifndef KERNEL_BASE_MESH_SUBCELLS_HPP
#define KERNEL_BASE_MESH_SUBCELLS_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh/vertex.hpp>
#include <kernel/base_mesh/cell.hpp>
#include <kernel/base_mesh/cell_1d_edge.hpp>
#include <kernel/base_mesh/cell_2d_quad.hpp>
#include <kernel/base_mesh/cell_2d_tri.hpp>
#include <kernel/graph.hpp>

/// FEAST namespace
namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
  namespace BaseMesh
  {
    // forward declaration
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class FileParser;

    /**
    * \brief class providing containers and functions for base mesh subcells
    *
    * To be specialised w.r.t. template parameter space_dim_.
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
    class Subcells
    {
    };


    /**
    * \brief class providing containers and functions for subcells of 1D cells, i.e. for vertices
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
    class Subcells<1, space_dim_, world_dim_>
    {
      /// shortcut to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;

      // Parsing of the mesh files is outsourced to class FileParser. Make this class friend such that it has
      // access to all members of this class.
      friend class FileParser<space_dim_, world_dim_>;


    private:

      /* *****************
      * member variables *
      *******************/
      /// array of vertices
      std::vector<Vertex_*> _vertices;


    protected:

      /* **********
      * functions *
      ************/
      /**
      * \brief templated function to remove items from the given vector
      *
      * The item is swapped to the end of the list and then deleted and removed.
      *
      * \tparam T_
      * class type of the item and the vector items
      *
      * \param[in] v
      * vector from which the item is to be removed
      *
      * \param[in] item
      * item to be removed from the vector
      */
//COMMENT_HILMAR: Move this code to some auxiliary class (it is also needed in bm.hpp)
      template<typename T_>
      inline void _remove_item(std::vector<T_>& v, T_ item)
      {
        CONTEXT("BaseMesh::Subcells::_remove_item()");
        assert(item->index() < v.size());
        v[item->index()] = v.back();
        v[item->index()]->set_index(item->index());
        v[v.size()-1] = item;
        item->set_index(v.size()-1);
        delete item;
        v.pop_back();
      }


      /* **********
      * functions *
      ************/
      /**
      * \brief adds given vertex to the base mesh's subcells and sets its index
      *
      * \param[in] v
      * vertex to be added to the vector #_vertices
      */
      inline void _add(Vertex_* v)
      {
        CONTEXT("BaseMesh::Subcells::_add()");
        _vertices.push_back(v);
        v->set_index(_vertices.size()-1);
      }


      /**
      * \brief deletes given vertex
      *
      * Deleting the given vertex is an O(1) operation since its position within the vector is given by its index.
      *
      * \param[in] v
      * vertex to be removed from the vector #_vertices
      */
      inline void _remove(Vertex_* v)
      {
        CONTEXT("BaseMesh::Subcells::_remove()");
        _remove_item<Vertex_*>(_vertices, v);
      }


    public:

      /* ***************************
      * constructors & destructors *
      *****************************/

      /// default CTOR
      Subcells()
        : _vertices(nullptr)
      {
        CONTEXT("BaseMesh::Subcells::Subcells()");
      }

      /// DTOR
      virtual ~Subcells()
      {
        CONTEXT("BaseMesh::Subcells::~Subcells()");
        // delete all vertices and their associated information
        // (pop_back calls destructor of the element being removed, so do not use an iterator because fiddling about
        // with the std::vector invalidates it. Also, pop_back calls default DTOR of Vertex_*, so we have to manually
        // call delete here)
        while (!_vertices.empty())
        {
          delete _vertices.back();
          _vertices.pop_back();
        }
      }


      /* *********************************
      * getters & setters & manipulators *
      ************************************/

      /// returns number of vertices in this mesh
      inline index_t_glob num_vertices() const
      {
        CONTEXT("BaseMesh::Subcells::num_vertices()");
        return _vertices.size();
      }


      /// returns vertex at given index
      inline Vertex_* vertex(index_t_glob const index)
      {
        CONTEXT("BaseMesh::Subcells::vertex()");
        assert(index < num_vertices());
        return _vertices[index];
      }


      /**
      * \brief adds vertex created during subdivision of the corresponding edge
      *
      * \param[in] subdiv_data
      * object containing pointer to the vertex to be added
      */
      inline void add_created_subcells(SubdivisionData<1, space_dim_, world_dim_>* subdiv_data)
      {
        CONTEXT("BaseMesh::Subcells::add_created_subcells()");
        _add(subdiv_data->created_vertex);
      }


      /**
      * \brief prints this subcells object to the given ostream.
      *
      * According to http://www.cplusplus.com/reference/iostream, ostream is the superclass for both stdout, stderr and
      * a (stream associated with) an arbitrary already opened ASCII file, so this routine can be used for logging
      * and and printf()-style debugging at the same time. Neat.
      *
      * \param[in,out] stream
      * stream to dump the subcells into
      */
      void print(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Subcells::print()");
        stream << _vertices.size() << " vertices" << std::endl;
        for(unsigned int ivert(0) ; ivert < _vertices.size() ; ++ivert)
        {
          _vertices[ivert]->print(stream);
          stream << std::endl;
        }
      }
    }; // class Subcells<1, space_dim_, world_dim_>




    /**
    * \brief class providing containers and functions for subcells of 2D cells, i.e. for vertices and edges
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
    class Subcells<2, space_dim_, world_dim_>
     : public Subcells<1, space_dim_, world_dim_>
    {
      /// shortcut to save typing of template parameters
      typedef Cell<1, space_dim_, world_dim_> Cell_1D_;

      /// shortcut to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;

      // Parsing of the mesh files is outsourced to class FileParser. Make this class friend such that it has
      // access to all members of this class.
      friend class FileParser<space_dim_, world_dim_>;


    private:

      /* *****************
      * member variables *
      *******************/
      /// array of edges
      std::vector<Cell_1D_*> _edges;


    protected:

      /* **********
      * functions *
      ************/
      /// adds given vertex to the base mesh's subcells and sets its index (wrapper calling function in parent class)
      inline void _add(Vertex_* v)
      {
        CONTEXT("BaseMesh::Subcells::_add()");
        Subcells<1, space_dim_, world_dim_>::_add(v);
      }


      /// deletes given vertex (wrapper calling function in parent class)
      inline void _remove(Vertex_* v)
      {
        CONTEXT("BaseMesh::Subcells::_remove()");
        Subcells<1, space_dim_, world_dim_>::_remove(v);
      }


      /**
      * \brief adds given edge to the base mesh's subcells and sets its index
      *
      * \param[in] e
      * edge to be added to the vector #_edges
      */
      inline void _add(Cell_1D_* e)
      {
        CONTEXT("BaseMesh::Subcells::_add()");
        _edges.push_back(e);
        e->set_index(_edges.size()-1);
      }


      /**
      * \brief deletes given edge
      *
      * Deleting the given edge is an O(1) operation since its position within the vector is given by its index.
      *
      * \param[in] e
      * edge to be removed from the vector #_edges
      */
      inline void _remove(Cell_1D_* e)
      {
        CONTEXT("BaseMesh::Subcells::_remove()");
        // the template keyword is necessary here, otherwise the compiler cannot parse the expression
        this->template _remove_item<Cell_1D_*>(_edges, e);
      }


    public:

      /* ***************************
      * constructors & destructors *
      *****************************/

      ///default CTOR
      Subcells()
        : Subcells<1, space_dim_, world_dim_>(),
        _edges(nullptr)
      {
        CONTEXT("BaseMesh::Subcells::Subcells()");
      }


      /// DTOR (automatically virtual since DTOR of base class Item is virtual)
      ~Subcells()
      {
        CONTEXT("BaseMesh::Subcells::~Subcells()");
        // delete all edges and their associated information
        // (pop_back calls destructor of the element being removed, so do not use an iterator because fiddling about
        // with the std::vector invalidates it. Also, pop_back calls default DTOR of Cell_1D_*, so we have to manually
        // call delete here)
        while (!_edges.empty())
        {
          delete _edges.back();
          _edges.pop_back();
        }
      }


      /* *********************************
      * getters & setters & manipulators *
      ************************************/

      /// returns number of edges in this mesh (including inactive ones)
      inline index_t_glob num_edges() const
      {
        CONTEXT("BaseMesh::Subcells::num_edges()");
        // TODO: potentiell falsch, auch Kanten koennen inaktiv sein und duerfen dann beim Transfer zu den
        // Rechenprozessen nicht mitgezaehlt werden!
        return _edges.size();
      }


      /// returns edge at given index
      inline Cell_1D_* edge(index_t_glob const index)
      {
        CONTEXT("BaseMesh::Subcells::edge()");
        assert(index < num_edges());
        return _edges[index];
      }


      /**
      * \brief adds subcells created during subdivision to the corresponding subcell vectors
      *
      * \param[in] subdiv_data
      * object containing pointers to the subcells to be added
      */
      inline void add_created_subcells(SubdivisionData<2, space_dim_, world_dim_>* subdiv_data)
      {
        CONTEXT("BaseMesh::Subcells::add_created_subcells()");
        for(unsigned int i(0) ; i < subdiv_data->created_vertices.size() ; ++i)
        {
          _add(subdiv_data->created_vertices[i]);
        }
        for(unsigned int i(0) ; i < subdiv_data->created_edges.size() ; ++i)
        {
          _add(subdiv_data->created_edges[i]);
        }
      }


      /**
      * \brief prints this subcells object to the given ostream.
      *
      * According to http://www.cplusplus.com/reference/iostream, ostream is the superclass for both stdout, stderr and
      * a (stream associated with) an arbitrary already opened ASCII file, so this routine can be used for logging
      * and and printf()-style debugging at the same time. Neat.
      *
      * \param[in,out] stream
      * stream to dump this object into
      */
      void print(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Subcells::print()");
        Subcells<1, space_dim_, world_dim_>::print(stream);
        stream << _edges.size() << " edges" << std::endl;
        for(unsigned int iedge(0) ; iedge < _edges.size() ; ++iedge)
        {
          _edges[iedge]->print(stream);
          stream << std::endl;
        }
      }
    }; // class Subcells<2, space_dim_, world_dim_>




    /**
    * \brief class providing containers and functions for subcells of 3D cells, i.e. for vertices, edges and faces
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
    class Subcells<3, space_dim_, world_dim_>
      : public Subcells<2, space_dim_, world_dim_>
    {
      /// shortcut to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;

      /// shortcut to save typing of template parameters
      typedef Cell<1, space_dim_, world_dim_> Cell_1D_;

      /// shortcut to save typing of template parameters
      typedef Cell<2, space_dim_, world_dim_> Cell_2D_;

      // Parsing of the mesh files is outsourced to class FileParser. Make this class friend such that it has
      // access to all members of this class.
      friend class FileParser<space_dim_, world_dim_>;


    private:

      /// array of faces
      std::vector<Cell_2D_*> _faces;


    protected:
      /* **********
      * functions *
      ************/
      /// adds given vertex to the base mesh's subcells and sets its index (wrapper calling function in parent class)
      inline void _add(Vertex_* v)
      {
        CONTEXT("BaseMesh::Subcells::_add()");
        Subcells<1, space_dim_, world_dim_>::_add(v);
      }


      /// deletes given vertex (wrapper calling function in parent class)
      inline void _remove(Vertex_* v)
      {
        CONTEXT("BaseMesh::Subcells::_remove()");
        Subcells<1, space_dim_, world_dim_>::_remove(v);
      }


      /// adds given edge to the base mesh's subcells and sets its index (wrapper calling function in parent class)
      inline void _add(Cell_1D_* e)
      {
        CONTEXT("BaseMesh::Subcells::_add()");
        Subcells<2, space_dim_, world_dim_>::_add(e);
      }


      /// deletes given edge (wrapper calling function in parent class)
      inline void _remove(Cell_1D_* e)
      {
        CONTEXT("BaseMesh::Subcells::_remove()");
        Subcells<2, space_dim_, world_dim_>::_remove(e);
      }


      /**
      * \brief adds given face to the base mesh's subcells and sets its index
      *
      * \param[in] f
      * face to be added to the vector #_faces
      */
      inline void _add(Cell_2D_* f)
      {
        CONTEXT("BaseMesh::Subcells::_add()");
        _faces.push_back(f);
        f->set_index(_faces.size()-1);
      }


      /**
      * \brief deletes given face
      *
      * Deleting the given face is an O(1) operation since its position within the vector is given by its index.
      *
      * \param[in] f
      * face to be removed from the vector #_faces
      */
      inline void _remove(Cell_2D_* f)
      {
        CONTEXT("BaseMesh::Subcells::_remove()");
        // the template keyword is necessary here, otherwise the compiler cannot parse the expression
        this->template _remove_item<Cell_2D_*>(_faces, f);
      }


    public:

      ///default CTOR
      Subcells()
        : Subcells<2, space_dim_, world_dim_>(),
          _faces(nullptr)
      {
        CONTEXT("BaseMesh::Subcells::Subcells()");
      }


      /// DTOR (automatically virtual since DTOR of base class Item is virtual)
      ~Subcells()
      {
        CONTEXT("BaseMesh::Subcells::~Subcells()");
//COMMENT_HILMAR: nochmal nachgucken, ob bzw. wie der DTOR von Subcells<2, ...> aufgerufen wird.
        // delete all faces and their associated information
        // (pop_back calls destructor of the element being removed, so do not use an iterator because fiddling about
        // with the std::vector invalidates it. Also, pop_back calls default DTOR of Cell_2D_*, so we have to manually
        // call delete here)
        while (!_faces.empty())
        {
          delete _faces.back();
          _faces.pop_back();
        }
      }


      /// returns number of faces in this mesh (including inactive ones)
      inline index_t_glob num_faces() const
      {
        CONTEXT("BaseMesh::Subcells::num_faces()");
        // TODO: potentiell falsch, auch Faces koennen inaktiv sein und duerfen dann beim Transfer zu den Rechenprozessen
        // nicht mitgezaehlt werden!
        return _faces.size();
      }


      /// returns face at given index
      inline Cell_2D_* face(index_t_glob const index)
      {
        CONTEXT("BaseMesh::Subcells::face()");
        assert(index < _faces.size());
        return _faces[index];
      }


      /**
      * \brief adds subcells created during subdivision to the corresponding subcell vectors
      *
      * \param[in] subdiv_data
      * object containing pointers to the subcells to be added
      */
      inline void add_created_subcells(SubdivisionData<3, space_dim_, world_dim_>* subdiv_data)
      {
        CONTEXT("BaseMesh::Subcells::add_created_subcells()");
        for(unsigned int i(0) ; i < subdiv_data->created_vertices.size() ; ++i)
        {
          _add(subdiv_data->created_vertices[i]);
        }
        for(unsigned int i(0) ; i < subdiv_data->created_edges.size() ; ++i)
        {
          _add(subdiv_data->created_edges[i]);
        }
        for(unsigned int i(0) ; i < subdiv_data->created_faces.size() ; ++i)
        {
          _add(subdiv_data->created_faces[i]);
        }
      }


      /**
      * \brief prints this subcells object to the given ostream.
      *
      * According to http://www.cplusplus.com/reference/iostream, ostream is the superclass for both stdout, stderr and
      * a (stream associated with) an arbitrary already opened ASCII file, so this routine can be used for logging
      * and and printf()-style debugging at the same time. Neat.
      *
      * \param[in,out] stream
      * stream to dump this object into
      */
      void print(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Subcells::print()");
        Subcells<2, space_dim_, world_dim_>::print(stream);
        stream << _faces.size() << " faces" << std::endl;
        for (unsigned int iface(0) ; iface < _faces.size() ; ++iface)
        {
          _faces[iface]->print(stream);
          stream << std::endl;
        }
      }
    }; // class Subcells<3, space_dim_, world_dim_>
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_SUBCELLS_HPP
