#pragma once
#ifndef KERNEL_BASE_MESH_1D_HPP
#define KERNEL_BASE_MESH_1D_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_vertex.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_cell_1d_edge.hpp>

namespace FEAST
{
  namespace BaseMesh
  {
    /**
    * \brief 1D base mesh
    *
    * \author Dominik Goeddeke
    * \author Hilmar Wobker
    */
// COMMENT_HILMAR: Um Code-Redundanz zu vermeiden, koennten wir ueberlegen, eine Klasse BaseMesh einzufuehren,
// die als Elternklasse fuer BaseMesh1D/2D/3D dient. Und/Oder (aehnlich wie bei Cell) dimension-abhaengige Interfaces
// definieren.
    template<unsigned char world_dim_>
    class BaseMesh1D
    {
      /// shortcuts various cell types to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;
      typedef Edge<1, world_dim_> Edge_;
      typedef Cell<1, 1, world_dim_> Cell_;


    private:

      /* *****************
      * member variables *
      *******************/
      /// array of vertices
      std::vector<Vertex_*> _vertices;

      /// array of cells
      std::vector<Cell_*> _cells;

      /* **********
      * functions *
      ************/
      /**
      * \brief templated function to remove vertices and cells from the corresponding vector
      *
      * The item is swapped to the end of the list and then deleted and removed.
      * TODO: This code is replicated in all base_mesh_xD classes since there is no generic class to inherit from
      */
      template<typename T_>
      inline void remove(std::vector<T_>& v, T_ item)
      {
        assert(item->index() < v.size());
        v[item->index()] = v.back();
        v[item->index()]->set_index(item->index());
        v[v.size()-1] = item;
        item->set_index(v.size()-1);
        delete item;
        v.pop_back();
      }

    public:

      /* ***************************
      * constructors & destructors *
      *****************************/
      /// default CTOR, currently generates a test mesh
      BaseMesh1D()
      {
        /* Base mesh example consisting of three quads and two tris:
        *    v0---e0---v1---e1---v2---e2---v3
        */

        // create the four vertices
        // v0
        Vertex_* v = new Vertex_();
        v->set_coord(0, 0.0);
        add(v);

        // v1
        v = new Vertex_();
        v->set_coord(0, 1.0);
        add(v);

        // v2
        v = new Vertex_();
        v->set_coord(0, 2.0);
        add(v);

        // v3
        v = new Vertex_();
        v->set_coord(0, 3.0);
        add(v);

        // create the three edges (=cells)
        // e0
        Edge_* e = new Edge_(_vertices[0], _vertices[1], 0);
        add(e);
        // e1
        e = new Edge_(_vertices[1], _vertices[2], 0);
        add(e);
        // e2
        e = new Edge_(_vertices[2], _vertices[3], 0);
        add(e);

        // set neighbourhood information (emulated file parser part 2)
        _cells[0]->add_neighbour(SDIM_VERTEX, 1, _cells[1]);
        _cells[1]->add_neighbour(SDIM_VERTEX, 0, _cells[0]);
        _cells[1]->add_neighbour(SDIM_VERTEX, 1, _cells[2]);
        _cells[2]->add_neighbour(SDIM_VERTEX, 0, _cells[1]);
      }

      /// default destructor
      ~BaseMesh1D()
      {
        // delete all cells and their associated information
        // (pop_back calls destructor of the element being removed, so do not use an iterator because fiddling about with
        // the std::vector invalidates it. Also, pop_back calls default DTOR of BMI*, so we have to manually call
        // delete here)
        while (!_cells.empty())
        {
          delete _cells.back();
          _cells.pop_back();
        }

        // delete all vertices and their associated information
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
      inline global_index_t num_vertices() const
      {
        return _vertices.size();
      }

      /// returns number of cells in this mesh (this is potentially expensive)
      inline global_index_t num_cells() const
      {
        // TODO: nochmal implementieren, wenn inaktive immer schoen ans Ende geschoben werden und es einen index gibt,
        // der die Position des letzten aktiven merkt
        global_index_t counter = 0;
        for (global_index_t i(0) ; i < _cells.size() ; ++i)
        {
          if (_cells[i]->active())
          {
            counter++;
          }
        }
        return counter;
      }

      /// returns vertex at given index
      inline Vertex_* vertex(global_index_t const index)
      {
        assert(index < _vertices.size());
        return _vertices[index];
      }

      /// returns cell at given index
      inline Cell_* cell(global_index_t const index) const
      {
        assert(index < _cells.size());
        return _cells[index];
      }

      /// adds given vertex to base mesh and sets its index
      inline void add(Vertex_* v)
      {
        _vertices.push_back(v);
        v->set_index(_vertices.size()-1);
      }

      /// adds given cell to base mesh and sets its index
      inline void add(Cell_* c)
      {
        _cells.push_back(c);
        c->set_index(_cells.size()-1);
      }

      /// deletes given vertex
      inline void remove(Vertex_* v)
      {
        remove<Vertex_*>(_vertices, v);
      }

      /// deletes given cell
      inline void remove(Cell_* c)
      {
        remove<Cell_*>(_cells, c);
      }


      inline void add_created_items(SubdivisionData<1, 1, world_dim_>* subdiv_data)
      {
        add(subdiv_data->created_vertex);
        for(unsigned int i(0) ; i < subdiv_data->created_cells.size() ; ++i)
        {
          add(subdiv_data->created_cells[i]);
        }
      }


      /// validates the base mesh and all of its cells
      void validate() const
      {
        std::cout << "Validating cells..." << std::endl;
// COMMENT_HILMAR: I thought, iterators were simple to use, i.e. like this:
//        for (std::vector<Cell_*>::iterator it = _cells.begin() ; it != _cells.end() ; ++it)
//        {
//          it->validate();
//        }
// But that seems not to be the case... :-(
// So, I use a standard counter for the for loop.
        for (unsigned int icell(0) ; icell < _cells.size() ; ++icell)
        {
          cell(icell)->validate();
        }
        std::cout << "...done!" << std::endl;

        // COMMENT_HILMAR: add further validations...
      }


      /**
      * \brief prints this BaseMesh to the given ostream.
      *
      * According to http://www.cplusplus.com/reference/iostream, ostream is the superclass for both stdout, stderr and
      * a (stream associated with) an arbitrary already opened ASCII file, so this routine can be used for logging
      * and and printf()-style debugging at the same time. Neat.
      *
      * \param[in] stream
      *            Stream to dump this BaseMesh into.
      */
      void print(std::ostream& stream)
      {
        stream << "---------------------------------------------------" << std::endl;
        stream << "|               DUMPING BASE MESH                  " << std::endl;
        stream << "---------------------------------------------------" << std::endl;
        stream << num_vertices() <<" Vertices:" << std::endl;
        for (unsigned int ivert(0) ; ivert < _vertices.size() ; ++ivert)
        {
          _vertices[ivert]->print(stream);
          stream << std::endl;
        }
        stream << num_cells() << " Cells" << std::endl;
        for (unsigned int icell(0) ; icell < _cells.size() ; ++icell)
        {
          _cells[icell]->print(stream);
          stream << std::endl;
        }
        stream << "---------------------------------------------------" << std::endl;
      }
    }; // class BaseMesh1D
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_1D_HPP
