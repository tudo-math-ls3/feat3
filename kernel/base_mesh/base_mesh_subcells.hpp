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

namespace FEAST
{
  namespace BaseMesh
  {
    /**
    * \brief class providing containers and functions for base mesh subcells
    *
    * To be specialised w.r.t. template parameter space_dim_.
    *
    * \author Hilmar Wobker
    */
    template<unsigned char space_dim_, unsigned char world_dim_>
    class Subcells
    {
    }

    /**
    * \brief class providing containers and functions for 2D base mesh subcells, i.e. for edges
    *
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class Subcells<2, world_dim_>
    {
      /// shortcut to save typing of template parameters
      typedef Cell<1, 2, world_dim_> Cell_1D_;


    private:

      /* *****************
      * member variables *
      *******************/
      /// array of edges
      std::vector<Cell_1D_*> _edges;


      /* **********
      * functions *
      ************/
      /**
      * \brief templated function to remove vertices, edges and cells from the corresponding vector
      *
      * The item is swapped to the end of the list and then deleted and removed.
      * TODO: This code is replicated in all base_mesh_xD classes since there is no generic class to inherit from
      */
      template<typename T_>
      inline void _remove(std::vector<T_>& v, T_ item)
      {
        assert(item->index() < v.size());
        v[item->index()] = v.back();
        v[item->index()]->set_index(item->index());
        v[v.size()-1] = item;
        item->set_index(v.size()-1);
        delete item;
        v.pop_back();
      }


      /// adds given edge to container and sets its index
      inline void _add(Cell_1D_* e)
      {
        _edges.push_back(e);
        e->set_index(_edges.size()-1);
      }


      /// deletes given edge
      inline void _remove(Cell_1D_* e)
      {
        _remove<Cell_1D_*>(_edges, e);
      }



    public:

      /* ***************************
      * constructors & destructors *
      *****************************/

      /**
      * \brief  default CTOR for a 2D base mesh
      *
      * Ceates a base mesh basing on the provided mesh file.
      *
      * \param[in] file_name
      * name of the mesh file
      */
      Subcells<2, world_dim_>()
        : _edges(nullptr)
      {
      }

      /// default destructor
      ~BaseMesh2D()
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

        // delete all edges and their associated information
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
      inline global_index_t num_edges() const
      {
        // TODO: potentiell falsch, auch Kanten koennen inaktiv sein und duerfen dann beim Transfer zu den
        // Rechenprozessen nicht mitgezaehlt werden!
        return _edges.size();
      }


      /// returns number of active cells in this mesh (this is potentially expensive)
      inline global_index_t num_active_cells() const
      {
        // TODO: nochmal implementieren, wenn inaktive immer schoen ans Ende geschoben werden und es einen index gibt,
        // der die Position des letzten aktiven merkt
        global_index_t counter = 0;
        for(global_index_t i(0) ; i < num_cells() ; ++i)
        {
          if(cell(i)->active())
          {
            counter++;
          }
        }
        return counter;
      }


      /// returns edge at given index
      inline Cell_1D_* edge(global_index_t const index)
      {
        assert(index < num_edges());
        return _edges[index];
      }


      /// returns graph pointer #_graph
      inline Graph* graph() const
      {
        return _graph;
      }

      inline void add_created_subcells(SubdivisionData<2, 2, world_dim_>* subdiv_data)
      {
        for(unsigned int i(0) ; i < subdiv_data->created_edges.size() ; ++i)
        {
          _add(subdiv_data->created_edges[i]);
        }
      }


      /// validates the subcells
      void validate() const
      {
//        std::cout << "Validating subcells..." << std::endl;
//        for(unsigned int icell(0) ; icell < _cells.size() ; ++icell)
//        {
//          cell(icell)->validate();
//        }
//        std::cout << "...done!" << std::endl;
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
        stream << _vertices.size() <<" vertices:" << std::endl;
        for(unsigned int ivert(0) ; ivert < _vertices.size() ; ++ivert)
        {
          _vertices[ivert]->print(stream);
          stream << std::endl;
        }
        stream << _edges.size() << " edges" << std::endl;
        for(unsigned int iedge(0) ; iedge < _edges.size() ; ++iedge)
        {
          _edges[iedge]->print(stream);
          stream << std::endl;
        }
      }
    }; // class Subcells
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_SUBCELLS_HPP
