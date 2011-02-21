#pragma once
#ifndef KERNEL_BASE_MESH_BM_HPP
#define KERNEL_BASE_MESH_BM_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh/vertex.hpp>
#include <kernel/base_mesh/cell.hpp>
#include <kernel/base_mesh/subcells.hpp>
#include <kernel/graph.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    // forward declaration
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class FileParser;

    /**
    * \brief base mesh
    *
    * long description missing
    *
    * \note This class is not called 'BaseMesh' since this would collide with the namespace.
    *
    * \author Dominik Goeddeke
    * \author Hilmar Wobker
    */
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class BM
    {
      /// shortcuts to save typing of template parameters
      typedef Cell<space_dim_, space_dim_, world_dim_> Cell_;

      // Parsing of the mesh files is outsourced to class FileParser. Make this class friend such that it has
      // access to all members of this class.
      friend class FileParser<space_dim_, world_dim_>;


    private:

      /* *****************
      * member variables *
      *******************/
      /// object containing the vectors of subcells
      Subcells<space_dim_, space_dim_, world_dim_> _subcells;

      /// vector of cells (of full dimension)
      std::vector<Cell_*> _cells;

      /// graph describing the connectivity of the base mesh
// COMMENT_HILMAR: probably not really needed... just temporarily placed here, later we need a graph structure for
// the connectivity of process patches.
      Graph* _graph;

      /* **********
      * functions *
      ************/
      /**
      * \brief templated function to remove items from the given vector
      *
      * The item is swapped to the end of the list and then deleted and removed.
      */
//COMMENT_HILMAR: Move this code to some parent class BaseMeshAux or similar.
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

      /// adds given cell to base mesh and sets its index
      inline void _add(Cell_* c)
      {
        _cells.push_back(c);
        c->set_index(_cells.size()-1);
      }

      /// deletes given cell
      inline void _remove(Cell_* c)
      {
        _remove<Cell_*>(_cells, c);
      }




    public:

      /* ***************************
      * constructors & destructors *
      *****************************/

      /**
      * \brief  default CTOR for a base mesh
      */
      BM()
        : _cells(nullptr),
          _graph(nullptr)
      {
      }

      /// default destructor
      ~BM()
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
      }


      /* *********************************
      * getters & setters & manipulators *
      ************************************/

      /// returns number of cells in this mesh (including inactive ones)
      inline global_index_t num_cells() const
      {
        return _cells.size();
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


      /// returns cell at given index
      inline Cell_* cell(global_index_t const index) const
      {
        assert(index < num_cells());
        return _cells[index];
      }


      /// returns graph pointer #_graph
      inline Graph* graph() const
      {
        return _graph;
      }


      /// adds (sub)cells created during subdivision to the corresponding (sub)cell vectors
      inline void add_created_items(SubdivisionData<space_dim_, space_dim_, world_dim_>* subdiv_data)
      {
        _subcells.add_created_subcells(subdiv_data);
        for(unsigned int i(0) ; i < subdiv_data->created_cells.size() ; ++i)
        {
          _add(subdiv_data->created_cells[i]);
        }
      }


      /**
      * \brief sets numbers (not indices!) in all active cells (has to be called after subdivision of a cell etc.)
      *
      * \note Stupidly runs linearly through the vector of cells and overwrites existing numbers.
      * A different numbering strategy (try to keep existing numbers, fill gaps, ...) may be more clever... this will
      * be modified later.
      */
      inline void set_cell_numbers() const
      {
        global_index_t counter = 0;
        for(global_index_t i(0) ; i < num_cells() ; ++i)
        {
          if(cell(i)->active())
          {
            cell(i)->set_number(counter);
            counter++;
          }
          else
          {
            cell(i)->unset_number();
          }
        }
      }


      /// creates connectivity graph from information stored in this base mesh
// COMMENT_HILMAR: This is just an intermediate solution to artificially connect the base mesh to the load balancer.
// I.e., we assume here that each process receives exactly one BMC and that the connectivity graph relevant for the
// load balancer actually is the connectivity graph of the base mesh. Later, there will be the matrix patch layer and
// the process patch layer, which both have their own connectivity structure. The load balancer then actually needs the
// connectivity graph of the process patch layer. We also do not distinguish between edge and vertex neighbours here.
      void create_graph()
      {
        global_index_t n_active_cells = num_active_cells();
        // allocate index array
        unsigned int* index = new unsigned int[n_active_cells + 1];

        // graph data structure is filled by two sweeps through the cell list
        // first sweep: count neighbours of each cell, and maintain running total to fill index array
        // treat last index entry separately because cell array has one less entry than index array
        unsigned int num_neighbours_so_far = 0;
        // counter for active cells
        global_index_t ipos(0);
        for (global_index_t icell=0 ; icell < num_cells() ; ++icell)
// TODO: wir brauchen einen iterator fuer aktive Zellen!
        {
          if(cell(icell)->active())
          {
            // set neighbours counted so far
            index[ipos] = num_neighbours_so_far;
//std::cout << "Setting index[" << ipos << "] = " << num_neighbours_so_far << std::endl;
            // count neighbours (here: edge neighbours and vertex neighbours)
            for (unsigned char sdim(0) ; sdim < space_dim_ ; ++sdim)
            {
              num_neighbours_so_far += cell(icell)->num_neighbours_subdim((subdim)sdim);
            }
            ++ipos;
          }
        }
        index[n_active_cells] = num_neighbours_so_far;
//std::cout << "Setting index[" << n_active_cells << "] = " << num_neighbours_so_far << std::endl;

        // second sweep through data structure
        // second sweep adds actual neighbour cell numbers in the appropriate places into array neighbours
        // again, treat last loop instance separately
        unsigned int* neighbours = new unsigned int[index[n_active_cells]];
        num_neighbours_so_far = 0;
        for (global_index_t icell=0 ; icell < n_active_cells ; icell++)
// TODO: wir brauchen einen iterator fuer aktive Zellen!
        {
          Cell_* c = cell(icell);
          if (c->active())
          {
            for (unsigned char sdim(0) ; sdim < space_dim_ ; ++sdim)
            {
              for(unsigned char item(0) ; item < c->num_subitems_per_subdim((subdim)sdim) ; ++item)
              {
                std::vector<Cell_*>& neigh_cells = c->neighbours_item((subdim)sdim, item);
                for(unsigned char k(0) ; k < neigh_cells.size() ; ++k)
                {
                  neighbours[num_neighbours_so_far] = neigh_cells[k]->number();
//std::cout << "neighbours[" << num_neighbours_so_far << "] = " << neighbours[num_neighbours_so_far] << std::endl;
                  ++num_neighbours_so_far;
                }
              }
            }
          }
        }

        // now, create graph object
        // temporarily, do not distinguish edge neighbours and diagonal neighbours
        if (_graph != nullptr)
        {
          delete _graph;
          _graph = nullptr;
        }
        _graph = new Graph(num_cells(), index, neighbours);
      }


      /// validates the base mesh and all of its cells
      void validate() const
      {
        std::cout << "Validating cells..." << std::endl;
        for(unsigned int icell(0) ; icell < _cells.size() ; ++icell)
        {
          cell(icell)->validate();
        }
        std::cout << "...done!" << std::endl;
//        _subcells.validate();
        // COMMENT_HILMAR: add further validations...
      }


      /**
      * \brief prints this base mesh to the given ostream.
      *
      * According to http://www.cplusplus.com/reference/iostream, ostream is the superclass for both stdout, stderr and
      * a (stream associated with) an arbitrary already opened ASCII file, so this routine can be used for logging
      * and and printf()-style debugging at the same time. Neat.
      *
      * \param[in] stream
      *            Stream to dump this base mesh into.
      */
      void print(std::ostream& stream)
      {
        stream << "---------------------------------------------------" << std::endl;
        stream << "|               DUMPING BASE MESH                  " << std::endl;
        stream << "---------------------------------------------------" << std::endl;
        _subcells.print(stream);
        stream << _cells.size() << " cells" << std::endl;
        for(unsigned int icell(0) ; icell < _cells.size() ; ++icell)
        {
          _cells[icell]->print(stream);
          stream << std::endl;
        }
        stream << "---------------------------------------------------" << std::endl;
      }
    }; // class BM
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_BASE_MESH_HPP
