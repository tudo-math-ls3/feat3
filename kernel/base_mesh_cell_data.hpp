#pragma once
#ifndef KERNEL_BASE_MESH_CELL_DATA_HPP
#define KERNEL_BASE_MESH_CELL_DATA_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_item.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    /**
    * \brief keywords for the subdimensions
    *
    * See CellData for explanation of the term subdimension.
    * The values serve as array indics, so must not be changed.
    */
    enum subdim
    {
      /// subdimension 0 = vertex
      SDIM_VERTEX = 0,
      /// subdimension 1 = edge
      SDIM_EDGE = 1,
      /// subdimension 2 = face
      SDIM_FACE = 2
    };

    /// forward declaration of class BaseMesh::Cell
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Cell;


    /**
    * \brief emtpy cell data class definition for cells with dimension smaller than space dimension
    *
    * Avoid instatiation of cell-specific data in shapes with dimension smaller than space dimension.
    * The class is only implemented for cell_dim_ = space_dim_ (see below).
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    * \author Peter Zajac
    */
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class CellData
      : public Item
    {
    private:

    public:
      /// dummy function called by cells with dimension smaller than space dimension
      inline void _init_neighbours(unsigned char array_size, unsigned char num_subcells_per_subdim[])
      {
        // do nothing here
      }

      /// dummy function called by cells with dimension smaller than space dimension
      inline void add_neighbour(
        subdim subdim,
        unsigned char item,
        Cell<cell_dim_, space_dim_, world_dim_>* neighbour)
      {
        // do nothing here
      }

      /// dummy function called by cells with dimension smaller than space dimension
      inline void check_neighbourhood() const
      {
      }

      /// dummy function called by cells with dimension smaller than space dimension
      inline void print(std::ostream& stream)
      {
        // do nothing here
      }

    };


    /**
    * \brief class storing cell-specific data like neighbourhood information
    *
    * Only implemented for cell_dim_ = space_dim_.
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    * \author Peter Zajac
    */
    template<
      unsigned char cell_space_dim_,
      unsigned char world_dim_>
    class CellData<cell_space_dim_, cell_space_dim_, world_dim_>
      : public Item
    {
    private:

      /**
      * \brief array of number of subcells per subdimension
      *
      * Definitions of two terms we use:
      * 1) "subdimension" = space dimension minus co-dimension where
      *      in 3D: codim 1 -> faces, codim 2 -> edges, codim 3 -> vertices
      *      in 2D: codim 1 -> edges, codim 2 -> vertices
      *      in 1D: codim 1 -> vertices
      *    Co-dimension 0 is not needed here. Hence, the number of subdimensions is equal to the space dimension.
      * 2) "subcells" = items of lower dimension in a cell, distinguished by subdimension. Examples:
      *      in 3D: a hexa has 6 subcells of subdimension 2 (=faces), 12 subcells of subdimension 1 (=edges)
      *             and 8 subcells of subdimension 0 (=vertices)
      *      in 2D: a tri has 3 subcells of subdimension 1 (=edges) and 3 subcells of subdimension 0 (=vertices)
      *      in 1D: an edge has 2 subcells of subdimension 0 (=vertices)
      */
      unsigned char _num_subcells_per_subdim[cell_space_dim_];

      /**
      * \brief  two-dimensional array of vectors of neighbour cells:
      *
      * First dimension: subdimension at which neighbours are regarded
      *    (0 = vertex neighbours, 1 = edge neighbours (only in 2D/3D), 2 = face neighbours (only in 3D))
      * Second dimension: index of the item in the cell (vertices, edges, faces). Examples:
      * _neighbours[SDIM_VERTEX][1]: vector of vertex neighbours at vertex 1  (SDIM_VERTEX = 0)
      * _neighbours[SDIM_EDGE][3]: vector of edge neighbours at edge 3  (SDIM_EDGE = 1)
      * _neighbours[SDIM_FACE][4]: vector of face neighbours at face 4  (SDIM_FACE = 2)
      */
      std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>* _neighbours[cell_space_dim_];


    protected:

      /// function for initialising the neighbour arrays/vectors
      void _init_neighbours(unsigned char array_size, unsigned char num_subcells_per_subdim[])
      {
        assert(array_size == cell_space_dim_);
        for(int i(0) ; i < cell_space_dim_ ; ++i)
        {
          _num_subcells_per_subdim[i] = num_subcells_per_subdim[i];
          _neighbours[i]
            = new std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>[_num_subcells_per_subdim[i]];
          // TODO: Sind die einzelen std::vector damit schon ordentlich initialisiert? Muesste eigentlich...
        }
      }


    public:

      ~CellData()
      {
        // TODO: deallocate _neighbours !!!
      }

      inline unsigned char num_subcells_per_subdim(subdim subdim)
      {
        assert(subdim < cell_space_dim_);
        return _num_subcells_per_subdim[subdim];
      }

      /// returns one specific neighbour for given subdim, item and index
      inline Cell<cell_space_dim_, cell_space_dim_, world_dim_>* neighbour(
        subdim subdim,
        unsigned char item,
        unsigned char index) const
      {
        assert(subdim < world_dim_);
        assert(item < _num_subcells_per_subdim[subdim]);
        assert(index < _neighbours[subdim][item].size());
        return _neighbours[subdim][item][index];
      }

      /// returns vector of neighbours for given subdimension and given item
      inline std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>& neighbours(
        subdim subdim,
        unsigned char item) const
      {
        return _neighbours[subdim][item];
      }

      /// returns array of vectors of neighbours for given subdimension
      inline std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>* neighbours_subdim(subdim subdim) const
      {
        return _neighbours[subdim];
      }

      /// add neighbour to the vector of neighbours of given subdim and item
      inline void add_neighbour(
        subdim subdim,
        unsigned char item,
        Cell<cell_space_dim_, cell_space_dim_, world_dim_>* neighbour)
      {
        _neighbours[subdim][item].push_back(neighbour);
      }

      /// checks whether neighbours are set correctly
      inline void check_neighbourhood() const
      {
//std::cout << "starting neighbourhood check" << std::endl;
        for(unsigned char sdim(0) ; sdim < cell_space_dim_; ++sdim)
        {
//std::cout << "sdim " <<  (int) sdim << std::endl;
          for(unsigned char item(0) ; item < _num_subcells_per_subdim[sdim] ; ++item)
          {
//std::cout << "item " <<  (int) item << std::endl;
            std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>&
              neighs = neighbours((subdim)sdim, item);
            for(unsigned int n = 0 ; n < neighs.size() ; n++)
            {
//std::cout << "n " <<  n << std::endl;
              if(!neighs[n]->active())
              {
                std::cerr << "Cell ";
                neighs[n]->print_index(std::cerr);
                std::cerr << ", being the subdim-" << (int)sdim << "-neighbour of cell ";
                this->print_index(std::cerr);
                std::cerr << " (at item " << (int)item << ", " << n << "-th pos.), has children"
                          << " which must not be the case!" << std::endl;
                std::cerr << "It seems the neighbours of cell ";
                this->print_index(std::cerr);
                std::cerr << " have not been updated correctly!" << std::endl;
              }
              else
              {
                std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>*
                  neighs_of_neigh = neighs[n]->neighbours_subdim((subdim)sdim);
                bool neighbour_found = false;
                for(unsigned char item_nn = 0 ; item_nn < neighs[n]->num_subcells_per_subdim((subdim)sdim) ; item_nn++)
                {
//std::cout << "item_nn " << item_nn << std::endl;
                  for(unsigned int nn = 0 ; nn < neighs_of_neigh[item_nn].size() ; nn++)
                  {
//std::cout << "nn " << nn << std::endl;
                    if(neighs_of_neigh[item_nn][nn] == this)
                    {
                      std::cout << "subdim-" << (int)sdim << "-neighbourhood between cell ";
                      this->print_index(std::cout);
                      std::cout << " (at item " << (int)item << ", " << n << "-th pos.) and cell ";
                      neighs[n]->print_index(std::cout);
                      std::cout << " (item = " << (int)item_nn << ", " << nn << "-th pos.)" << std::endl;
                      neighbour_found = true;
                      break;
                    }
                  }  // for 0 <= nn < neighs_of_neigh[item_nn].size()
                  if(neighbour_found)
                  {
                    break;
                  }
                } // for 0 <= item_nn < neighs[n]->num_subcells_per_subdim((subdim)sdim)
                if(!neighbour_found)
                {
                  std::cerr << "No neighbour found!" << std::endl;
                  std::cerr << "subdim-" << (int)sdim << "-neighbourhood between cell ";
                  this->print_index(std::cerr);
                  std::cerr << " (at item " << (int)item << ", " << n << "-th pos.) and cell ";
                  neighs[n]->print_index(std::cerr);
                  std::cerr << std::endl;
                }
              }
            } // for 0 <= n < neighs.size()
          } // for 0 <= item < _num_subcells_per_subdim[sdim]
        } // for 0 <= sdim < cell_space_dim
      }

      /// print neighbourhood information
      inline void print(std::ostream& stream)
      {
        // print neighbourhood information into the next line
        stream << std::endl << "    [N:  ";
        for(unsigned char subdim(0) ; subdim < cell_space_dim_ ; ++subdim)
        {
          if(subdim == 0)
          {
            stream << "V( ";
          }
          else if(subdim == 1)
          {
            stream << ", E( ";
          }
          else if(subdim == 2)
          {
            stream << ", F( ";
          }
          else
          {
            stream << ", X( ";
          }
          for(unsigned char item(0) ; item < _num_subcells_per_subdim[subdim] ; ++item)
          {
            if (_neighbours[subdim][item].size() > 0)
            {
              _neighbours[subdim][item][0]->print_index(stream);
              for(unsigned char k(1) ; k < _neighbours[subdim][item].size() ; ++k)
              {
                stream << ", ";
                _neighbours[subdim][item][k]->print_index(stream);
              }
              if(item < _num_subcells_per_subdim[subdim]-1)
              {
                stream << " | ";
              }
              else
              {
                stream << ")";
              }
            }
            else
            {
              if(item < _num_subcells_per_subdim[subdim]-1)
              {
                stream << "- | ";
              }
              else
              {
                stream << "- )";
              }
            }
          }
        }
        stream << "]";
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_DATA_HPP
