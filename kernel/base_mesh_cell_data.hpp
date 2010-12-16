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
      /// dummy function called by constructors of cells with dimension smaller than space dimension
      inline void _init_neighbours(unsigned char array_size, unsigned char num_subcells_per_subdim[])
      {
        // do nothing here
      }
      /// dummy print function
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
      *             and 8 subcells of subdimension 1 (=vertices)
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
      * _neighbours[0][1]: vector of vertex neighbours at the vertex 1
      * _neighbours[1][3]: vector of edge neighbours at the edge 3
      * _neighbours[2][4]: vector of face neighbours at the face 4
      */
      std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>* _neighbours[cell_space_dim_];

    protected:
      void _init_neighbours(unsigned char array_size, unsigned char num_subcells_per_subdim[])
      {
        assert(array_size == cell_space_dim_);
        for(int i(0) ; i < cell_space_dim_ ; ++i)
        {
          _num_subcells_per_subdim[i] = num_subcells_per_subdim[i];
          _neighbours[i] = new std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>[_num_subcells_per_subdim[i]];
          // TODO: Sind die einzelen std::vector damit schon ordentlich initialisiert? Muesste eigentlich...
        }
      }

    public:

//      CellData(unsigned char num_subcells_per_subdim[])
//      {
//        for(int i(0) ; i < cell_space_dim_ ; ++i)
//        {
//          _neighbours[i] = new std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>[_num_subcells_per_subdim[i]];
//        }
//      }

      ~CellData()
      {
        // TODO: deallocate _neighbours !!!
      }


      inline Cell<cell_space_dim_, cell_space_dim_, world_dim_>* neighbour(
        unsigned char subdim,
        unsigned char item,
        unsigned char index) const
      {
        assert(subdim < world_dim_);
        assert(item < _num_subcells_per_subdim[subdim]);
        assert(index < _neighbours[subdim][item].size());
        return _neighbours[subdim][item][index];
      }


      inline void add_neighbour(
        unsigned char subdim,
        unsigned char item,
        Cell<cell_space_dim_, cell_space_dim_, world_dim_>* neighbour)
      {
        _neighbours[subdim][item].push_back(neighbour);
      }


      inline void print(std::ostream& stream)
      {
        stream << "[ N:  ";
        for(unsigned char subdim(0) ; subdim < cell_space_dim_ ; ++subdim)
        {
          if(subdim == 0)
          {
            stream << "V( ";
          }
          else if(subdim == 1)
          {
            stream << "E( ";
          }
          else if(subdim == 2)
          {
            stream << "F( ";
          }
          else
          {
            stream << "X( ";
          }
          for(unsigned char item(0) ; item < _num_subcells_per_subdim[subdim] ; ++item)
          {
            if (_neighbours[subdim][item].size() > 0)
            {
              stream << _neighbours[subdim][item][0]->index();
              for(unsigned char k(1) ; k < _neighbours[subdim][item].size() ; ++k)
              {
                stream << ", " << _neighbours[subdim][item][k]->index();
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
          stream << " ";
        }
        stream << "]";
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_DATA_HPP
