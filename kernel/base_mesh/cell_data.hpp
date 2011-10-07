/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR below.
 * Generally, Peter wanted to take a deeper look at the base mesh implementation.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_BASE_MESH_CELL_DATA_HPP
#define KERNEL_BASE_MESH_CELL_DATA_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/base_mesh/item.hpp>

// includes, system
#include <iostream> // for std::ostream
#include <vector>   // for std::vector

namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
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



    /**
    * \brief general cell information which is needed by Cell and CellData class
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
// COMMENT_HILMAR: Das erscheint irgendwie maechtig umstaendlich, ist aber leider aufgrund der Besonderheit der
//   "dazwischen geschobenen" Klasse CellData<...> (die fuer cell_dim_ != space_dim_ leer ist) noetig.

// COMMENT_HILMAR: Im Moment gibt es hier nur das array _num_subitems_per_subdim[], eventuell kommt aber noch mehr
//   dazu. Sollte sich herausstellen, dass doch nix mehr dazu kommt, dann kann man diese Klasse hier auch wieder
//   abschaffen und dieses einzige array per Code-Verdopplung in die beiden Varianten der Klasse CellData verschieben.
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class CellInfo
      : public Item
    {

    private:

      /**
      * \brief array of number of subitems per subdimension
      *
      * Definitions of two terms we use:
      * 1) "subdimension" = space dimension minus co-dimension where
      *      in 3D: codim 1 -> faces, codim 2 -> edges, codim 3 -> vertices
      *      in 2D: codim 1 -> edges, codim 2 -> vertices
      *      in 1D: codim 1 -> vertices
      *    Co-dimension 0 is not needed here. Hence, the number of subdimensions is equal to the space dimension.
      * 2) "subitems" = items of lower dimension in a cell, distinguished by subdimension. Examples:
      *      in 3D: a hexa has 6 subitems of subdimension 2 (=faces), 12 subitems of subdimension 1 (=edges)
      *             and 8 subitems of subdimension 0 (=vertices)
      *      in 2D: a tri has 3 subitems of subdimension 1 (=edges) and 3 subitems of subdimension 0 (=vertices)
      *      in 1D: an edge has 2 subitems of subdimension 0 (=vertices)
      */
      unsigned char _num_subitems_per_subdim[cell_dim_];
// COMMENT_HILMAR: Eigentlich ganz schoener overkill, das fuer jede Zelle extra abzuspeichern. Besser waere evtl.,
//   die Informationen irgendwo statisch pro Zelltyp 1x abzuspeichern und dann mit einer Art Zelltyp-ID abzufragen.


    protected:

      /**
      * \brief setter for the number of subitems per subdimension
      *
      * \param[in] array_size
      * size of the array num_subitems_per_subdim
      *
      * \param[in] num_subitems_per_subdim
      * number of subitems per subdimension
      */
      void _set_num_subitems_per_subdim(
        unsigned char array_size,
        unsigned char num_subitems_per_subdim[])
      {
        CONTEXT("BaseMesh::CellInfo::_set_num_subitems_per_subdim()");
        ASSERT(array_size == cell_dim_, "Array size " + stringify(array_size) + " must be equal to cell dimension"
               + stringify(cell_dim_) + ".");
        for(int i(0) ; i < array_size ; ++i)
        {
          _num_subitems_per_subdim[i] = num_subitems_per_subdim[i];
        }
      }


    public:

      /**
      * \brief getter for the number of subitems of the given subdimension
      *
      * \param[in] subdim
      * desired subdimension
      *
      * \return number of subitems of the given subdimension
      */
      inline unsigned char num_subitems_per_subdim(subdim subdim) const
      {
        CONTEXT("BaseMesh::CellInfo::num_subitems_per_subdim()");
        ASSERT(subdim < cell_dim_, "Subdimension " + stringify(subdim) + " must be smaller than cell dimension"
               + stringify(cell_dim_) + ".");
        return _num_subitems_per_subdim[subdim];
      }


    };


    /// forward declaration of class Cell
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Cell;


    /**
    * \brief emtpy cell data class definition for cells with dimension smaller than space dimension
    *
    * Avoid instatiation of cell-specific data in cells with dimension smaller than space dimension.
    * The class is only implemented for cell_dim_ = space_dim_ (see below).
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
    class CellData
      : public CellInfo<cell_dim_, space_dim_, world_dim_>
    {

    private:


    public:

      /// dummy function called by cells with dimension smaller than space dimension
      inline void _init_neighbours()
      {
        // do nothing here
      }

      /// dummy function called by cells with dimension smaller than space dimension
      inline void add_neighbour(
        subdim,
        unsigned char,
        Cell<cell_dim_, space_dim_, world_dim_>*)
      {
        // do nothing here
      }

      /// dummy function called by cells with dimension smaller than space dimension
      inline void print(std::ostream&) const
      {
        // do nothing here
      }

    };


    /**
    * \brief class storing cell-specific data like neighbourhood information
    *
    * This data is only stored for cells of full dimension, i.e., for those having space dimension.
    *
    * \tparam cell_space_dim_
    * cell and space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in
    * a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_space_dim_,
      unsigned char world_dim_>
    class CellData<cell_space_dim_, cell_space_dim_, world_dim_>
      : public CellInfo<cell_space_dim_, cell_space_dim_, world_dim_>
    {

    private:

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
      void _init_neighbours()
      {
        CONTEXT("BaseMesh::CellData::_init_neighbours()");
        for(int sdim(0) ; sdim < cell_space_dim_ ; ++sdim)
        {
          _neighbours[sdim]
            = new std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>
                    [this->num_subitems_per_subdim((subdim)sdim)];
        }
      }


    public:

      /// DTOR (automatically virtual since DTOR of base class Item is virtual)
      ~CellData()
      {
        CONTEXT("BaseMesh::CellData::~CellData()");
        for(int sdim(0) ; sdim < cell_space_dim_ ; ++sdim)
        {
          for(unsigned int item(0) ; item < this->num_subitems_per_subdim((subdim)sdim) ; ++item)
          {
            _neighbours[sdim][item].clear();
          }
          delete [] _neighbours[sdim];
        }
      }


      /// returns number of neighbours for given subdimension and given item
      inline unsigned int num_neighbours_item(
        subdim subdim,
        unsigned char item) const
      {
        CONTEXT("BaseMesh::CellData::num_neighbours_item()");
        ASSERT(subdim < cell_space_dim_, "Subdimension " + stringify(subdim)
               + " must not exceed cell-space dimension" + stringify(cell_space_dim_) + ".");
        ASSERT(item < this->num_subitems_per_subdim(subdim), "Item index " + stringify(item)
               + " must not exceed number of subitems " + stringify(this->num_subitems_per_subdim(subdim))
               + " of subdimension " + stringify(subdim) + ".");
        return _neighbours[subdim][item].size();
      }


      /// returns number of neighbours summed over all items of given subdimension
      inline unsigned int num_neighbours_subdim(subdim subdim) const
      {
        CONTEXT("BaseMesh::CellData::num_neighbours_subdim()");
        ASSERT(subdim < cell_space_dim_, "Subdimension " + stringify(subdim)
               + " must be smaller than cell-space dimension" + stringify(cell_space_dim_) + ".");
        unsigned int num_neighbours(0);
        for(unsigned int item(0) ; item < this->num_subitems_per_subdim(subdim) ; ++item)
        {
          num_neighbours += (unsigned int)(_neighbours[subdim][item].size());
        }
        return num_neighbours;
      }


      /// returns one specific neighbour for given subdim, item and index
      inline Cell<cell_space_dim_, cell_space_dim_, world_dim_>* neighbour(
        subdim subdim,
        unsigned char item,
        unsigned char index) const
      {
        CONTEXT("BaseMesh::CellData::neighbour()");
        ASSERT(subdim < cell_space_dim_, "Subdimension " + stringify(subdim)
               + " must not exceed cell-space dimension" + stringify(cell_space_dim_) + ".");
        ASSERT(item < this->num_subitems_per_subdim(subdim), "Item index " + stringify(item)
               + " must not exceed number of subitems " + stringify(this->num_subitems_per_subdim(subdim))
               + " of subdimension " + stringify(subdim) + ".");
        ASSERT(index < _neighbours[subdim][item].size(), "Index " + stringify(index)
               + " must not exceed neighbour vector size " + stringify(_neighbours[subdim][item].size())
               + " for item " + stringify(item) + " of subdimension " + stringify(subdim) + ".");
        return _neighbours[subdim][item][index];
      }


      /// returns vector of neighbours for given subdimension and given item
      inline std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>& neighbours_item(
        subdim subdim,
        unsigned char item) const
      {
        CONTEXT("BaseMesh::CellData::neighbours_item()");
        ASSERT(subdim < cell_space_dim_, "Subdimension " + stringify(subdim)
               + " must not exceed cell-space dimension" + stringify(cell_space_dim_) + ".");
        ASSERT(item < this->num_subitems_per_subdim(subdim), "Item index " + stringify(item)
               + " must not exceed number of subitems " + stringify(this->num_subitems_per_subdim(subdim))
               + " of subdimension " + stringify(subdim) + ".");
        return _neighbours[subdim][item];
      }


      /// returns array of vectors of neighbours for given subdimension
      inline std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>* neighbours_subdim(subdim subdim) const
      {
        CONTEXT("BaseMesh::CellData::neighbours_subdim()");
        ASSERT(subdim < cell_space_dim_, "Subdimension " + stringify(subdim)
               + " must not exceed cell-space dimension" + stringify(cell_space_dim_) + ".");
        return _neighbours[subdim];
      }


      /// add neighbour to the vector of neighbours of given subdim and item
      inline void add_neighbour(
        subdim subdim,
        unsigned char item,
        Cell<cell_space_dim_, cell_space_dim_, world_dim_>* neighbour)
      {
        CONTEXT("BaseMesh::CellData::add_neighbour()");
        ASSERT(subdim < cell_space_dim_, "Subdimension " + stringify(subdim)
               + " must not exceed cell-space dimension" + stringify(cell_space_dim_) + ".");
        ASSERT(item < this->num_subitems_per_subdim(subdim), "Item index " + stringify(item)
               + " must not exceed number of subitems " + stringify(this->num_subitems_per_subdim(subdim))
               + " of subdimension " + stringify(subdim) + ".");
        _neighbours[subdim][item].push_back(neighbour);
      }


      /// print neighbourhood information
      inline void print(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::CellData::print()");
        // print neighbourhood information into the next line
        stream << std::endl << "    [N:  ";
        for(unsigned char sdim(0) ; sdim < cell_space_dim_ ; ++sdim)
        {
          if(sdim == 0)
          {
            stream << "V( ";
          }
          else if(sdim == 1)
          {
            stream << ", E( ";
          }
          else if(sdim == 2)
          {
            stream << ", F( ";
          }
          else
          {
            stream << ", X( ";
          }
          for(unsigned char item(0) ; item < this->num_subitems_per_subdim((subdim)sdim) ; ++item)
          {
            if (_neighbours[sdim][item].size() > 0)
            {
              _neighbours[sdim][item][0]->print_index(stream);
              for(unsigned char k(1) ; k < _neighbours[sdim][item].size() ; ++k)
              {
                stream << ", ";
                _neighbours[sdim][item][k]->print_index(stream);
              }
              if(item < this->num_subitems_per_subdim((subdim)sdim)-1)
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
              if(item < this->num_subitems_per_subdim((subdim)sdim)-1)
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
      } // print()
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // KERNEL_BASE_MESH_CELL_DATA_HPP
