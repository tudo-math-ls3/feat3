#pragma once
#ifndef KERNEL_FOUNDATION_NEIGHBOURHOOD_HPP
#define KERNEL_FOUNDATION_NEIGHBOURHOOD_HPP

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

// includes, system

namespace FEAST
{
  /**
   * \brief class providing a neighbourhood data structure for defining connectivity of process patches
   *
   * This data structure represents the part of the global graph this process is associated with. It basically stores
   * the number of neighbours and the ranks of the neighbours. The data structure can be constructed from a global
   * graph with the help of the MPI functions
   *   MPI_Dist_graph_neighbors_count(...)
   * to get the number of neighbours and
   *   MPI_Dist_graph_neighbors(...)
   * to get the ranks of the neighbours. With the help of this data structure the global MPI topology graph can be
   * created via the function MPI_Dist_graph_create(...).
   *
   * \todo This is only a very rough first version, which will be surely adapted to our needs...
   *
   * \author Hilmar Wobker
   * \author Peter Zajac
   */
  class Neighbourhood
  {
  private:
    /// number of neighbours
    int _num_neighbours;

    /**
     * \brief neighbour ranks
     * Dimension: #_num_neighbours
     */
    int* _neighbours;

  public:
    /// CTOR
    Neighbourhood(
      int num_neighbours,
      const int* neighbours) :
      _num_neighbours(num_neighbours),
      _neighbours(nullptr)
    {
      CONTEXT("Neighbourhood::Neighbourhood()");
      ASSERT_(num_neighbours >= 0);
      ASSERT_(neighbours != nullptr);

      _neighbours = new int[num_neighbours];
      for(int i(0); i < num_neighbours; ++i)
      {
        _neighbours[i] = neighbours[i];
      }
    }

    /// DTOR
    virtual ~Neighbourhood()
    {
      if(_neighbours != nullptr)
      {
        delete [] _neighbours;
        _neighbours = nullptr;
      }
    }

    /**
     * \brief Returns the number of neighbours.
     * \returns
     * The number of neighbours.
     */
    inline int num_neighbours() const
    {
      return _num_neighbours;
    }

    /**
     * \brief Returns the neighbour array.
     * \returns
     * An array containing the ranks of the neighbours.
     */
    inline const int* neighbours() const
    {
      return _neighbours;
    }

    /**
     * \brief Dumps the neighbourhood into a stream.
     *
     * This function dumps the neighbourhood into an output stream.
     *
     * \param[in,out] stream
     * The stream where to dump to.
     */
    template<typename Stream_>
    void dump(Stream_& stream) const
    {
      CONTEXT("Neighbourhood::dump()");
      stream << "neighbourhood: ";
      for(int i(0) ; i < _num_neighbours ; ++i)
      {
        stream << _neighbours[i] << " ";
      }
      stream << std::endl;
    }

    /**
     * \brief Dumps the neighbourhood into a string.
     *
     * \returns
     * The neighbourhood dump as a string.
     */
    inline String dump() const
    {
      CONTEXT("Neighbourhood::dump()");
      std::ostringstream oss;
      dump(oss);
      return oss.str();
    }
  }; // class Neighbourhood
} // namespace FEAST

#endif // KERNEL_FOUNDATION_NEIGHBOURHOOD_HPP
