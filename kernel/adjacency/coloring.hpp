// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ADJACENCY_COLOURING_HPP
#define KERNEL_ADJACENCY_COLOURING_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Adjacency
  {
    // forward declaration
    class Graph;

    /**
     * \brief Coloring object implementation
     *
     * \todo detailed description
     *
     * \author Constantin Christof
     */
    class Coloring
    {
    protected:
      using IndexVector = std::vector < Index>;

      /// total number of colors used
      Index _num_colors;

      /**
       * \brief coloring vector
       *
       * Dimension: #_num_nodes
       * The coloring vector is defined as the second row of the following chart:
       *    node number :         0                   1           ...         num_nodes-1
       *    color number: color of node 0    color of node 1    ...  color of node num_nodes-1
       *
       */
      IndexVector _coloring;

    public:

      /**
       * \brief Default constructor.
       *
       * This constructor creates a new empty coloring object, but does not allocate any vectors.
       */
      Coloring();

      /**
       * \brief Allocation Constructor.
       *
       * This constructor creates a new coloring object and allocates
       * a coloring vector of length num_nodes.
       *
       * \note This constructor does not initialize the allocated vector -- it has to be initialized by the user
       * after construction.
       *
       * \param[in] num_nodes
       * The total number of nodes.
       *
       * \param[in] num_colors
       * The number of colors.
       *
       */
      Coloring(
        Index num_nodes,
        Index num_colors);

      /**
       * \brief Array Constructor.
       *
       * This constructor creates a new coloring object from a given array.
       *
       * \param[in] num_nodes
       * The total number of nodes.
       *
       * \param[in] coloring
       * The coloring array
       *
       */
      Coloring(
        Index num_nodes,
        Index* coloring);

      /**
      * \brief Vector Constructor.
      *
      * This constructor creates a new coloring object from a given vector.
      *
      * \param[in] num_colors
      * The total number of colors.
      *
      * \param[in] coloring
      * The coloring vector
      *
      */
      Coloring(
        Index num_colors,
        const IndexVector& coloring);

      /**
       * \brief Creation out of a given Graph
       *
       * This constructor creates a coloring object out of a graph. The returned coloring
       * satisfies the condition that adjacent nodes do not have the same color.
       *
       * \param[in] graph
       */
      explicit Coloring(const Graph& graph);

      /**
       * \brief Creation out of a given Graph with a prescribed order
       *
       * This constructor creates a coloring object out of a graph. The returned coloring
       * satisfies the condition that adjacent nodes do not have the same color. The algorithm
       * proceeds through the nodes as given by the parameter order.
       *
       * \param[in] graph
       *
       * \param[in] order
       * Permutation array that describes how the algorithm is supposed to proceed through the nodes.
       */
      explicit Coloring(const Graph& graph, const Index* order);

      /// move ctor
      Coloring(Coloring&& other);

      /// move-assign operator
      Coloring& operator=(Coloring&& other);

      /// virtual destructor
      virtual ~Coloring();

      /**
       * \brief Clones this coloring.
       *
       * \returns A deep-copy of this coloring.
       */
      Coloring clone() const
      {
        if (!_coloring.empty())
          //return Coloring(Index(_coloring.size()), _coloring.data());
          return Coloring(_num_colors, _coloring);
        else
          return Coloring();
      }

      /**
       * \brief Creates a color partition graph.
       *
       * This function creates a graph out of this coloring object.
       * The graph's image is the set of nodes, the domain is the set of colors.
       */
      Graph create_partition_graph() const;

      /**
       * \brief Returns the coloring array.
       * \returns The coloring array.
       */
      Index* get_coloring()
      {
        return _coloring.data();
      }

      /** \copydoc get_coloring() */
      const Index* get_coloring() const
      {
        return _coloring.data();
      }

      /**
       * \brief Returns the total number of nodes.
       *
       * \returns The total number of nodes.
       */
      Index get_num_nodes() const
      {
        return Index(_coloring.size());
      }

      /**
       * \brief Returns the maximum color index.
       *
       * \returns The maximum color index.
       */
      Index get_max_color() const
      {
        return _num_colors - 1;
      }
    }; // class Coloring
  } // namespace Adjacency
} // namespace FEAT

#endif // KERNEL_ADJACENCY_COLOURING_HPP
