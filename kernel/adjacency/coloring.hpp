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
      /// total number of nodes
      Index _num_nodes;

      /// total number of colors used
      Index _num_colors;

      /**
       * \brief coloring array
       *
       * Dimension: #_num_nodes
       * The coloring array is defined as the second row of the following chart:
       *    node number :         0                   1           ...         num_nodes-1
       *    color number: color of node 0    color of node 1    ...  color of node num_nodes-1
       *
       */
      Index* _coloring;

    public:

      /**
       * \brief Default constructor.
       *
       * This constructor creates a new empty coloring object, but does not allocate any arrays.
       */
      Coloring();

      /**
       * \brief Allocation Constructor.
       *
       * This constructor creates a new coloring object and allocates
       * a coloring array of length num_nodes.
       *
       * \note This constructor does not initialize the allocated array -- it has to be initialized by the user
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
        if(_coloring != nullptr)
          return Coloring(_num_nodes, _coloring);
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
        return _coloring;
      }

      /** \copydoc get_coloring() */
      const Index* get_coloring() const
      {
        return _coloring;
      }

      /**
       * \brief Returns the total number of nodes.
       *
       * \returns The total number of nodes.
       */
      Index get_num_nodes() const
      {
        return _num_nodes;
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
