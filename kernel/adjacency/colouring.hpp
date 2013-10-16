#pragma once
#ifndef KERNEL_ADJACENCY_COLOURING_HPP
#define KERNEL_ADJACENCY_COLOURING_HPP 1

// includes, FEAST
#include <kernel/util/assertion.hpp>

namespace FEAST
{
  namespace Adjacency
  {
    // forward declaration
    class Graph;

    /**
     * \brief Colouring object implementation
     *
     * \todo detailed description
     *
     * \author Constantin Christof
     */
    class Colouring
    {
    protected:
      /// total number of nodes
      Index _num_nodes;

      /// total number of colours used
      Index _num_colours;

      /**
       * \brief colouring array
       *
       * Dimension: #_num_nodes
       * The colouring array is defined as the second row of the following chart:
       *    node number :         0                   1           ...         num_nodes-1
       *    colour number: colour of node 0    colour of node 1    ...  colour of node num_nodes-1
       *
       */
      Index* _colouring;

    public:

      /**
       * \brief Default constructor.
       *
       * This constructor creates a new empty colouring object, but does not allocate any arrays.
       */
      Colouring();

      /**
       * \brief Allocation Constructor.
       *
       * This constructor creates a new colouring object and allocates
       * a colouring array of length num_nodes.
       *
       * \note This constructor does not initialise the allocated array -- it has to be initialised by the user
       * after construction.
       *
       * \param[in] num_nodes
       * The total number of nodes.
       *
       * \param[in] num_colours
       * The number of colours.
       *
       */
      Colouring(
        Index num_nodes,
        Index num_colours);

      /**
       * \brief Array Constructor.
       *
       * This constructor creates a new colouring object from a given array.
       *
       * \param[in] num_nodes
       * The total number of nodes.
       *
       * \param[in] colouring
       * The colouring array
       *
       */
      Colouring(
        Index num_nodes,
        Index* colouring);

      /**
       * \brief Creation out of a given Graph
       *
       * This constructor creates a colouring object out of a graph. The returned colouring
       * satisfies the condition that adjacent nodes do not have the same colour.
       *
       * \param[in] graph
       */
      explicit Colouring(const Graph& graph);

      /**
       * \brief Creation out of a given Graph with a prescribed order
       *
       * This constructor creates a colouring object out of a graph. The returned colouring
       * satisfies the condition that adjacent nodes do not have the same colour. The algorithm
       * proceeds through the nodes as given by the parameter order.
       *
       * \param[in] graph
       *
       * \param[in] order
       * Permutation array that describes how the algorithm is supposed to proceed through the nodes.
       */
      explicit Colouring(const Graph& graph, const Index* order);

      /// virtual destructor
      virtual ~Colouring();

      /**
       * \brief Returns the colouring array.
       * \returns The colouring array.
       */
      Index* get_colouring()
      {
        CONTEXT("Colouring::get_colouring()");
        return _colouring;
      }

      /** \copydoc get_colouring() */
      const Index* get_colouring() const
      {
        CONTEXT("Colouring::get_colouring()");
        return _colouring;
      }

      /**
       * \brief Returns the total number of nodes.
       *
       * \returns The total number of nodes.
       */
      Index get_num_nodes() const
      {
        CONTEXT("Colouring::get_num_nodes()");
        return _num_nodes;
      }

      /**
       * \brief Returns the maximum colour index.
       *
       * \returns The maximum colour index.
       */
      Index get_max_colour() const
      {
        CONTEXT("Colouring::get_max_colours()");
        return _num_colours - 1;
      }
    }; // class Colouring
  } // namespace Adjacency
} // namespace FEAST

#endif // KERNEL_ADJACENCY_COLOURING_HPP
