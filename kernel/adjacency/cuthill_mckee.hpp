#pragma once
#ifndef KERNEL_ADJACENCY_CUTHILLMCKEE_HPP
#define KERNEL_ADJACENCY_CUTHILLMCKEE_HPP 1

// includes, FEAST
#include <kernel/adjacency/permutation.hpp>
#include <kernel/adjacency/graph.hpp>

namespace FEAST
{
  namespace Adjacency
  {
    /**
     * \brief Cuthill McKee implementation
     *
     * \todo detailed description
     *
     * \author Constantin Christof
     */
    class CuthillMcKee
    {
    public:

     /**
      * \brief Root type enumeration
      *
      * This enumeration specifies how the root node of the Cuthill-McKee algorithm is chosen.
      */
     enum RootType
      {
        /**
         * \brief Default root
         *
         * In this mode, the first admissible node is chosen.
         */
        root_default = 0,

        /**
         * \brief Minimum-degree root
         *
         * In this mode, a node of minimum degree is chosen as the root node.
         */
        root_minimum_degree = 1,

        /**
         * \brief Maximum-degree root
         *
         * In this mode, a node of maximum degree is chosen as the root node.
         */
        root_maximum_degree = 2
      };

      /**
       * \brief Sort type enumeration
       *
       * This enumeration specifies the order, that is used to arrange the nodes of each level.
       */
      enum SortType
      {
        /**
         * \brief Default sorting type
         *
         * In this mode, the nodes are not sorted at all.
         */
        sort_default = 0,

        /**
         * \brief Ascending order
         *
         * In this mode, the nodes are arranged in an ascending order corresponding to their degrees.
         */
        sort_asc = 1,

        /**
         * \brief Descending order
         *
         * In this mode, the nodes are arranged in a descending order corresponding to their degrees.
         */
        sort_desc = 2
      };

      /**
       * \brief Cuthill-McKee permutation computation functon
       *
       * This function creates a Cuthill-McKee permutation of the given graph.
       *
       * \param[in] graph
       * The graph, the permutation is calculated for.
       *
       * \param[in] reverse
       * This bool determines, if the reverse Cuthill-McKee permutation should be calculated.
       * If \c true, then the reversed permutation is used.
       *
       * \param[in] root_type
       * This parameter determines the way, the root nodes are chosen.
       *
       * \param[in] sort_type
       * This parameter determines, which sorting is used in each level of the Cuthill-McKee algorithm.
       *
       * \returns
       * The Cuthill-McKee permutation.
       */
      static Permutation compute(
        const Graph& graph,
        bool reverse = false,
        CuthillMcKee::RootType root_type = root_default,
        CuthillMcKee::SortType sort_type = sort_default);
    }; // class CuthillMcKee
  } // namespace Adjacency
} // namespace FEAST

#endif // KERNEL_ADJACENCY_CUTHILLMCKEE_HPP
