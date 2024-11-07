// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/adjacency/permutation.hpp>
#include <kernel/adjacency/graph.hpp>

// includes, system
#include <vector>

namespace FEAT
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
     enum class RootType
      {
        /**
         * \brief Default root
         *
         * In this mode, the first admissible node is chosen.
         */
        standard = 0,

        /**
         * \brief Minimum-degree root
         *
         * In this mode, a node of minimum degree is chosen as the root node.
         */
        minimum_degree = 1,

        /**
         * \brief Maximum-degree root
         *
         * In this mode, a node of maximum degree is chosen as the root node.
         */
        maximum_degree = 2
      };

      /**
       * \brief Sort type enumeration
       *
       * This enumeration specifies the order, that is used to arrange the nodes of each level.
       */
      enum class SortType
      {
        /**
         * \brief Default sorting type
         *
         * In this mode, the nodes are not sorted at all.
         */
        standard = 0,

        /**
         * \brief Ascending order
         *
         * In this mode, the nodes are arranged in an ascending order corresponding to their degrees.
         */
        asc = 1,

        /**
         * \brief Descending order
         *
         * In this mode, the nodes are arranged in a descending order corresponding to their degrees.
         */
        desc = 2
      };

      /**
       * \brief Cuthill-McKee permutation computation function
       *
       * This function creates a Cuthill-McKee permutation of the given graph.
       *
       * \param[in] graph
       * The \transient graph that the permutation is calculated for.
       *
       * \param[in] reverse
       * This bool determines, if the reverse Cuthill-McKee permutation should be calculated.
       * If \c true, then the reversed permutation is used.
       *
       * \param[in] RootType::type
       * This parameter determines the way, the root nodes are chosen.
       *
       * \param[in] SortType::type
       * This parameter determines, which sorting is used in each level of the Cuthill-McKee algorithm.
       *
       * \returns
       * The Cuthill-McKee permutation.
       */
      static Permutation compute(
        std::vector<Index>& layers,
        const Graph& graph,
        bool reverse = false,
        CuthillMcKee::RootType r_type = RootType::standard,
        CuthillMcKee::SortType t_type = SortType::standard);

      /**
       * \brief Cuthill-McKee permutation computation function
       *
       * This function creates a Cuthill-McKee permutation of the given graph.
       *
       * \param[in] graph
       * The \transient graph that the permutation is calculated for.
       *
       * \param[in] reverse
       * This bool determines, if the reverse Cuthill-McKee permutation should be calculated.
       * If \c true, then the reversed permutation is used.
       *
       * \param[in] RootType::type
       * This parameter determines the way, the root nodes are chosen.
       *
       * \param[in] SortType::type
       * This parameter determines, which sorting is used in each level of the Cuthill-McKee algorithm.
       *
       * \returns
       * The Cuthill-McKee permutation.
       */
      static Permutation compute(
        const Graph& graph,
        bool reverse = false,
        CuthillMcKee::RootType R_type = RootType::standard,
        CuthillMcKee::SortType s_type = SortType::standard);
    }; // class CuthillMcKee
  } // namespace Adjacency
} // namespace FEAT
