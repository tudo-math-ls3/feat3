// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Adjacency
  {
    Permutation CuthillMcKee::compute(
      const Graph& graph,
      bool reverse,
      CuthillMcKee::RootType root_type,
      CuthillMcKee::SortType sort_type)
    {
      // get the number of nodes
      Index num_nodes = graph.get_num_nodes_domain();

      // get the graph's data
      const Index* domain_ptr = graph.get_domain_ptr();
      const Index* image_idx = graph.get_image_idx();

      // a vector which keeps track of which nodes have already been processed
      std::vector<bool> node_mask(num_nodes, false);

      // a vector storing the degree of each node
      std::vector<Index> node_degree(num_nodes);
      for(Index j(0); j < num_nodes; ++j)
      {
        node_degree[j] = graph.degree(j);
      }

      // auxiliary variables
      Index lvl1 = 0;
      Index lvl2 = 0;
      Index lvl3 = 0;

      // create permutation
      Permutation perm(graph.get_num_nodes_domain());
      Index* permutation = perm.get_perm_pos();

      while(lvl2 < num_nodes)
      {
        // initialize root
        Index root = num_nodes + 1;

        // get the root node
        switch(root_type)
        {

        // default: choose the first node possible
        case root_default:
          for(Index j(0); j < num_nodes; ++j)
          {
            // check if the node is unmarked
            if(!node_mask[j])
            {
              root = j;
              break;
            }
          }
          break;

        // choose the node of minimum degree
        case root_minimum_degree:
          {
            Index min = num_nodes + 1;
            for(Index j(0); j < num_nodes; ++j)
            {
              if((node_degree[j] < min) && !node_mask[j])
              {
                root = j;
                min = node_degree[j];
              }
            }
          }
          break;

        // choose the node of maximum degree
        case root_maximum_degree:
          {
            Index max = 0;
            for(Index j(0); j < num_nodes; ++j)
            {
              if((node_degree[j] > max) && !node_mask[j])
              {
                root = j;
                max = node_degree[j];
              }
            }
            break;
          }

        // else: error
        default:
          XABORTM("Invalid root_type parameter!");
        }

        // if something very odd has happened
        XASSERTM(root < num_nodes, "No root node found!");

        // initialize the first level of the root node
        ++lvl1;
        lvl2 = lvl1;
        lvl3 = lvl1;

        // add the root into the permutation array
        Index lvl_root = lvl1;
        permutation[lvl1 - 1] = root;
        node_mask[root] = 1;

        // loop through the adjancy levels of the root
        while(lvl2 < num_nodes)
        {

          // loop through all nodes in the current level
          for(Index i(lvl1 - 1); i < lvl2 ; ++i)
          {
            // get the node's index
            Index n = permutation[i];

            // Go through all nodes which are adjacent to node n
            for(Index j(domain_ptr[n]); j < domain_ptr[n+1] ; ++j)
            {
              // Get the index of the adjacent node
              Index k = image_idx[j];

              // has this node been processed?
              if(!node_mask[k])
              {
                ++lvl3;
                permutation[lvl3 - 1] = k;
                node_mask[k] = true;
              }
            } //j loop
          } // i loop

          if(lvl3 <= lvl2)
          {
            break;
          }

          // sorting
          switch(sort_type)
          {
          // default: do nothing
          case sort_default:
            break;

          // descending order
          case sort_desc:

            // if there is nothing to do
            if(lvl2 + 1 >= lvl3)
            {
              break;
            }

            // sorting algorithm: linear insertion
            for(Index i(lvl2); i < lvl3; ++i)
            {
              Index x = node_degree[permutation[i]];
              Index y = permutation[i];
              Index k(i);
              for(; (k > lvl2) && (x > node_degree[permutation[k-1]]); --k)
              {
                permutation[k] = permutation[k-1];
              }
              permutation[k] = y;
            }
            break;

          // ascending order
          case sort_asc:

            // if there is nothing to do
            if(lvl2 + 1 >= lvl3)
            {
              break;
            }

            // sorting algorithm: linear insertion
            for(Index i(lvl2); i < lvl3; ++i)
            {
              Index x = node_degree[permutation[i]];
              Index y = permutation[i];
              Index k(i);
              for(; (k > lvl2) && (x < node_degree[permutation[k-1]]); --k)
              {
                permutation[k] = permutation[k-1];
              }
              permutation[k] = y;
            }
            break;

          // else: error
          default:
            XABORTM("Invalid run_type parameter!");
          }

          // get to the next adjacency level
          lvl1 = lvl2+1;
          lvl2 = lvl3;

        } //while(lvl2 < num_nodes)

        // reverse?
        if(reverse)
        {
          for(Index ind(lvl_root - 1), jnd(lvl2 - 1); ind < jnd; ++ind, --jnd)
          {
            Index n = permutation[ind];
            permutation[ind] = permutation[jnd];
            permutation[jnd] = n;
          }
        }

      } //while(lvl2 < num_nodes) (the separability loop)

      // compute swap array
      perm.calc_swap_from_perm();

      // return permutation
      return perm;
    }
  } // namespace Adjacency
} // namespace FEAT
