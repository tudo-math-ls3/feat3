// includes, FEAST
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <kernel/adjacency/graph.hpp>

namespace FEAST
{
  namespace Adjacency
  {
    CuthillMcKee::CuthillMcKee(
      const Graph& graph,
      bool reverse,
      CuthillMcKee::RootType root_type,
      CuthillMcKee::SortType sort_type)
        :
      _perm(graph.get_num_nodes_domain())
    {
      CONTEXT("CuthillMcKee::CuthillMcKee() [graph]");

      // get the number of nodes
      Index num_nodes = graph.get_num_nodes_domain();

      // get the graph's data
      const Index* domain_ptr = graph.get_domain_ptr();
      const Index* image_idx = graph.get_image_idx();

      // auxiliary array storing the node degrees
      int* node_degree = new int[num_nodes];
      for(Index j(0); j < num_nodes; ++j)
      {
        node_degree[j] = graph.degree(j);
      }

      // auxiliary variables
      Index lvl1 = 0;
      Index lvl2 = 0;
      Index lvl3, lvl_root, n, k, ind, jnd;
      int root, min, max;

      // permutation array
      Index* permutation = _perm.get_perm_pos();

      while(lvl2 < num_nodes)
      {
        // initialise root
        root = -1;

        // get the root node
        switch(root_type)
        {

        // default: choose the first node possible
        case root_default:
          for(Index j(0); j < num_nodes; ++j)
          {
            if(node_degree[j] > 0)
            {
              root = j;
              break;
            }
          }
          break;

        // choose the node of minimum degree
        case root_minimum_degree:
          min = num_nodes + 1;
          for(Index j(0); j < num_nodes; ++j)
          {
            if(node_degree[j] < min && node_degree[j] > 0)
            {
              root = j;
              min = node_degree[j];
            }
          }
          break;

        // choose the node of maximum degree
        case root_maximum_degree:
          max = 0;
          for(Index j(0); j < num_nodes; ++j)
          {
            if(node_degree[j] > max)
            {
              root = j;
              max = node_degree[j];
            }
          }
          break;

        // else: error
        default:
          throw InternalError("Invalid root_type parameter!");
        }

        // if something very odd has happend
        if(root < 0)
        {
          throw InternalError("No root found!");
        }

        // initialise the first level of the root node
        ++lvl1;
        lvl2 = lvl1;
        lvl3 = lvl1;

        // add the root into the permutation array
        lvl_root = lvl1;
        permutation[lvl1 - 1] = root;
        node_degree[root] = -node_degree[root];

        // loop through the adjancy levels of the root
        while(lvl2 < num_nodes)
        {

          // loop through all nodes in the current level
          for(Index i(lvl1 - 1); i < lvl2 ; ++i)
          {
            // get the node's index
            n = permutation[i];

            // Go through all nodes which are adjacent to node n
            for(Index j(domain_ptr[n]); j < domain_ptr[n+1] ; ++j)
              {
                // Get the index of the adjacent node
                k = image_idx[j];

                // has this node been processed?
                if(node_degree[k] > 0)
                {
                  ++lvl3;
                  permutation[lvl3 - 1] = k;
                  node_degree[k] = -node_degree[k];
                } //if

              }//j loop
          }// i loop

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
              int x = node_degree[permutation[i]];
              Index y = permutation[i];
              Index k(i);
              for(; (k > lvl2) && (x > node_degree[permutation[k-1]]); --k)
              {
                permutation[k] = permutation[k-1];
              }
              permutation[k] = y;
            }
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
              int x = node_degree[permutation[i]];
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
            throw InternalError("Invalid run_type parameter!");
          }

          // get to the next adjacency level
          lvl1 = lvl2+1;
          lvl2 = lvl3;

        } //while(lvl2 < num_nodes)

        // reverse?
        if(reverse)
        {
          ind = lvl_root - 1;
          jnd = lvl2 - 1;
          while(ind < jnd)
          {
            n = permutation[ind];
            permutation[ind] = permutation[jnd];
            permutation[jnd] = n;
            ++ind;
            --jnd;
          }
        }

      } //while(lvl2 < num_nodes) (the separability loop)

      // create Cuthill-McKee permutation
      _perm.calc_swap_from_perm();

      // deleting the auxiliary arrays
      delete [] node_degree;
    }
  } // namespace Adjacency
} // namespace FEAST
