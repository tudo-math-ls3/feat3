// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/adjacency/adjactor.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/util/dist.hpp>

#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>

namespace FEAT::Global
{
  /// Distributed flood fill
  class FloodFill
  {
  public:
    /**
     * \brief Distributed flood fill
     *
     * \param[in] first_node_of_rank Vector that contains the first node of each rank
     * \param[in] graph Subgraph of current rank
     * \param[in] comm Communicator that graph is distributed over
     *
     * \returns A Colors x Nodes graph which assigns a color to each local node of the input graph
     *
     * \note This method assumes that the input graph on rank `i` contains all nodes from `first_node_on_rank[i]` to
     * `first_node_on_rank[i + 1] - 1`.
     */
    static Adjacency::Graph
    run(const std::vector<Index>& first_node_of_rank, const Adjacency::Graph& graph, Dist::Comm& comm)
    {
      // This distributed flood fill works essentially like a normal flood fill with two exceptions:
      // - When starting a new "flood" the seed node is selected by the smallest rank with uncolored nodes left
      // - When encountering a node of another rank, that node is added to that ranks stack
      // For each seed the flood fill is iterated until all ranks stacks are empty.
      // The flood fill is done if there are no more uncolored nodes on any rank.

      // NOTE(mmuegge): There is probably an algorithm that scales better,
      // which performs a normal flood fill on each rank at the same time
      // and the reconciles the different colors used afterwards.

      // The starting rank for each flood is selected by a bidding proccess,
      // where each ranks either submits its own rank or max_bid.
      // The rank matching the smallest bid starts the next flood.
      static constexpr int max_bid = std::numeric_limits<int>::max();

      // Sentinel value for uncolored nodes.
      static constexpr Index uncolored(~Index(0));

      const int rank = comm.rank();
      const Index num_nodes = graph.get_num_nodes_domain();

      // Stack of nodes to be processed on this rank
      std::vector<Index> stack;

      // Colors for this ranks nodes
      std::vector<Index> colors(num_nodes, uncolored);

      // Color of current flood.
      Index current_color = 0;

      while(true)
      {
        // Find rank that starts next flood fill
        // If no such rank exists, break
        auto first_uncolored = std::find_if(colors.begin(), colors.end(), [](Index c) { return c == uncolored; });
        bool all_colored = first_uncolored == colors.end();

        // Submit either our rank or max_bid as bid.
        // Smallest non-max_bid bid starts the next flood.
        const int bid = all_colored ? max_bid : rank;
        int winner(0);
        comm.allreduce(&bid, &winner, 1, Dist::op_min);

        if(winner == max_bid)
        {
          // All ranks submitted max_bid.
          // There is thus no rank remaining that could start a flood fill.
          // We are done.
          break;
        }

        if(winner == rank)
        {
          // We won the bid
          // Initialize our stack with first uncolored node
          const auto idx = (Index)std::distance(colors.begin(), first_uncolored);
          stack.push_back(idx);
          colors[idx] = current_color;
        }

        while(true)
        {
          // Buffers for communicating foreign nodes
          std::vector<int> ranks;
          std::vector<Index> nodes;

          // Process local stack
          while(!stack.empty())
          {
            const Index current_node = stack.back();
            stack.pop_back();

            for(const Index neighbor : Adjacency::ImageRange(graph, current_node))
            {
              const int node_rank = get_rank_of_node(first_node_of_rank, neighbor);

              // We assume that the graph gives us "global" indices
              // We have to convert these to rank-local indices first before we
              // can do anything useful with them
              const Index local_neighbor = neighbor - first_node_of_rank[std::size_t(node_rank)];

              if(node_rank == rank)
              {
                // We own this node
                // If it has not already been visited add it to the stack
                if(colors[local_neighbor] == uncolored)
                {
                  stack.push_back(local_neighbor);
                  colors[local_neighbor] = current_color;
                }
              }
              else
              {
                // Some other rank owns this node.
                // Communicate that connection to the other rank later
                ranks.push_back(node_rank);
                nodes.push_back(local_neighbor);
              }
            }
          }

          // Send and receive from other ranks
          std::vector<Index> visited_by_others = exchange_nodes(comm, ranks, nodes);

          // Reset buffers for next iteration
          ranks.clear();
          nodes.clear();

          // Add newly visited nodes to our stack
          for(Index node : visited_by_others)
          {
            if(colors[node] == uncolored)
            {
              // Node has not been visited on this rank
              stack.push_back(node);
              colors[node] = current_color;
            }
          }

          // Determine if all ranks are done
          // If so, break

          const int done = stack.empty() ? 1 : 0;
          int all_done(0);
          comm.allreduce(&done, &all_done, 1, Dist::op_min);

          if(all_done == 1)
          {
            break;
          }
        }

        // Go to next color
        current_color++;
      }

      // Everything is colored.
      // Construct output color -> nodes graph

      // current_color is the last color used when all ranks determined they were done
      const Index num_colors = current_color;

      Adjacency::Graph result(num_colors, num_nodes, num_nodes);

      // Count occurrences of each color to determine non-zeros for each color
      std::vector<Index> color_counts(num_colors, 0);
      for(Index color : colors)
      {
        color_counts[color]++;
      }

      // Write domain ptr using color_counts
      Index* dom_ptr = result.get_domain_ptr();
      dom_ptr[0] = 0;
      for(Index i(1); i < num_colors; i++)
      {
        dom_ptr[i] = color_counts[i - 1];
      }

      // Fill in nodes
      Index* idx_ptr = result.get_image_idx();
      for(Index color(0); color < num_colors; color++)
      {
        for(Index node(0); node < num_nodes; node++)
        {
          if(colors[node] == color)
          {
            *idx_ptr = node;
            idx_ptr++;
          }
        }
      }

      return result;
    }

  private:
    /// alltoallv for exchanging nodes that have been visited in the last flood fill iteration
    static std::vector<Index>
    exchange_nodes(Dist::Comm& comm, const std::vector<int>& ranks, const std::vector<Index>& nodes)
    {
      // Build sendcounts and receivecounts
      std::vector<int> sendcounts(std::size_t(comm.size()), 0);
      std::vector<int> receivecounts(std::size_t(comm.size()), 0);

      for(int rank : ranks)
      {
        sendcounts[std::size_t(rank)]++;
      }

      comm.alltoall(sendcounts.data(), 1, receivecounts.data(), 1);

      // Build displacements
      std::vector<int> senddisplacements;
      std::vector<int> receivedisplacements;
      std::exclusive_scan(sendcounts.begin(), sendcounts.end(), std::back_inserter(senddisplacements), 0);
      std::exclusive_scan(receivecounts.begin(), receivecounts.end(), std::back_inserter(receivedisplacements), 0);

      // Build sendbuffer
      const std::size_t sendbuffer_size(ranks.size());
      std::vector<Index> sendbuffer(sendbuffer_size, 0);

      std::vector<int> write_offsets = senddisplacements;
      for(Index i(0); i < ranks.size(); i++)
      {
        const auto rank = std::size_t(ranks[i]);
        sendbuffer[std::size_t(write_offsets[rank])] = nodes[i];
        write_offsets[rank]++;
      }

      // Allocate receivebuffer
      const std::size_t receivebuffer_size = std::size_t(receivedisplacements.back() + receivecounts.back());
      std::vector<Index> receivebuffer(receivebuffer_size, 0);

      // Exchange nodes between all ranks
      comm.alltoallv(
        sendbuffer.data(),
        sendcounts.data(),
        senddisplacements.data(),
        receivebuffer.data(),
        receivecounts.data(),
        receivedisplacements.data());

      return receivebuffer;
    }

    static int get_rank_of_node(const std::vector<Index>& first_of_rank, Index node)
    {
      int rank = (int)first_of_rank.size() - 1;
      for(auto it = first_of_rank.rbegin(); it != first_of_rank.rend(); ++it, --rank)
      {
        if(*it <= node)
        {
          return rank;
        }
      }
      XABORTM("FloodFill: Could not find rank for node");
    }
  };

} // namespace FEAT::Global
