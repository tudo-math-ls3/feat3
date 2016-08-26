#pragma once
#ifndef KERNEL_ADJACENCY_UNI_EDGE_SET_HPP
#define KERNEL_ADJACENCY_UNI_EDGE_SET_HPP 1

#include <kernel/adjacency/graph.hpp>

#include <vector>

namespace FEAT
{
  namespace Adjacency
  {
    class UniEdgeSet
    {
    protected:
      Graph _graph;

    public:
      UniEdgeSet()
      {
      }

      explicit UniEdgeSet(const Graph& graph)
      {
        XASSERT(graph.get_num_nodes_domain() == graph.get_num_nodes_image());

        const Index n = graph.get_num_nodes_domain();
        const Index* gptr = graph.get_domain_ptr();
        const Index* gidx = graph.get_image_idx();

        // count the number of edges
        Index m = Index(0);
        for(Index i(0); i < n; ++i)
        {
          for(Index j(gptr[i]); j < gptr[i+1]; ++j)
            m += Index(gidx[j] > i ? 1 : 0);
        }

        // allocate our graph arrays
        _graph = Adjacency::Graph(n, n, m);
        Index* xptr = _graph.get_domain_ptr();
        Index* xidx = _graph.get_image_idx();

        // build edge array
        xptr[0] = Index(0);
        for(Index i(0); i < n; ++i)
        {
          Index k(xptr[i]);
          for(Index j(gptr[i]); j < gptr[i+1]; ++j)
          {
            if(gidx[j] > i)
            {
              xidx[k] = gidx[j];
              ++k;
            }
          }
          xptr[i+1] = k;
        }
      }

      UniEdgeSet(UniEdgeSet&& other) :
        _graph(std::forward<Graph>(other._graph))
      {
      }

      virtual ~UniEdgeSet()
      {
      }

      Index find_edge(Index v0, Index v1) const
      {
        Index i = (v0 < v1 ? v0 : v1);
        Index j = (v0 < v1 ? v1 : v0);

        const Index* xptr = _graph.get_domain_ptr();
        const Index* xidx = _graph.get_image_idx();

        for(Index k(xptr[i]); k < xptr[i+1]; ++k)
        {
          if(xidx[k] == j)
            return k;
        }
        return ~Index(0);
      }
    };
  } // namespace Adjacency
} // namespace FEAT

#endif // KERNEL_ADJACENCY_UNI_EDGE_SET_HPP
