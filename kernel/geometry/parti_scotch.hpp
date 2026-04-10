// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>

#if defined(FEAT_HAVE_SCOTCH) || defined(DOXYGEN)

#include <kernel/adjacency/coloring.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/util/dist.hpp>

namespace FEAT::Geometry
{
  class PartiScotch
  {
    /// Main communicator
    const Dist::Comm& _comm;

    /// Communicator for performing partitioning
    Dist::Comm _work_comm;

    /// Desired number of partitions
    Index _num_parts{0};

    /// Global number of elements
    Index _global_num_elements{0};

    /// Local number of elements
    Index _local_num_elements{0};

    /// Mininum number of elements on each rank
    Index _min_elems;

    /// Maximum number of ranks to use _for partitioning_
    int _max_work_procs;

    /// First element of this process
    Index _first_elem{0};

    /// Subgraph for this rank
    Adjacency::Graph _subgraph;

    /// Element partitioning
    std::vector<Index> _parts;

    /// Elements-at-rank coloring
    Adjacency::Coloring _coloring;


  public:
    explicit PartiScotch(const Dist::Comm& comm, Index min_elems = 1000U, int max_work_procs = 1000) :
      _comm(comm),
      _work_comm(Dist::Comm::null()),
      _min_elems(min_elems),
      _max_work_procs(max_work_procs)
    {
      XASSERT(max_work_procs >= 0);
      XASSERT(min_elems > 0U);
    }

    template<int nf_>
    bool execute(
      const IndexSet<nf_>& faces_at_elem,
      Index num_parts,
      const std::vector<Real>& weights)
    {
      Adjacency::Graph graph(Adjacency::RenderType::as_is, faces_at_elem);
      _create_subgraph(graph, num_parts);
      return _execute(weights);
    }

    Adjacency::Graph build_elems_at_rank() const
    {
      return _coloring.create_partition_graph();
    }

  protected:
    void _create_subgraph(const Adjacency::Graph& faces_at_elem, Index num_parts);
    bool _execute(const std::vector<Real>& weights);
    bool _gather_coloring();
    bool _broadcast_coloring(bool ok);
  }; // class PartiScotch
} // namespace FEAT::Geometry

#endif
