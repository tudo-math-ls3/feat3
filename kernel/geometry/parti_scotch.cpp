// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/parti_scotch.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/math.hpp>
#include <numeric>

#if defined(FEAT_HAVE_SCOTCH) || defined(DOXYGEN)

#include <scotch.h>
#include <ptscotch.h>

namespace FEAT::Geometry
{
  /// RAII-wrapper for Scotch graphs
  struct ScotchGraph
  {
    SCOTCH_Dgraph* inner;
    std::vector<SCOTCH_Num> vertloctab;
    std::vector<SCOTCH_Num> edgeloctab;
    std::vector<SCOTCH_Num> veloloctab;

    explicit ScotchGraph(const Dist::Comm& comm) : inner(SCOTCH_dgraphAlloc())
    {
      XASSERT(SCOTCH_dgraphInit(inner, comm.mpi_comm()) == 0);
    }

    ScotchGraph(ScotchGraph&) = delete;
    ScotchGraph(ScotchGraph&&) = default;

    ScotchGraph& operator=(const ScotchGraph&) = delete;
    ScotchGraph& operator=(ScotchGraph&&) = default;

    ~ScotchGraph()
    {
      SCOTCH_dgraphExit(inner);
      SCOTCH_dgraphFree(inner);
    }
  };

  namespace
  {
    // Convert FEAT graph to SCOTCH graph
    [[nodiscard]] bool
    graph_to_scotch(ScotchGraph& scotch, const Adjacency::Graph& graph, Index first_elem, const std::vector<Real>& weights)
    {
      // NOTE: On SCOTCH naming:
      // glb - global, i.e. across all ranks
      // loc - local, i.e. on this rank
      // gst - ghost

      // Base value for all array indexing
      // 0 for C, 1 for Fortran
      static constexpr SCOTCH_Num baseval(0);

      // Array of start indices in edgeloctab and edgegsttab of vertex adjacency sub-arrays
      scotch.vertloctab.clear();

      // Array holding the adjacency array of every local vertex
      scotch.edgeloctab.clear();

      {
        const Index* dom_ptr = graph.get_domain_ptr();
        const Index* img_idx = graph.get_image_idx();
        for(Index i(0); i < graph.get_num_nodes_domain(); ++i)
        {
          scotch.vertloctab.push_back(SCOTCH_Num(scotch.edgeloctab.size()));
          for(Index k(dom_ptr[i]); k < dom_ptr[i+1]; ++k)
          {
            if(img_idx[k] != first_elem + i)
            {
              scotch.edgeloctab.push_back(SCOTCH_Num(img_idx[k]));
            }
          }
        }
        scotch.vertloctab.push_back(SCOTCH_Num(scotch.edgeloctab.size()));
      }

      // Array of after-last indices in edgeloctab and edgegsttab of vertex adjacency sub-arrays.
      // Because our edgeloctab is stored contiguously, we do not have to allocate memory for this,
      // and can instead reuse vertloctab with an offset of one.
      SCOTCH_Num* const vendloctab = scotch.vertloctab.data() + 1;

      // Number of local vertices borne by this process
      const auto vertlocnbr = SCOTCH_Num(graph.get_num_nodes_domain());

      XASSERT(vertlocnbr == SCOTCH_Num(scotch.vertloctab.size() - 1));

      // Maximum number of local vertices to be created on this process.
      // Must be set to vertlocnbr for graphs without holes in their global
      // numbering.
      const SCOTCH_Num vertlocmax = vertlocnbr;

      // Local number of arcs (directed edges), including arcs to local vertices as well as to ghost vertices.
      const auto edgelocnbr = SCOTCH_Num(scotch.edgeloctab.size());

      XASSERT(scotch.vertloctab.back() == edgelocnbr);

      // Lower bounded by the minimum size of teh edge array required to encompass all used adjacency values;
      // it can be set to edgelocnbr if the edge array is compact (which our edge array is)
      const SCOTCH_Num edgelocsiz = edgelocnbr;

      XASSERT(scotch.vertloctab.front() == baseval);

      // Local vertex load array, holding integer load associated with every vertex
      scotch.veloloctab.clear();

      // set up weights if desired
      if(!weights.empty())
      {
        XASSERT(weights.size() == std::size_t(vertlocnbr));

        scotch.veloloctab.resize(std::size_t(vertlocnbr));

        // SCOTCH expects integer weights, but we have float weights, so we first have to compute
        // the total sum of all weights and then choose a scaling factor so that the sum of the
        // scaled weights is a large integer not greater than 10^9; then we can simply multiply
        // the weights by that scaling factor and convert to integers.

        // compute the total sum of all weights
        Real wsum = std::accumulate(weights.begin(), weights.end(), Real(0));
        XASSERTM(wsum > 1e-3, "weight sum appears to be zero");

        // compute integer weights for METIS
        Real scale = 1E+9 / wsum;
        for(std::size_t i(0); i < std::size_t(vertlocnbr); ++i)
        {
          scotch.veloloctab[i] = SCOTCH_Num(scale * weights[i]);
        }
      }

      // Optional array of local vertex labels
      SCOTCH_Num* const vlblloctab = nullptr;

      // Optional array holding the local and ghost indices of neighbors of local vertices.
      SCOTCH_Num* const edgegsttab = nullptr;

      // Optional array holding the integer load associated with every arc.
      // Matching arcs should always have identical loads. (You hear that, ParMETIS?)
      SCOTCH_Num* const edloloctab = nullptr;

      // NOTE: For PT-SCOTCH to operate properly, either all the arrays of a kind must
      // be set to null on all processes, or else all of them must be non null.

      return SCOTCH_dgraphBuild(
               scotch.inner,
               baseval,
               vertlocnbr,
               vertlocmax,
               scotch.vertloctab.data(),
               vendloctab,
               scotch.veloloctab.empty() ? nullptr : scotch.veloloctab.data(),
               vlblloctab,
               edgelocnbr,
               edgelocsiz,
               scotch.edgeloctab.data(),
               edgegsttab,
               edloloctab) == 0;
    }
  } // namespace

  void PartiScotch::_create_subgraph(const Adjacency::Graph& faces_at_elem, const Index num_parts)
  {
    // set number of desired partitions and total number of elements
    this->_num_parts = num_parts;
    this->_global_num_elements = faces_at_elem.get_num_nodes_domain();

    // first of all, determine how many processes we want to use for Zoltan
    int num_procs = Math::max(1, int(faces_at_elem.get_num_nodes_domain() / _min_elems));
    if(_max_work_procs > 0)
    {
      Math::mini(num_procs, _max_work_procs);
    }
    Math::mini(num_procs, int(num_parts));
    Math::mini(num_procs, _comm.size());

    // do we need to create a sub-communicator for ParMETIS?
    if(num_procs < _comm.size())
    {
      _work_comm = _comm.comm_create_range_incl(num_procs);
    }
    else
    {
      _work_comm = _comm.comm_dup();
    }

    // get my rank and size
    const int z_rank = _work_comm.rank();
    const int z_size = _work_comm.size();

    // make sure that the number of desired parts is not smaller than the size of the communicator
    XASSERTM(Index(z_size) <= num_parts, "thou shall not create less partitions than thou hast ranks");

    // does this process not participate in the sub-communicator?
    if(z_size == 0)
    {
      this->_first_elem = this->_local_num_elements = Index(0);
      return;
    }

    // get number of elements and faces
    const Index num_elems = faces_at_elem.get_num_nodes_domain();
    const Index num_faces = faces_at_elem.get_num_nodes_image();

    // determine the first and the last element of this rank
    Index first_elem = (Index(z_rank) * num_elems) / Index(z_size);
    Index end_elem = (Index(z_rank + 1) * num_elems) / Index(z_size);
    Index my_num_elems = end_elem - first_elem;

    const Index* dom_ptr = faces_at_elem.get_domain_ptr();
    const Index* img_idx = faces_at_elem.get_image_idx();

    // create a sub-graph for this rank
    const Index first_ptr = dom_ptr[first_elem];
    const Index my_num_idx = dom_ptr[end_elem] - first_ptr;
    Adjacency::Graph faces_at_my_elem(my_num_elems, num_faces, my_num_idx);

    Index* my_dom_ptr = faces_at_my_elem.get_domain_ptr();
    Index* my_img_idx = faces_at_my_elem.get_image_idx();
    for(Index i(0); i <= my_num_elems; ++i)
    {
      my_dom_ptr[i] = dom_ptr[first_elem + i] - first_ptr;
    }
    for(Index i(0); i < my_num_idx; ++i)
    {
      my_img_idx[i] = img_idx[first_ptr + i];
    }

    // transpose graph for this rank
    Adjacency::Graph elems_at_faces(Adjacency::RenderType::transpose, faces_at_elem);

    // combine graph into our sub-hypergraph
    this->_first_elem = first_elem;
    this->_local_num_elements = my_num_elems;
    this->_subgraph = Adjacency::Graph(Adjacency::RenderType::injectify_sorted, faces_at_my_elem, elems_at_faces);
    this->_parts.resize(std::size_t(my_num_elems));
  }

  bool PartiScotch::_execute(const std::vector<Real>& weights)
  {
    // If the size of the sub-communicator is zero, then this process does not participate
    // in the actual partitioning process and it only needs to participate in the broadcast
    if(this->_work_comm.size() == 0)
    {
      return this->_broadcast_coloring(true);
    }

    // Create SCOTCH graph
    ScotchGraph scotch_graph(_work_comm);
    if(!graph_to_scotch(scotch_graph, _subgraph, _first_elem,  weights))
    {
      return this->_broadcast_coloring(false);
    }

    // Check graph
    if(SCOTCH_dgraphCheck(scotch_graph.inner) != 0)
    {
      return this->_broadcast_coloring(false);
    }

    // Create strategy for partitioning
    SCOTCH_Strat partitioning_strategy;

    if(SCOTCH_stratInit(&partitioning_strategy) != 0)
    {
      return this->_broadcast_coloring(false);
    }

    // NOTE: Initialized strategy object will be filled with default strategy once
    // passed to any other operation.

    // Partition distributed graph
    std::vector<SCOTCH_Num> partloctab(_subgraph.get_num_nodes_domain());
    if(SCOTCH_dgraphPart(scotch_graph.inner, SCOTCH_Num(_num_parts), &partitioning_strategy, partloctab.data()) != 0)
    {
      SCOTCH_stratFree(&partitioning_strategy);
      return this->_broadcast_coloring(false);
    }

    SCOTCH_stratFree(&partitioning_strategy);

    for(Index i(0); i < _local_num_elements; ++i)
    {
      _parts[std::size_t(i)] = Index(partloctab[std::size_t(i)]);
    }

    // Broadcast coloring if gathering was ok
    return this->_broadcast_coloring(this->_gather_coloring());
  }

  // NOTE: Copied from the ParMETIS partitioner
  bool PartiScotch::_gather_coloring()
  {
    // get my rank and size
    const int scotch_rank = _work_comm.rank();
    const auto scotch_size = Index(_work_comm.size());

    // determine maximum export size
    Index max_local_elems(0);
    _work_comm.allreduce(&this->_local_num_elements, &max_local_elems, 1u, Dist::op_max);
    max_local_elems += 1; // for actual size

    // allocate send buffer
    auto send_buf = std::vector<Index>(std::size_t(max_local_elems));
    send_buf[0] = _local_num_elements;
    for(Index i(0); i < _local_num_elements; ++i)
    {
      send_buf[std::size_t(i + 1)] = _parts[std::size_t(i)];
    }

    // allocate receive buffer and gather coloring at rank 0
    std::vector<Index> recv_buf(scotch_rank == 0 ? std::size_t(max_local_elems * scotch_size) : std::size_t(0));
    _work_comm.gather(send_buf.data(), max_local_elems, recv_buf.data(), max_local_elems, 0);

    // all non-root ranks are done here
    if(scotch_rank > 0)
    {
      return true;
    }

    // allocate coloring on root
    this->_coloring = Adjacency::Coloring(this->_global_num_elements, this->_num_parts);
    Index* col = this->_coloring.get_coloring();

    // build coloring vector
    for(Index i(0), k(0); i < scotch_size; ++i)
    {
      const Index off = i * max_local_elems;
      const Index n = recv_buf[std::size_t(off)];
      for(Index j(1); j <= n; ++j, ++k)
      {
        col[k] = recv_buf[std::size_t(off + j)];
      }
    }

    // make sure we have at least 1 element per color/partition
    std::vector<int> color_counts(std::size_t(this->_num_parts), 0);
    for(Index i(0); i < this->_global_num_elements; ++i)
    {
      XASSERT(col[i] < this->_num_parts);
      ++color_counts[std::size_t(col[i])];
    }
    for(Index i(0); i < this->_num_parts; ++i)
    {
      if(color_counts[std::size_t(i)] == 0)
      {
        return false;
      }
    }

    return true;
  }

  // NOTE: Copied from the ParMETIS partitioner
  bool PartiScotch::_broadcast_coloring(bool ok)
  {
    int was_ok(ok ? 0 : 1);
    int was_all_ok(0);

    this->_comm.allreduce(&was_ok, &was_all_ok, std::size_t(1), Dist::op_sum);
    if(was_all_ok > 0)
    {
      return false;
    }

    if(this->_comm.rank() > 0)
    {
      this->_coloring = Adjacency::Coloring(this->_global_num_elements, this->_num_parts);
    }

    this->_comm.bcast(this->_coloring.get_coloring(), this->_global_num_elements, 0);
    return true;
  }
} // namespace FEAT::Geometry

#else

// dummy function to suppress linker warnings
void scotch_not_linked() {}

#endif
