// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2022 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/parti_parmetis.hpp>

#ifdef FEAT_HAVE_PARMETIS
#include <parmetis.h>

using namespace FEAT;
using namespace Geometry;

namespace FEAT
{
  namespace Geometry
  {
    PartiParMETIS::PartiParMETIS(const Dist::Comm& comm) :
      _comm(comm),
      _sub_comm(Dist::Comm::null()),
      _num_elems(0u),
      _num_parts(0u),
      _first_elem(0u),
      _num_local_elems(0u),
      _subgraph(),
      _midpoints(),
      _parts(),
      _coloring()
    {
    }

    PartiParMETIS::~PartiParMETIS()
    {
    }

    void PartiParMETIS::_create_subgraph(const Adjacency::Graph& faces_at_elem, const Index num_parts)
    {
      // set number of desired partitions and total number of elements
      this->_num_parts = num_parts;
      this->_num_elems = faces_at_elem.get_num_nodes_domain();

      // first of all, determine how many processes we want to use for Zoltan
      int num_procs = Math::max(1, int(faces_at_elem.get_num_nodes_domain() / min_elems));
      Math::mini(num_procs, max_procs);
      Math::mini(num_procs, _comm.size());

      // do we need to create a sub-communicator for ParMETIS?
      if(num_procs < _comm.size())
        _sub_comm = _comm.comm_create_range_incl(num_procs);
      else
        _sub_comm = _comm.comm_dup();

      // get my rank and size
      const int z_rank = _sub_comm.rank();
      const int z_size = _sub_comm.size();

      // make sure that the number of desired parts is not smaller than the size of the communicator
      XASSERTM(Index(z_size) <= num_parts, "thou shall not create less partitions than thou hast ranks");

      // does this process not participate in the sub-communicator?
      if(z_size == 0)
      {
        this->_first_elem = this->_num_local_elems = Index(0);
        return;
      }

      // get number of elements and faces
      const Index num_elems = faces_at_elem.get_num_nodes_domain();
      const Index num_faces = faces_at_elem.get_num_nodes_image();

      // determine the first and the last element of this rank
      Index first_elem = (Index(z_rank) * num_elems) / Index(z_size);
      Index end_elem = (Index(z_rank+1) * num_elems) / Index(z_size);
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
        my_dom_ptr[i] = dom_ptr[first_elem + i] - first_ptr;
      for(Index i(0); i < my_num_idx; ++i)
        my_img_idx[i] = img_idx[first_ptr + i];

      // transpose graph for this rank
      Adjacency::Graph elems_at_faces(Adjacency::RenderType::transpose, faces_at_elem);

      // combine graph into our sub-hypergraph
      this->_first_elem = first_elem;
      this->_num_local_elems = my_num_elems;
      this->_subgraph = Adjacency::Graph(Adjacency::RenderType::injectify_sorted, faces_at_my_elem, elems_at_faces);
      this->_parts.resize(my_num_elems);
    }

    bool PartiParMETIS::_execute()
    {
      // if the size of the sub-communicator is zero, then this process does not participate
      // in the actual partitioning process and it only needs to participate in the broadcast
      bool is_ok = true;
      if(this->_sub_comm.size() > 0)
      {
        // apply Zoltan partitioner
        is_ok = this->_apply_parmetis();
      }

      // broadcast coloring from root rank and return
      return this->_broadcast_coloring(is_ok);
    }

    bool PartiParMETIS::_apply_parmetis()
    {
      // convert all inputs to ParMETIS formats
      idx_t wgtflag(0), numflag(0), ncon(1), edgecut(0);
      idx_t nparts = idx_t(_num_parts);
      idx_t ndims = idx_t(_midpoints.size() / _num_local_elems);
      std::vector<idx_t> options(4, 0);
      std::vector<real_t> tpwgts(_num_parts, real_t(1) / real_t(_num_parts));
      std::vector<real_t> ubvec(1, real_t(1.05));
      std::vector<idx_t> part(_num_local_elems, 0);

      // convert midpoints to xyz
      std::vector<real_t> xyz;
      if(!_midpoints.empty())
      {
        xyz.resize(_midpoints.size());
        for(std::size_t i(0); i < _midpoints.size(); ++i)
          xyz[i] = real_t(_midpoints[i]);
      }

      // set up vertex distribution array
      const int z_size = this->_sub_comm.size();
      std::vector<idx_t> vtxdist;
      for(int i(0); i <= z_size; ++i)
        vtxdist.push_back(idx_t(_num_elems*i) / idx_t(z_size));

      // create adjacency arrays from our sub-graph without the self-adjacency
      std::vector<idx_t> xadj, adjncy;
      const Index nv = _subgraph.get_num_nodes_domain();
      const Index ne = _subgraph.get_num_indices();
      xadj.reserve(nv + 1u);
      adjncy.reserve(ne - nv);
      const Index* dom_ptr = _subgraph.get_domain_ptr();
      const Index* img_idx = _subgraph.get_image_idx();
      for(Index i(0); i < _subgraph.get_num_nodes_domain(); ++i)
      {
        xadj.push_back(idx_t(adjncy.size()));
        for(Index k(dom_ptr[i]); k < dom_ptr[i+1]; ++k)
        {
          if(img_idx[k] != _first_elem + i)
            adjncy.push_back(idx_t(img_idx[k]));
        }
      }
      xadj.push_back(idx_t(adjncy.size()));

      // call ParMETIS
      int rtn = ParMETIS_V3_PartGeomKway(
        vtxdist.data(),   // in: first dof for each processor, on all processes, length P+1
        xadj.data(),      // in: dom_ptr of adjacency Graph (w/o self-adjacency)
        adjncy.data(),    // in: img_idx of adjacency Graph (w/o self-adjacency)
        nullptr,          // in: vertex weights (may be set to nullptr)
        nullptr,          // in: edge weights (may be set to nullptr)
        &wgtflag,         // in: 0 for no weights
        &numflag,         // in: 0 for C numbering
        &ndims,           // in: dimension of mesh
        xyz.data(),       // in: vertex coords (element midpoints), process-local only
        &ncon,            // in: =1
        &nparts,          // in: number of desired partitions
        tpwgts.data(),    // in: array of 1/nparts
        ubvec.data(),     // in: array of 1.05
        options.data(),   // parameter array
        &edgecut,         // out: number of cut edges
        part.data(),      // out: partition for each local node
        &_sub_comm.mpi_comm() // in: self explanatory
      );

      // convert parts
      _parts.resize(_num_local_elems);
      for(Index i(0); i < _num_local_elems; ++i)
        _parts[i] = Index(part[i]);

      // gather the coloring
      return _gather_coloring() && (rtn == METIS_OK);
    }

    bool PartiParMETIS::_gather_coloring()
    {
      // get my rank and size
      const int z_rank = _sub_comm.rank();
      const Index z_size = Index(_sub_comm.size());

      // determine maximum export size
      Index max_local_elems(0);
      _sub_comm.allreduce(&this->_num_local_elems, &max_local_elems, 1u, Dist::op_max);
      max_local_elems += 1; // for actual size

      // allocate send buffer
      std::vector<Index> send_buf(max_local_elems);
      send_buf[0] = _num_local_elems;
      for(Index i(0); i < _num_local_elems; ++i)
        send_buf[i+1] = this->_parts[i];

      // allocate receive buffer and gather coloring at rank 0
      std::vector<Index> recv_buf(z_rank == 0 ? max_local_elems * z_size : std::size_t(0));
      _sub_comm.gather(send_buf.data(), max_local_elems, recv_buf.data(), max_local_elems, 0);

      // all ranks > 0 are done here
      if(z_rank > 0)
        return true;

      // compute total element count
      Index num_elems(0);
      for(Index i(0); i < z_size; ++i)
        num_elems += recv_buf[i*max_local_elems];

      // allocate ranks vector
      this->_coloring = Adjacency::Coloring(this->_num_elems, this->_num_parts);
      Index* col = this->_coloring.get_coloring();

      // build ranks vector
      for(Index i(0), k(0); i < z_size; ++i)
      {
        Index off = i*max_local_elems;
        Index n = recv_buf[off];
        for(Index  j(1); j <= n; ++j, ++k)
          col[k] = recv_buf[off+j];
      }

      // build color counts
      std::vector<int> color_counts(this->_num_parts, 0);
      for(Index i(0); i < this->_num_elems; ++i)
        ++color_counts[col[i]];

      // make sure we have at least 1 element per color/partition
      for(Index i(0); i < this->_num_parts; ++i)
      {
        if(color_counts[i] == 0)
          return false;
      }

      // okay
      return true;
    }

    bool PartiParMETIS::_broadcast_coloring(bool metis_ok)
    {
      // first of all, let's compare our ParMETIS return codes;
      // all processes which did not participate in the partitioning have a 'true' status
      int was_ok(metis_ok ? 0 : 1),  was_all_ok(0);
      this->_comm.allreduce(&was_ok, &was_all_ok, std::size_t(1), Dist::op_sum);
      if(was_all_ok > 0)
        return false; // at least 1 process failed

      // allocate coloring if not on root
      if(this->_comm.rank() > 0)
        this->_coloring = Adjacency::Coloring(this->_num_elems, this->_num_parts);

      // broadcast coloring
      this->_comm.bcast(this->_coloring.get_coloring(), this->_num_elems, 0);

      // okay
      return true;
    }

  } // namespace Geometry
} // namespace FEAT

#else // no FEAT_HAVE_PARMETIS

// dummy function to suppress linker warnings
void parmetis_not_linked() {}

#endif // FEAT_HAVE_PARMETIS
