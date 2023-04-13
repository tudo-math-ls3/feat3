// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/parti_zoltan.hpp>

#ifdef FEAT_HAVE_ZOLTAN
#include <zoltan.h>

using namespace FEAT;
using namespace Geometry;

/**
 * \brief Callback function for Zoltan's <c>ZOLTAN_NUM_OBJ_FN</c> interface
 *
 * This functions returns the number of elements that have been assigned to this process.
 *
 * \param[inout] data
 * A pointer to the PartiZoltan::Hypergraph object containing the hypergraph for this process.
 *
 * \param[out] ierr
 * A pointer to the object that receives the error code
 *
 * \returns The number of elements currently assigned to this process.
 */
extern "C" static int feat_zoltan_num_elems(void* data, int* ierr)
{
  *ierr = ZOLTAN_OK;
  const Adjacency::Graph& grp = reinterpret_cast<const PartiZoltan::Hypergraph*>(data)->sub_graph;
  return int(grp.get_num_nodes_domain());
}

/**
 * \brief Callback function for Zoltan's <c>ZOLTAN_OBJ_LIST_FN </c> interface
 *
 * This functions specifies the global IDs of the elements assigned to this process.
 *
 * \param[inout] data
 * A pointer to the PartiZoltan::Hypergraph object containing the hypergraph for this process.
 *
 * \param[out] ierr
 * A pointer to the object that receives the error code
 */
extern "C" static void feat_zoltan_elem_ids(void* data, int size_gid, int size_lid,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR DOXY(local_id), int wgt_dim, float* DOXY(obj_wgts), int* ierr)
{
  XASSERT(size_gid == 1);
  XASSERT(size_lid == 0);
  XASSERT(wgt_dim == 0);
  const PartiZoltan::Hypergraph& hyg = *reinterpret_cast<const PartiZoltan::Hypergraph*>(data);
  const Index n = hyg.sub_graph.get_num_nodes_domain();
  for(Index i(0); i < n; ++i)
    global_id[i] = ZOLTAN_ID_TYPE(hyg.first_elem + i);
  *ierr = ZOLTAN_OK;
}

/**
 * \brief Callback function for Zoltan's <c>ZOLTAN_HG_SIZE_CS_FN</c> interface
 *
 * \param[inout] data
 * A pointer to the PartiZoltan::Hypergraph object containing the hypergraph for this process.
 *
 * \param[out] ierr
 * A pointer to the object that receives the error code
 */
extern "C" static void feat_zoltan_hypergraph_size(void* data, int* num_lists, int* num_nonzeroes, int* format, int* ierr)
{
  *ierr = ZOLTAN_OK;
  const Adjacency::Graph& grp = reinterpret_cast<const PartiZoltan::Hypergraph*>(data)->sub_graph;
  *num_lists = int(grp.get_num_nodes_domain());
  *num_nonzeroes = int(grp.get_num_indices());
  *format = ZOLTAN_COMPRESSED_EDGE;
  *ierr = ZOLTAN_OK;
}

/**
 * \brief Callback function for Zoltan's <c>ZOLTAN_HG_CS_FN</c> interface
 *
 * \param[inout] data
 * A pointer to the PartiZoltan::Hypergraph object containing the hypergraph for this process.
 *
 * \param[out] ierr
 * A pointer to the object that receives the error code
 */
extern "C" static void feat_zoltan_hypergraph_data(void* data, int size_gid, int num_edges, int num_nonzeroes,
  int format, ZOLTAN_ID_PTR edgeGID, int* vtxPtr, ZOLTAN_ID_PTR vtxGID, int* ierr)
{
  const Adjacency::Graph& grp = reinterpret_cast<const PartiZoltan::Hypergraph*>(data)->sub_graph;
  XASSERT(size_gid == 1);
  XASSERT(num_edges == int(grp.get_num_nodes_domain()));
  XASSERT(num_nonzeroes == int(grp.get_num_indices()));
  XASSERT(format == ZOLTAN_COMPRESSED_EDGE);

  const Index* dom_ptr = grp.get_domain_ptr();
  const Index* img_idx = grp.get_image_idx();

  // copy domain ptr = hyperedge ptr
  for(int i(0); i < num_edges; ++i)
  {
    edgeGID[i] = ZOLTAN_ID_TYPE(i);
    vtxPtr[i] = int(dom_ptr[i]);
  }

  // copy image idx = vertex idx
  for(int i(0); i < num_nonzeroes; ++i)
  {
    vtxGID[i] = ZOLTAN_ID_TYPE(img_idx[i]);
  }
  *ierr = ZOLTAN_OK;
}

namespace FEAT
{
  namespace Geometry
  {
    PartiZoltan::PartiZoltan(const Dist::Comm& comm) :
      _comm(comm),
      _zoltan_comm(Dist::Comm::null()),
      _hypergraph(),
      _num_elems(0u),
      _num_parts(0u),
      _coloring(),
      _zz(nullptr)
    {
    }

    PartiZoltan::~PartiZoltan()
    {
      if(_zz != nullptr)
      {
        Zoltan_Destroy(reinterpret_cast<Zoltan_Struct**>(&_zz));
        _zz = nullptr;
      }
    }

    bool PartiZoltan::execute(const Adjacency::Graph& faces_at_elem, const Index num_parts)
    {
      // set number of desired partitions and total number of elements
      this->_num_parts = num_parts;
      this->_num_elems = faces_at_elem.get_num_nodes_domain();

      // first of all, determine how many processes we want to use for Zoltan
      int num_procs = Math::max(1, int(faces_at_elem.get_num_nodes_domain() / min_elems));
      Math::mini(num_procs, max_procs);
      Math::mini(num_procs, _comm.size());

      // do we need to create a sub-communicator for Zoltan?
      if(num_procs < _comm.size())
        _zoltan_comm = _comm.comm_create_range_incl(num_procs);
      else
        _zoltan_comm = _comm.comm_dup();

      // get my rank and size
      //const int z_rank = _zoltan_comm.rank();
      const int z_size = _zoltan_comm.size();

      // make sure that the number of desired parts is not smaller than the size of the communicator
      XASSERTM(Index(z_size) <= num_parts, "thou shall not create less partitions than thou hast ranks");

      // if the size of the Zoltan communicator is zero, then this process does not participate
      // in the actual partitioning process and it only needs to participate in the broadcast
      bool is_ok = true;
      if(z_size > 0)
      {
        // create our sub-hypergraph
        this->_create_hypergraph(faces_at_elem);

        // apply Zoltan partitioner
        is_ok = this->_apply_zoltan();
      }

      // broadcast coloring from root rank and return
      return this->_broadcast_coloring(is_ok);
    }

    void PartiZoltan::_create_hypergraph(const Adjacency::Graph& faces_at_elem)
    {
      // get number of elements and faces
      const Index num_elems = faces_at_elem.get_num_nodes_domain();
      const Index num_faces = faces_at_elem.get_num_nodes_image();

      // get my rank and size
      const int z_rank = _zoltan_comm.rank();
      const int z_size = _zoltan_comm.size();

      // this should not happen
      XASSERT(z_size > 0);

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
      this->_hypergraph.first_elem = first_elem;
      this->_hypergraph.sub_graph = Adjacency::Graph(Adjacency::RenderType::injectify_sorted, faces_at_my_elem, elems_at_faces);
    }

    bool PartiZoltan::_apply_zoltan()
    {
      // create Zoltan structure
      Zoltan_Struct* zz = Zoltan_Create(_zoltan_comm.mpi_comm());
      this->_zz = zz;

      // define Zoltan parameters
      Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
      Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");
      Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG");
      Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
      Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
      Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
      Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");
      Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", stringify(this->_num_parts).c_str());
      Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
      Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");
      //Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
      Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", "0.5");

      // set hypergraph interface callbacks
      Zoltan_Set_Num_Obj_Fn(zz, feat_zoltan_num_elems, &this->_hypergraph);
      Zoltan_Set_Obj_List_Fn(zz, feat_zoltan_elem_ids, &this->_hypergraph);
      Zoltan_Set_HG_Size_CS_Fn(zz, feat_zoltan_hypergraph_size, &this->_hypergraph);
      Zoltan_Set_HG_CS_Fn(zz, feat_zoltan_hypergraph_data, &this->_hypergraph);

      int changes(0), num_gid_entries(0), num_lid_entries(0), num_import(0), num_export(0);
      ZOLTAN_ID_PTR import_global_ids(nullptr);
      ZOLTAN_ID_PTR import_local_ids(nullptr);
      int* import_procs(nullptr);
      int* import_parts(nullptr);
      ZOLTAN_ID_PTR export_global_ids(nullptr);
      ZOLTAN_ID_PTR export_local_ids(nullptr);
      int* export_procs(nullptr);
      int* export_parts(nullptr);

      int rtn = Zoltan_LB_Partition(zz, &changes, &num_gid_entries, &num_lid_entries,
        &num_import, &import_global_ids, &import_local_ids, &import_procs, &import_parts,
        &num_export, &export_global_ids, &export_local_ids, &export_procs, &export_parts);

      // create coloring from Zoltan partitioning
      bool gather_ok = true;
      if(rtn == ZOLTAN_OK)
        gather_ok = this->_gather_coloring(num_export, export_parts);

      Zoltan_LB_Free_Part(&import_global_ids, &import_local_ids, &import_procs, &import_parts);
      Zoltan_LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs, &export_parts);

      Zoltan_Destroy(&zz);
      this->_zz = nullptr;

      return gather_ok && (rtn == ZOLTAN_OK);
    }

    bool PartiZoltan::_gather_coloring(const int num_export, const int* export_parts)
    {
      // get my rank and size
      const int z_rank = _zoltan_comm.rank();
      const int z_size = _zoltan_comm.size();

      // determine maximum export size
      int max_export(0);
      _zoltan_comm.allreduce(&num_export, &max_export, 1u, Dist::op_max);
      max_export += 1; // for actual size

      // allocate send buffer
      std::vector<int> send_buf(max_export);
      send_buf[0] = num_export;
      for(int i(0); i < num_export; ++i)
        send_buf[i+1] = export_parts[i];

      // allocate receive buffer and gather coloring at rank 0
      std::vector<int> recv_buf(z_rank == 0 ? max_export * z_size : std::size_t(0));
      _zoltan_comm.gather(send_buf.data(), max_export, recv_buf.data(), max_export, 0);

      // all ranks > 0 are done here
      if(z_rank > 0)
        return true;

      // compute total element count
      int num_elems(0);
      for(int i(0); i < z_size; ++i)
        num_elems += recv_buf[i*max_export];

      // allocate ranks vector
      this->_coloring = Adjacency::Coloring(num_elems, this->_num_parts);
      Index* col = this->_coloring.get_coloring();

      // build ranks vector
      for(int i(0), k(0); i < z_size; ++i)
      {
        int off = i*max_export;
        int n = recv_buf[off];
        for(int j(1); j <= n; ++j, ++k)
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

    bool PartiZoltan::_broadcast_coloring(bool zoltan_ok)
    {
      // first of all, let's compare our Zoltan return codes;
      // all processes which did not participate in the partitioning have a 'true' status
      int was_ok(zoltan_ok ? 0 : 1),  was_all_ok(0);
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

#else // no FEAT_HAVE_ZOLTAN

// dummy function to suppress linker warnings
void zoltan_did_not_make_it_to_the_feat_parti() {}

#endif // FEAT_HAVE_ZOLTAN
