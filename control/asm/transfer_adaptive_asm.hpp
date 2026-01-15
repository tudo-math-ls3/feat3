// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include "kernel/adjacency/adjactor.hpp"
#include "kernel/adjacency/base.hpp"
#include "kernel/adjacency/graph.hpp"
#include "kernel/geometry/adaptive_mesh.hpp"
#include "kernel/geometry/adaptive_mesh_layer.hpp"
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>

namespace FEAT::Control::Asm
{
  /**
  * \brief Assembles a scalar transfer for a pair of system levels
  *
  * \param[in] virt_lvl_f
  * A \transient reference to the virtual fine domain level for which the transfer is to be assembled
  *
  * \param[in] virt_lvl_c
  * A \transient reference to the virtual coarse domain level for which the transfer is to be assembled
  *
  * \param[in] cubature
  * The name of the cubature rule to be used for integration.
  *
  * \param[in] trunc
  * Specifies whether the truncation operator is also to be assembled along with the prolongation and restriction
  * operators.
  *
  * \param[in] shrink
  * Specifies whether the transfer matrices are to be shrunk, i.e. whether all entries whose absolute value is
  * smaller than 1e-3 times the largest absolute value are to be removed from the matrices.
  *
  * \param[in] space_lambda
  * A lambda expression that takes a const reference to a DomainLevel_ and returns a const reference to the
  * finite element space on that level.
  *
  * \param[out] transfer
  * A \transient reference to the local transfer object that is to be assembled.
  *
  * \param[in] muxer
  * A \transient reference to the coarse muxer on the fine level.
  *
  * \param[in] gate_f
  * A \transient reference to the fine level gate.
  *
  * \param[in] gate_c
  * A \transient reference to the coarse level gate, if this process participates in the coarse level,
  * otherwise \transient a reference to the same object as gate_f, thus indicating that this process does not
  * participate in the parent layer.
  */
  template<typename DomainLevel_, typename SpaceLambda_, typename Transfer_, typename Muxer_, typename Gate_>
  void asm_transfer_adaptive_scalar(
  const Domain::VirtualLevel<DomainLevel_>& virt_lvl_f,
  const Domain::VirtualLevel<DomainLevel_>& virt_lvl_c,
  const String& cubature, bool trunc, bool shrink,
  SpaceLambda_&& space_lambda, Transfer_& transfer, const Muxer_& muxer, const Gate_& gate_f, const Gate_& gate_c)
  {
  using DataType = typename Transfer_::DataType;
  using IndexType = typename Transfer_::IndexType;
  using ScalarVectorType = LAFEM::DenseVector<DataType, IndexType>;

  // get fine and coarse domain levels
  const DomainLevel_& level_f = *virt_lvl_f;
  const DomainLevel_& level_c = virt_lvl_c.is_child() ? virt_lvl_c.level_c() : *virt_lvl_c;

  // get fine and coarse spaces
  const auto& space_f = *space_lambda(level_f);
  const auto& space_c = *space_lambda(level_c);

  // get local transfer matrices
  auto& loc_prol = transfer.get_mat_prol();
  auto& loc_rest = transfer.get_mat_rest();
  auto& loc_trunc = transfer.get_mat_trunc();

  // assemble structure?
  if (loc_prol.empty())
  {
    if(level_c.mesh_layer && level_f.mesh_layer)
    {
      // Both levels are local refinements
      Geometry::AdaptiveChildMapping<typename DomainLevel_::AdaptiveMeshLayerType> c2f_mapping(*level_c.mesh_layer, *level_f.mesh_layer);
      Assembly::SymbolicAssembler::assemble_matrix_intermesh(loc_prol, space_f, space_c, c2f_mapping);
    }
    else if(!level_c.mesh_layer && level_f.mesh_layer)
    {
      // Coarse level is foundation mesh
      Geometry::FoundationAdaptiveAdjactor<typename DomainLevel_::AdaptiveMeshType> foundation_to_adaptive(level_f.mesh_layer->adaptive_mesh());
      typename DomainLevel_::AdaptiveMeshLayerType layer_zero(level_f.mesh_layer->adaptive_mesh_ptr(), Geometry::Layer{0});
      Geometry::AdaptiveChildMapping<typename DomainLevel_::AdaptiveMeshLayerType> zero_to_one(layer_zero, *level_f.mesh_layer);

      Adjacency::CompositeAdjactor adj(foundation_to_adaptive, zero_to_one);

      Assembly::SymbolicAssembler::assemble_matrix_intermesh(loc_prol, space_f, space_c, adj);
    }
    else
    {
      XASSERT(!level_c.mesh_layer);
      XASSERT(!level_f.mesh_layer);
      // Both level are regular refinements
      Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_f, space_c);
    }
  }

  // create local weight vectors
  auto loc_vec_weight_p = loc_prol.create_vector_l();
  auto loc_vec_weight_t = loc_prol.create_vector_r();

  // assemble prolongation matrix
  loc_prol.format();
  loc_vec_weight_p.format();
  loc_vec_weight_t.format();

  // if we need a truncation matrix, then compute its sparsity pattern now
  if(trunc)
      loc_trunc.transpose(loc_prol);

  // assemble prolongation matrix
  if(level_c.mesh_layer && level_f.mesh_layer)
  {
    // Both levels are local refinements
    Geometry::AdaptiveChildMapping<typename DomainLevel_::AdaptiveMeshLayerType> c2f_mapping(*level_c.mesh_layer, *level_f.mesh_layer);
    Assembly::SymbolicAssembler::assemble_matrix_intermesh(loc_prol, space_f, space_c, c2f_mapping);

    Adjacency::Graph f2c_mapping(Adjacency::RenderType::transpose, c2f_mapping);

    Assembly::GridTransfer::assemble_intermesh_transfer(loc_prol, loc_vec_weight_p, space_f, space_c, f2c_mapping, cubature);

    if(trunc)
      Assembly::GridTransfer::assemble_intermesh_transfer(loc_trunc, loc_vec_weight_t, space_f, space_c, f2c_mapping, cubature);
  }
  else if(!level_c.mesh_layer && level_f.mesh_layer)
  {
    // Coarse level is foundation mesh
    Geometry::FoundationAdaptiveAdjactor<typename DomainLevel_::AdaptiveMeshType> foundation_to_adaptive(level_f.mesh_layer->adaptive_mesh());
    typename DomainLevel_::AdaptiveMeshLayerType layer_zero(level_f.mesh_layer->adaptive_mesh_ptr(), Geometry::Layer{0});
    Geometry::AdaptiveChildMapping<typename DomainLevel_::AdaptiveMeshLayerType> zero_to_one(layer_zero, *level_f.mesh_layer);

    Adjacency::CompositeAdjactor adj(foundation_to_adaptive, zero_to_one);
    Adjacency::Graph adj_transpose(Adjacency::RenderType::transpose, adj);

    Assembly::GridTransfer::assemble_intermesh_transfer(loc_prol, loc_vec_weight_p, space_f, space_c, adj_transpose, cubature);

    if(trunc)
      Assembly::GridTransfer::assemble_intermesh_transfer(loc_trunc, loc_vec_weight_t, space_f, space_c, adj_transpose, cubature);
  }
  else
  {
    // Both level are regular refinements
    Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight_p, space_f, space_c, cubature);
    if(trunc)
        Assembly::GridTransfer::assemble_truncation(loc_trunc, loc_vec_weight_t, space_f, space_c, cubature);
  }

  // synchronize weight vector using the fine level gate
  gate_f.sync_0(loc_vec_weight_p);

  // invert components
  loc_vec_weight_p.component_invert(loc_vec_weight_p);

  // scale prolongation matrix
  loc_prol.scale_rows(loc_prol, loc_vec_weight_p);

  // shrink prolongation matrix if desired
  if(shrink)
      loc_prol.shrink(DataType(1E-3) * loc_prol.max_abs_element());

  // copy and transpose
  loc_rest = loc_prol.transpose();

  // do we need a truncation operator?
  if(!trunc)
      return;

  // We now need to synchronize the weight vector in analogy to the prolongation matrix assembly.
  // Note that the weight vector is now a coarse-level vector, so we need to synchronize over
  // the coarse-level gate. This may be a bit more complicated if the coarse level is a ghost
  // level, as in this case we have to join/split around the sync operation.

  // synchronize weight vector using the muxer/gate
  if(!virt_lvl_c.is_child())
  {
      // The coarse level is a simple (non-child) level that exists on all processes,
      // so simply sync over the coarse-level gate:
      gate_c.sync_0(loc_vec_weight_t);
  }
  else if(virt_lvl_c.is_parent())
  {
      // The coarse level is a child level and this is one of the parent processes which contain
      // the coarse-level gate. So we first need to join the weight vector onto the parent processes
      // first, then sync that joined vector over the parent gate and finally split the result
      // to all child processes -- this "emulates" a sync over the (non-existent) coarse-level
      // child gate, which is what we actually require here...

      // create temporary vector on parent partitioning
      ScalarVectorType loc_tmp(gate_c.get_freqs().size());

      // join child weights over muxer
      muxer.join(loc_vec_weight_t, loc_tmp);

      // sync over coarse gate
      gate_c.sync_0(loc_tmp);

      // split over muxer
      muxer.split(loc_vec_weight_t, loc_tmp);
  }
  else // ghost
  {
      // The coarse level is a ghost level, i.e. a child but not a parent. In this case, we
      // only have to participate in the join/send operations of the muxer, which are part
      // of the operations executed on the parents handled by the else-if case above.

      muxer.join_send(loc_vec_weight_t);

      // parent performs sync over its gate here (see above else-if)

      muxer.split_recv(loc_vec_weight_t);
  }

  // invert components
  loc_vec_weight_t.component_invert(loc_vec_weight_t);

  // scale prolongation matrix
  loc_trunc.scale_rows(loc_trunc, loc_vec_weight_t);

  // shrink prolongation matrix if desired
  if(shrink)
      loc_trunc.shrink(DataType(1E-3) * loc_trunc.max_abs_element());
  }
} // namespace FEAT::Control::Asm
