// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_ASSEMBLY_TRANSFER_ASM_HPP
#define CONTROL_ASSEMBLY_TRANSFER_ASM_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>

#include <control/domain/domain_control.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Asm
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
      void asm_transfer_scalar(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_f,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_c,
        const String& cubature, bool trunc, bool shrink,
        SpaceLambda_&& space_lambda, Transfer_& transfer, const Muxer_& muxer, const Gate_& gate_f, const Gate_& gate_c)
      {
        typedef typename Transfer_::DataType DataType;
        typedef typename Transfer_::IndexType IndexType;
        typedef LAFEM::DenseVector<DataType, IndexType> ScalarVectorType;

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
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_f, space_c);

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
        Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight_p, space_f, space_c, cubature);

        // assemble truncation matrix if desired
        if(trunc)
          Assembly::GridTransfer::assemble_truncation(loc_trunc, loc_vec_weight_t, space_f, space_c, cubature);

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

      /**
       * \brief Assembles a blocked transfer for a pair of system levels
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
      void asm_transfer_blocked(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_f,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_c,
        const String& cubature, bool trunc, bool shrink,
        SpaceLambda_&& space_lambda, Transfer_& transfer, const Muxer_& muxer, const Gate_& gate_f, const Gate_& gate_c)
      {
        typedef typename Transfer_::DataType DataType;
        typedef typename Transfer_::IndexType IndexType;
        typedef typename Gate_::MirrorType MirrorType;
        typedef LAFEM::DenseVector<DataType, IndexType> ScalarVectorType;

        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_f;
        const DomainLevel_& level_c = virt_lvl_c.is_child() ? virt_lvl_c.level_c() : *virt_lvl_c;

        // get fine and coarse spaces
        const auto& space_f = *space_lambda(level_f);
        const auto& space_c = *space_lambda(level_c);

        // get local transfer matrices
        auto& loc_prol_wrapped = transfer.get_mat_prol();
        auto& loc_rest_wrapped = transfer.get_mat_rest();
        auto& loc_trunc_wrapped = transfer.get_mat_trunc();

        // get the unwrapped types
        auto& loc_prol = loc_prol_wrapped.unwrap();
        auto& loc_rest = loc_rest_wrapped.unwrap();
        auto& loc_trunc = loc_trunc_wrapped.unwrap();

        // assemble structure?
        if (loc_prol.empty())
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_f, space_c);

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
        Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight_p, space_f, space_c, cubature);

        // assemble truncation matrix if desired
        if(trunc)
          Assembly::GridTransfer::assemble_truncation(loc_trunc, loc_vec_weight_t, space_f, space_c, cubature);

        // create a temporary scalar gate
        Global::Gate<ScalarVectorType, MirrorType> scalar_gate;
        scalar_gate.convert(gate_f, loc_vec_weight_p.clone(LAFEM::CloneMode::Allocate));

        // synchronize weight vector using the scalar gate
        scalar_gate.sync_0(loc_vec_weight_p);

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
          XASSERT(&gate_c != &gate_f);

          // The coarse level is a simple (non-child) level that exists on all processes,
          // so simply synchronize over the coarse-level gate:

          // create a temporary scalar gate on the coarse level
          Global::Gate<ScalarVectorType, MirrorType> scalar_gate_c;
          scalar_gate_c.convert(gate_c, ScalarVectorType(space_c.get_num_dofs()));

          // sync over scalar coarse gate
          scalar_gate_c.sync_0(loc_vec_weight_t);
        }
        else if(virt_lvl_c.is_parent())
        {
          XASSERT(&gate_c != &gate_f);

          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then sync that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a sync over the (non-existent) coarse-level
          // child gate, which is what we actually require here...

          // get the number of the coarse level dofs from the gate's frequency vector, which is different
          // than space_c.get_num_dofs(), because the first one is defined on the parent patch, whereas
          // the latter one is defined on the child patch
          const Index num_dofs_c = gate_c.get_freqs().size();

          // create a temporary scalar gate on the coarse level
          Global::Gate<ScalarVectorType, MirrorType> scalar_gate_c;
          scalar_gate_c.convert(gate_c, ScalarVectorType(num_dofs_c));

          // create a temporary scalar muxer on the fine level
          Global::Muxer<ScalarVectorType, MirrorType> scalar_muxer_f;
          scalar_muxer_f.convert(muxer, ScalarVectorType(space_f.get_num_dofs()));

          // create temporary vector on parent partitioning
          ScalarVectorType loc_tmp(num_dofs_c);

          // join child weights over muxer
          scalar_muxer_f.join(loc_vec_weight_t, loc_tmp);

          // sync over scalar coarse gate
          scalar_gate_c.sync_0(loc_tmp);

          // split over muxer
          scalar_muxer_f.split(loc_vec_weight_t, loc_tmp);
        }
        else // ghost
        {
          // The coarse level is a ghost level, i.e. a child but not a parent. In this case, we
          // only have to participate in the join/send operations of the muxer, which are part
          // of the operations executed on the parents handled by the else-if case above.

          // create a temporary scalar muxer on the fine level
          Global::Muxer<ScalarVectorType, MirrorType> scalar_muxer_f;
          scalar_muxer_f.convert(muxer, ScalarVectorType(space_f.get_num_dofs()));

          scalar_muxer_f.join_send(loc_vec_weight_t);

          // parent performs sync over its gate here (see above else-if)

          scalar_muxer_f.split_recv(loc_vec_weight_t);
        }

        // invert components
        loc_vec_weight_t.component_invert(loc_vec_weight_t);

        // scale prolongation matrix
        loc_trunc.scale_rows(loc_trunc, loc_vec_weight_t);

        // shrink prolongation matrix if desired
        if(shrink)
          loc_trunc.shrink(DataType(1E-3) * loc_trunc.max_abs_element());
      }
    } // namespace Asm
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_ASSEMBLY_TRANSFER_ASM_HPP
