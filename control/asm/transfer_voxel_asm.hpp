// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>

#include <control/domain/voxel_domain_control.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Asm
    {
      namespace VoxelAux
      {
        /// Deslags a vector, i.e. removes all entries whose indices are not present in the given mirror.
        template<typename Vector_, typename Mirror_>
        static Vector_ deslag_vector(const Vector_& vector, const Mirror_& mirror)
        {
          const Index num_idx = mirror.num_indices();
          Vector_ vector_x(num_idx);
          const auto* v = vector.elements();
          auto* x = vector_x.elements();
          const auto* idx_f = mirror.indices();
          for(Index i(0); i < num_idx; ++i)
            x[i] = v[idx_f[i]];
          return vector_x;
        }

        /// Deslags a matrix by rows, i.e. removes all rows whose row indices are not present in the given mirror.
        template<typename Matrix_, typename Mirror_>
        static Matrix_ deslag_matrix_rows(const Matrix_& matrix, const Mirror_& row_mirror)
        {
          const Index num_idx = row_mirror.num_indices();
          const auto* idx_f = row_mirror.indices();

          const auto* row_ptr_s = matrix.row_ptr();
          const auto* col_idx_s = matrix.col_ind();
          Index num_nze = 0u;
          for(Index i(0); i < num_idx; ++i)
            num_nze += row_ptr_s[idx_f[i] + 1u] - row_ptr_s[idx_f[i]];

          Matrix_ matrix_x(num_idx, matrix.columns(), num_nze);
          auto* row_ptr_x = matrix_x.row_ptr();
          auto* col_idx_x = matrix_x.col_ind();
          auto* val_x = matrix_x.val();
          const auto* val_s = matrix.val();
          row_ptr_x[0] = 0u;
          for(Index i(0); i < num_idx; ++i)
          {
            Index row_s = idx_f[i];
            Index k(row_ptr_x[i]);
            for(auto j(row_ptr_s[row_s]); j < row_ptr_s[row_s+1u]; ++j, ++k)
            {
              col_idx_x[k] = col_idx_s[j];
              val_x[k] = val_s[j];
            }
            row_ptr_x[i+1u] = k;
          }
          return matrix_x;
        }

        /// Deslags a matrix by columns, i.e. removes all columns whose row indices are not present in the given mirror.
        template<typename Matrix_, typename Mirror_>
        static Matrix_ deslag_matrix_cols(const Matrix_& matrix, const Mirror_& col_mirror)
        {
          const Index num_idx = col_mirror.num_indices();
          const auto* idx_f = col_mirror.indices();

          std::vector<Index> col_map(matrix.columns(), ~Index(0));
          for(Index i(0); i < num_idx; ++i)
            col_map[idx_f[i]] = i;

          const Index num_rows = matrix.rows();
          const auto* row_ptr_s = matrix.row_ptr();
          const auto* col_idx_s = matrix.col_ind();
          Index used_elems = matrix.used_elements();
          Index num_nze = 0u;
          for(Index i(0); i < used_elems; ++i)
            num_nze += (col_map[col_idx_s[i]] != ~Index(0));

          Matrix_ matrix_x(num_rows, num_idx, num_nze);
          auto* row_ptr_x = matrix_x.row_ptr();
          auto* col_idx_x = matrix_x.col_ind();
          auto* val_x = matrix_x.val();
          const auto* val_s = matrix.val();
          row_ptr_x[0] = 0u;
          for(Index i(0); i < num_rows; ++i)
          {
            Index k(row_ptr_x[i]);
            for(auto j(row_ptr_s[i]); j < row_ptr_s[i+1u]; ++j)
            {
              if(col_map[col_idx_s[j]] != ~Index(0))
              {
                col_idx_x[k] = col_map[col_idx_s[j]];
                val_x[k] = val_s[j];
                ++k;
              }
            }
            row_ptr_x[i+1u] = k;
          }
          return matrix_x;
        }
      } // namespace VoxelAux

      /**
       * \brief Assembles a scalar transfer operator for a system level defined on a voxel domain
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
      void asm_transfer_voxel_scalar(
        const Control::Domain::VirtualLevel<Control::Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_f,
        const Control::Domain::VirtualLevel<Control::Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_c,
        const String& cubature, bool trunc, bool shrink,
        SpaceLambda_&& space_lambda, Transfer_& transfer, const Muxer_& muxer, const Gate_& gate_f, const Gate_& gate_c)
      {
        typedef typename Transfer_::DataType DataType;
        typedef typename Transfer_::IndexType IndexType;
        typedef typename Gate_::MirrorType MirrorType;
        typedef LAFEM::DenseVector<DataType, IndexType> ScalarVectorType;

        // do we have a slag level?
        const bool have_slag(virt_lvl_f->slag_level);

        // we have: C -> S >> F, where '->' is the prolongation and '>>' is the deslagging operator
        const DomainLevel_& level_f = have_slag ? *virt_lvl_f->slag_level : static_cast<const DomainLevel_&>(*virt_lvl_f);
        const DomainLevel_& level_c = virt_lvl_c.is_child() ? virt_lvl_c.level_c() : *virt_lvl_c;

        //const auto& space_f = level_f.space;
        const auto& space_f = *space_lambda(level_f);
        const auto& space_c = *space_lambda(level_c);

        // get local transfer matrices
        auto& loc_prol = transfer.get_mat_prol();
        auto& loc_rest = transfer.get_mat_rest();
        auto& loc_trunc = transfer.get_mat_trunc();

        // assemble structure?
        if (loc_prol.empty())
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_f, space_c);

          // if we need a truncation matrix, then compute its sparsity pattern now
        if(trunc)
          loc_trunc.transpose(loc_prol);

        // create local weight vectors
        auto loc_vec_weight_p = loc_prol.create_vector_l();
        auto loc_vec_weight_t = loc_prol.create_vector_r();

        // assemble prolongation matrix
        loc_prol.format();
        loc_vec_weight_p.format();
        loc_vec_weight_t.format();

        // assemble prolongation matrix
        Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight_p, space_f, space_c, cubature);

        // assemble truncation matrix if desired
        if(trunc)
          Assembly::GridTransfer::assemble_truncation(loc_trunc, loc_vec_weight_t, space_f, space_c, cubature);

        // do we have a slag level?
        if(have_slag)
        {
          // assemble a slag mirror
          const auto* slag_part = level_f.find_patch_part(-1);
          XASSERTM(slag_part != nullptr, "slag patch part not found!");
          MirrorType slag_mirror;
          Assembly::MirrorAssembler::assemble_mirror(slag_mirror, space_f, *slag_part);

          loc_vec_weight_p = VoxelAux::deslag_vector(loc_vec_weight_p, slag_mirror);
          loc_prol = VoxelAux::deslag_matrix_rows(loc_prol, slag_mirror);
          if(trunc)
          {
            loc_trunc = VoxelAux::deslag_matrix_cols(loc_trunc, slag_mirror);
          }
        }

        // synchronize weight vector using the gate
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
        // level, as in this case we have to join/split around the synchronize operation.

        // synchronize weight vector using the muxer/gate
        if(!virt_lvl_c.is_child())
        {
          XASSERT(&gate_f != &gate_c);

          // The coarse level is a simple (non-child) level that exists on all processes,
          // so simply synchronize over the coarse-level gate:
          gate_c.sync_0(loc_vec_weight_t);
        }
        else if(virt_lvl_c.is_parent())
        {
          XASSERT(&gate_f != &gate_c);

          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then synch that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a synch over the (non-existent) coarse-level
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
       * \brief Assembles a blocked transfer operator for a system level defined on a voxel domain
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
      static void asm_transfer_voxel_blocked(
        const Control::Domain::VirtualLevel<Control::Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_f,
        const Control::Domain::VirtualLevel<Control::Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_c,
        const String& cubature, bool trunc, bool shrink,
        SpaceLambda_&& space_lambda, Transfer_& transfer, const Muxer_& muxer, const Gate_& gate_f, const Gate_& gate_c)
      {
        typedef typename Transfer_::DataType DataType;
        typedef typename Transfer_::IndexType IndexType;
        typedef typename Gate_::MirrorType MirrorType;
        typedef LAFEM::DenseVector<DataType, IndexType> ScalarVectorType;

        // do we have a slag level?
        const bool have_slag(virt_lvl_f->slag_level);

        // we have: C -> S >> F, where '->' is the prolongation and '>>' is the deslagging operator
        const DomainLevel_& level_f = have_slag ? *virt_lvl_f->slag_level : static_cast<const DomainLevel_&>(*virt_lvl_f);
        const DomainLevel_& level_c = virt_lvl_c.is_child() ? virt_lvl_c.level_c() : *virt_lvl_c;

        //const auto& space_f = level_f.space;
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

        // if we need a truncation matrix, then compute its sparsity pattern now
        if(trunc)
          loc_trunc.transpose(loc_prol);

        // create local weight vector
        auto loc_vec_weight_p = loc_prol.create_vector_l();
        auto loc_vec_weight_t = loc_prol.create_vector_r();

        // assemble prolongation matrix
        loc_prol.format();
        loc_vec_weight_p.format();
        loc_vec_weight_t.format();

        // assemble prolongation matrix
        Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight_p, space_f, space_c, cubature);

        // assemble truncation matrix if desired
        if(trunc)
          Assembly::GridTransfer::assemble_truncation(loc_trunc, loc_vec_weight_t, space_f, space_c, cubature);

        // do we have a slag level?
        if(have_slag)
        {
          // assemble a slag mirror
          const auto* slag_part = level_f.find_patch_part(-1);
          XASSERTM(slag_part != nullptr, "slag patch part not found!");
          MirrorType slag_mirror;
          Assembly::MirrorAssembler::assemble_mirror(slag_mirror, space_f, *slag_part);

          loc_vec_weight_p = VoxelAux::deslag_vector(loc_vec_weight_p, slag_mirror);
          loc_prol = VoxelAux::deslag_matrix_rows(loc_prol, slag_mirror);
          if(trunc)
          {
            loc_trunc = VoxelAux::deslag_matrix_cols(loc_trunc, slag_mirror);
          }
        }

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
        // level, as in this case we have to join/split around the synchronize operation.

        // synchronize weight vector using the muxer/gate
        if(!virt_lvl_c.is_child())
        {
          XASSERT(&gate_f != &gate_c);

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
          XASSERT(&gate_f != &gate_c);

          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then synch that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a synch over the (non-existent) coarse-level
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
