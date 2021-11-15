// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_BLOCKED_BASIC_HPP
#define CONTROL_BLOCKED_BASIC_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/global/muxer.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>

namespace FEAT
{
  namespace Control
  {
    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename BlockedMatrix_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename TransferMatrix_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>
      >
    class BlockedBasicSystemLevel
    {
    public:
      static_assert(std::is_same<MemType_, typename BlockedMatrix_::MemType>::value, "MemType mismatch!");
      static_assert(std::is_same<DataType_, typename BlockedMatrix_::DataType>::value, "DataType mismatch!");
      static_assert(std::is_same<IndexType_, typename BlockedMatrix_::IndexType>::value, "IndexType mismatch!");

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;
      static constexpr int dim = dim_;

      /// define local blocked matrix type
      typedef BlockedMatrix_ LocalBlockedMatrix;

      /// define local system matrix type
      typedef BlockedMatrix_ LocalSystemMatrix;

      /// define local transfer matrix type
      typedef TransferMatrix_ LocalSystemTransferMatrix;

      /// define local system vector type
      typedef typename LocalSystemMatrix::VectorTypeR LocalSystemVector;

      /// define local system transfer operator type
      typedef typename LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define system gate
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;

      // define system muxer
      typedef Global::Muxer<LocalSystemVector, SystemMirror> SystemMuxer;

      /// define global system vector type
      typedef Global::Vector<LocalSystemVector, SystemMirror> GlobalSystemVector;

      /// define global system matrix type
      typedef Global::Matrix<LocalSystemMatrix, SystemMirror, SystemMirror> GlobalSystemMatrix;

      /// define global system operator type
      typedef Global::Transfer<LocalSystemTransfer, SystemMirror> GlobalSystemTransfer;

      /* ***************************************************************************************** */

      /// our system gate
      SystemGate gate_sys;

      /// our coarse-level system muxer
      SystemMuxer coarse_muxer_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;

      /// our global transfer operator
      GlobalSystemTransfer transfer_sys;

      /// CTOR
      BlockedBasicSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys),
        transfer_sys(&coarse_muxer_sys)
      {
      }

      virtual ~BlockedBasicSystemLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->gate_sys.bytes() + this->coarse_muxer_sys.bytes()
          + this->transfer_sys.bytes() + this->matrix_sys.local().bytes();
      }

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const BlockedBasicSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        gate_sys.convert(other.gate_sys);
        coarse_muxer_sys.convert(other.coarse_muxer_sys);
        matrix_sys.convert(&gate_sys, &gate_sys, other.matrix_sys);
        transfer_sys.convert(&coarse_muxer_sys, other.transfer_sys);
      }

      template<typename DomainLevel_>
      void assemble_gate(const Domain::VirtualLevel<DomainLevel_>& virt_dom_lvl)
      {
        const auto& dom_level = virt_dom_lvl.level();
        const auto& dom_layer = virt_dom_lvl.layer();
        const auto& space = dom_level.space;

        // set the gate comm
        this->gate_sys.set_comm(dom_layer.comm_ptr());

        // loop over all ranks
        for(Index i(0); i < dom_layer.neighbor_count(); ++i)
        {
          int rank = dom_layer.neighbor_rank(i);

          // try to find our halo
          auto* halo = dom_level.find_halo_part(rank);
          XASSERT(halo != nullptr);

          // assemble the mirror
          SystemMirror mirror_sys;
          Assembly::MirrorAssembler::assemble_mirror(mirror_sys, space, *halo);

          // push mirror into gate
          this->gate_sys.push(rank, std::move(mirror_sys));
        }

        // create local template vector
        LocalSystemVector tmpl_s(space.get_num_dofs());

        // compile gate
        this->gate_sys.compile(std::move(tmpl_s));
      }

      template<typename DomainLevel_>
      void assemble_coarse_muxer(const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse)
      {
        // assemble muxer parent
        if(virt_lvl_coarse.is_parent())
        {
          XASSERT(virt_lvl_coarse.is_child());

          const auto& layer_c = virt_lvl_coarse.layer_c();
          const DomainLevel_& level_p = virt_lvl_coarse.level_p();

          // loop over all children
          for(Index i(0); i < layer_c.child_count(); ++i)
          {
            const auto* child = level_p.find_patch_part(int(i));
            XASSERT(child != nullptr);
            SystemMirror child_mirror;
            Assembly::MirrorAssembler::assemble_mirror(child_mirror, level_p.space, *child);
            this->coarse_muxer_sys.push_child(std::move(child_mirror));
          }
        }

        // assemble muxer child
        if(virt_lvl_coarse.is_child())
        {
          const auto& layer_c = virt_lvl_coarse.layer_c();
          const DomainLevel_& level_c = virt_lvl_coarse.level_c();

          // manually set up an identity gather/scatter matrix
          Index n = level_c.space.get_num_dofs();
          SystemMirror parent_mirror(n, n);
          auto* idx = parent_mirror.indices();
          for(Index i(0); i < n; ++i)
            idx[i] = i;

          // set parent and sibling comms
          this->coarse_muxer_sys.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            std::move(parent_mirror)
          );

          // compile muxer
          LocalSystemVector vec_tmp(level_c.space.get_num_dofs());
          this->coarse_muxer_sys.compile(vec_tmp);
        }
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_transfer(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space;
        const auto& space_c = level_c.space;

        // get local transfer operator
        LocalSystemTransfer& loc_trans = this->transfer_sys.local();

        // get local transfer matrices
        LocalSystemTransferMatrix& loc_prol_wrapped = loc_trans.get_mat_prol();
        LocalSystemTransferMatrix& loc_rest_wrapped = loc_trans.get_mat_rest();

        // get the unwrapped types
        auto& loc_prol = loc_prol_wrapped.unwrap();
        auto& loc_rest = loc_rest_wrapped.unwrap();

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_f, space_c);
        }

        // create local blocked weight vector
        LocalSystemVector loc_vec_weight = loc_prol_wrapped.create_vector_l();

        // create local weight vector
        auto loc_scal_vec_weight = loc_prol.create_vector_l();

        // get the data arrays of the weight vectors
        auto* v_wb = loc_vec_weight.elements();
        auto* v_ws = loc_scal_vec_weight.elements();

        // assemble prolongation matrix
        loc_prol.format();
        loc_scal_vec_weight.format();

        // assemble prolongation matrix
        Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_scal_vec_weight, space_f, space_c, cubature);

        // copy weights from scalar to blocked
        for(Index i(0); i < loc_prol.rows(); ++i)
          v_wb[i] = v_ws[i];

        // synchronize blocked weight vector
        this->gate_sys.sync_0(loc_vec_weight);

        // copy weights from blocked to scalar
        for(Index i(0); i < loc_prol.rows(); ++i)
          v_ws[i] = v_wb[i][0];

        // invert components
        loc_scal_vec_weight.component_invert(loc_scal_vec_weight);

        // scale prolongation matrix
        loc_prol.scale_rows(loc_prol, loc_scal_vec_weight);

        // copy and transpose
        loc_rest = loc_prol.transpose();

        // compile global transfer
        transfer_sys.compile();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_truncation(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const BlockedBasicSystemLevel* sys_lvl_coarse = nullptr)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT((sys_lvl_coarse != nullptr) || !virt_lvl_coarse.is_parent());

        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space;
        const auto& space_c = level_c.space;

        // get local transfer operator
        LocalSystemTransfer& loc_trans = this->transfer_sys.local();

        // get local transfer matrices
        const LocalSystemTransferMatrix& loc_rest_wrapped = loc_trans.get_mat_rest();
        LocalSystemTransferMatrix& loc_trunc_wrapped = loc_trans.get_mat_trunc();

        const auto& loc_rest = loc_rest_wrapped.unwrap();
        auto& loc_trunc = loc_trunc_wrapped.unwrap();

        // restriction matrix must be already assembled
        XASSERTM(loc_rest.size() > Index(0), "you need to call 'assemble_transfer' first");

        // clone sparsity pattern of restriction matrix
        loc_trunc = loc_rest.clone(LAFEM::CloneMode::Layout);

        // create local weight vector
        LocalSystemVector loc_vec_weight = loc_trunc_wrapped.create_vector_l();

        // create a scalar weight vector for the assembly
        auto loc_scal_vec_weight = loc_trunc.create_vector_l();

        // get the data arrays of the weight vectors
        auto* v_wb = loc_vec_weight.elements();
        auto* v_ws = loc_scal_vec_weight.elements();

        // format
        loc_trunc.format();
        loc_scal_vec_weight.format();

        // assemble truncation matrix
        Assembly::GridTransfer::assemble_truncation(loc_trunc, loc_scal_vec_weight, space_f, space_c, cubature);

        // copy weights from scalar to blocked
        for(Index i(0); i < loc_trunc.rows(); ++i)
          v_wb[i] = v_ws[i];

        // We now need to synchronize the weight vector in analogy to the prolongation matrix assembly.
        // Note that the weight vector is now a coarse-level vector, so we need to synchronize over
        // the coarse-level gate. This may be a bit more complicated if the coarse level is a ghost
        // level, as in this case we have to join/split around the synch operation.

        // synchronize weight vector using the muxer/gate
        if(!virt_lvl_coarse.is_child())
        {
          // The coarse level is a simple (non-child) level that exists on all processes,
          // so simply synch over the coarse-level gate:
          sys_lvl_coarse->gate_sys.sync_0(loc_vec_weight);
        }
        else if(virt_lvl_coarse.is_parent())
        {
          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then synch that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a synch over the (non-existent) coarse-level
          // child gate, which is what we actually require here...

          // create temporary vector on parent partitioning
          LocalSystemVector loc_tmp = sys_lvl_coarse->gate_sys._freqs.clone(LAFEM::CloneMode::Allocate);

          // join child weights over muxer
          this->coarse_muxer_sys.join(loc_vec_weight, loc_tmp);

          // sync over coarse gate
          sys_lvl_coarse->gate_sys.sync_0(loc_tmp);

          // split over muxer
          this->coarse_muxer_sys.split(loc_vec_weight, loc_tmp);
        }
        else // ghost
        {
          // The coarse level is a ghost level, i.e. a child but not a parent. In this case, we
          // only have to participate in the join/send operations of the muxer, which are part
          // of the operations executed on the parents handled by the else-if case above.

          this->coarse_muxer_sys.join_send(loc_vec_weight);

          // parent performs sync over its gate here (see above else-if)

          this->coarse_muxer_sys.split_recv(loc_vec_weight);
        }

        // copy weights from blocked to scalar
        for(Index i(0); i < loc_trunc.rows(); ++i)
          v_ws[i] = v_wb[i][0];

        // invert components
        loc_scal_vec_weight.component_invert(loc_scal_vec_weight);

        // scale reduction matrix
        loc_trunc.scale_rows(loc_trunc, loc_scal_vec_weight);
      }
    }; // class BlockedBasicSystemLevel<...>


    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename BlockedMatrix_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename TransferMatrix_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>>
    class BlockedUnitFilterSystemLevel :
      public BlockedBasicSystemLevel<dim_, MemType_, DataType_, IndexType_, BlockedMatrix_, TransferMatrix_>
    {
    public:
      typedef BlockedBasicSystemLevel<dim_, MemType_, DataType_, IndexType_, BlockedMatrix_, TransferMatrix_> BaseClass;

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define local filter type
      typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename Mem2_, typename DT2_, typename IT2_, typename BlockedMatrix2_>
      using BaseType = BlockedUnitFilterSystemLevel<dim_, Mem2_, DT2_, IT2_, BlockedMatrix2_>;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      BlockedUnitFilterSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return BaseClass::bytes () + filter_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source BlockedUnitFilterSystemLevel content as content of current BlockedUnitFilterSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const BlockedUnitFilterSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }


      template<typename DomainLevel_, typename Space_>
      void assemble_homogeneous_unit_filter(const DomainLevel_& dom_level, const Space_& space)
      {
        auto& loc_filter = this->filter_sys.local();

        // create unit-filter assembler
        Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> unit_asm;

        std::deque<String> part_names = dom_level.get_mesh_node()->get_mesh_part_names(true);
        for(const auto& name : part_names)
        {
          auto* mesh_part_node = dom_level.get_mesh_node()->find_mesh_part_node(name);
          XASSERT(mesh_part_node != nullptr);

          // let's see if we have that mesh part
          // if it is nullptr, then our patch is not adjacent to that boundary part
          auto* mesh_part = mesh_part_node->get_mesh();
          if (mesh_part != nullptr)
          {
            // add to boundary assembler
            unit_asm.add_mesh_part(*mesh_part);
          }
        }

        // finally, assemble the filter
        unit_asm.assemble(loc_filter, space);
      }
    }; // class BlockedUnitFilterSystemLevel<...>
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_BLOCKED_BASIC_HPP
