// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_STOKES_BLOCKED_HPP
#define CONTROL_STOKES_BLOCKED_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/tuple_diag_matrix.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/domain_assembler_helpers.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/gpdv_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/muxer.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/mean_filter.hpp>
#include <kernel/global/transfer.hpp>

namespace FEAT
{
  namespace Control
  {
    template
    <
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
    >
    class StokesBlockedSystemLevel
    {
    public:
      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;
      static constexpr int dim = dim_;

      // scalar types
      typedef ScalarMatrix_ LocalScalarMatrix;
      typedef typename LocalScalarMatrix::VectorTypeL LocalScalarVector;

      // define local blocked matrix type
      typedef MatrixBlockA_ LocalMatrixBlockA;
      typedef MatrixBlockB_ LocalMatrixBlockB;
      typedef MatrixBlockD_ LocalMatrixBlockD;
      typedef LocalScalarMatrix LocalSchurMatrix;
      typedef LAFEM::SaddlePointMatrix<LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockD> LocalSystemMatrix;

      // define local vector types
      typedef typename LocalMatrixBlockB::VectorTypeL LocalVeloVector;
      typedef typename LocalMatrixBlockD::VectorTypeL LocalPresVector;
      typedef LAFEM::TupleVector<LocalVeloVector, LocalPresVector> LocalSystemVector;

      // define local transfer matrix types
      typedef TransferMatrixV_ LocalVeloTransferMatrix;
      typedef TransferMatrixP_ LocalPresTransferMatrix;
      typedef LAFEM::TupleDiagMatrix<TransferMatrixV_, TransferMatrixP_> LocalSystemTransferMatrix;

      // define local transfer operators
      typedef LAFEM::Transfer<LocalVeloTransferMatrix> LocalVeloTransfer;
      typedef LAFEM::Transfer<LocalPresTransferMatrix> LocalPresTransfer;
      typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

      // define mirror types
      typedef LAFEM::VectorMirror<MemType, DataType, IndexType> ScalarMirror;
      typedef ScalarMirror VeloMirror;
      typedef ScalarMirror PresMirror;
      typedef LAFEM::TupleMirror<VeloMirror, PresMirror> SystemMirror;

      // define gates
      typedef Global::Gate<LocalVeloVector, VeloMirror> VeloGate;
      typedef Global::Gate<LocalPresVector, PresMirror> PresGate;
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;

      // define muxers
      typedef Global::Muxer<LocalVeloVector, VeloMirror> VeloMuxer;
      typedef Global::Muxer<LocalPresVector, PresMirror> PresMuxer;
      typedef Global::Muxer<LocalSystemVector, SystemMirror> SystemMuxer;

      // define global vector types
      typedef Global::Vector<LocalVeloVector, VeloMirror> GlobalVeloVector;
      typedef Global::Vector<LocalPresVector, PresMirror> GlobalPresVector;
      typedef Global::Vector<LocalSystemVector, SystemMirror> GlobalSystemVector;

      // define global matrix types
      typedef Global::Matrix<LocalMatrixBlockA, VeloMirror, VeloMirror> GlobalMatrixBlockA;
      typedef Global::Matrix<LocalMatrixBlockB, VeloMirror, PresMirror> GlobalMatrixBlockB;
      typedef Global::Matrix<LocalMatrixBlockD, PresMirror, VeloMirror> GlobalMatrixBlockD;
      typedef Global::Matrix<LocalScalarMatrix, PresMirror, PresMirror> GlobalSchurMatrix;
      typedef Global::Matrix<LocalSystemMatrix, SystemMirror, SystemMirror> GlobalSystemMatrix;

      // define global transfer types
      typedef Global::Transfer<LocalVeloTransfer, VeloMirror> GlobalVeloTransfer;
      typedef Global::Transfer<LocalPresTransfer, PresMirror> GlobalPresTransfer;
      typedef Global::Transfer<LocalSystemTransfer, SystemMirror> GlobalSystemTransfer;

      /* ***************************************************************************************** */

      /// our system gate
      VeloGate gate_velo;
      PresGate gate_pres;
      SystemGate gate_sys;

      /// our coarse-level system muxer
      VeloMuxer coarse_muxer_velo;
      PresMuxer coarse_muxer_pres;
      SystemMuxer coarse_muxer_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;
      GlobalMatrixBlockA matrix_a;
      GlobalMatrixBlockB matrix_b;
      GlobalMatrixBlockD matrix_d;
      GlobalSchurMatrix matrix_s;

      /// our global transfer operator
      GlobalVeloTransfer transfer_velo;
      GlobalPresTransfer transfer_pres;
      GlobalSystemTransfer transfer_sys;

      /// CTOR
      StokesBlockedSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys),
        matrix_a(&gate_velo, &gate_velo),
        matrix_b(&gate_velo, &gate_pres),
        matrix_d(&gate_pres, &gate_velo),
        matrix_s(&gate_pres, &gate_pres),
        transfer_velo(&coarse_muxer_velo),
        transfer_pres(&coarse_muxer_pres),
        transfer_sys(&coarse_muxer_sys)
      {
      }

      virtual ~StokesBlockedSystemLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->matrix_sys.local().bytes () + this->matrix_s.local().bytes() + transfer_sys.bytes();
      }

      void compile_system_matrix()
      {
        matrix_sys.local().block_a() = matrix_a.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().block_b() = matrix_b.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().block_d() = matrix_d.local().clone(LAFEM::CloneMode::Shallow);
      }

      void compile_system_transfer()
      {
        transfer_sys.get_mat_prol().template at<0,0>() = transfer_velo.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<0,0>() = transfer_velo.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_prol().template at<1,1>() = transfer_pres.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<1,1>() = transfer_pres.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_trunc().template at<0,0>() = transfer_velo.get_mat_trunc().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_trunc().template at<1,1>() = transfer_pres.get_mat_trunc().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.compile();
      }

      template<typename M_, typename D_, typename I_, typename SMA_, typename SMB_, typename SMD_, typename SM_, typename TV_, typename TP_>
      void convert(const StokesBlockedSystemLevel<dim_, M_, D_, I_, SMA_, SMB_, SMD_, SM_, TV_, TP_> & other)
      {
        gate_velo.convert(other.gate_velo);
        gate_pres.convert(other.gate_pres);
        gate_sys.convert(other.gate_sys);

        coarse_muxer_velo.convert(other.coarse_muxer_velo);
        coarse_muxer_pres.convert(other.coarse_muxer_pres);
        coarse_muxer_sys.convert(other.coarse_muxer_sys);

        matrix_a.convert(&gate_velo, &gate_velo, other.matrix_a);
        matrix_b.convert(&gate_velo, &gate_pres, other.matrix_b);
        matrix_d.convert(&gate_pres, &gate_velo, other.matrix_d);
        matrix_s.convert(&gate_pres, &gate_pres, other.matrix_s);

        transfer_velo.convert(&coarse_muxer_velo, other.transfer_velo);
        transfer_pres.convert(&coarse_muxer_pres, other.transfer_pres);

        this->compile_system_matrix();
        this->compile_system_transfer();
      }

      template<typename DomainLevel_>
      void assemble_gates(const Domain::VirtualLevel<DomainLevel_>& virt_dom_lvl)
      {
        const auto& dom_level = virt_dom_lvl.level();
        const auto& dom_layer = virt_dom_lvl.layer();
        const auto& space_velo = dom_level.space_velo;
        const auto& space_pres = dom_level.space_pres;

        // set the gate comm
        this->gate_velo.set_comm(dom_layer.comm_ptr());
        this->gate_pres.set_comm(dom_layer.comm_ptr());
        this->gate_sys.set_comm(dom_layer.comm_ptr());

        // loop over all ranks
        for(Index i(0); i < dom_layer.neighbor_count(); ++i)
        {
          int rank = dom_layer.neighbor_rank(i);

          // try to find our halo
          auto* halo = dom_level.find_halo_part(rank);
          XASSERT(halo != nullptr);

          // create (empty) velocity mirror
          VeloMirror mirror_velo;
          Assembly::MirrorAssembler::assemble_mirror(mirror_velo, space_velo, *halo);

          // create (empty) pressure mirror
          PresMirror mirror_pres;
          Assembly::MirrorAssembler::assemble_mirror(mirror_pres, space_pres, *halo);

          // create a system mirror
          SystemMirror mirror_sys(mirror_velo.clone(), mirror_pres.clone());

          // push mirrors into gates
          if(!mirror_velo.empty())
            this->gate_velo.push(rank, std::move(mirror_velo));
          if(!mirror_pres.empty())
            this->gate_pres.push(rank, std::move(mirror_pres));
          if(!mirror_sys.empty())
            this->gate_sys.push(rank, std::move(mirror_sys));
        }

        // create local template vectors
        LocalVeloVector tmpl_v(space_velo.get_num_dofs());
        LocalPresVector tmpl_p(space_pres.get_num_dofs());
        LocalSystemVector tmpl_s(tmpl_v.clone(), tmpl_p.clone());

        // compile gates
        this->gate_velo.compile(std::move(tmpl_v));
        this->gate_pres.compile(std::move(tmpl_p));
        this->gate_sys.compile(std::move(tmpl_s));
      }

      template<typename DomainLevel_>
      void assemble_coarse_muxers(const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse)
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
            SystemMirror child_mirror_sys;
            VeloMirror& child_mirror_v = child_mirror_sys.template at<0>();
            PresMirror& child_mirror_p = child_mirror_sys.template at<1>();
            Assembly::MirrorAssembler::assemble_mirror(child_mirror_v, level_p.space_velo, *child);
            Assembly::MirrorAssembler::assemble_mirror(child_mirror_p, level_p.space_pres, *child);
            this->coarse_muxer_velo.push_child(child_mirror_v.clone(LAFEM::CloneMode::Shallow));
            this->coarse_muxer_pres.push_child(child_mirror_p.clone(LAFEM::CloneMode::Shallow));
            this->coarse_muxer_sys.push_child(std::move(child_mirror_sys));
          }
        }

        // assemble muxer child
        if(virt_lvl_coarse.is_child())
        {
          const auto& layer_c = virt_lvl_coarse.layer_c();
          const DomainLevel_& level_c = virt_lvl_coarse.level_c();

          SystemMirror parent_mirror_sys;
          VeloMirror& parent_mirror_v = parent_mirror_sys.template at<0>();
          PresMirror& parent_mirror_p = parent_mirror_sys.template at<1>();

          // manually set up an identity gather/scatter matrix
          {
            Index n = level_c.space_velo.get_num_dofs();
            parent_mirror_v = ScalarMirror(n, n);
            auto* idx = parent_mirror_v.indices();
            for(Index i(0); i < n; ++i)
              idx[i] = i;
          }
          {
            Index n = level_c.space_pres.get_num_dofs();
            parent_mirror_p = ScalarMirror(n, n);
            auto* idx = parent_mirror_p.indices();
            for(Index i(0); i < n; ++i)
              idx[i] = i;
          }

          // set parent and sibling comms
          this->coarse_muxer_velo.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            parent_mirror_v.clone(LAFEM::CloneMode::Shallow)
          );
          this->coarse_muxer_pres.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            parent_mirror_p.clone(LAFEM::CloneMode::Shallow)
          );
          this->coarse_muxer_sys.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            std::move(parent_mirror_sys)
          );

          // compile muxer
          LocalVeloVector tmpl_v(level_c.space_velo.get_num_dofs());
          LocalPresVector tmpl_p(level_c.space_pres.get_num_dofs());
          LocalSystemVector tmpl_s(tmpl_v.clone(), tmpl_p.clone());
          this->coarse_muxer_velo.compile(tmpl_v);
          this->coarse_muxer_pres.compile(tmpl_p);
          this->coarse_muxer_sys.compile(tmpl_s);
        }
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_velocity_transfer(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_velo;
        const auto& space_c = level_c.space_velo;

        // get local transfer operator
        LocalVeloTransfer& loc_trans = this->transfer_velo.local();

        // get local transfer matrices
        LocalVeloTransferMatrix& loc_prol_wrapped = loc_trans.get_mat_prol();
        LocalVeloTransferMatrix& loc_rest_wrapped = loc_trans.get_mat_rest();

        // get the unwrapped types
        typename LocalVeloTransferMatrix::BaseClass& loc_prol = loc_prol_wrapped;
        typename LocalVeloTransferMatrix::BaseClass& loc_rest = loc_rest_wrapped;

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_f, space_c);
        }

        // create a local weight vector
        LocalVeloVector loc_vec_weight = loc_prol_wrapped.create_vector_l();

        // create a scalar weight vector for the assembly
        LocalScalarVector loc_scal_vec_weight = loc_prol.create_vector_l();

        // get the data arrays of the weight vectors
        auto* v_wb = loc_vec_weight.elements();
        auto* v_ws = loc_scal_vec_weight.elements();

        // assemble prolongation matrix
        {
          loc_prol.format();
          loc_scal_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_scal_vec_weight,
            space_f, space_c, cubature);

          // copy weights from scalar to blocked
          for(Index i(0); i < loc_prol.rows(); ++i)
            v_wb[i] = v_ws[i];

          // synchronize blocked weight vector
          this->gate_velo.sync_0(loc_vec_weight);

          // copy weights from blocked to scalar
          for(Index i(0); i < loc_prol.rows(); ++i)
            v_ws[i] = v_wb[i][0];

          // invert weight components
          loc_scal_vec_weight.component_invert(loc_scal_vec_weight);

          // scale prolongation matrix
          loc_prol.scale_rows(loc_prol, loc_scal_vec_weight);

          // copy and transpose
          loc_rest = loc_prol.transpose();
        }

        // compile velocity transfer
        this->transfer_velo.compile();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_pressure_transfer(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_pres;
        const auto& space_c = level_c.space_pres;

        // get local transfer operator
        LocalPresTransfer& loc_trans = this->transfer_pres.local();

        // get local transfer matrices
        LocalPresTransferMatrix& loc_prol = loc_trans.get_mat_prol();
        LocalPresTransferMatrix& loc_rest = loc_trans.get_mat_rest();

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_f, space_c);
        }

        // get local pressure weight vector
        LocalPresVector loc_vec_weight = loc_prol.create_vector_l();

        // assemble prolongation matrix
        {
          loc_prol.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight,
            space_f, space_c, cubature);

          // synchronize weight vector
          this->gate_pres.sync_0(loc_vec_weight);

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol.scale_rows(loc_prol, loc_vec_weight);

          // copy and transpose
          loc_rest = loc_prol.transpose();
        }

        // compile pressure transfer
        this->transfer_pres.compile();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_transfers(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        this->assemble_velocity_transfer(virt_lvl_fine, virt_lvl_coarse, cubature);
        this->assemble_pressure_transfer(virt_lvl_fine, virt_lvl_coarse, cubature);

        this->compile_system_transfer();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_velocity_truncation(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const StokesBlockedSystemLevel* sys_lvl_coarse = nullptr)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT((sys_lvl_coarse != nullptr) || !virt_lvl_coarse.is_parent());

        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_velo;
        const auto& space_c = level_c.space_velo;

        // get local transfer operator
        LocalVeloTransfer& loc_trans = this->transfer_velo.local();

        // get local transfer matrices
        const LocalVeloTransferMatrix& loc_rest_wrapped = loc_trans.get_mat_rest();
        LocalVeloTransferMatrix& loc_trunc_wrapped = loc_trans.get_mat_trunc();

        // get the matrix blocks
        const typename LocalVeloTransferMatrix::BaseClass& loc_rest = loc_rest_wrapped;
        typename LocalVeloTransferMatrix::BaseClass& loc_trunc = loc_trunc_wrapped;

        // restriction matrix must be already assembled
        XASSERTM(loc_rest.size() > Index(0), "you need to call 'assemble_velocity_transfer' first");

        // clone sparsity pattern of restriction matrix
        loc_trunc = loc_rest.clone(LAFEM::CloneMode::Layout);

        // create local weight vector
        LocalVeloVector loc_vec_weight = loc_rest_wrapped.create_vector_l();

        // create a scalar weight vector for the assembly
        LocalScalarVector loc_scal_vec_weight = loc_rest.create_vector_l();

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
          sys_lvl_coarse->gate_velo.sync_0(loc_vec_weight);
        }
        else if(virt_lvl_coarse.is_parent())
        {
          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then synch that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a synch over the (non-existent) coarse-level
          // child gate, which is what we actually require here...

          // create temporary vector on parent partitioning
          LocalVeloVector loc_tmp = sys_lvl_coarse->gate_velo._freqs.clone(LAFEM::CloneMode::Allocate);

          // join child weights over muxer
          this->coarse_muxer_velo.join(loc_vec_weight, loc_tmp);

          // sync over coarse gate
          sys_lvl_coarse->gate_velo.sync_0(loc_tmp);

          // split over muxer
          this->coarse_muxer_velo.split(loc_vec_weight, loc_tmp);
        }
        else // ghost
        {
          // The coarse level is a ghost level, i.e. a child but not a parent. In this case, we
          // only have to participate in the join/send operations of the muxer, which are part
          // of the operations executed on the parents handled by the else-if case above.

          this->coarse_muxer_velo.join_send(loc_vec_weight);

          // parent performs sync over its gate here (see above else-if)

          this->coarse_muxer_velo.split_recv(loc_vec_weight);
        }

        // copy weights from blocked to scalar
        for(Index i(0); i < loc_trunc.rows(); ++i)
          v_ws[i] = v_wb[i][0];

        // invert components
        loc_scal_vec_weight.component_invert(loc_scal_vec_weight);

        // scale reduction matrix
        loc_trunc.scale_rows(loc_trunc, loc_scal_vec_weight);
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_pressure_truncation(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const StokesBlockedSystemLevel* sys_lvl_coarse = nullptr)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT((sys_lvl_coarse != nullptr) || !virt_lvl_coarse.is_parent());

        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_pres;
        const auto& space_c = level_c.space_pres;

        // get local transfer operator
        LocalPresTransfer& loc_trans = this->transfer_pres.local();

        // get local transfer matrices
        const LocalPresTransferMatrix& loc_rest = loc_trans.get_mat_rest();
        LocalPresTransferMatrix& loc_trunc = loc_trans.get_mat_trunc();

        // restriction matrix must be already assembled
        XASSERTM(loc_rest.size() > Index(0), "you need to call 'assemble_prescity_transfer' first");

        // clone sparsity pattern of restriction matrix
        loc_trunc = loc_rest.clone(LAFEM::CloneMode::Layout);

        // create local weight vector
        LocalPresVector loc_vec_weight = loc_trunc.create_vector_l();

        // format
        loc_trunc.format();
        loc_vec_weight.format();

        // assemble truncation matrix
        Assembly::GridTransfer::assemble_truncation(loc_trunc, loc_vec_weight, space_f, space_c, cubature);

        // We now need to synchronize the weight vector in analogy to the prolongation matrix assembly.
        // Note that the weight vector is now a coarse-level vector, so we need to synchronize over
        // the coarse-level gate. This may be a bit more complicated if the coarse level is a ghost
        // level, as in this case we have to join/split around the synch operation.

        // synchronize weight vector using the muxer/gate
        if(!virt_lvl_coarse.is_child())
        {
          // The coarse level is a simple (non-child) level that exists on all processes,
          // so simply synch over the coarse-level gate:
          sys_lvl_coarse->gate_pres.sync_0(loc_vec_weight);
        }
        else if(virt_lvl_coarse.is_parent())
        {
          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then synch that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a synch over the (non-existent) coarse-level
          // child gate, which is what we actually require here...

          // create temporary vector on parent partitioning
          LocalPresVector loc_tmp = sys_lvl_coarse->gate_pres._freqs.clone(LAFEM::CloneMode::Allocate);

          // join child weights over muxer
          this->coarse_muxer_pres.join(loc_vec_weight, loc_tmp);

          // sync over coarse gate
          sys_lvl_coarse->gate_pres.sync_0(loc_tmp);

          // split over muxer
          this->coarse_muxer_pres.split(loc_vec_weight, loc_tmp);
        }
        else // ghost
        {
          // The coarse level is a ghost level, i.e. a child but not a parent. In this case, we
          // only have to participate in the join/send operations of the muxer, which are part
          // of the operations executed on the parents handled by the else-if case above.

          this->coarse_muxer_pres.join_send(loc_vec_weight);

          // parent performs sync over its gate here (see above else-if)

          this->coarse_muxer_pres.split_recv(loc_vec_weight);
        }

        // invert components
        loc_vec_weight.component_invert(loc_vec_weight);

        // scale reduction matrix
        loc_trunc.scale_rows(loc_trunc, loc_vec_weight);
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_truncations(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const StokesBlockedSystemLevel* sys_lvl_coarse = nullptr)
      {
        this->assemble_velocity_truncation(virt_lvl_fine, virt_lvl_coarse, cubature, sys_lvl_coarse);
        this->assemble_pressure_truncation(virt_lvl_fine, virt_lvl_coarse, cubature, sys_lvl_coarse);
        this->compile_system_transfer();
      }

      template<typename SpaceVelo_, typename SpacePres_, typename Cubature_>
      void assemble_grad_div_matrices(const SpaceVelo_& space_velo, const SpacePres_& space_pres, const Cubature_& cubature)
      {
        Assembly::GradPresDivVeloAssembler::assemble(this->matrix_b.local(), this->matrix_d.local(),
        space_velo, space_pres, cubature);
      }

      template<typename Trafo_, typename SpaceVelo_, typename SpacePres_>
      void assemble_grad_div_matrices(Assembly::DomainAssembler<Trafo_>& dom_asm,
        const SpaceVelo_& space_velo, const SpacePres_& space_pres, const String& cubature)
      {
        // assemble matrix structure of B
        if(this->matrix_b.local().empty())
          Assembly::SymbolicAssembler::assemble_matrix_std2(this->matrix_b.local(), space_velo, space_pres);

        // assemble matrix B
        Assembly::Common::GradientTestOperatorBlocked<dim_> grad_op;
        Assembly::assemble_bilinear_operator_matrix_2(
          dom_asm, this->matrix_b.local(), grad_op, space_velo, space_pres, cubature, -DataType(1));

        // transpose to obtain matrix D
        this->matrix_d.local().transpose(this->matrix_b.local());
      }

      template<typename SpaceVelo_>
      void assemble_velo_struct(const SpaceVelo_& space_velo)
      {
        // assemble matrix structure
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_a.local(), space_velo);
      }

      template<typename SpacePres_>
      void assemble_pres_struct(const SpacePres_& space_pres)
      {
        // assemble matrix structure
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_s.local(), space_pres);
      }
    }; // class StokesBlockedSystemLevel<...>

    template
    <
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
    >
    class StokesBlockedUnitVeloNonePresSystemLevel :
      public StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
      MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalVeloFilter;
      typedef LAFEM::NoneFilter<MemType_, DataType_, IndexType_> LocalPresFilter;
      typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresFilter> LocalSystemFilter;

      // define global filter types
      typedef Global::Filter<LocalVeloFilter, typename BaseClass::VeloMirror> GlobalVeloFilter;
      typedef Global::Filter<LocalPresFilter, typename BaseClass::PresMirror> GlobalPresFilter;
      typedef Global::Filter<LocalSystemFilter, typename BaseClass::SystemMirror> GlobalSystemFilter;

      // (global) filters
      GlobalSystemFilter filter_sys;
      GlobalVeloFilter filter_velo;
      GlobalPresFilter filter_pres;

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->filter_sys.bytes() + BaseClass::bytes();
      }

      void compile_system_filter()
      {
        filter_sys.local().template at<0>() = filter_velo.local().clone(LAFEM::CloneMode::Shallow);
        filter_sys.local().template at<1>() = filter_pres.local().clone(LAFEM::CloneMode::Shallow);
      }
    };

    template
    <
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
    >
    class StokesBlockedSlipUnitVeloNonePresSystemLevel :
      public StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
      MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::SlipFilter<MemType_, DataType_, IndexType_, dim_> LocalVeloSlipFilter;
      typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalVeloUnitFilter;
      typedef LAFEM::FilterChain<LocalVeloSlipFilter, LocalVeloUnitFilter> LocalVeloFilter;
      typedef LAFEM::NoneFilter<MemType_, DataType_, IndexType_> LocalPresFilter;
      typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresFilter> LocalSystemFilter;

      // define global filter types
      typedef Global::Filter<LocalVeloFilter, typename BaseClass::VeloMirror> GlobalVeloFilter;
      typedef Global::Filter<LocalPresFilter, typename BaseClass::PresMirror> GlobalPresFilter;
      typedef Global::Filter<LocalSystemFilter, typename BaseClass::SystemMirror> GlobalSystemFilter;

      // (global) filters
      GlobalSystemFilter filter_sys;
      GlobalVeloFilter filter_velo;
      GlobalPresFilter filter_pres;

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->filter_sys.bytes() + BaseClass::bytes();
      }

      void compile_system_filter()
      {
        filter_sys.local().template at<0>() = filter_velo.local().clone(LAFEM::CloneMode::Shallow);
        filter_sys.local().template at<1>() = filter_pres.local().clone(LAFEM::CloneMode::Shallow);
      }
    };

    /**
     * \brief System level using a MeanFilter for the pressure
     *
     * This is necessary when there are only Dirichlet BCs for the velocity
     */
    template
    <
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
    >
    struct StokesBlockedUnitVeloMeanPresSystemLevel :
      public StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
      typedef StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
      MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalVeloFilter;
      typedef Global::MeanFilter<MemType_, DataType_, IndexType_> LocalPresFilter;
      typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresFilter> LocalSystemFilter;

      // define global filter types
      typedef Global::Filter<LocalVeloFilter, typename BaseClass::VeloMirror> GlobalVeloFilter;
      typedef Global::Filter<LocalPresFilter, typename BaseClass::PresMirror> GlobalPresFilter;
      typedef Global::Filter<LocalSystemFilter, typename BaseClass::SystemMirror> GlobalSystemFilter;

      // (global) filters
      GlobalSystemFilter filter_sys;
      GlobalVeloFilter filter_velo;
      GlobalPresFilter filter_pres;

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->matrix_sys.local().bytes () + this->matrix_s.local().bytes() + filter_sys.local().bytes();
      }

      void compile_system_filter()
      {
        filter_sys.local().template at<0>() = filter_velo.local().clone(LAFEM::CloneMode::Shallow);
        filter_sys.local().template at<1>() = filter_pres.local().clone(LAFEM::CloneMode::Shallow);
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source StokesUnitVeloMeanPresSystemLevel content as content of current StokesUnitVeloMeanPresSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const StokesBlockedUnitVeloMeanPresSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_velo.convert(other.filter_velo);
        filter_pres.convert(other.filter_pres);

        compile_system_filter();
      }

      template<typename SpacePres_, typename Cubature_>
      void assemble_pressure_mean_filter(const SpacePres_& space_pres, const Cubature_& cubature)
      {
        // get our local pressure filter
        LocalPresFilter& fil_loc_p = this->filter_pres.local();

        // create two global vectors
        typename BaseClass::GlobalPresVector vec_glob_v(&this->gate_pres), vec_glob_w(&this->gate_pres);

        // fetch the local vectors
        typename BaseClass::LocalPresVector& vec_loc_v = vec_glob_v.local();
        typename BaseClass::LocalPresVector& vec_loc_w = vec_glob_w.local();

        // fetch the frequency vector of the pressure gate
        typename BaseClass::LocalPresVector& vec_loc_f = this->gate_pres._freqs;

        // assemble the mean filter
        Assembly::MeanFilterAssembler::assemble(vec_loc_v, vec_loc_w, space_pres, cubature);

        // synchronize the vectors
        vec_glob_v.sync_1();
        vec_glob_w.sync_0();

        // build the mean filter
        fil_loc_p = LocalPresFilter(vec_loc_v.clone(), vec_loc_w.clone(), vec_loc_f.clone(), this->gate_pres.get_comm());
      }
    }; // struct StokesBlockedUnitVeloMeanPresSystemLevel<...>

    template
    <
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
    >
    class StokesBlockedSlipUnitVeloMeanPresSystemLevel :
      public StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
      MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::SlipFilter<MemType_, DataType_, IndexType_, dim_> LocalVeloSlipFilter;
      typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalVeloUnitFilter;
      typedef LAFEM::FilterChain<LocalVeloSlipFilter, LocalVeloUnitFilter> LocalVeloFilter;
      typedef Global::MeanFilter<MemType_, DataType_, IndexType_> LocalPresFilter;
      typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresFilter> LocalSystemFilter;

      // define global filter types
      typedef Global::Filter<LocalVeloFilter, typename BaseClass::VeloMirror> GlobalVeloFilter;
      typedef Global::Filter<LocalPresFilter, typename BaseClass::PresMirror> GlobalPresFilter;
      typedef Global::Filter<LocalSystemFilter, typename BaseClass::SystemMirror> GlobalSystemFilter;

      // (global) filters
      GlobalSystemFilter filter_sys;
      GlobalVeloFilter filter_velo;
      GlobalPresFilter filter_pres;

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->filter_sys.bytes() + BaseClass::bytes();
      }

      void compile_system_filter()
      {
        filter_sys.local().template at<0>() = filter_velo.local().clone(LAFEM::CloneMode::Shallow);
        filter_sys.local().template at<1>() = filter_pres.local().clone(LAFEM::CloneMode::Shallow);
      }

      /**
        *
        * \brief Conversion method
        *
        * Use source StokesUnitVeloMeanPresSystemLevel content as content of current StokesUnitVeloMeanPresSystemLevel.
        *
        */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const StokesBlockedUnitVeloMeanPresSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_velo.convert(other.filter_velo);
        filter_pres.convert(other.filter_pres);

        compile_system_filter();
      }

      template<typename SpacePres_, typename Cubature_>
      void assemble_global_filters(const SpacePres_& space_pres, const Cubature_& cubature)
      {
        // get our local pressure filter
        LocalPresFilter& fil_loc_p = this->filter_pres.local();

        // create two global vectors
        typename BaseClass::GlobalPresVector vec_glob_v(&this->gate_pres), vec_glob_w(&this->gate_pres);

        // fetch the local vectors
        typename BaseClass::LocalPresVector& vec_loc_v = vec_glob_v.local();
        typename BaseClass::LocalPresVector& vec_loc_w = vec_glob_w.local();

        // fetch the frequency vector of the pressure gate
        typename BaseClass::LocalPresVector& vec_loc_f = this->gate_pres._freqs;

        // assemble the mean filter
        Assembly::MeanFilterAssembler::assemble(vec_loc_v, vec_loc_w, space_pres, cubature);

        // synchronize the vectors
        vec_glob_v.sync_1();
        vec_glob_w.sync_0();

        // build the mean filter
        fil_loc_p = LocalPresFilter(vec_loc_v.clone(), vec_loc_w.clone(), vec_loc_f.clone(), this->gate_pres.get_comm());

        // Sync the filter vector in the SlipFilter
        const typename BaseClass::VeloGate& my_col_gate(this->gate_velo);

        // For all slip filters...
        //for(auto& it : filter_sys.local().template at<0>())
        //{

        auto& it = filter_velo.local().template at<0>();
        // get the filter vector
        auto& slip_filter_vector = it.get_filter_vector();

        if(slip_filter_vector.used_elements() > 0)
        {
          // Temporary DenseVector for syncing
          typename BaseClass::LocalVeloVector tmp(slip_filter_vector.size(), DataType_(0));

          auto* tmp_elements = tmp.template elements<LAFEM::Perspective::native>();
          auto* sfv_elements = slip_filter_vector.template elements<LAFEM::Perspective::native>();

          // Copy sparse filter vector contents to DenseVector
          for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
          {
            Index idense(slip_filter_vector.indices()[isparse]);
            tmp_elements[idense] = sfv_elements[isparse];
          }

          my_col_gate.sync_0(tmp);
          // Copy sparse filter vector contents to DenseVector
          for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
          {
            Index idense(slip_filter_vector.indices()[isparse]);
            tmp_elements[idense].normalize();
            sfv_elements[isparse] = tmp_elements[idense];

          }
        }
        else
        {
          // Temporary DenseVector for syncing
          typename BaseClass::LocalVeloVector tmp(slip_filter_vector.size(), DataType_(0));
          my_col_gate.sync_0(tmp);
        }
      }
    };
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STOKES_BLOCKED_HPP
