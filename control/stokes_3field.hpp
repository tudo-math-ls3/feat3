#pragma once
#ifndef CONTROL_STOKES_3FIELD_HPP
#define CONTROL_STOKES_3FIELD_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>
#include <kernel/lafem/null_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
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
    /**
     *
     *  /  A  B  R  \   / v \
     *  |  D  .  .  | * | p |
     *  \  K  .  M  /   \ s /
     *
     */
    template
    <
      int dim_,
      int nsc_ = (dim_*(dim_+1))/2,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
      typename MatrixBlockM_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, nsc_, nsc_>,
      typename MatrixBlockK_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, nsc_, dim_>,
      typename MatrixBlockL_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, nsc_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixS_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, nsc_>
    >
    class Stokes3FieldSystemLevel
    {
    public:
      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;
      static constexpr int dim = dim_;
      static constexpr int nsc = nsc_;

      // scalar types
      typedef ScalarMatrix_ LocalScalarMatrix;
      typedef typename LocalScalarMatrix::VectorTypeL LocalScalarVector;

      // define local blocked matrix types
      typedef MatrixBlockA_ LocalMatrixBlockA;
      typedef MatrixBlockB_ LocalMatrixBlockB;
      typedef MatrixBlockD_ LocalMatrixBlockD;
      typedef MatrixBlockM_ LocalMatrixBlockM;
      typedef MatrixBlockK_ LocalMatrixBlockK;
      typedef MatrixBlockL_ LocalMatrixBlockL;
      typedef LAFEM::NullMatrix<MemType_, DataType_, IndexType_,   1,   1> NullMatrixBlockPP;
      typedef LAFEM::NullMatrix<MemType_, DataType_, IndexType_, nsc,   1> NullMatrixBlockSP;
      typedef LAFEM::NullMatrix<MemType_, DataType_, IndexType_,   1, nsc> NullMatrixBlockPS;
      typedef LAFEM::TupleMatrix<
        LAFEM::TupleMatrixRow<LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockL>,
        LAFEM::TupleMatrixRow<LocalMatrixBlockD, NullMatrixBlockPP, NullMatrixBlockPS>,
        LAFEM::TupleMatrixRow<LocalMatrixBlockK, NullMatrixBlockSP, LocalMatrixBlockM>
      >  LocalSystemMatrix;

      // define local vector types
      typedef typename LocalMatrixBlockB::VectorTypeL LocalVeloVector;
      typedef typename LocalMatrixBlockD::VectorTypeL LocalPresVector;
      typedef typename LocalMatrixBlockM::VectorTypeL LocalStressVector;
      typedef LAFEM::TupleVector<LocalVeloVector, LocalPresVector, LocalStressVector> LocalSystemVector;

      // define local transfer matrix types
      typedef TransferMatrixV_ LocalVeloTransferMatrix;
      typedef TransferMatrixP_ LocalPresTransferMatrix;
      typedef TransferMatrixS_ LocalStressTransferMatrix;
      typedef LAFEM::TupleDiagMatrix<TransferMatrixV_, TransferMatrixP_, TransferMatrixS_> LocalSystemTransferMatrix;

      // define local transfer operators
      typedef LAFEM::Transfer<LocalVeloTransferMatrix> LocalVeloTransfer;
      typedef LAFEM::Transfer<LocalPresTransferMatrix> LocalPresTransfer;
      typedef LAFEM::Transfer<LocalStressTransferMatrix> LocalStressTransfer;
      typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

      // define mirror types
      typedef LAFEM::VectorMirror<MemType, DataType, IndexType> ScalarMirror;
      typedef ScalarMirror VeloMirror;
      typedef ScalarMirror PresMirror;
      typedef ScalarMirror StressMirror;
      typedef LAFEM::TupleMirror<VeloMirror, PresMirror, StressMirror> SystemMirror;

      // define gates
      typedef Global::Gate<LocalVeloVector, VeloMirror> VeloGate;
      typedef Global::Gate<LocalPresVector, PresMirror> PresGate;
      typedef Global::Gate<LocalStressVector, StressMirror> StressGate;
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;

      // define muxers
      typedef Global::Muxer<LocalVeloVector, VeloMirror> VeloMuxer;
      typedef Global::Muxer<LocalPresVector, PresMirror> PresMuxer;
      typedef Global::Muxer<LocalStressVector, StressMirror> StressMuxer;
      typedef Global::Muxer<LocalSystemVector, SystemMirror> SystemMuxer;

      // define global vector types
      typedef Global::Vector<LocalVeloVector, VeloMirror> GlobalVeloVector;
      typedef Global::Vector<LocalPresVector, PresMirror> GlobalPresVector;
      typedef Global::Vector<LocalStressVector, StressMirror> GlobalStressVector;
      typedef Global::Vector<LocalSystemVector, SystemMirror> GlobalSystemVector;

      // define global matrix types
      typedef Global::Matrix<LocalMatrixBlockA, VeloMirror, VeloMirror> GlobalMatrixBlockA;
      typedef Global::Matrix<LocalMatrixBlockB, VeloMirror, PresMirror> GlobalMatrixBlockB;
      typedef Global::Matrix<LocalMatrixBlockD, PresMirror, VeloMirror> GlobalMatrixBlockD;
      typedef Global::Matrix<LocalMatrixBlockK, StressMirror, VeloMirror> GlobalMatrixBlockK;
      typedef Global::Matrix<LocalMatrixBlockL, VeloMirror, StressMirror> GlobalMatrixBlockL;
      typedef Global::Matrix<LocalMatrixBlockM, StressMirror, StressMirror> GlobalMatrixBlockM;
      typedef Global::Matrix<LocalSystemMatrix, SystemMirror, SystemMirror> GlobalSystemMatrix;

      // define global transfer types
      typedef Global::Transfer<LocalVeloTransfer, VeloMirror> GlobalVeloTransfer;
      typedef Global::Transfer<LocalPresTransfer, PresMirror> GlobalPresTransfer;
      typedef Global::Transfer<LocalStressTransfer, StressMirror> GlobalStressTransfer;
      typedef Global::Transfer<LocalSystemTransfer, SystemMirror> GlobalSystemTransfer;

      /* ***************************************************************************************** */

      /// our system gate
      VeloGate gate_velo;
      PresGate gate_pres;
      StressGate gate_stress;
      SystemGate gate_sys;

      /// our coarse-level system muxer
      VeloMuxer coarse_muxer_velo;
      PresMuxer coarse_muxer_pres;
      StressMuxer coarse_muxer_stress;
      SystemMuxer coarse_muxer_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;
      GlobalMatrixBlockA matrix_a;
      GlobalMatrixBlockB matrix_b;
      GlobalMatrixBlockD matrix_d;
      GlobalMatrixBlockM matrix_m;
      GlobalMatrixBlockK matrix_k;
      GlobalMatrixBlockL matrix_l;

      /// null-matrix blocks
      NullMatrixBlockPP null_matrix_pp;
      NullMatrixBlockSP null_matrix_sp;
      NullMatrixBlockPS null_matrix_ps;

      /// our global transfer operator
      GlobalVeloTransfer transfer_velo;
      GlobalPresTransfer transfer_pres;
      GlobalStressTransfer transfer_stress;
      GlobalSystemTransfer transfer_sys;

      /// CTOR
      Stokes3FieldSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys),
        matrix_a(&gate_velo, &gate_velo),
        matrix_b(&gate_velo, &gate_pres),
        matrix_d(&gate_pres, &gate_velo),
        matrix_m(&gate_stress, &gate_stress),
        matrix_k(&gate_stress, &gate_velo),
        matrix_l(&gate_velo, &gate_stress),
        transfer_velo(&coarse_muxer_velo),
        transfer_pres(&coarse_muxer_pres),
        transfer_stress(&coarse_muxer_stress),
        transfer_sys(&coarse_muxer_sys)
      {
      }

      virtual ~Stokes3FieldSystemLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return (*this->matrix_sys).bytes () + transfer_sys.bytes();
      }

      void compile_system_matrix()
      {
        matrix_sys.local().template at<0,0>() = matrix_a.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().template at<0,1>() = matrix_b.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().template at<0,2>() = matrix_l.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().template at<1,0>() = matrix_d.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().template at<1,1>() = null_matrix_pp.clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().template at<1,2>() = null_matrix_ps.clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().template at<2,0>() = matrix_k.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().template at<2,1>() = null_matrix_sp.clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().template at<2,2>() = matrix_m.local().clone(LAFEM::CloneMode::Shallow);
      }

      void compile_system_transfer()
      {
        transfer_sys.get_mat_prol().template at<0,0>() = transfer_velo.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<0,0>() = transfer_velo.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_trunc().template at<0,0>() = transfer_velo.get_mat_trunc().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_prol().template at<1,1>() = transfer_pres.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<1,1>() = transfer_pres.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_trunc().template at<1,1>() = transfer_pres.get_mat_trunc().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_prol().template at<2,2>() = transfer_stress.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<2,2>() = transfer_stress.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_trunc().template at<2,2>() = transfer_stress.get_mat_trunc().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.compile();
      }

      template<typename DomainLevel_>
      void assemble_gates(const Domain::VirtualLevel<DomainLevel_>& virt_dom_lvl)
      {
        const auto& dom_level = virt_dom_lvl.level();
        const auto& dom_layer = virt_dom_lvl.layer();
        const auto& space_velo = dom_level.space_velo;
        const auto& space_pres = dom_level.space_pres;
        const auto& space_stress = dom_level.space_stress;

        // set the gate comm
        this->gate_velo.set_comm(dom_layer.comm_ptr());
        this->gate_pres.set_comm(dom_layer.comm_ptr());
        this->gate_stress.set_comm(dom_layer.comm_ptr());
        this->gate_sys.set_comm(dom_layer.comm_ptr());

        // loop over all ranks
        for(Index i(0); i < dom_layer.neighbour_count(); ++i)
        {
          int rank = dom_layer.neighbour_rank(i);

          // try to find our halo
          auto* halo = dom_level.find_halo_part(rank);
          XASSERT(halo != nullptr);

          // create (empty) velocity mirror
          VeloMirror mirror_velo;
          Assembly::MirrorAssembler::assemble_mirror(mirror_velo, space_velo, *halo);

          // create (empty) pressure mirror
          PresMirror mirror_pres;
          Assembly::MirrorAssembler::assemble_mirror(mirror_pres, space_pres, *halo);

          // create (empty) stress mirror
          StressMirror mirror_stress;
          Assembly::MirrorAssembler::assemble_mirror(mirror_stress, space_stress, *halo);

          // create a system mirror
          SystemMirror mirror_sys(mirror_velo.clone(), mirror_pres.clone(), mirror_stress.clone());

          // push mirrors into gates
          if(!mirror_velo.empty())
            this->gate_velo.push(rank, std::move(mirror_velo));
          if(!mirror_pres.empty())
            this->gate_pres.push(rank, std::move(mirror_pres));
          if(!mirror_stress.empty())
            this->gate_stress.push(rank, std::move(mirror_stress));
          if(!mirror_sys.empty())
            this->gate_sys.push(rank, std::move(mirror_sys));
        }

        // create local template vectors
        LocalVeloVector tmpl_v(space_velo.get_num_dofs());
        LocalPresVector tmpl_p(space_pres.get_num_dofs());
        LocalStressVector tmpl_s(space_stress.get_num_dofs());
        LocalSystemVector tmpl_sys(tmpl_v.clone(), tmpl_p.clone(), tmpl_s.clone());

        // compile gates
        this->gate_velo.compile(std::move(tmpl_v));
        this->gate_pres.compile(std::move(tmpl_p));
        this->gate_stress.compile(std::move(tmpl_s));
        this->gate_sys.compile(std::move(tmpl_sys));
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
            StressMirror& child_mirror_s = child_mirror_sys.template at<2>();
            Assembly::MirrorAssembler::assemble_mirror(child_mirror_v, level_p.space_velo, *child);
            Assembly::MirrorAssembler::assemble_mirror(child_mirror_p, level_p.space_pres, *child);
            Assembly::MirrorAssembler::assemble_mirror(child_mirror_s, level_p.space_stress, *child);
            this->coarse_muxer_velo.push_child(child_mirror_v.clone(LAFEM::CloneMode::Shallow));
            this->coarse_muxer_pres.push_child(child_mirror_p.clone(LAFEM::CloneMode::Shallow));
            this->coarse_muxer_stress.push_child(child_mirror_s.clone(LAFEM::CloneMode::Shallow));
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
          StressMirror& parent_mirror_s = parent_mirror_sys.template at<2>();

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
          {
            Index n = level_c.space_stress.get_num_dofs();
            parent_mirror_s = ScalarMirror(n, n);
            auto* idx = parent_mirror_s.indices();
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
          this->coarse_muxer_stress.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            parent_mirror_s.clone(LAFEM::CloneMode::Shallow)
          );
          this->coarse_muxer_sys.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            std::move(parent_mirror_sys)
          );

          // compile muxer
          LocalVeloVector tmpl_v(level_c.space_velo.get_num_dofs());
          LocalPresVector tmpl_p(level_c.space_pres.get_num_dofs());
          LocalStressVector tmpl_s(level_c.space_stress.get_num_dofs());
          LocalSystemVector tmpl_sys(tmpl_v.clone(), tmpl_p.clone(), tmpl_s.clone());
          this->coarse_muxer_velo.compile(tmpl_v);
          this->coarse_muxer_pres.compile(tmpl_p);
          this->coarse_muxer_stress.compile(tmpl_s);
          this->coarse_muxer_sys.compile(tmpl_sys);
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

          // synchronise blocked weight vector
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

          // synchronise weight vector
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
      void assemble_stress_transfer(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_stress;
        const auto& space_c = level_c.space_stress;

        // get local transfer operator
        LocalStressTransfer& loc_trans = this->transfer_stress.local();

        // get local transfer matrices
        LocalStressTransferMatrix& loc_prol_wrapped = loc_trans.get_mat_prol();
        LocalStressTransferMatrix& loc_rest_wrapped = loc_trans.get_mat_rest();

        // get the unwrapped types
        typename LocalStressTransferMatrix::BaseClass& loc_prol = loc_prol_wrapped;
        typename LocalStressTransferMatrix::BaseClass& loc_rest = loc_rest_wrapped;

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_f, space_c);
        }

        // create a local weight vector
        LocalStressVector loc_vec_weight = loc_prol_wrapped.create_vector_l();

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

          // synchronise blocked weight vector
          this->gate_stress.sync_0(loc_vec_weight);

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

        // compile stresscity transfer
        this->transfer_stress.compile();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_transfers(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        this->assemble_velocity_transfer(virt_lvl_fine, virt_lvl_coarse, cubature);
        this->assemble_pressure_transfer(virt_lvl_fine, virt_lvl_coarse, cubature);
        this->assemble_stress_transfer(virt_lvl_fine, virt_lvl_coarse, cubature);

        this->compile_system_transfer();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_velocity_truncation(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const Stokes3FieldSystemLevel* sys_lvl_coarse = nullptr)
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

        // We now need to synchronise the weight vector in analogy to the prolongation matrix assembly.
        // Note that the weight vector is now a coarse-level vector, so we need to synchronise over
        // the coarse-level gate. This may be a bit more complicated if the coarse level is a ghost
        // level, as in this case we have to join/split around the synch operation.

        // synchronise weight vector using the muxer/gate
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
        const Stokes3FieldSystemLevel* sys_lvl_coarse = nullptr)
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

        // We now need to synchronise the weight vector in analogy to the prolongation matrix assembly.
        // Note that the weight vector is now a coarse-level vector, so we need to synchronise over
        // the coarse-level gate. This may be a bit more complicated if the coarse level is a ghost
        // level, as in this case we have to join/split around the synch operation.

        // synchronise weight vector using the muxer/gate
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
      void assemble_stress_truncation(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const Stokes3FieldSystemLevel* sys_lvl_coarse = nullptr)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT((sys_lvl_coarse != nullptr) || !virt_lvl_coarse.is_parent());

        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_stress;
        const auto& space_c = level_c.space_stress;

        // get local transfer operator
        LocalStressTransfer& loc_trans = this->transfer_stress.local();

        // get local transfer matrices
        const LocalStressTransferMatrix& loc_rest_wrapped = loc_trans.get_mat_rest();
        LocalStressTransferMatrix& loc_trunc_wrapped = loc_trans.get_mat_trunc();

        // get the matrix blocks
        const typename LocalStressTransferMatrix::BaseClass& loc_rest = loc_rest_wrapped;
        typename LocalStressTransferMatrix::BaseClass& loc_trunc = loc_trunc_wrapped;

        // restriction matrix must be already assembled
        XASSERTM(loc_rest.size() > Index(0), "you need to call 'assemble_stresscity_transfer' first");

        // clone sparsity pattern of restriction matrix
        loc_trunc = loc_rest.clone(LAFEM::CloneMode::Layout);

        // create local weight vector
        LocalStressVector loc_vec_weight = loc_rest_wrapped.create_vector_l();

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

        // We now need to synchronise the weight vector in analogy to the prolongation matrix assembly.
        // Note that the weight vector is now a coarse-level vector, so we need to synchronise over
        // the coarse-level gate. This may be a bit more complicated if the coarse level is a ghost
        // level, as in this case we have to join/split around the synch operation.

        // synchronise weight vector using the muxer/gate
        if(!virt_lvl_coarse.is_child())
        {
          // The coarse level is a simple (non-child) level that exists on all processes,
          // so simply synch over the coarse-level gate:
          sys_lvl_coarse->gate_stress.sync_0(loc_vec_weight);
        }
        else if(virt_lvl_coarse.is_parent())
        {
          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then synch that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a synch over the (non-existent) coarse-level
          // child gate, which is what we actually require here...

          // create temporary vector on parent partitioning
          LocalStressVector loc_tmp = sys_lvl_coarse->gate_stress._freqs.clone(LAFEM::CloneMode::Allocate);

          // join child weights over muxer
          this->coarse_muxer_stress.join(loc_vec_weight, loc_tmp);

          // sync over coarse gate
          sys_lvl_coarse->gate_stress.sync_0(loc_tmp);

          // split over muxer
          this->coarse_muxer_stress.split(loc_vec_weight, loc_tmp);
        }
        else // ghost
        {
          // The coarse level is a ghost level, i.e. a child but not a parent. In this case, we
          // only have to participate in the join/send operations of the muxer, which are part
          // of the operations executed on the parents handled by the else-if case above.

          this->coarse_muxer_stress.join_send(loc_vec_weight);

          // parent performs sync over its gate here (see above else-if)

          this->coarse_muxer_stress.split_recv(loc_vec_weight);
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
      void assemble_truncations(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const Stokes3FieldSystemLevel* sys_lvl_coarse = nullptr)
      {
        this->assemble_velocity_truncation(virt_lvl_fine, virt_lvl_coarse, cubature, sys_lvl_coarse);
        this->assemble_pressure_truncation(virt_lvl_fine, virt_lvl_coarse, cubature, sys_lvl_coarse);
        this->assemble_stress_truncation(virt_lvl_fine, virt_lvl_coarse, cubature, sys_lvl_coarse);
        this->compile_system_transfer();
      }

      template<typename SpaceVelo_, typename SpacePres_, typename Cubature_>
      void assemble_grad_div_matrices(const SpaceVelo_& space_velo, const SpacePres_& space_pres, const Cubature_& cubature)
      {
        Assembly::GradPresDivVeloAssembler::assemble(this->matrix_b.local(), this->matrix_d.local(),
        space_velo, space_pres, cubature);
      }

      template<typename SpaceVelo_, typename SpacePres_, typename SpaceStress_>
      void assemble_structs(const SpaceVelo_& space_velo, const SpacePres_& space_pres, const SpaceStress_& space_stress)
      {
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_a.local(), space_velo);
        Assembly::SymbolicAssembler::assemble_matrix_std2(this->matrix_b.local(), space_velo, space_pres);
        Assembly::SymbolicAssembler::assemble_matrix_std2(this->matrix_d.local(), space_pres, space_velo);
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_m.local(), space_stress);
        Assembly::SymbolicAssembler::assemble_matrix_std2(this->matrix_l.local(), space_velo, space_stress);
        Assembly::SymbolicAssembler::assemble_matrix_std2(this->matrix_k.local(), space_stress, space_velo);
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->null_matrix_pp, space_pres);
        Assembly::SymbolicAssembler::assemble_matrix_std2(this->null_matrix_ps, space_pres, space_stress);
        Assembly::SymbolicAssembler::assemble_matrix_std2(this->null_matrix_sp, space_stress, space_pres);
      }
    }; // class Stokes3FieldSystemLevel<...>
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STOKES_3FIELD_HPP
