// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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

#include <control/domain/domain_control.hpp>
#include <control/asm/gate_asm.hpp>
#include <control/asm/muxer_asm.hpp>
#include <control/asm/splitter_asm.hpp>
#include <control/asm/transfer_asm.hpp>
#include <control/asm/transfer_voxel_asm.hpp>
#include <control/asm/unit_filter_asm.hpp>
#include <control/asm/mean_filter_asm.hpp>
#include <control/asm/slip_filter_asm.hpp>

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
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
      typename MatrixBlockM_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, nsc_, nsc_>,
      typename MatrixBlockK_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, nsc_, dim_>,
      typename MatrixBlockL_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, nsc_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixS_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, nsc_>
    >
    class Stokes3FieldSystemLevel
    {
    public:
      // basic types
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
      typedef LAFEM::NullMatrix<DataType_, IndexType_,   1,   1> NullMatrixBlockPP;
      typedef LAFEM::NullMatrix<DataType_, IndexType_, nsc,   1> NullMatrixBlockSP;
      typedef LAFEM::NullMatrix<DataType_, IndexType_,   1, nsc> NullMatrixBlockPS;
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
      typedef LAFEM::VectorMirror<DataType, IndexType> ScalarMirror;
      typedef ScalarMirror VeloMirror;
      typedef ScalarMirror PresMirror;
      typedef ScalarMirror StressMirror;
      typedef LAFEM::TupleMirror<VeloMirror, PresMirror, StressMirror> SystemMirror;

      // define gates
      typedef Global::Gate<LocalVeloVector, VeloMirror> VeloGate;
      typedef Global::Gate<LocalPresVector, PresMirror> PresGate;
      typedef Global::Gate<LocalStressVector, StressMirror> StressGate;
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;
      typedef Global::Gate<LocalScalarVector, ScalarMirror> ScalarGate;

      // define muxers
      typedef Global::Muxer<LocalVeloVector, VeloMirror> VeloMuxer;
      typedef Global::Muxer<LocalPresVector, PresMirror> PresMuxer;
      typedef Global::Muxer<LocalStressVector, StressMirror> StressMuxer;
      typedef Global::Muxer<LocalSystemVector, SystemMirror> SystemMuxer;
      typedef Global::Muxer<LocalScalarVector, ScalarMirror> ScalarMuxer;

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
      ScalarGate gate_scalar_velo;
      ScalarGate gate_scalar_stress;

      /// our coarse-level system muxer
      VeloMuxer coarse_muxer_velo;
      PresMuxer coarse_muxer_pres;
      StressMuxer coarse_muxer_stress;
      SystemMuxer coarse_muxer_sys;
      ScalarMuxer coarse_muxer_scalar_velo;
      ScalarMuxer coarse_muxer_scalar_stress;

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

      // no copies, no problems
      Stokes3FieldSystemLevel(const Stokes3FieldSystemLevel&) = delete;
      Stokes3FieldSystemLevel& operator=(const Stokes3FieldSystemLevel&) = delete;

      virtual ~Stokes3FieldSystemLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->matrix_sys.local().bytes () + transfer_sys.bytes();
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
        Asm::asm_gate(virt_dom_lvl, virt_dom_lvl->space_velo, this->gate_velo, true);
        Asm::asm_gate(virt_dom_lvl, virt_dom_lvl->space_pres, this->gate_pres, true);
        Asm::asm_gate(virt_dom_lvl, virt_dom_lvl->space_stress, this->gate_stress, true);
        Asm::build_gate_tuple(this->gate_sys, this->gate_velo, this->gate_pres, this->gate_stress);
        this->gate_scalar_velo.convert(this->gate_velo, LocalScalarVector(virt_dom_lvl->space_velo.get_num_dofs()));
        this->gate_scalar_stress.convert(this->gate_stress, LocalScalarVector(virt_dom_lvl->space_stress.get_num_dofs()));
      }

      template<typename DomainLevel_>
      void assemble_coarse_muxers(const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse)
      {
        Asm::asm_muxer(virt_lvl_coarse, [](const DomainLevel_& dl){return &dl.space_velo;}, this->coarse_muxer_velo);
        Asm::asm_muxer(virt_lvl_coarse, [](const DomainLevel_& dl){return &dl.space_pres;}, this->coarse_muxer_pres);
        Asm::asm_muxer(virt_lvl_coarse, [](const DomainLevel_& dl){return &dl.space_stress;}, this->coarse_muxer_stress);
        Asm::build_muxer_tuple(this->coarse_muxer_sys, this->gate_sys.get_freqs(), this->coarse_muxer_velo, this->coarse_muxer_pres, this->coarse_muxer_stress);
        this->coarse_muxer_scalar_velo.convert(this->coarse_muxer_velo, this->gate_scalar_velo.get_freqs().clone(LAFEM::CloneMode::Shallow));
        this->coarse_muxer_scalar_stress.convert(this->coarse_muxer_stress, this->gate_scalar_stress.get_freqs().clone(LAFEM::CloneMode::Shallow));
      }

      template<typename DomainLevel_>
      void assemble_transfers(
        const Stokes3FieldSystemLevel& sys_lvl_coarse,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const String& cubature, bool trunc_v = false, bool trunc_p = false, bool trunc_s = false, bool shrink = true)
      {
        Asm::asm_transfer_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_v, shrink,
          [](const DomainLevel_& dl) {return &dl.space_velo;},
          this->transfer_velo.local(), this->coarse_muxer_velo, this->gate_velo, sys_lvl_coarse.gate_velo);
        Asm::asm_transfer_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_p, shrink,
          [](const DomainLevel_& dl) {return &dl.space_pres;},
          this->transfer_pres.local(), this->coarse_muxer_pres, this->gate_pres, sys_lvl_coarse.gate_pres);
        Asm::asm_transfer_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_s, shrink,
          [](const DomainLevel_& dl) {return &dl.space_stress;},
          this->transfer_stress.local(), this->coarse_muxer_stress, this->gate_stress, sys_lvl_coarse.gate_stress);

        this->transfer_velo.compile();
        this->transfer_pres.compile();
        this->transfer_stress.compile();
        this->compile_system_transfer();
      }

      template<typename DomainLevel_>
      void assemble_transfers(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const String& cubature, bool trunc_v = false, bool trunc_p = false, bool trunc_s = false, bool shrink = true)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT(!virt_lvl_coarse.is_parent());

        Asm::asm_transfer_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_v, shrink,
          [](const DomainLevel_& dl) {return &dl.space_velo;},
          this->transfer_velo.local(), this->coarse_muxer_velo, this->gate_velo, this->gate_velo);
        Asm::asm_transfer_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_p, shrink,
          [](const DomainLevel_& dl) {return &dl.space_pres;},
          this->transfer_pres.local(), this->coarse_muxer_pres, this->gate_pres, this->gate_pres);
        Asm::asm_transfer_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_s, shrink,
          [](const DomainLevel_& dl) {return &dl.space_stress;},
          this->transfer_stress.local(), this->coarse_muxer_stress, this->gate_stress, this->gate_stress);

        this->transfer_velo.compile();
        this->transfer_pres.compile();
        this->transfer_stress.compile();
        this->compile_system_transfer();
      }

      template<typename SpaceVelo_, typename SpacePres_, typename Cubature_>
      void assemble_grad_div_matrices(const SpaceVelo_& space_velo, const SpacePres_& space_pres, const Cubature_& cubature)
      {
        Assembly::GradPresDivVeloAssembler::assemble(this->matrix_b.local(), this->matrix_d.local(), space_velo, space_pres, cubature);
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
