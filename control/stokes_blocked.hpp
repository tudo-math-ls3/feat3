// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/filter_sequence.hpp>
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
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
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
#include <kernel/global/splitter.hpp>

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
    template
    <
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixS_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>
    >
    class StokesBlockedSystemLevel
    {
    public:
      // basic types
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
      typedef TransferMatrixS_ LocalScalarTransferMatrix;
      typedef LAFEM::TupleDiagMatrix<TransferMatrixV_, TransferMatrixP_> LocalSystemTransferMatrix;

      // define local transfer operators
      typedef LAFEM::Transfer<LocalVeloTransferMatrix> LocalVeloTransfer;
      typedef LAFEM::Transfer<LocalPresTransferMatrix> LocalPresTransfer;
      typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;
      typedef LAFEM::Transfer<LocalScalarTransferMatrix> LocalScalarTransfer;

      // define mirror types
      typedef LAFEM::VectorMirror<DataType, IndexType> ScalarMirror;
      typedef ScalarMirror VeloMirror;
      typedef ScalarMirror PresMirror;
      typedef LAFEM::TupleMirror<VeloMirror, PresMirror> SystemMirror;

      // define gates
      typedef Global::Gate<LocalVeloVector, VeloMirror> VeloGate;
      typedef Global::Gate<LocalPresVector, PresMirror> PresGate;
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;
      typedef Global::Gate<LocalScalarVector, ScalarMirror> ScalarGate;

      // define muxers
      typedef Global::Muxer<LocalVeloVector, VeloMirror> VeloMuxer;
      typedef Global::Muxer<LocalPresVector, PresMirror> PresMuxer;
      typedef Global::Muxer<LocalSystemVector, SystemMirror> SystemMuxer;
      typedef Global::Muxer<LocalScalarVector, ScalarMirror> ScalarMuxer;

      // define splitters
      typedef Global::Splitter<LocalVeloVector, VeloMirror> VeloSplitter;
      typedef Global::Splitter<LocalPresVector, PresMirror> PresSplitter;
      typedef Global::Splitter<LocalSystemVector, SystemMirror> SystemSplitter;
      typedef Global::Splitter<LocalScalarVector, ScalarMirror> ScalarSplitter;

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
      typedef Global::Transfer<LocalScalarTransfer, ScalarMirror> GlobalScalarTransfer;

      /* ***************************************************************************************** */

      /// our system gate
      VeloGate gate_velo;
      PresGate gate_pres;
      SystemGate gate_sys;
      ScalarGate gate_scalar_velo;

      /// our coarse-level system muxer
      VeloMuxer coarse_muxer_velo;
      PresMuxer coarse_muxer_pres;
      SystemMuxer coarse_muxer_sys;
      ScalarMuxer coarse_muxer_scalar_velo;

      /// our base-mesh multiplexer
      VeloSplitter base_splitter_velo;
      PresSplitter base_splitter_pres;
      SystemSplitter base_splitter_sys;
      ScalarSplitter base_splitter_scalar_velo;

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
      GlobalScalarTransfer transfer_scalar_velo;

      /// local type-1 system matrix
      LocalSystemMatrix local_matrix_sys_type1;

      /// CTOR
      StokesBlockedSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys),
        matrix_a(&gate_velo, &gate_velo),
        matrix_b(&gate_velo, &gate_pres),
        matrix_d(&gate_pres, &gate_velo),
        matrix_s(&gate_pres, &gate_pres),
        transfer_velo(&coarse_muxer_velo),
        transfer_pres(&coarse_muxer_pres),
        transfer_sys(&coarse_muxer_sys),
        transfer_scalar_velo(&coarse_muxer_scalar_velo)
      {
      }

      // no copies, no problems
      StokesBlockedSystemLevel(const StokesBlockedSystemLevel&) = delete;
      StokesBlockedSystemLevel& operator=(const StokesBlockedSystemLevel&) = delete;

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

      void compile_scalar_transfer()
      {
        transfer_scalar_velo.get_mat_prol().clone(this->transfer_velo.get_mat_prol().unwrap(), LAFEM::CloneMode::Shallow);
        transfer_scalar_velo.get_mat_rest().clone(this->transfer_velo.get_mat_rest().unwrap(), LAFEM::CloneMode::Shallow);
        transfer_scalar_velo.get_mat_trunc().clone(this->transfer_velo.get_mat_trunc().unwrap(), LAFEM::CloneMode::Shallow);
        transfer_scalar_velo.compile();
      }

      template<typename D_, typename I_, typename SMA_, typename SMB_, typename SMD_, typename SM_, typename TV_, typename TP_>
      void convert(const StokesBlockedSystemLevel<dim_, D_, I_, SMA_, SMB_, SMD_, SM_, TV_, TP_> & other)
      {
        gate_velo.convert(other.gate_velo);
        gate_pres.convert(other.gate_pres);
        Asm::build_gate_tuple(this->gate_sys, this->gate_velo, this->gate_pres);
        this->gate_scalar_velo.convert(this->gate_velo, LocalScalarVector(this->gate_velo.get_freqs().template size<LAFEM::Perspective::native>()));

        coarse_muxer_velo.convert(other.coarse_muxer_velo);
        coarse_muxer_pres.convert(other.coarse_muxer_pres);
        Asm::build_muxer_tuple(this->coarse_muxer_sys, this->gate_sys.get_freqs(), this->coarse_muxer_velo, this->coarse_muxer_pres);
        this->coarse_muxer_scalar_velo.convert(this->coarse_muxer_velo, this->gate_scalar_velo.get_freqs().clone(LAFEM::CloneMode::Shallow));

        base_splitter_velo.convert(other.base_splitter_velo);
        base_splitter_pres.convert(other.base_splitter_pres);
        Asm::build_splitter_tuple(this->base_splitter_sys, this->gate_sys.get_freqs(), this->base_splitter_velo, this->base_splitter_pres);
        this->base_splitter_scalar_velo.convert(this->base_splitter_velo, this->gate_scalar_velo.get_freqs().clone(LAFEM::CloneMode::Shallow));

        transfer_velo.convert(&coarse_muxer_velo, other.transfer_velo);
        transfer_pres.convert(&coarse_muxer_pres, other.transfer_pres);
        this->compile_system_transfer();
        this->compile_scalar_transfer();

        matrix_a.convert(&gate_velo, &gate_velo, other.matrix_a);
        matrix_b.convert(&gate_velo, &gate_pres, other.matrix_b);
        matrix_d.convert(&gate_pres, &gate_velo, other.matrix_d);
        matrix_s.convert(&gate_pres, &gate_pres, other.matrix_s);
        this->compile_system_matrix();
      }

      GlobalSystemVector create_global_vector_sys() const
      {
        return GlobalSystemVector(&this->gate_sys, this->gate_sys.get_freqs().clone(LAFEM::CloneMode::Layout));
      }

      GlobalVeloVector create_global_vector_velo() const
      {
        return GlobalVeloVector(&this->gate_velo, this->gate_velo.get_freqs().clone(LAFEM::CloneMode::Layout));
      }

      GlobalPresVector create_global_vector_pres() const
      {
        return GlobalPresVector(&this->gate_pres, this->gate_pres.get_freqs().clone(LAFEM::CloneMode::Layout));
      }

      template<typename DomainLevel_>
      void assemble_gates(const Domain::VirtualLevel<DomainLevel_>& virt_dom_lvl)
      {
        Asm::asm_gate(virt_dom_lvl, virt_dom_lvl->space_velo, this->gate_velo, true);
        Asm::asm_gate(virt_dom_lvl, virt_dom_lvl->space_pres, this->gate_pres, true);
        Asm::build_gate_tuple(this->gate_sys, this->gate_velo, this->gate_pres);
        this->gate_scalar_velo.convert(this->gate_velo, LocalScalarVector(virt_dom_lvl->space_velo.get_num_dofs()));
      }

      template<typename DomainLevel_>
      void assemble_coarse_muxers(const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse)
      {
        Asm::asm_muxer(virt_lvl_coarse, [](const DomainLevel_& dl){return &dl.space_velo;}, this->coarse_muxer_velo);
        Asm::asm_muxer(virt_lvl_coarse, [](const DomainLevel_& dl){return &dl.space_pres;}, this->coarse_muxer_pres);
        Asm::build_muxer_tuple(this->coarse_muxer_sys, this->gate_sys.get_freqs(), this->coarse_muxer_velo, this->coarse_muxer_pres);
        this->coarse_muxer_scalar_velo.convert(this->coarse_muxer_velo, this->gate_scalar_velo.get_freqs().clone(LAFEM::CloneMode::Shallow));
      }

      template<typename DomainLevel_>
      void assemble_base_splitters(const Domain::VirtualLevel<DomainLevel_>& virt_lvl)
      {
        Asm::asm_splitter(virt_lvl, [](const DomainLevel_& dl){return &dl.space_velo;}, this->base_splitter_velo);
        Asm::asm_splitter(virt_lvl, [](const DomainLevel_& dl){return &dl.space_pres;}, this->base_splitter_pres);
        Asm::build_splitter_tuple(this->base_splitter_sys, this->gate_sys.get_freqs(), this->base_splitter_velo, this->base_splitter_pres);
        this->base_splitter_scalar_velo.convert(this->base_splitter_velo, this->gate_scalar_velo.get_freqs().clone(LAFEM::CloneMode::Shallow));
      }

      template<typename DomainLevel_>
      void assemble_transfers(
        const StokesBlockedSystemLevel& sys_lvl_coarse,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const String& cubature, bool trunc_v = false, bool trunc_p = false, bool shrink = true)
      {
        Asm::asm_transfer_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_v, shrink,
          [](const DomainLevel_& dl) {return &dl.space_velo;},
          this->transfer_velo.local(), this->coarse_muxer_velo, this->gate_velo, sys_lvl_coarse.gate_velo);
        Asm::asm_transfer_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_p, shrink,
          [](const DomainLevel_& dl) {return &dl.space_pres;},
          this->transfer_pres.local(), this->coarse_muxer_pres, this->gate_pres, sys_lvl_coarse.gate_pres);

        this->transfer_velo.compile();
        this->transfer_pres.compile();
        this->compile_system_transfer();
        this->compile_scalar_transfer();
      }

      template<typename DomainLevel_>
      void assemble_transfers(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const String& cubature, bool trunc_v = false, bool trunc_p = false, bool shrink = true)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT(!virt_lvl_coarse.is_parent());

        Asm::asm_transfer_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_v, shrink,
          [](const DomainLevel_& dl) {return &dl.space_velo;},
          this->transfer_velo.local(), this->coarse_muxer_velo, this->gate_velo, this->gate_velo);
        Asm::asm_transfer_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_p, shrink,
          [](const DomainLevel_& dl) {return &dl.space_pres;},
          this->transfer_pres.local(), this->coarse_muxer_pres, this->gate_pres, this->gate_pres);

        this->transfer_velo.compile();
        this->transfer_pres.compile();
        this->compile_system_transfer();
        this->compile_scalar_transfer();
      }

      template<typename DomainLevel_>
      void assemble_transfers_voxel(
        const StokesBlockedSystemLevel& sys_lvl_coarse,
        const Domain::VirtualLevel<Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_fine,
        const Domain::VirtualLevel<Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_coarse,
        const String& cubature, bool trunc_v = false, bool trunc_p = false, bool shrink = true)
      {
        Asm::asm_transfer_voxel_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_v, shrink,
          [](const DomainLevel_& dl) {return &dl.space_velo;},
          this->transfer_velo.local(), this->coarse_muxer_velo, this->gate_velo, sys_lvl_coarse.gate_velo);
        Asm::asm_transfer_voxel_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_p, shrink,
          [](const DomainLevel_& dl) {return &dl.space_pres;},
          this->transfer_pres.local(), this->coarse_muxer_pres, this->gate_pres, sys_lvl_coarse.gate_pres);

        this->transfer_velo.compile();
        this->transfer_pres.compile();
        this->compile_system_transfer();
        this->compile_scalar_transfer();
      }

      template<typename DomainLevel_>
      void assemble_transfers_voxel(
        const Domain::VirtualLevel<Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_fine,
        const Domain::VirtualLevel<Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_coarse,
        const String& cubature, bool trunc_v = false, bool trunc_p = false, bool shrink = true)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT(!virt_lvl_coarse.is_parent());

        Asm::asm_transfer_voxel_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_v, shrink,
          [](const DomainLevel_& dl) {return &dl.space_velo;},
          this->transfer_velo.local(), this->coarse_muxer_velo, this->gate_velo, this->gate_velo);
        Asm::asm_transfer_voxel_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc_p, shrink,
          [](const DomainLevel_& dl) {return &dl.space_pres;},
          this->transfer_pres.local(), this->coarse_muxer_pres, this->gate_pres, this->gate_pres);

        this->transfer_velo.compile();
        this->transfer_pres.compile();
        this->compile_system_transfer();
        this->compile_scalar_transfer();
      }

      template<typename SpaceVelo_, typename SpacePres_, typename Cubature_>
      void assemble_grad_div_matrices(const SpaceVelo_& space_velo, const SpacePres_& space_pres, const Cubature_& cubature)
      {
        Assembly::GradPresDivVeloAssembler::assemble(this->matrix_b.local(), this->matrix_d.local(), space_velo, space_pres, cubature);
      }

      template<typename Trafo_, typename SpaceVelo_, typename SpacePres_>
      void assemble_grad_div_matrices(Assembly::DomainAssembler<Trafo_>& dom_asm,
        const SpaceVelo_& space_velo, const SpacePres_& space_pres, const String& cubature)
      {
        // assemble matrix structure of B
        if(this->matrix_b.local().empty())
          Assembly::SymbolicAssembler::assemble_matrix_std2(this->matrix_b.local(), space_velo, space_pres);

        // assemble matrix B
        this->matrix_b.local().format();
        Assembly::Common::GradientTestOperatorBlocked<dim_> grad_op;
        Assembly::assemble_bilinear_operator_matrix_2(
          dom_asm, this->matrix_b.local(), grad_op, space_velo, space_pres, cubature, -DataType(1));

        // transpose to obtain matrix D
        this->matrix_d.local().transpose(this->matrix_b.local());
      }

      template<typename Trafo_, typename SpaceVelo_, typename SpacePres_, typename AsmDT_>
      void assemble_grad_div_matrices_high_prec(Assembly::DomainAssembler<Trafo_>& dom_asm,
        const SpaceVelo_& space_velo, const SpacePres_& space_pres, const String& cubature, const AsmDT_ asm_weight)
      {
        // assemble matrix structure of B
        if(this->matrix_b.local().empty())
          Assembly::SymbolicAssembler::assemble_matrix_std2(this->matrix_b.local(), space_velo, space_pres);

        // assemble matrix B
        this->matrix_b.local().format();
        {
          // create local matrix with asm datatype
          typename LocalMatrixBlockB::template ContainerTypeByDI<AsmDT_, IndexType_> tmp_mat_b;
          tmp_mat_b.convert(this->matrix_b.local());
          tmp_mat_b.format();

          Assembly::Common::GradientTestOperatorBlocked<dim_> grad_op;
          Assembly::assemble_bilinear_operator_matrix_2(
            dom_asm, tmp_mat_b, grad_op, space_velo, space_pres, cubature, -asm_weight);

          // convert matrix
          this->matrix_b.local().convert(tmp_mat_b);
        }

        // transpose to obtain matrix D
        this->matrix_d.local().transpose(this->matrix_b.local());
      }

      template<typename SpaceVelo_>
      void assemble_velo_struct(const SpaceVelo_& space_velo)
      {
        // assemble matrix structure
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_a.local(), space_velo);
        this->matrix_a.local().format();
      }

      template<typename SpacePres_>
      void assemble_pres_struct(const SpacePres_& space_pres)
      {
        // assemble matrix structure
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_s.local(), space_pres);
        this->matrix_s.local().format();
      }

      void compile_local_matrix_sys_type1()
      {
        this->local_matrix_sys_type1.block_a() = matrix_a.convert_to_1();
        this->local_matrix_sys_type1.block_b() = matrix_b.local().clone(LAFEM::CloneMode::Weak);
        this->local_matrix_sys_type1.block_d() = matrix_d.local().clone(LAFEM::CloneMode::Weak);
      }
    }; // class StokesBlockedSystemLevel<...>

    template
    <
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>
    >
    class StokesBlockedUnitVeloNonePresSystemLevel :
      public StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::UnitFilterBlocked<DataType_, IndexType_, dim_> LocalVeloFilter;
      typedef LAFEM::NoneFilter<DataType_, IndexType_> LocalPresFilter;
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

      void compile_local_matrix_sys_type1()
      {
        BaseClass::compile_local_matrix_sys_type1();

        // apply filter to A and B
        this->filter_velo.local().filter_mat(this->local_matrix_sys_type1.block_a());
        this->filter_velo.local().filter_offdiag_row_mat(this->local_matrix_sys_type1.block_b());
      }
    }; // class StokesBlockedUnitVeloNonePresSystemLevel

    template
    <
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>
    >
    class StokesBlockedSlipUnitVeloNonePresSystemLevel :
      public StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
      MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::SlipFilter<DataType_, IndexType_, dim_> LocalVeloSlipFilter;
      typedef LAFEM::UnitFilterBlocked<DataType_, IndexType_, dim_> LocalVeloUnitFilter;
      typedef LAFEM::FilterChain<LocalVeloSlipFilter, LocalVeloUnitFilter> LocalVeloFilter;
      typedef LAFEM::NoneFilter<DataType_, IndexType_> LocalPresFilter;
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

      void compile_local_matrix_sys_type1()
      {
        BaseClass::compile_local_matrix_sys_type1();

        // apply filter to A and B
        this->filter_velo.local().template at<1>().filter_mat(this->local_matrix_sys_type1.block_a());
        this->filter_velo.local().template at<1>().filter_offdiag_row_mat(this->local_matrix_sys_type1.block_b());
      }
    }; // class StokesBlockedSlipUnitVeloNonePresSystemLevel

    /**
     * \brief System level using a MeanFilter for the pressure
     *
     * This is necessary when there are only Dirichlet BCs for the velocity
     */
    template
    <
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>
    >
    struct StokesBlockedUnitVeloMeanPresSystemLevel :
      public StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::UnitFilterBlocked<DataType_, IndexType_, dim_> LocalVeloFilter;
      typedef Global::MeanFilter<DataType_, IndexType_> LocalPresFilter;
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

      template<typename D_, typename I_, typename SM_>
      void convert(const StokesBlockedUnitVeloMeanPresSystemLevel<dim_, D_, I_, SM_> & other)
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

      void compile_local_matrix_sys_type1()
      {
        BaseClass::compile_local_matrix_sys_type1();

        // apply filter to A and B
        this->filter_velo.local().filter_mat(this->local_matrix_sys_type1.block_a());
        this->filter_velo.local().filter_offdiag_row_mat(this->local_matrix_sys_type1.block_b());
      }
    }; // struct StokesBlockedUnitVeloMeanPresSystemLevel<...>

    template
    <
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>
    >
    class StokesBlockedSlipUnitVeloMeanPresSystemLevel :
      public StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::SlipFilter<DataType_, IndexType_, dim_> LocalVeloSlipFilter;
      typedef LAFEM::UnitFilterBlocked<DataType_, IndexType_, dim_> LocalVeloUnitFilter;
      typedef LAFEM::FilterChain<LocalVeloSlipFilter, LocalVeloUnitFilter> LocalVeloFilter;
      typedef Global::MeanFilter<DataType_, IndexType_> LocalPresFilter;
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

      template<typename D_, typename I_, typename SM_>
      void convert(const StokesBlockedUnitVeloMeanPresSystemLevel<dim_, D_, I_, SM_> & other)
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

        // synchronize the slip filter
        Asm::sync_slip_filter(this->gate_velo, filter_velo.local().template at<0>());
      }

      void compile_local_matrix_sys_type1()
      {
        BaseClass::compile_local_matrix_sys_type1();

        // apply filter to A and B
        this->filter_velo.local().template at<1>().filter_mat(this->local_matrix_sys_type1.block_a());
        this->filter_velo.local().template at<1>().filter_offdiag_row_mat(this->local_matrix_sys_type1.block_b());
      }
    }; // class StokesBlockedSlipUnitVeloMeanPresSystemLevel

    /**
     * \brief Stokes blocked System level combining all supported types of boundary conditions
     *
     * \author Peter Zajac
     */
    template
    <
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>
    >
    class StokesBlockedCombinedSystemLevel :
      public StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local velocity filter chain
      typedef LAFEM::SlipFilter<DataType_, IndexType_, dim_> LocalVeloSlipFilter;
      typedef LAFEM::UnitFilterBlocked<DataType_, IndexType_, dim_> LocalVeloUnitFilter;
      typedef LAFEM::FilterSequence<LocalVeloSlipFilter> LocalVeloSlipFilterSeq;
      typedef LAFEM::FilterSequence<LocalVeloUnitFilter> LocalVeloUnitFilterSeq;
      typedef LAFEM::FilterChain<LocalVeloSlipFilterSeq, LocalVeloUnitFilterSeq> LocalVeloFilter;

      // define local pressure filter chain
      typedef Global::MeanFilter<DataType_, IndexType_> LocalPresMeanFilter;
      typedef LAFEM::UnitFilter<DataType_, IndexType_> LocalPresUnitFilter;
      typedef LAFEM::FilterChain<LocalPresMeanFilter, LocalPresUnitFilter> LocalPresFilter;

      // define local system filter
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

      template<typename D_, typename I_, typename SM_>
      void convert(const StokesBlockedCombinedSystemLevel<dim_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_velo.convert(other.filter_velo);
        filter_pres.convert(other.filter_pres);

        compile_system_filter();
      }

      LocalVeloSlipFilterSeq& get_local_velo_slip_filter_seq()
      {
        return this->filter_velo.local().template at<0>();
      }

      LocalVeloUnitFilterSeq& get_local_velo_unit_filter_seq()
      {
        return this->filter_velo.local().template at<1>();
      }

      LocalPresMeanFilter& get_local_pres_mean_filter()
      {
        return this->filter_pres.local().template at<0>();
      }

      LocalPresUnitFilter& get_local_pres_unit_filter()
      {
        return this->filter_pres.local().template at<1>();
      }

      template<typename DomainLevel_, typename Space_>
      void assemble_velocity_unit_filter(const DomainLevel_& dom_level, const Space_& space, const String& name, const String& mesh_parts)
      {
        auto& loc_filter = get_local_velo_unit_filter_seq().find_or_add(name);
        loc_filter.clear();
        Asm::asm_unit_filter_blocked_homogeneous(loc_filter, dom_level, space, mesh_parts);
      }

      template<typename DomainLevel_, typename Space_, typename Function_>
      void assemble_velocity_unit_filter(const DomainLevel_& dom_level, const Space_& space, const String& name, const String& mesh_parts, const Function_& function)
      {
        auto& loc_filter = get_local_velo_unit_filter_seq().find_or_add(name);
        loc_filter.clear();
        Asm::asm_unit_filter_blocked(loc_filter, dom_level, space, mesh_parts, function);
      }

      template<typename DomainLevel_, typename Space_>
      void assemble_velocity_slip_filter(const DomainLevel_& dom_level, const Space_& space, const String& name, const String& mesh_parts)
      {
        auto& loc_filter = get_local_velo_slip_filter_seq().find_or_add(name);
        loc_filter.clear();
        Asm::asm_slip_filter(loc_filter, dom_level, space, mesh_parts);
        Asm::sync_slip_filter(this->gate_velo, loc_filter);
      }

      template<typename SpacePres_>
      void assemble_pressure_mean_filter(const SpacePres_& space_pres, const String& cubature_name)
      {
        this->get_local_pres_mean_filter() = Asm::asm_mean_filter(this->gate_pres, space_pres, cubature_name);
      }

      void sync_velocity_slip_filters()
      {
        for(auto& it : get_local_velo_slip_filter_seq())
          Asm::sync_slip_filter(this->gate_velo, it.second);
      }

      void compile_local_matrix_sys_type1()
      {
        BaseClass::compile_local_matrix_sys_type1();
        for(const auto& loc_fil : get_local_velo_unit_filter_seq())
        {
          loc_fil.second.filter_mat(this->local_matrix_sys_type1.block_a());
          loc_fil.second.filter_offdiag_row_mat(this->local_matrix_sys_type1.block_b());
        }
      }
    }; // class StokesBlockedCombinedSystemLevel<...>
  } // namespace Control
} // namespace FEAT
