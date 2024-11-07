// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
#include <kernel/lafem/saddle_point_matrix.hpp>
//#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/tuple_diag_matrix.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/domain_assembler_helpers.hpp>
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
#include <kernel/global/symmetric_lumped_schur_matrix.hpp>
#include <kernel/global/transfer.hpp>

#include <control/domain/domain_control.hpp>
#include <control/asm/gate_asm.hpp>
#include <control/asm/muxer_asm.hpp>
#include <control/asm/splitter_asm.hpp>
#include <control/asm/transfer_asm.hpp>

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
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>
    >
    class ScalarMixedSystemLevel
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
      typedef LAFEM::SaddlePointMatrix<LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockD> LocalSystemMatrix;

      // define local vector types
      typedef typename LocalMatrixBlockB::VectorTypeL LocalVeloVector;
      typedef typename LocalMatrixBlockD::VectorTypeL LocalPresVector;
      typedef LAFEM::TupleVector<LocalVeloVector, LocalPresVector> LocalSystemVector;

      // define local filter types
      typedef LAFEM::SlipFilter<DataType_, IndexType_, dim_> LocalVeloSlipFilter;
      typedef LAFEM::FilterSequence<LocalVeloSlipFilter> LocalVeloFilter;
      typedef LAFEM::NoneFilter<DataType_, IndexType_> LocalPresFilter;
      typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresFilter> LocalSystemFilter;

      // define local transfer matrix types
      typedef TransferMatrixV_ LocalVeloTransferMatrix;
      typedef TransferMatrixP_ LocalPresTransferMatrix;
      typedef LAFEM::TupleDiagMatrix<TransferMatrixV_, TransferMatrixP_> LocalSystemTransferMatrix;

      // define local transfer operators
      typedef LAFEM::Transfer<LocalVeloTransferMatrix> LocalVeloTransfer;
      typedef LAFEM::Transfer<LocalPresTransferMatrix> LocalPresTransfer;
      typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

      // define mirror types
      typedef LAFEM::VectorMirror<DataType, IndexType> ScalarMirror;
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

      // define global filter types
      typedef Global::Filter<LocalVeloFilter, VeloMirror> GlobalVeloFilter;
      typedef Global::Filter<LocalPresFilter, PresMirror> GlobalPresFilter;
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      // define global matrix types
      typedef Global::Matrix<LocalMatrixBlockA, VeloMirror, VeloMirror> GlobalMatrixBlockA;
      typedef Global::Matrix<LocalMatrixBlockB, VeloMirror, PresMirror> GlobalMatrixBlockB;
      typedef Global::Matrix<LocalMatrixBlockD, PresMirror, VeloMirror> GlobalMatrixBlockD;
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
      GlobalMatrixBlockA matrix_a;
      GlobalVeloVector lumped_matrix_a;
      GlobalMatrixBlockB matrix_b;
      GlobalMatrixBlockD matrix_d;
      GlobalSystemMatrix matrix_sys;

      /// our global transfer operator
      GlobalVeloTransfer transfer_velo;
      GlobalPresTransfer transfer_pres;
      GlobalSystemTransfer transfer_sys;

      // (global) filters
      GlobalSystemFilter filter_sys;
      GlobalVeloFilter filter_velo;
      GlobalPresFilter filter_pres;

      /**
       * \brief Empty standard constructor
       */
      ScalarMixedSystemLevel() :
        matrix_a(&gate_velo, &gate_velo),
        lumped_matrix_a(&gate_velo),
        matrix_b(&gate_velo, &gate_pres),
        matrix_d(&gate_pres, &gate_velo),
        matrix_sys(&gate_sys, &gate_sys),
        transfer_velo(&coarse_muxer_velo),
        transfer_pres(&coarse_muxer_pres),
        transfer_sys(&coarse_muxer_sys)
      {
      }

      /**
       * \brief Constructor for using essential boundary conditions
       *
       * \param[in] neumann_list
       * List of MeshPart names where essential boundary conditions (in this case: Neumann) are to be enforced.
       *
       */
      explicit ScalarMixedSystemLevel(const std::deque<String>& neumann_list) :
        matrix_a(&gate_velo, &gate_velo),
        lumped_matrix_a(&gate_velo),
        matrix_b(&gate_velo, &gate_pres),
        matrix_d(&gate_pres, &gate_velo),
        transfer_velo(&coarse_muxer_velo),
        transfer_pres(&coarse_muxer_pres),
        transfer_sys(&coarse_muxer_sys),
        filter_velo(neumann_list)
      {
      }

      // no copies, no problems
      ScalarMixedSystemLevel(const ScalarMixedSystemLevel&) = delete;
      ScalarMixedSystemLevel& operator=(const ScalarMixedSystemLevel&) = delete;

      virtual ~ScalarMixedSystemLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        //return (*this->matrix_sys).bytes () + (*this->matrix_s).bytes() + transfer_sys.bytes();
        return matrix_sys.bytes() + filter_sys.bytes() + transfer_sys.bytes();
      }

      void compile_system_matrix()
      {
        if(lumped_matrix_a.local().size() == Index(0))
        {
          lumped_matrix_a.local() = matrix_a.local().create_vector_l();
        }

        matrix_a.lump_rows(lumped_matrix_a);

        matrix_sys.local().block_a() = matrix_a.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().block_b() = matrix_b.local().clone(LAFEM::CloneMode::Shallow);
        matrix_sys.local().block_d() = matrix_d.local().clone(LAFEM::CloneMode::Shallow);
      }

      void compile_system_transfer()
      {
        transfer_sys.get_mat_prol().template at<0,0>()
          = transfer_velo.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<0,0>()
          = transfer_velo.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_prol().template at<1,1>()
          = transfer_pres.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<1,1>()
          = transfer_pres.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
      }

      template<typename D_, typename I_, typename SMA_, typename SMB_, typename SMD_, typename SM_, typename TV_, typename TP_>
      void convert(const ScalarMixedSystemLevel<dim_, D_, I_, SMA_, SMB_, SMD_, SM_, TV_, TP_> & other)
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

        transfer_velo.convert(&coarse_muxer_velo, other.transfer_velo);
        transfer_pres.convert(&coarse_muxer_pres, other.transfer_pres);

        this->compile_system_matrix();
        this->compile_system_transfer();

        filter_velo.convert(other.filter_velo);
        filter_pres.convert(other.filter_pres);

        this->compile_system_filter();
      }

      template<typename DomainLevel_>
      void assemble_gates(const Domain::VirtualLevel<DomainLevel_>& virt_dom_lvl)
      {
        Asm::asm_gate(virt_dom_lvl, virt_dom_lvl->space_velo, this->gate_velo, true);
        Asm::asm_gate(virt_dom_lvl, virt_dom_lvl->space_pres, this->gate_pres, true);
        Asm::build_gate_tuple(this->gate_sys, this->gate_velo, this->gate_pres);
        //this->gate_scalar_velo.convert(this->gate_velo, LocalScalarVector(virt_dom_lvl->space_velo.get_num_dofs()));
      }

      template<typename DomainLevel_>
      void assemble_coarse_muxers(const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse)
      {
        Asm::asm_muxer(virt_lvl_coarse, [](const DomainLevel_& dl){return &dl.space_velo;}, this->coarse_muxer_velo);
        Asm::asm_muxer(virt_lvl_coarse, [](const DomainLevel_& dl){return &dl.space_pres;}, this->coarse_muxer_pres);
        Asm::build_muxer_tuple(this->coarse_muxer_sys, this->gate_sys.get_freqs(), this->coarse_muxer_velo, this->coarse_muxer_pres);
        //this->coarse_muxer_scalar_velo.convert(this->coarse_muxer_velo, this->gate_scalar_velo.get_freqs().clone(LAFEM::CloneMode::Shallow));
      }

      template<typename DomainLevel_>
      void assemble_base_splitters(const Domain::VirtualLevel<DomainLevel_>& virt_lvl)
      {
        Asm::asm_splitter(virt_lvl, [](const DomainLevel_& dl){return &dl.space_velo;}, this->base_splitter_velo);
        Asm::asm_splitter(virt_lvl, [](const DomainLevel_& dl){return &dl.space_pres;}, this->base_splitter_pres);
        Asm::build_splitter_tuple(this->base_splitter_sys, this->gate_sys.get_freqs(), this->base_splitter_velo, this->base_splitter_pres);
        //this->base_splitter_scalar_velo.convert(this->base_splitter_velo, this->gate_scalar_velo.get_freqs().clone(LAFEM::CloneMode::Shallow));
      }

      template<typename DomainLevel_>
      void assemble_transfers(
        const ScalarMixedSystemLevel& sys_lvl_coarse,
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
      }

      template<typename SpaceVelo_, typename SpacePres_, typename Cubature_>
      void assemble_grad_div_matrices(const SpaceVelo_& space_velo, const SpacePres_& space_pres, const Cubature_& cubature)
      {
        Assembly::GradPresDivVeloAssembler::assemble(this->matrix_b.local(), this->matrix_d.local(),
        space_velo, space_pres, cubature, DataType(1), DataType(1));
      }

      template<typename SpaceVelo_>
      void assemble_velo_struct(const SpaceVelo_& space_velo)
      {
        // assemble matrix structure
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_a.local(), space_velo);
      }

      void compile_system_filter()
      {
        filter_sys.local().template at<0>() = filter_velo.local().clone(LAFEM::CloneMode::Shallow);
        filter_sys.local().template at<1>() = filter_pres.local().clone(LAFEM::CloneMode::Shallow);
      }

      void assemble_global_filters()
      {
        // Sync the filter vector in the SlipFilter
        const VeloGate& my_col_gate(this->gate_velo);

        // For all slip filters...
        for(auto& it : filter_sys.local().template at<0>())
        {

          // Get the filter vector
          auto& slip_filter_vector = it.second.get_filter_vector();

          // If the filter  is not empty (meaning the MeshPart is on our patch):
          if(slip_filter_vector.used_elements() > 0)
          {
            // Temporary DenseVector for syncing
            LocalVeloVector tmp(slip_filter_vector.size(), DataType_(0));

            auto* tmp_elements = tmp.template elements<LAFEM::Perspective::native>();
            auto* sfv_elements = slip_filter_vector.template elements<LAFEM::Perspective::native>();

            // Copy sparse filter vector contents to DenseVector
            for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
            {
              Index idense(slip_filter_vector.indices()[isparse]);
              tmp_elements[idense] = sfv_elements[isparse];
            }

            // Synchronize the temporary DenseVector
            my_col_gate.sync_0(tmp);

            // Copy sparse filter vector contents to DenseVector
            for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
            {
              Index idense(slip_filter_vector.indices()[isparse]);
              tmp_elements[idense].normalize();
              sfv_elements[isparse] = tmp_elements[idense];

            }
          }
          // This happens if the non-empty MeshPart belonging to this filter is not present on this patch
          else
          {
            // Temporary DenseVector for syncing
            LocalVeloVector tmp(slip_filter_vector.size(), DataType_(0));
            my_col_gate.sync_0(tmp);
          }
        }
      }
    }; // class ScalarMixedSystemLevel<...>
  } // namespace Control
} // namespace FEAT
