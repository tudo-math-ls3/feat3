// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/global/muxer.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/domain_assembler_helpers.hpp>

#include <control/domain/domain_control.hpp>
#include <control/asm/gate_asm.hpp>
#include <control/asm/muxer_asm.hpp>
#include <control/asm/splitter_asm.hpp>
#include <control/asm/transfer_asm.hpp>
#include <control/asm/mean_filter_asm.hpp>
#include <control/asm/unit_filter_asm.hpp>
#include <control/asm/slip_filter_asm.hpp>

namespace FEAT
{
  namespace Control
  {
    template<
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename BlockedMatrix_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename TransferMatrix_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>
      >
    class BlockedBasicSystemLevel
    {
    public:
      static_assert(std::is_same<DataType_, typename BlockedMatrix_::DataType>::value, "DataType mismatch!");
      static_assert(std::is_same<IndexType_, typename BlockedMatrix_::IndexType>::value, "IndexType mismatch!");

      // basic types
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
      typedef LAFEM::VectorMirror<DataType, IndexType> SystemMirror;

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

      // no copies, no problems
      BlockedBasicSystemLevel(const BlockedBasicSystemLevel&) = delete;
      BlockedBasicSystemLevel& operator=(const BlockedBasicSystemLevel&) = delete;

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
        Asm::asm_gate(virt_dom_lvl, virt_dom_lvl->space, this->gate_sys, true);
      }

      template<typename DomainLevel_>
      void assemble_coarse_muxer(const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse)
      {
        Asm::asm_muxer(virt_lvl_coarse, [](const DomainLevel_& dl){return &dl.space;}, this->coarse_muxer_sys);
      }

      template<typename DomainLevel_>
      void assemble_transfer(
        const BlockedBasicSystemLevel& sys_lvl_coarse,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const String& cubature, bool trunc = false, bool shrink = true)
      {
        Asm::asm_transfer_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc, shrink,
          [](const DomainLevel_& dl) {return &dl.space;},
          this->transfer_sys.local(), this->coarse_muxer_sys, this->gate_sys, sys_lvl_coarse.gate_sys);

        this->transfer_sys.compile();
      }

      template<typename DomainLevel_>
      void assemble_transfer(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const String& cubature, bool trunc = false, bool shrink = true)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT(!virt_lvl_coarse.is_parent());

        Asm::asm_transfer_blocked(virt_lvl_fine, virt_lvl_coarse, cubature, trunc, shrink,
          [](const DomainLevel_& dl) {return &dl.space;},
          this->transfer_sys.local(), this->coarse_muxer_sys, this->gate_sys, this->gate_sys);

        this->transfer_sys.compile();
      }

      template<typename Trafo_, typename Space_>
      void assemble_laplace_matrix(Assembly::DomainAssembler<Trafo_>& dom_asm, const Space_& space, const String& cubature, const DataType nu = DataType(1))
      {
        // get local matrix
        auto& loc_matrix = this->matrix_sys.local();

        // assemble structure?
        if(loc_matrix.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_std1(loc_matrix, space);
        }

        // format and assemble Laplace
        loc_matrix.format();
        Assembly::Common::LaplaceOperatorBlocked<dim_> laplace_op;
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_matrix, laplace_op, space, cubature, nu);
      }
    }; // class BlockedBasicSystemLevel<...>


    template<
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename BlockedMatrix_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename TransferMatrix_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>>
    class BlockedUnitFilterSystemLevel :
      public BlockedBasicSystemLevel<dim_, DataType_, IndexType_, BlockedMatrix_, TransferMatrix_>
    {
    public:
      typedef BlockedBasicSystemLevel<dim_, DataType_, IndexType_, BlockedMatrix_, TransferMatrix_> BaseClass;

      // basic types
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<DataType, IndexType> SystemMirror;

      /// define local filter type
      typedef LAFEM::UnitFilterBlocked<DataType_, IndexType_, dim_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename DT2_, typename IT2_, typename BlockedMatrix2_>
      using BaseType = BlockedUnitFilterSystemLevel<dim_, DT2_, IT2_, BlockedMatrix2_>;

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

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const BlockedUnitFilterSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }

      template<typename DomainLevel_, typename Space_>
      void assemble_homogeneous_unit_filter(const DomainLevel_& dom_level, const Space_& space)
      {
        Asm::asm_unit_filter_blocked_homogeneous(this->filter_sys.local(), dom_level, space, "*");
      }
    }; // class BlockedUnitFilterSystemLevel<...>

    template<
      int dim_,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename BlockedMatrix_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
      typename TransferMatrix_ = LAFEM::SparseMatrixBWrappedCSR<DataType_, IndexType_, dim_>>
    class BlockedCombinedSystemLevel :
      public BlockedBasicSystemLevel<dim_, DataType_, IndexType_, BlockedMatrix_, TransferMatrix_>
    {
    public:
      typedef BlockedBasicSystemLevel<dim_, DataType_, IndexType_, BlockedMatrix_, TransferMatrix_> BaseClass;

      // basic types
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<DataType, IndexType> SystemMirror;

      /// define local filter types
      typedef LAFEM::UnitFilterBlocked<DataType_, IndexType_, dim_> LocalUnitFilter;
      typedef LAFEM::SlipFilter<DataType_, IndexType_, dim_> LocalSlipFilter;
      //typedef LAFEM::MeanFilterBlocked<DataType_, IndexType_, dim_> LocalMeanFilter;
      typedef LAFEM::FilterSequence<LocalUnitFilter> LocalUnitFilterSeq;
      typedef LAFEM::FilterSequence<LocalSlipFilter> LocalSlipFilterSeq;
      typedef LAFEM::FilterChain<LocalUnitFilterSeq, LocalSlipFilterSeq> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename DT2_, typename IT2_, typename BlockedMatrix2_>
      using BaseType = BlockedUnitFilterSystemLevel<dim_, DT2_, IT2_, BlockedMatrix2_>;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      BlockedCombinedSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return BaseClass::bytes () + filter_sys.bytes();
      }

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const BlockedCombinedSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }

      LocalUnitFilterSeq& get_local_unit_filter_seq()
      {
        return this->filter_sys.local().template at<0>();
      }

      LocalUnitFilter& get_local_unit_filter(const String& name)
      {
        return get_local_unit_filter_seq().find_or_add(name);
      }

      LocalSlipFilterSeq& get_local_slip_filter_seq()
      {
        return this->filter_sys.local().template at<1>();
      }

      LocalSlipFilter& get_local_slip_filter(const String& name)
      {
        return get_local_slip_filter_seq().find_or_add(name);
      }

      //LocalMeanFilter& get_local_mean_filter()
      //{
        //return this->filter_sys.local().template at<1>();
      //}

      template<typename DomainLevel_, typename Space_>
      void assemble_homogeneous_unit_filter(const DomainLevel_& dom_level, const Space_& space, const String& name, const String& mesh_parts)
      {
        auto& loc_filter = get_local_unit_filter(name);
        loc_filter.clear();
        Asm::asm_unit_filter_blocked_homogeneous(loc_filter, dom_level, space, mesh_parts);
      }

      template<typename DomainLevel_, typename Space_, typename Function_>
      void assemble_unit_filter(const DomainLevel_& dom_level, const Space_& space, const String& name, const String& mesh_parts, const Function_& function)
      {
        auto& loc_filter = get_local_unit_filter(name);
        loc_filter.clear();
        Asm::asm_unit_filter_blocked(loc_filter, dom_level, space, mesh_parts, function);
      }

      template<typename DomainLevel_, typename Space_>
      void assemble_slip_filter(const DomainLevel_& dom_level, const Space_& space, const String& name, const String& mesh_parts)
      {
        auto& loc_filter = get_local_slip_filter(name);
        loc_filter.clear();
        Asm::asm_slip_filter(loc_filter, dom_level, space, mesh_parts);
        Asm::sync_slip_filter(this->gate_sys, loc_filter);
      }

      //template<typename Space_>
      //void assemble_mean_filter(const Space_& space, const String& cubature)
      //{
      //  get_local_mean_filter() = Asm::asm_mean_filter(this->gate_sys, space, cubature);
      //}
    }; // class BlockedCombinedSystemLevel<...>
  } // namespace Control
} // namespace FEAT
