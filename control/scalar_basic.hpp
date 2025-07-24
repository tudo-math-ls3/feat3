// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/muxer.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/global/splitter.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/mean_filter.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/voxel_assembly/poisson_assembler.hpp>

#include <control/domain/domain_control.hpp>
#include <control/asm/gate_asm.hpp>
#include <control/asm/muxer_asm.hpp>
#include <control/asm/splitter_asm.hpp>
#include <control/asm/transfer_asm.hpp>
#include <control/asm/transfer_voxel_asm.hpp>
#include <control/asm/mean_filter_asm.hpp>
#include <control/asm/unit_filter_asm.hpp>

#include <deque>

namespace FEAT
{
  namespace Control
  {
    template<
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>,
      typename TransferMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_> >
    class ScalarBasicSystemLevel
    {
    public:
      static_assert(std::is_same<DataType_, typename ScalarMatrix_::DataType>::value, "DataType mismatch!");
      static_assert(std::is_same<IndexType_, typename ScalarMatrix_::IndexType>::value, "IndexType mismatch!");

      // basic types
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define local scalar matrix type
      typedef ScalarMatrix_ LocalScalarMatrix;

      /// define local system matrix type
      typedef ScalarMatrix_ LocalSystemMatrix;

      /// define local transfer matrix type
      typedef TransferMatrix_ LocalSystemTransferMatrix;

      /// define local system vector type
      typedef typename LocalSystemMatrix::VectorTypeR LocalSystemVector;

      /// define local system transfer operator type
      typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

      /// define system mirror type
      typedef LAFEM::VectorMirror<DataType, IndexType> SystemMirror;

      /// define system gate
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;

      /// define system muxer
      typedef Global::Muxer<LocalSystemVector, SystemMirror> SystemMuxer;

      /// define system splitter
      typedef Global::Splitter<LocalSystemVector, SystemMirror> SystemSplitter;

      /// define global system vector type
      typedef Global::Vector<LocalSystemVector, SystemMirror> GlobalSystemVector;

      /// define global system matrix type
      typedef Global::Matrix<LocalSystemMatrix, SystemMirror, SystemMirror> GlobalSystemMatrix;

      /// define global system transfer operator type
      typedef Global::Transfer<LocalSystemTransfer, SystemMirror> GlobalSystemTransfer;

      /* ***************************************************************************************** */

      /// our system gate
      SystemGate gate_sys;

      /// our coarse-level system muxer
      SystemMuxer coarse_muxer_sys;

      /// our base-mesh multiplexer
      SystemSplitter base_splitter_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;

      /// our global transfer operator
      GlobalSystemTransfer transfer_sys;

      /// CTOR
      ScalarBasicSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys),
        transfer_sys(&coarse_muxer_sys)
      {
      }

      // no copies, no problems
      ScalarBasicSystemLevel(const ScalarBasicSystemLevel&) = delete;
      ScalarBasicSystemLevel& operator=(const ScalarBasicSystemLevel&) = delete;

      virtual ~ScalarBasicSystemLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->gate_sys.bytes() + this->coarse_muxer_sys.bytes() + this->base_splitter_sys.bytes()
          + this->transfer_sys.bytes() + this->matrix_sys.local().bytes();
      }

      template<typename D_, typename I_, typename SM_>
      void convert(const ScalarBasicSystemLevel<D_, I_, SM_> & other)
      {
        gate_sys.convert(other.gate_sys);
        coarse_muxer_sys.convert(other.coarse_muxer_sys);
        base_splitter_sys.convert(other.base_splitter_sys);
        matrix_sys.convert(&gate_sys, &gate_sys, other.matrix_sys);
        transfer_sys.convert(&coarse_muxer_sys, other.transfer_sys);
      }

      template<typename DomainLevel_, typename SpaceType_>
      void assemble_gate(const Domain::VirtualLevel<DomainLevel_>& virt_dom_lvl, const SpaceType_& space)
      {
        Asm::asm_gate(virt_dom_lvl, space, this->gate_sys, true);
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
      void assemble_base_splitter(const Domain::VirtualLevel<DomainLevel_>& virt_lvl)
      {
        Asm::asm_splitter(virt_lvl, [](const DomainLevel_& dl){return &dl.space;}, this->base_splitter_sys);
      }

      template<typename DomainLevel_>
      void assemble_transfer(
        const ScalarBasicSystemLevel& sys_lvl_coarse,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const String& cubature, bool trunc = false, bool shrink = true)
      {
        Asm::asm_transfer_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc, shrink,
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

        Asm::asm_transfer_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc, shrink,
          [](const DomainLevel_& dl) {return &dl.space;},
          this->transfer_sys.local(), this->coarse_muxer_sys, this->gate_sys, this->gate_sys);

        this->transfer_sys.compile();
      }

      template<typename DomainLevel_>
      void assemble_transfer_voxel(
        const ScalarBasicSystemLevel& sys_lvl_coarse,
        const Domain::VirtualLevel<Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_fine,
        const Domain::VirtualLevel<Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_coarse,
        const String& cubature, bool trunc = false, bool shrink = true)
      {
        Asm::asm_transfer_voxel_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc, shrink,
          [](const DomainLevel_& dl) {return &dl.space;},
          this->transfer_sys.local(), this->coarse_muxer_sys, this->gate_sys, sys_lvl_coarse.gate_sys);

        this->transfer_sys.compile();
      }

      template<typename DomainLevel_>
      void assemble_transfer_voxel(
        const Domain::VirtualLevel<Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_fine,
        const Domain::VirtualLevel<Domain::VoxelDomainLevelWrapper<DomainLevel_>>& virt_lvl_coarse,
        const String& cubature, bool trunc = false, bool shrink = true)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT(!virt_lvl_coarse.is_parent());

        Asm::asm_transfer_voxel_scalar(virt_lvl_fine, virt_lvl_coarse, cubature, trunc, shrink,
          [](const DomainLevel_& dl) {return &dl.space;},
          this->transfer_sys.local(), this->coarse_muxer_sys, this->gate_sys, this->gate_sys);

        this->transfer_sys.compile();
      }

      template<typename Space_, typename Cubature_>
      void assemble_laplace_matrix(const Space_& space, const Cubature_& cubature, const DataType nu = DataType(1))
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
        Assembly::Common::LaplaceOperator laplace_op;
        Assembly::BilinearOperatorAssembler::assemble_matrix1(loc_matrix, laplace_op, space, cubature, nu);
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
        Assembly::Common::LaplaceOperator laplace_op;
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_matrix, laplace_op, space, cubature, nu);
      }

      template<typename Space_>
      void assemble_laplace_voxel_based(const Adjacency::Coloring& coloring, const Space_& space, const String& cubature, const DataType nu = DataType(1))
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
        VoxelAssembly::VoxelPoissonAssembler<Space_, DataType, IndexType> voxel_assembler(space, coloring);
        voxel_assembler.assemble_matrix1(loc_matrix, space, Cubature::DynamicFactory(cubature), nu);
      }

      template<typename Space_>
      void symbolic_assembly_std1(const Space_& space)
      {
        // get local matrix
        auto& loc_matrix = this->matrix_sys.local();

        ASSERTM(loc_matrix.empty(), "Called symbolic assembly on non empty matrix");
        // assemble structure?
        Assembly::SymbolicAssembler::assemble_matrix_std1(loc_matrix, space);

      }
    }; // class ScalarBasicSystemLevel<...>

    template<
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_> >
    class ScalarUnitFilterSystemLevel :
      public ScalarBasicSystemLevel<DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<DataType_, IndexType_, ScalarMatrix_> BaseClass;

      // basic types
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<DataType, IndexType> SystemMirror;

      /// define local filter type
      typedef LAFEM::UnitFilter<DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename DT2_, typename IT2_, typename ScalarMatrix2_>
      using BaseType = ScalarUnitFilterSystemLevel<DT2_, IT2_, ScalarMatrix2_>;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      ScalarUnitFilterSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return BaseClass::bytes() + this->filter_sys.local().bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarUnitFilterSystemLevel content as content of current ScalarUnitFilterSystemLevel.
       *
       */
      template<typename D_, typename I_, typename SM_>
      void convert(const ScalarUnitFilterSystemLevel<D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }

      template<typename DomainLevel_, typename Space_>
      void assemble_homogeneous_unit_filter(const DomainLevel_& dom_level, const Space_& space)
      {
        Asm::asm_unit_filter_scalar_homogeneous(this->filter_sys.local(), dom_level, space, "*");
      }

    }; // class ScalarUnitFilterSystemLevel<...>

    template<
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_> >
    class ScalarMeanFilterSystemLevel :
      public ScalarBasicSystemLevel<DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<DataType_, IndexType_, ScalarMatrix_> BaseClass;

      // basic types
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<DataType, IndexType> SystemMirror;

      /// define local filter type
      typedef Global::MeanFilter<DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename DT2_, typename IT2_, typename ScalarMatrix2_>
      using BaseType = ScalarMeanFilterSystemLevel<DT2_, IT2_, ScalarMatrix2_>;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      ScalarMeanFilterSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return BaseClass::bytes() + this->filter_sys.local().bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarMeanFilterSystemLevel content as content of current ScalarMeanFilterSystemLevel.
       *
       */
      template<typename D_, typename I_, typename SM_>
      void convert(const ScalarMeanFilterSystemLevel<D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }

      template<typename Space_, typename Cubature_>
      void assemble_mean_filter(const Space_& space, const Cubature_& cubature)
      {
        this->filter_sys.local() = Asm::asm_mean_filter(this->gate_sys, space, cubature);
      }
    }; // class ScalarMeanFilterSystemLevel<...>

    template<
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_> >
    class ScalarCombinedSystemLevel :
      public ScalarBasicSystemLevel<DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<DataType_, IndexType_, ScalarMatrix_> BaseClass;

      // basic types
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<DataType, IndexType> SystemMirror;

      /// define local filter types
      typedef LAFEM::UnitFilter<DataType, IndexType> LocalUnitFilter;
      typedef Global::MeanFilter<DataType, IndexType> LocalMeanFilter;
      typedef LAFEM::FilterSequence<LocalUnitFilter> LocalUnitFilterSeq;
      typedef LAFEM::FilterChain<LocalUnitFilterSeq, LocalMeanFilter> LocalSystemFilter;

      /// define global filter types
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename DT2_, typename IT2_, typename ScalarMatrix2_>
      using BaseType = ScalarCombinedSystemLevel<DT2_, IT2_, ScalarMatrix2_>;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      ScalarCombinedSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return BaseClass::bytes() + this->filter_sys.local().bytes();
      }

      template<typename D_, typename I_, typename SM_>
      void convert(const ScalarMeanFilterSystemLevel<D_, I_, SM_> & other)
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

      LocalMeanFilter& get_local_mean_filter()
      {
        return this->filter_sys.local().template at<1>();
      }

      template<typename DomainLevel_, typename Space_>
      void assemble_unit_filter(const DomainLevel_& dom_level, const Space_& space, const String& name, const String& mesh_part_names)
      {
        auto& loc_filter = get_local_unit_filter(name);
        loc_filter.clear();
        Asm::asm_unit_filter_scalar_homogeneous(loc_filter, dom_level, space, mesh_part_names);
      }

      template<typename DomainLevel_, typename Space_, typename Function_>
      void assemble_unit_filter(const DomainLevel_& dom_level, const Space_& space, const String& name, const String& mesh_part_names, const Function_& function)
      {
        auto& loc_filter = get_local_unit_filter(name);
        loc_filter.clear();
        Asm::asm_unit_filter_scalar(loc_filter, dom_level, space, mesh_part_names, function);
      }

      template<typename Space_>
      void assemble_mean_filter(const Space_& space, const String& cubature)
      {
        get_local_mean_filter() = Asm::asm_mean_filter(this->gate_sys, space, cubature);
      }
    }; // class ScalarCombinedSystemLevel<...>
  } // namespace Control
} // namespace FEAT
