#pragma once
#ifndef CONTROL_SCALAR_BASIC_HPP
#define CONTROL_SCALAR_BASIC_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/comm_base.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/foundation_gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/mean_filter.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>

#include <control/domain/domain_control.hpp>

#include <deque>

namespace FEAST
{
  namespace Control
  {
    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarBasicSystemLevel
    {
    public:
      static_assert(std::is_same<MemType_, typename ScalarMatrix_::MemType>::value, "MemType missmatch!");
      static_assert(std::is_same<DataType_, typename ScalarMatrix_::DataType>::value, "DataType missmatch!");
      static_assert(std::is_same<IndexType_, typename ScalarMatrix_::IndexType>::value, "IndexType missmatch!");

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define local scalar matrix type
      typedef ScalarMatrix_ LocalScalarMatrix;

      /// define local system matrix type
      typedef ScalarMatrix_ LocalSystemMatrix;

      /// define local system vector type
      typedef typename LocalSystemMatrix::VectorTypeR LocalSystemVector;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define system gate
      typedef Global::FoundationGate<LocalSystemVector, SystemMirror> SystemGate;

      /// define global system vector type
      typedef Global::Vector<LocalSystemVector> GlobalSystemVector;

      /// define global system matrix type
      typedef Global::Matrix<LocalSystemMatrix> GlobalSystemMatrix;

      /* ***************************************************************************************** */

      /// our system gate
      SystemGate gate_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;

      /// CTOR
      ScalarBasicSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys)
      {
      }

      virtual ~ScalarBasicSystemLevel()
      {
      }

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarBasicSystemLevel<M_, D_, I_, SM_> & other)
      {
        gate_sys.convert(other.gate_sys);
        matrix_sys.convert(&gate_sys, &gate_sys, other.matrix_sys);
      }
    }; // class ScalarBasicSystemLevel<...>

    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarUnitFilterSystemLevel :
      public ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      /// define local filter type
      typedef LAFEM::UnitFilter<MemType_, DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter> GlobalSystemFilter;

      /// Our class base type
      template <typename Mem2_, typename DT2_, typename IT2_, typename ScalarMatrix2_>
      using BaseType = class ScalarUnitFilterSystemLevel<Mem2_, DT2_, IT2_, ScalarMatrix2_>;

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
        return this->matrix_sys.bytes () + filter_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarUnitFilterSystemLevel content as content of current ScalarUnitFilterSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarUnitFilterSystemLevel<M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }
    }; // class ScalarUnitFilterSystemLevel<...>

    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarMeanFilterSystemLevel :
      public ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      /// define local filter type
      typedef Global::MeanFilter<MemType_, DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter> GlobalSystemFilter;

      /// Our class base type
      template <typename Mem2_, typename DT2_, typename IT2_, typename ScalarMatrix2_>
      using BaseType = class ScalarMeanFilterSystemLevel<Mem2_, DT2_, IT2_, ScalarMatrix2_>;

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
        return this->matrix_sys.bytes () + filter_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarMeanFilterSystemLevel content as content of current ScalarMeanFilterSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarMeanFilterSystemLevel<M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }
    }; // class ScalarMeanFilterSystemLevel<...>

    template<typename SystemLevel_, typename ScalarMatrix_ = typename SystemLevel_::LocalScalarMatrix>
    class ScalarBasicTransferLevel
    {
    public:
      /// our local transfer matrix type
      typedef ScalarMatrix_ LocalSystemTransferMatrix;

      /// our global transfer matrix type
      typedef Global::Matrix<LocalSystemTransferMatrix> GlobalSystemTransferMatrix;

      /// Our class base type
      template <typename SystemLevel2_>
      using BaseType = class ScalarBasicTransferLevel<SystemLevel2_>;

      /// our global transfer matrices
      GlobalSystemTransferMatrix prol_sys;

      /// \copydoc ScalarBasicTransferLevel::prol_sys
      GlobalSystemTransferMatrix rest_sys;

      ScalarBasicTransferLevel()
      {
      }

      /// CTOR
      explicit ScalarBasicTransferLevel(SystemLevel_& lvl_coarse, SystemLevel_& lvl_fine) :
        prol_sys(&lvl_fine.gate_sys, &lvl_coarse.gate_sys),
        rest_sys(&lvl_coarse.gate_sys, &lvl_fine.gate_sys)
      {
      }

      virtual ~ScalarBasicTransferLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return prol_sys.bytes() + rest_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarBasicTransferLevel content as content of current ScalarBasicTransferLevel.
       *
       * \warning The provided SystemLevels must already be converted to the matching
       * configuration, as they contain the used gateways.
       *
       */
      template <typename SL_, typename SM_>
      void convert(SystemLevel_ & lvl_coarse , SystemLevel_ & lvl_fine, const ScalarBasicTransferLevel<SL_, SM_> & other)
      {
        prol_sys.convert(&lvl_fine.gate_sys, &lvl_coarse.gate_sys, other.prol_sys);
        rest_sys.convert(&lvl_coarse.gate_sys, &lvl_fine.gate_sys, other.rest_sys);
      }
    }; // class ScalarBasicTransferLevel<...>

    template<typename Space_>
    class ScalarBasicAssemblerLevel
    {
    public:
      typedef Space_ SpaceType;
      typedef typename SpaceType::TrafoType TrafoType;
      typedef typename TrafoType::MeshType MeshType;
      typedef Control::Domain::DomainLevel<MeshType> DomainLevelType;
      typedef Control::Domain::DomainLayer<MeshType> DomainLayerType;

      DomainLevelType& domain_level;
      MeshType& mesh;
      TrafoType trafo;
      SpaceType space;
      Cubature::DynamicFactory cubature;

      explicit ScalarBasicAssemblerLevel(DomainLevelType& dom_lvl) :
        domain_level(dom_lvl),
        mesh(domain_level.get_mesh()),
        trafo(mesh),
        space(trafo),
        cubature("auto-degree:" + stringify(Math::sqr(SpaceType::local_degree)+2))
      {
      }

      virtual ~ScalarBasicAssemblerLevel()
      {
      }

      //template<typename MemType_, typename DataType_, typename IndexType_, template<typename, typename, typename> class ScalarMatrix_>
      //void assemble_gates(const DomainLayerType& dom_layer, ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>& sys_level)
      template<typename SystemLevel_>
      void assemble_gates(const DomainLayerType& dom_layer, SystemLevel_& sys_level)
      {
        // get our gate
        typename SystemLevel_::SystemGate& gate_sys = sys_level.gate_sys;

        // loop over all ranks
        for(Index i(0); i < dom_layer.size(); ++i)
        {
          Index rank = dom_layer.get_rank(i);
          Index ctag = dom_layer.get_ctag(i);

          // try to find our halo
          auto* halo = domain_level.find_halo_part(rank);
          if (halo == nullptr)
            throw InternalError("ERROR: Halo not found!");

          // assemble the mirror
          typename SystemLevel_::SystemMirror mirror_sys;
          Assembly::MirrorAssembler::assemble_mirror(mirror_sys, space, *halo);

          // push mirror into gate
          gate_sys.push(rank,ctag, std::move(mirror_sys));
        }

        // create local template vector
        typename SystemLevel_::LocalSystemVector tmpl_s(space.get_num_dofs());

        // compile gate
        gate_sys.compile(std::move(tmpl_s));
      }

      template<typename TransferLevel_>
      void assemble_system_transfer(TransferLevel_& trans_level, ScalarBasicAssemblerLevel& level_coarse)
      {
        // get global transfer matrices
        typename TransferLevel_::GlobalSystemTransferMatrix& glob_prol = trans_level.prol_sys;
        typename TransferLevel_::GlobalSystemTransferMatrix& glob_rest = trans_level.rest_sys;

        // get local transfer matrices
        typename TransferLevel_::LocalSystemTransferMatrix& loc_prol = (*glob_prol);
        typename TransferLevel_::LocalSystemTransferMatrix& loc_rest = (*glob_rest);

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, this->space, level_coarse.space);
        }

        // create a global pressure weight vector
        auto glob_vec_weight = glob_prol.create_vector_l();

        // get local pressure weight vector
        auto& loc_vec_weight = (*glob_vec_weight);

        // assemble prolongation matrix
        {
          loc_prol.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight,
            this->space, level_coarse.space, this->cubature);

          // synchronise weight vector
          glob_vec_weight.sync_0();

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol.scale_rows(loc_prol, loc_vec_weight);

          // copy and transpose
          loc_rest = loc_prol.transpose();
        }
      }
    }; // class ScalarBasicAssemblerLevel<...>

  } // namespace Control
} // namespace FEAST

#endif // CONTROL_SCALAR_BASIC_HPP
