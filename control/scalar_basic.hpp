#pragma once
#ifndef CONTROL_SCALAR_BASIC_HPP
#define CONTROL_SCALAR_BASIC_HPP 1

#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/geometry/domain_control.hpp>
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

namespace FEAST
{
  namespace Control
  {
    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      template<typename,typename,typename> class ScalarMatrix_ = LAFEM::SparseMatrixCSR>
    class ScalarBasicSystemLevel
    {
    public:
      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define local scalar matrix type
      typedef ScalarMatrix_<MemType, DataType, IndexType> LocalScalarMatrix;

      /// define local system matrix type
      typedef ScalarMatrix_<MemType, DataType, IndexType> LocalSystemMatrix;

      /// define local system vector type
      typedef typename LocalSystemMatrix::VectorTypeR LocalSystemVector;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType, DataType, IndexType> SystemMirror;

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

    public:
      ScalarBasicSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys)
      {
      }

      virtual ~ScalarBasicSystemLevel()
      {
      }
    }; // class ScalarBasicSystemLevel<...>

    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      template<typename, typename, typename> class ScalarMatrix_ = LAFEM::SparseMatrixCSR>
    class ScalarUnitFilterSystemLevel :
      public ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      /// define local filter type
      typedef LAFEM::UnitFilter<MemType_, DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter> GlobalSystemFilter;

      /// our global system filter
      GlobalSystemFilter filter_sys;

    public:
      ScalarUnitFilterSystemLevel() :
        BaseClass()
      {
      }
    }; // class ScalarUnitFilterSystemLevel<...>

    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      template<typename, typename, typename> class ScalarMatrix_ = LAFEM::SparseMatrixCSR>
    class ScalarMeanFilterSystemLevel :
      public ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      /// define local filter type
      typedef Global::MeanFilter<MemType_, DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter> GlobalSystemFilter;

      /// our global system filter
      GlobalSystemFilter filter_sys;

    public:
      ScalarMeanFilterSystemLevel() :
        BaseClass()
      {
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

      /// our two neighbour system levels
      SystemLevel_& level_coarse;
      SystemLevel_& level_fine;

      /// our global transfer matrices
      GlobalSystemTransferMatrix prol_sys, rest_sys;

    public:
      explicit ScalarBasicTransferLevel(SystemLevel_& lvl_coarse, SystemLevel_& lvl_fine) :
        level_coarse(lvl_coarse),
        level_fine(lvl_fine),
        prol_sys(&level_fine.gate_sys, &level_coarse.gate_sys),
        rest_sys(&level_coarse.gate_sys, &level_fine.gate_sys)
      {
      }

      virtual ~ScalarBasicTransferLevel()
      {
      }
    };

    template<typename Space_>
    class ScalarBasicAssemblerLevel
    {
    public:
      typedef Space_ SpaceType;
      typedef typename SpaceType::TrafoType TrafoType;
      typedef typename TrafoType::MeshType MeshType;
      typedef Geometry::DomainLevel<MeshType> DomainLevelType;
      typedef Geometry::DomainLayer<MeshType> DomainLayerType;

    public:
      DomainLevelType& domain_level;
      MeshType& mesh;
      TrafoType trafo;
      SpaceType space;
      Cubature::DynamicFactory cubature;

    public:
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
          Assembly::SymbolicMatrixAssembler<Assembly::Stencil::StandardRefinement>::assemble(
            loc_prol, this->space, level_coarse.space);
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
