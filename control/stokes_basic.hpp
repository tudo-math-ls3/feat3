#pragma once
#ifndef CONTROL_STOKES_BASIC_HPP
#define CONTROL_STOKES_BASIC_HPP 1

#include <kernel/geometry/export_vtk.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/power_mirror.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/power_filter.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_diag_matrix.hpp>
#include <kernel/lafem/power_diag_matrix.hpp>
#include <kernel/lafem/power_col_matrix.hpp>
#include <kernel/lafem/power_row_matrix.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/global/foundation_gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/mean_filter.hpp>

#include <control/domain/domain_control.hpp>

namespace FEAST
{
  namespace Control
  {
    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    struct StokesBasicSystemLevel
    {
      static_assert(std::is_same<MemType_, typename ScalarMatrix_::MemType>::value, "MemType missmatch!");
      static_assert(std::is_same<DataType_, typename ScalarMatrix_::DataType>::value, "DataType missmatch!");
      static_assert(std::is_same<IndexType_, typename ScalarMatrix_::IndexType>::value, "IndexType missmatch!");

      // basic types
      static constexpr int dim = dim_;
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      // define local matrix types
      typedef ScalarMatrix_ LocalScalarMatrix;
      typedef LAFEM::PowerDiagMatrix<LocalScalarMatrix, dim> LocalMatrixBlockA;
      typedef LAFEM::PowerColMatrix<LocalScalarMatrix, dim> LocalMatrixBlockB;
      typedef LAFEM::PowerRowMatrix<LocalScalarMatrix, dim> LocalMatrixBlockD;
      typedef LAFEM::SaddlePointMatrix<LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockD> LocalSystemMatrix;

      // define local vector types
      typedef typename LocalScalarMatrix::VectorTypeR LocalScalarVector;
      typedef LAFEM::PowerVector<LocalScalarVector, dim> LocalVeloVector;
      typedef LocalScalarVector LocalPresVector;
      typedef LAFEM::TupleVector<LocalVeloVector, LocalPresVector> LocalSystemVector;

      // define mirror types
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> ScalarMirror;
      typedef LAFEM::PowerMirror<ScalarMirror, dim> VeloMirror;
      typedef ScalarMirror PresMirror;
      typedef LAFEM::TupleMirror<VeloMirror, PresMirror> SystemMirror;

      // define gates
      typedef Global::FoundationGate<LocalVeloVector, VeloMirror> VeloGate;
      typedef Global::FoundationGate<LocalPresVector, PresMirror> PresGate;
      typedef Global::FoundationGate<LocalSystemVector, SystemMirror> SystemGate;

      // define global vector types
      typedef Global::Vector<LocalVeloVector> GlobalVeloVector;
      typedef Global::Vector<LocalPresVector> GlobalPresVector;
      typedef Global::Vector<LocalSystemVector> GlobalSystemVector;

      // define global matrix types
      typedef Global::Matrix<LocalMatrixBlockA> GlobalMatrixBlockA;
      typedef Global::Matrix<LocalMatrixBlockB> GlobalMatrixBlockB;
      typedef Global::Matrix<LocalMatrixBlockD> GlobalMatrixBlockD;
      typedef Global::Matrix<LocalScalarMatrix> GlobalScalarMatrix;
      typedef Global::Matrix<LocalSystemMatrix> GlobalSystemMatrix;

      /* ***************************************************************************************** */

      // gates
      VeloGate gate_velo;
      PresGate gate_pres;
      SystemGate gate_sys;

      // (global) matrices
      GlobalSystemMatrix matrix_sys;
      GlobalMatrixBlockA matrix_a;
      GlobalMatrixBlockB matrix_b;
      GlobalMatrixBlockD matrix_d;
      GlobalScalarMatrix matrix_s;

      StokesBasicSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys),
        matrix_a(&gate_velo, &gate_velo),
        matrix_b(&gate_velo, &gate_pres),
        matrix_d(&gate_pres, &gate_velo),
        matrix_s(&gate_pres, &gate_pres)
      {
      }

      virtual ~StokesBasicSystemLevel()
      {
      }

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const StokesBasicSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        gate_velo.convert(other.gate_velo);
        gate_pres.convert(other.gate_pres);
        gate_sys.convert(other.gate_sys);

        matrix_a.convert(&gate_velo, &gate_velo, other.matrix_a);
        matrix_b.convert(&gate_velo, &gate_pres, other.matrix_b);
        matrix_d.convert(&gate_pres, &gate_velo, other.matrix_d);
        matrix_s.convert(&gate_pres, &gate_pres, other.matrix_s);

        (*matrix_sys).block_a() = (*matrix_a).clone(LAFEM::CloneMode::Shallow);
        (*matrix_sys).block_b() = (*matrix_b).clone(LAFEM::CloneMode::Shallow);
        (*matrix_sys).block_d() = (*matrix_d).clone(LAFEM::CloneMode::Shallow);
      }
    }; // struct StokesBasicSystemLevel<...>

    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    struct StokesUnitVeloNonePresSystemLevel :
      public StokesBasicSystemLevel<dim_, MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
      typedef StokesBasicSystemLevel<dim_, MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      // define local filter types
      typedef LAFEM::UnitFilter<MemType_, DataType_, IndexType_> UnitVeloFilter;
      typedef LAFEM::PowerFilter<UnitVeloFilter, dim_> LocalVeloFilter;
      typedef LAFEM::NoneFilter<MemType_, DataType_, IndexType_> LocalPresFilter;
      typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresFilter> LocalSystemFilter;

      // define global filter types
      typedef Global::Filter<LocalVeloFilter> GlobalVeloFilter;
      typedef Global::Filter<LocalPresFilter> GlobalPresFilter;
      typedef Global::Filter<LocalSystemFilter> GlobalSystemFilter;

      // (global) filters
      GlobalSystemFilter filter_sys;
      GlobalVeloFilter filter_velo;
      GlobalPresFilter filter_pres;

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->matrix_sys.bytes () + this->matrix_s.bytes() + filter_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source StokesUnitVeloNonePresSystemLevel content as content of current StokesUnitVeloNonePresSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const StokesUnitVeloNonePresSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_velo.convert(other.filter_velo);
        filter_pres.convert(other.filter_pres);

        (*filter_sys).template at<0>() = (*filter_velo).clone(LAFEM::CloneMode::Shallow);
        (*filter_sys).template at<1>() = (*filter_pres).clone(LAFEM::CloneMode::Shallow);
      }
    }; // struct StokesUnitVeloNonePresSystemLevel<...>


    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    struct StokesUnitVeloMeanPresSystemLevel :
      public StokesBasicSystemLevel<dim_, MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
      typedef StokesBasicSystemLevel<dim_, MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      // define local filter types
      typedef LAFEM::UnitFilter<MemType_, DataType_, IndexType_> UnitVeloFilter;
      typedef LAFEM::PowerFilter<UnitVeloFilter, dim_> LocalVeloFilter;
      typedef Global::MeanFilter<MemType_, DataType_, IndexType_> LocalPresFilter;
      //typedef LAFEM::NoneFilter<MemType_, DataType_, IndexType_> LocalPresFilter;
      typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresFilter> LocalSystemFilter;

      // define global filter types
      typedef Global::Filter<LocalVeloFilter> GlobalVeloFilter;
      typedef Global::Filter<LocalPresFilter> GlobalPresFilter;
      typedef Global::Filter<LocalSystemFilter> GlobalSystemFilter;

      // (global) filters
      GlobalSystemFilter filter_sys;
      GlobalVeloFilter filter_velo;
      GlobalPresFilter filter_pres;

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->matrix_sys.bytes () + this->matrix_s.bytes() + filter_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source StokesUnitVeloMeanPresSystemLevel content as content of current StokesUnitVeloMeanPresSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const StokesUnitVeloMeanPresSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_velo.convert(other.filter_velo);
        filter_pres.convert(other.filter_pres);

        (*filter_sys).template at<0>() = (*filter_velo).clone(LAFEM::CloneMode::Shallow);
        (*filter_sys).template at<1>() = (*filter_pres).clone(LAFEM::CloneMode::Shallow);
      }
    }; // struct StokesUnitVeloNonePresSystemLevel<...>

    template<typename SystemLevel_>
    class StokesBasicTransferLevel
    {
    public:
      typedef SystemLevel_ SystemLevel;
      static constexpr int dim = SystemLevel_::dim;

      typedef typename SystemLevel_::LocalScalarVector LocalScalarVector;
      typedef typename SystemLevel_::LocalScalarMatrix LocalScalarMatrix;
      typedef LAFEM::PowerDiagMatrix<LocalScalarMatrix, dim> LocalVeloTransferMatrix;
      typedef LocalScalarMatrix LocalPresTransferMatrix;
      typedef LAFEM::TupleDiagMatrix<LocalVeloTransferMatrix, LocalPresTransferMatrix> LocalSystemTransferMatrix;

      typedef Global::Matrix<LocalVeloTransferMatrix> GlobalVeloTransferMatrix;
      typedef Global::Matrix<LocalPresTransferMatrix> GlobalPresTransferMatrix;
      typedef Global::Matrix<LocalSystemTransferMatrix> GlobalSystemTransferMatrix;

      // (global) transfer matrices
      GlobalVeloTransferMatrix prol_velo, rest_velo;
      GlobalPresTransferMatrix prol_pres, rest_pres;
      GlobalSystemTransferMatrix prol_sys, rest_sys;

      explicit StokesBasicTransferLevel()
      {
      }

      explicit StokesBasicTransferLevel(SystemLevel_& lvl_coarse, SystemLevel_& lvl_fine) :
        prol_velo(&lvl_fine.gate_velo, &lvl_coarse.gate_velo),
        rest_velo(&lvl_coarse.gate_velo, &lvl_fine.gate_velo),
        prol_pres(&lvl_fine.gate_pres, &lvl_coarse.gate_pres),
        rest_pres(&lvl_coarse.gate_pres, &lvl_fine.gate_pres),
        prol_sys(&lvl_fine.gate_sys, &lvl_coarse.gate_sys),
        rest_sys(&lvl_coarse.gate_sys, &lvl_fine.gate_sys)
      {
      }

      virtual ~StokesBasicTransferLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return prol_sys.bytes () + rest_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source StokesBasicTransferLevel content as content of current StokesBasicTransferLevel.
       *
       */
      template <typename SL_>
      void convert(SystemLevel_ & lvl_coarse, SystemLevel_ & lvl_fine, const StokesBasicTransferLevel<SL_> & other)
      {
        prol_velo.convert(&lvl_fine.gate_velo, &lvl_coarse.gate_velo, other.prol_velo);
        rest_velo.convert(&lvl_coarse.gate_velo, &lvl_fine.gate_velo, other.rest_velo);
        prol_pres.convert(&lvl_fine.gate_pres, &lvl_coarse.gate_pres, other.prol_pres);
        rest_pres.convert(&lvl_coarse.gate_pres, &lvl_fine.gate_pres, other.rest_pres);

        (*prol_sys).template at<0,0>() = (*prol_velo).clone(LAFEM::CloneMode::Shallow);
        (*prol_sys).template at<1,1>() = (*prol_pres).clone(LAFEM::CloneMode::Shallow);
        (*rest_sys).template at<0,0>() = (*rest_velo).clone(LAFEM::CloneMode::Shallow);
        (*rest_sys).template at<1,1>() = (*rest_pres).clone(LAFEM::CloneMode::Shallow);
      }
    }; // struct StokesBasicTransferLevel<...>

    template<
      typename SpaceVelo_,
      typename SpacePres_>
    class StokesBasicAssemblerLevel
    {
    public:
      typedef SpaceVelo_ SpaceVeloType;
      typedef SpacePres_ SpacePresType;
      typedef typename SpaceVelo_::TrafoType TrafoType;
      typedef typename TrafoType::MeshType MeshType;
      typedef Control::Domain::DomainLevel<MeshType> DomainLevelType;
      typedef Control::Domain::DomainLayer<MeshType> DomainLayerType;

    public:
      DomainLevelType& domain_level;
      MeshType& mesh;
      TrafoType trafo;
      SpaceVeloType space_velo;
      SpacePresType space_pres;
      Cubature::DynamicFactory cubature;

    public:
      explicit StokesBasicAssemblerLevel(DomainLevelType& dom_lvl) :
        domain_level(dom_lvl),
        mesh(domain_level.get_mesh()),
        trafo(mesh),
        space_velo(trafo),
        space_pres(trafo),
        cubature("auto-degree:" + stringify(Math::sqr(SpaceVeloType::local_degree)+2))
      {
      }

      virtual ~StokesBasicAssemblerLevel()
      {
      }

      template<typename SystemLevel_>
      void assemble_gates(const DomainLayerType& dom_layer, SystemLevel_& sys_level)
      {
        // create our gates
        typename SystemLevel_::VeloGate& gate_velo = sys_level.gate_velo;
        typename SystemLevel_::PresGate& gate_pres = sys_level.gate_pres;
        typename SystemLevel_::SystemGate& gate_sys  = sys_level.gate_sys;

        // loop over all ranks
        for(Index i(0); i < dom_layer.size(); ++i)
        {
          Index rank = dom_layer.get_rank(i);
          Index ctag = dom_layer.get_ctag(i);

          // try to find our halo
          auto* halo = domain_level.find_halo_part(rank);
          if(halo == nullptr)
            throw InternalError("ERROR: Halo not found!");

          // assemble the velocity components mirror
          typename SystemLevel_::ScalarMirror mirror_velo_comp;
          Assembly::MirrorAssembler::assemble_mirror(mirror_velo_comp, this->space_velo, *halo);

          // create a velocity mirror
          typename SystemLevel_::VeloMirror mirror_velo(mirror_velo_comp.clone());

          // create (empty) pressure mirror
          typename SystemLevel_::PresMirror mirror_pres;
          //Assembly::MirrorAssembler::assemble_mirror(mirror_pres, this->space_pres, *halo);

          // create a system mirror
          typename SystemLevel_::SystemMirror mirror_sys(mirror_velo.clone(), mirror_pres.clone());

          // push mirror into gates
          gate_velo.push(rank, ctag, std::move(mirror_velo));
          //gate_pres.push(rank, ctag, std::move(mirror_pres));
          gate_sys.push(rank, ctag, std::move(mirror_sys));
        }

        // create local template vectors
        typename SystemLevel_::LocalVeloVector tmpl_v(space_velo.get_num_dofs());
        typename SystemLevel_::LocalPresVector tmpl_p(space_pres.get_num_dofs());
        typename SystemLevel_::LocalSystemVector tmpl_s(tmpl_v.clone(), tmpl_p.clone());

        // compile gates
        gate_velo.compile(std::move(tmpl_v));
        gate_pres.compile(std::move(tmpl_p));
        gate_sys.compile(std::move(tmpl_s));
      }

      template<typename SystemLevel_>
      void assemble_velocity_matrices(SystemLevel_& sys_level)
      {
        // get the global matrix A
        typename SystemLevel_::GlobalMatrixBlockA& mat_glob_a = sys_level.matrix_a;

        // get the local matrix A
        typename SystemLevel_::LocalMatrixBlockA& mat_loc_a = *mat_glob_a;

        // get the diagonal blocks
        typename SystemLevel_::LocalScalarMatrix& mat_loc_a1 = mat_loc_a.template at<0,0>();
        typename SystemLevel_::LocalScalarMatrix& mat_loc_a2 = mat_loc_a.template at<1,1>();

        // assemble matrix structure?
        if(mat_loc_a1.empty())
        {
          Assembly::SymbolicMatrixAssembler<>::assemble1(mat_loc_a1, space_velo);
          mat_loc_a2 = mat_loc_a1.clone(LAFEM::CloneMode::Layout);
        }

        // assemble velocity laplace matrix
        {
          mat_loc_a1.format();
          Assembly::Common::LaplaceOperator laplace_op;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_loc_a1, laplace_op, space_velo, cubature);
        }

        // copy into A2
        mat_loc_a2.copy(mat_loc_a1);
      }

      template<typename SystemLevel_>
      void assemble_grad_div_matrices(SystemLevel_& sys_level)
      {
        typedef typename SystemLevel_::DataType DataType;

        // get the global matrices B and D
        typename SystemLevel_::GlobalMatrixBlockB& mat_glob_b = sys_level.matrix_b;
        typename SystemLevel_::GlobalMatrixBlockD& mat_glob_d = sys_level.matrix_d;

        // get the local matrix B and D
        typename SystemLevel_::LocalMatrixBlockB& mat_loc_b = *mat_glob_b;
        typename SystemLevel_::LocalMatrixBlockD& mat_loc_d = *mat_glob_d;

        // get the matrix blocks
        typename SystemLevel_::LocalScalarMatrix& mat_loc_b1 = mat_loc_b.template at<0,0>();
        typename SystemLevel_::LocalScalarMatrix& mat_loc_b2 = mat_loc_b.template at<1,0>();
        typename SystemLevel_::LocalScalarMatrix& mat_loc_d1 = mat_loc_d.template at<0,0>();
        typename SystemLevel_::LocalScalarMatrix& mat_loc_d2 = mat_loc_d.template at<0,1>();

        // assemble matrix structure?
        if(mat_loc_b1.empty())
        {
          Assembly::SymbolicMatrixAssembler<>::assemble2(mat_loc_b1, space_velo, space_pres);
          mat_loc_b2 = mat_loc_b1.clone(LAFEM::CloneMode::Layout);
        }

        // assemble pressure gradient matrices
        {
          mat_loc_b1.format();
          mat_loc_b2.format();
          Assembly::Common::TestDerivativeOperator<0> der_x;
          Assembly::Common::TestDerivativeOperator<1> der_y;
          Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_loc_b1, der_x, space_velo, space_pres, cubature, -DataType(1));
          Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_loc_b2, der_y, space_velo, space_pres, cubature, -DataType(1));
        }

        // assemble velocity divergence matrices
        {
          // "assemble" by transposing B
          mat_loc_d1 = mat_loc_b1.transpose();
          mat_loc_d2 = mat_loc_b2.transpose();
        }
      }

      template<typename SystemLevel_>
      void assemble_schur_matrix(SystemLevel_& sys_level)
      {
        typedef typename SystemLevel_::DataType DataType;

        // get the global S matrix
        typename SystemLevel_::GlobalScalarMatrix& mat_glob_s = sys_level.matrix_s;

        // get the local matrix S
        typename SystemLevel_::LocalScalarMatrix& mat_loc_s = *mat_glob_s;

        // assemble matrix structure?
        if(mat_loc_s.empty())
        {
          Assembly::SymbolicMatrixAssembler<>::assemble1(mat_loc_s, space_pres);
        }

        // assemble schur matrix
        {
          mat_loc_s.format();
          Assembly::Common::IdentityOperator id_op;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_loc_s, id_op, space_pres, cubature, -DataType(1));
        }
      }

      template<typename SystemLevel_>
      void assemble_system_matrix(SystemLevel_& sys_level)
      {
        // assemble velocity matrices
        assemble_velocity_matrices(sys_level);
        assemble_grad_div_matrices(sys_level);

        // clone content into our global system matrix
        (*sys_level.matrix_sys).block_a() = (*sys_level.matrix_a).clone(LAFEM::CloneMode::Shallow);
        (*sys_level.matrix_sys).block_b() = (*sys_level.matrix_b).clone(LAFEM::CloneMode::Shallow);
        (*sys_level.matrix_sys).block_d() = (*sys_level.matrix_d).clone(LAFEM::CloneMode::Shallow);
      }

      template<typename TransferLevel_>
      void assemble_velocity_transfer(TransferLevel_& trans_level, StokesBasicAssemblerLevel& level_coarse)
      {
        // get global transfer matrices
        typename TransferLevel_::GlobalVeloTransferMatrix& glob_prol_v = trans_level.prol_velo;
        typename TransferLevel_::GlobalVeloTransferMatrix& glob_rest_v = trans_level.rest_velo;

        // get local transfer matrices
        typename TransferLevel_::LocalVeloTransferMatrix& loc_prol_v = (*glob_prol_v);
        typename TransferLevel_::LocalVeloTransferMatrix& loc_rest_v = (*glob_rest_v);

        // get the matrix blocks
        typename TransferLevel_::LocalScalarMatrix& loc_prol_vx = loc_prol_v.template at<0,0>();
        typename TransferLevel_::LocalScalarMatrix& loc_prol_vy = loc_prol_v.template at<1,1>();
        typename TransferLevel_::LocalScalarMatrix& loc_rest_vx = loc_rest_v.template at<0,0>();
        typename TransferLevel_::LocalScalarMatrix& loc_rest_vy = loc_rest_v.template at<1,1>();

        // assemble structure?
        if(loc_prol_vx.empty())
        {
          Assembly::SymbolicMatrixAssembler<Assembly::Stencil::StandardRefinement>::assemble(
            loc_prol_vx, this->space_velo, level_coarse.space_velo);

          loc_prol_vy = loc_prol_vx.clone(LAFEM::CloneMode::Layout);
        }

        // create a global velocity weight vector
        auto glob_vec_weight = glob_prol_v.create_vector_l();

        // get local velocity weight vector
        auto& loc_vec_weight = (*glob_vec_weight);

        // get local weight vector components
        auto& loc_vec_wx = loc_vec_weight.template at<0>();
        auto& loc_vec_wy = loc_vec_weight.template at<1>();

        // assemble prolongation matrix
        {
          loc_prol_vx.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol_vx, loc_vec_wx,
            this->space_velo, level_coarse.space_velo, this->cubature);

          // synchronise weight vector
          loc_vec_wy.copy(loc_vec_wx);
          glob_vec_weight.sync_0();

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol_vx.scale_rows(loc_prol_vx, loc_vec_wx);

          // copy and transpose
          loc_prol_vy = loc_prol_vx.clone(LAFEM::CloneMode::Shallow);
          loc_rest_vx = loc_prol_vx.transpose();
          loc_rest_vy = loc_rest_vx.clone(LAFEM::CloneMode::Shallow);
        }
      }

      template<typename TransferLevel_>
      void assemble_pressure_transfer(TransferLevel_& trans_level, StokesBasicAssemblerLevel& level_coarse)
      {
        // get global transfer matrices
        typename TransferLevel_::GlobalPresTransferMatrix& glob_prol_p = trans_level.prol_pres;
        typename TransferLevel_::GlobalPresTransferMatrix& glob_rest_p = trans_level.rest_pres;

        // get local transfer matrices
        typename TransferLevel_::LocalPresTransferMatrix& loc_prol_p = (*glob_prol_p);
        typename TransferLevel_::LocalPresTransferMatrix& loc_rest_p = (*glob_rest_p);

        // assemble structure?
        if(loc_prol_p.empty())
        {
          Assembly::SymbolicMatrixAssembler<Assembly::Stencil::StandardRefinement>::assemble(
            loc_prol_p, this->space_pres, level_coarse.space_pres);
        }

        // create a global pressure weight vector
        auto glob_vec_weight = glob_prol_p.create_vector_l();

        // get local pressure weight vector
        auto& loc_vec_weight = (*glob_vec_weight);

        // assemble prolongation matrix
        {
          loc_prol_p.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol_p, loc_vec_weight,
            this->space_pres, level_coarse.space_pres, this->cubature);

          // synchronise weight vector
          glob_vec_weight.sync_0();

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol_p.scale_rows(loc_prol_p, loc_vec_weight);

          // copy and transpose
          loc_rest_p = loc_prol_p.transpose();
        }
      }

      template<typename TransferLevel_>
      void assemble_system_transfer(TransferLevel_& trans_level, StokesBasicAssemblerLevel& level_coarse)
      {
        assemble_velocity_transfer(trans_level, level_coarse);
        assemble_pressure_transfer(trans_level, level_coarse);

        // clone content into our global transfer matrix
        (*trans_level.prol_sys).template at<0,0>() = (*trans_level.prol_velo).clone(LAFEM::CloneMode::Shallow);
        (*trans_level.prol_sys).template at<1,1>() = (*trans_level.prol_pres).clone(LAFEM::CloneMode::Shallow);
        (*trans_level.rest_sys).template at<0,0>() = (*trans_level.rest_velo).clone(LAFEM::CloneMode::Shallow);
        (*trans_level.rest_sys).template at<1,1>() = (*trans_level.rest_pres).clone(LAFEM::CloneMode::Shallow);
      }

      template<typename SolVector_>
      void write_vtk(const String& vtk_name, const SolVector_& vector, int rank, int nprocs) const
      {
        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(this->mesh);

        // project velocity and pressure
        LAFEM::DenseVector<Mem::Main, double, Index> vtx_vx, vtx_vy, vtx_p;
        Assembly::DiscreteVertexProjector::project(vtx_vx, vector.template at<0>().template at<0>(), this->space_velo);
        Assembly::DiscreteVertexProjector::project(vtx_vy, vector.template at<0>().template at<1>(), this->space_velo);
        Assembly::DiscreteCellProjector::project(vtx_p, vector.template at<1>(), this->space_pres, this->cubature);

        // write velocity
        exporter.add_field_vertex("velocity", vtx_vx.elements(), vtx_vy.elements());;

        // write pressure
        exporter.add_scalar_cell("pressure", vtx_p.elements());

        // finally, write the VTK file
        exporter.write(vtk_name, rank, nprocs);
      }
    };
  } // namespace Control
} // namespace FEAST

#endif // CONTROL_STOKES_BASIC_HPP
