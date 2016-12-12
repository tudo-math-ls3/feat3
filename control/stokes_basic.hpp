#pragma once
#ifndef CONTROL_STOKES_BASIC_HPP
#define CONTROL_STOKES_BASIC_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
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
#include <kernel/lafem/transfer.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/discrete_projector.hpp>
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

namespace FEAT
{
  namespace Control
  {
    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    struct StokesBasicSystemLevel
    {
      static_assert(std::is_same<MemType_, typename ScalarMatrix_::MemType>::value, "MemType mismatch!");
      static_assert(std::is_same<DataType_, typename ScalarMatrix_::DataType>::value, "DataType mismatch!");
      static_assert(std::is_same<IndexType_, typename ScalarMatrix_::IndexType>::value, "IndexType mismatch!");

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

      // define local transfer matrix types
      typedef TransferMatrix_ LocalScalarTransferMatrix;
      typedef LAFEM::PowerDiagMatrix<TransferMatrix_, dim_> LocalVeloTransferMatrix;
      typedef TransferMatrix_ LocalPresTransferMatrix;
      typedef LAFEM::TupleDiagMatrix<LocalVeloTransferMatrix, LocalPresTransferMatrix> LocalSystemTransferMatrix;

      // define local transfer types
      typedef LAFEM::Transfer<LocalVeloTransferMatrix> LocalVeloTransfer;
      typedef LAFEM::Transfer<LocalPresTransferMatrix> LocalPresTransfer;
      typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

      // define mirror types
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> ScalarMirror;
      typedef LAFEM::PowerMirror<ScalarMirror, dim> VeloMirror;
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

      // gates
      VeloGate gate_velo;
      PresGate gate_pres;
      SystemGate gate_sys;

      /// our coarse-level system muxer
      VeloMuxer coarse_muxer_velo;
      PresMuxer coarse_muxer_pres;
      SystemMuxer coarse_muxer_sys;

      // (global) matrices
      GlobalSystemMatrix matrix_sys;
      GlobalMatrixBlockA matrix_a;
      GlobalMatrixBlockB matrix_b;
      GlobalMatrixBlockD matrix_d;
      GlobalSchurMatrix matrix_s;

      /// our global transfer operator
      GlobalVeloTransfer transfer_velo;
      GlobalPresTransfer transfer_pres;
      GlobalSystemTransfer transfer_sys;

      StokesBasicSystemLevel() :
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

      virtual ~StokesBasicSystemLevel()
      {
      }

      void compile_system_transfer()
      {
        // clone content into our global transfer matrix
        transfer_sys.get_mat_prol().template at<0,0>().clone(transfer_velo.get_mat_prol(), LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<0,0>().clone(transfer_velo.get_mat_rest(), LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_prol().template at<1,1>().clone(transfer_pres.get_mat_prol(), LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<1,1>().clone(transfer_pres.get_mat_rest(), LAFEM::CloneMode::Shallow);
      }

      void compile_system_matrix()
      {
        (*matrix_sys).block_a() = (*matrix_a).clone(LAFEM::CloneMode::Shallow);
        (*matrix_sys).block_b() = (*matrix_b).clone(LAFEM::CloneMode::Shallow);
        (*matrix_sys).block_d() = (*matrix_d).clone(LAFEM::CloneMode::Shallow);
      }

      template<typename M_, typename D_, typename I_, typename SM_, typename TM_>
      void convert(const StokesBasicSystemLevel<dim_, M_, D_, I_, SM_, TM_> & other)
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

        compile_system_matrix();

        // clone content into our global transfer matrix
        compile_system_transfer();
      }

      /// \todo find out what to do for disc/cont pressure spaces here...
      template<typename DomainLayer_, typename DomainLevel_, typename SpaceVelo_, typename SpacePres_>
      void assemble_gates(const DomainLayer_& dom_layer, const DomainLevel_& dom_level,
        const SpaceVelo_& space_velo, const SpacePres_& space_pres)
      {
        // set the gate comm
        this->gate_velo.set_comm(dom_layer.get_comm());
        this->gate_pres.set_comm(dom_layer.get_comm());
        this->gate_sys.set_comm(dom_layer.get_comm());

        // loop over all ranks
        for(Index i(0); i < dom_layer.size(); ++i)
        {
          Index rank = dom_layer.get_rank(i);

          // try to find our halo
          auto* halo = dom_level.find_halo_part(rank);
          XASSERT(halo != nullptr);

          // assemble the velocity components mirror
          ScalarMirror mirror_velo_comp;
          Assembly::MirrorAssembler::assemble_mirror(mirror_velo_comp, space_velo, *halo);

          // create a velocity mirror
          VeloMirror mirror_velo(mirror_velo_comp.clone());

          // create (empty) pressure mirror
          PresMirror mirror_pres;
          Assembly::MirrorAssembler::assemble_mirror(mirror_pres, space_pres, *halo);

          // create a system mirror
          SystemMirror mirror_sys(mirror_velo.clone(), mirror_pres.clone());

          // push mirror into gates
          this->gate_velo.push(int(rank), std::move(mirror_velo));
          if(!mirror_pres.get_gather().empty())
            this->gate_pres.push(int(rank), std::move(mirror_pres));
          this->gate_sys.push(int(rank), std::move(mirror_sys));
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

      template<typename SpaceVelo_, typename Cubature_>
      void assemble_velocity_transfer(const SpaceVelo_& space_velo_fine, const SpaceVelo_& space_velo_coarse, const Cubature_& cubature)
      {
        // get local transfer operator
        LocalVeloTransfer& loc_trans = this->transfer_velo.local();

        // get local transfer matrices
        LocalVeloTransferMatrix& loc_prol_v = loc_trans.get_mat_prol();
        LocalVeloTransferMatrix& loc_rest_v = loc_trans.get_mat_rest();

        // get the matrix blocks
        LocalScalarTransferMatrix& loc_prol_vx = loc_prol_v.get(0,0);
        LocalScalarTransferMatrix& loc_rest_vx = loc_rest_v.get(0,0);

        // assemble structure?
        if(loc_prol_vx.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol_vx, space_velo_fine, space_velo_coarse);

          for(int i(1); i < loc_prol_v.num_row_blocks; ++i)
            loc_prol_v.get(i,i) = loc_prol_vx.clone(LAFEM::CloneMode::Layout);
        }

        // get local velocity weight vector
        LocalVeloVector loc_vec_weight = loc_prol_v.create_vector_l();

        // get local weight vector components
        auto& loc_vec_wx = loc_vec_weight.get(0);

        // assemble prolongation matrix
        {
          loc_prol_vx.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol_vx, loc_vec_wx,
            space_velo_fine, space_velo_coarse, cubature);

          // synchronise weight vector
          for(int i(1); i < loc_vec_weight.num_blocks; ++i)
            loc_vec_weight.get(i).copy(loc_vec_wx);

          this->gate_velo.sync_0(loc_vec_weight);

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol_vx.scale_rows(loc_prol_vx, loc_vec_wx);

          // copy and transpose
          loc_rest_vx = loc_prol_vx.transpose();
          for(int i(1); i < loc_prol_v.num_row_blocks; ++i)
          {
            loc_prol_v.get(i,i) = loc_prol_vx.clone(LAFEM::CloneMode::Shallow);
            loc_rest_v.get(i,i) = loc_rest_vx.clone(LAFEM::CloneMode::Shallow);
          }
        }
      }

      template<typename SpacePres_, typename Cubature_>
      void assemble_pressure_transfer(const SpacePres_& space_pres_fine, const SpacePres_& space_pres_coarse, const Cubature_& cubature)
      {
        // get local transfer operator
        LocalPresTransfer& loc_trans = this->transfer_pres.local();

        // get local transfer matrices
        LocalPresTransferMatrix& loc_prol_p = loc_trans.get_mat_prol();
        LocalPresTransferMatrix& loc_rest_p = loc_trans.get_mat_rest();

        // assemble structure?
        if(loc_prol_p.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol_p, space_pres_fine, space_pres_coarse);
        }

        LocalPresVector loc_vec_weight = loc_prol_p.create_vector_l();

        // assemble prolongation matrix
        {
          loc_prol_p.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol_p, loc_vec_weight,
            space_pres_fine, space_pres_coarse, cubature);

          // synchronise weight vector
          this->gate_pres.sync_0(loc_vec_weight);

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol_p.scale_rows(loc_prol_p, loc_vec_weight);

          // copy and transpose
          loc_rest_p = loc_prol_p.transpose();
        }
      }

      template<typename SpaceVelo_, typename SpacePres_, typename Cubature_>
      void assemble_transfers(
        const SpaceVelo_& space_velo_fine, const SpacePres_& space_pres_fine,
        const SpaceVelo_& space_velo_coarse, const SpacePres_& space_pres_coarse, const Cubature_& cubature)
      {
        this->assemble_velocity_transfer(space_velo_fine, space_velo_coarse, cubature);
        this->assemble_pressure_transfer(space_pres_fine, space_pres_coarse, cubature);
        this->compile_system_transfer();
      }

      template<typename SpaceVelo_, typename Cubature_>
      void assemble_velocity_laplace_matrix(const SpaceVelo_& space_velo, const Cubature_& cubature, const DataType nu = DataType(1))
      {
        // get the local matrix A
        LocalMatrixBlockA& mat_loc_a = this->matrix_a.local();

        // get the diagonal blocks
        LocalScalarMatrix& mat_loc_a1 = mat_loc_a.get(0,0);

        // assemble matrix structure?
        if(mat_loc_a1.empty())
        {
          // assemble matrix structure for A11
          Assembly::SymbolicAssembler::assemble_matrix_std1(mat_loc_a1, space_velo);

          // clone layout for A22,...,Ann
          for(int i(1); i < mat_loc_a.num_row_blocks; ++i)
            mat_loc_a.get(i,i) = mat_loc_a1.clone(LAFEM::CloneMode::Layout);
        }

        // assemble velocity laplace matrix
        {
          mat_loc_a1.format();
          Assembly::Common::LaplaceOperator laplace_op;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_loc_a1, laplace_op, space_velo, cubature, nu);
        }

        // copy data into A22,...,Ann
        for(int i(1); i < mat_loc_a.num_row_blocks; ++i)
          mat_loc_a.get(i,i).copy(mat_loc_a1);
      }

      template<typename SpaceVelo_, typename SpacePres_, typename Cubature_>
      void assemble_grad_div_matrices(const SpaceVelo_& space_velo, const SpacePres_& space_pres, const Cubature_& cubature)
      {
        // get the local matrix B and D
        LocalMatrixBlockB& mat_loc_b = this->matrix_b.local();
        LocalMatrixBlockD& mat_loc_d = this->matrix_d.local();

        // get the matrix blocks
        LocalScalarMatrix& mat_loc_b1 = mat_loc_b.get(0,0);

        // assemble matrix structure?
        if(mat_loc_b1.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_std2(mat_loc_b1, space_velo, space_pres);
          for(int i(1); i < mat_loc_b.num_row_blocks; ++i)
            mat_loc_b.get(i,0) = mat_loc_b1.clone(LAFEM::CloneMode::Layout);
        }

        // assemble pressure gradient matrices
        for(int ider(0); ider < mat_loc_b.num_row_blocks; ++ider)
        {
          LocalScalarMatrix& mat_loc_bi = mat_loc_b.get(ider,0);
          mat_loc_bi.format();
          Assembly::Common::TestDerivativeOperator deri(ider);
          Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_loc_bi, deri, space_velo, space_pres, cubature, -DataType(1));
        }

        // assemble velocity divergence matrices
        /// \todo share matrix structures of D_i
        for(int i(0); i < mat_loc_b.num_row_blocks; ++i)
          mat_loc_d.get(0,i) = mat_loc_b.get(i,0).transpose();
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
        return (*this->matrix_sys).bytes () + (*this->matrix_s).bytes() + (*filter_sys).bytes();
      }

      void compile_system_filter()
      {
        (*filter_sys).template at<0>() = (*filter_velo).clone(LAFEM::CloneMode::Shallow);
        (*filter_sys).template at<1>() = (*filter_pres).clone(LAFEM::CloneMode::Shallow);
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

        compile_system_filter();
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
        return (*this->matrix_sys).bytes () + (*this->matrix_s).bytes() + (*filter_sys).bytes();
      }

      void compile_system_filter()
      {
        (*filter_sys).template at<0>() = (*filter_velo).clone(LAFEM::CloneMode::Shallow);
        (*filter_sys).template at<1>() = (*filter_pres).clone(LAFEM::CloneMode::Shallow);
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
        typename BaseClass::LocalPresVector& vec_loc_v = *vec_glob_v;
        typename BaseClass::LocalPresVector& vec_loc_w = *vec_glob_w;

        // fetch the frequency vector of the pressure gate
        typename BaseClass::LocalPresVector& vec_loc_f = this->gate_pres._freqs;

        // assemble the mean filter
        Assembly::MeanFilterAssembler::assemble(vec_loc_v, vec_loc_w, space_pres, cubature);

        // synchronise the vectors
        vec_glob_v.sync_1();
        vec_glob_w.sync_0();

        // build the mean filter
        fil_loc_p = LocalPresFilter(vec_loc_v.clone(), vec_loc_w.clone(), vec_loc_f.clone(), this->gate_pres.get_comm());
      }

    }; // struct StokesUnitVeloNonePresSystemLevel<...>

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

    public:
      explicit StokesBasicAssemblerLevel(DomainLevelType& dom_lvl) :
        domain_level(dom_lvl),
        mesh(domain_level.get_mesh()),
        trafo(mesh),
        space_velo(trafo),
        space_pres(trafo)
      {
      }

      virtual ~StokesBasicAssemblerLevel()
      {
      }

      template<typename SolVector_>
      void write_vtk(const String& vtk_name, const SolVector_& vector, const Dist::Comm& comm) const
      {
        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(this->mesh);

        // project velocity and pressure
        LAFEM::DenseVector<Mem::Main, double, Index> vtx_vx, vtx_vy, vtx_vz;
        Assembly::DiscreteVertexProjector::project(vtx_vx, vector.template at<0>().get(0), this->space_velo);
        Assembly::DiscreteVertexProjector::project(vtx_vy, vector.template at<0>().get(1), this->space_velo);
        if(vector.num_blocks > 2)
        {
          Assembly::DiscreteVertexProjector::project(vtx_vz, vector.template at<0>().get(2), this->space_velo);
          // write 3D velocity
          exporter.add_vertex_vector("velocity", vtx_vx.elements(), vtx_vy.elements(), vtx_vz.elements());
        }
        else
        {
          // write 2D velocity
          exporter.add_vertex_vector("velocity", vtx_vx.elements(), vtx_vy.elements());
        }

        // project pressure
        Cubature::DynamicFactory cubature("auto-degree:2");
        LAFEM::DenseVector<Mem::Main, double, Index> vtx_p;
        Assembly::DiscreteCellProjector::project(vtx_p, vector.template at<1>(), this->space_pres, cubature);

        // write pressure
        exporter.add_cell_scalar("pressure", vtx_p.elements());

        // finally, write the VTK file
        exporter.write(vtk_name, comm);
      }
    };
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STOKES_BASIC_HPP
