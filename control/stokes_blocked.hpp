#pragma once
#ifndef CONTROL_STOKES_BLOCKED_HPP
#define CONTROL_STOKES_BLOCKED_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/mean_filter.hpp>
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

#include <control/blocked_basic.hpp>
#include <control/scalar_basic.hpp>

namespace FEAT
{
  namespace Control
  {
    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
      >
    class StokesBlockedSystemLevel
    {
    public:
      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;
      static constexpr int dim = dim_;

      /// our compatible velocity sub-system type
      typedef BlockedBasicSystemLevel<dim, MemType_, DataType_, IndexType_, MatrixBlockA_, TransferMatrixV_> VeloSystemLevel;
      /// our compatible pressure sub-system type
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_, TransferMatrixV_> PresSystemLevel;

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
      typedef LAFEM::TupleDiagMatrix<TransferMatrixV_, TransferMatrixP_> LocalSystemTransferMatrix;

      // define local transfer operators
      typedef LAFEM::Transfer<LocalVeloTransferMatrix> LocalVeloTransfer;
      typedef LAFEM::Transfer<LocalPresTransferMatrix> LocalPresTransfer;
      typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

      // define mirror types
      typedef LAFEM::VectorMirror<MemType, DataType, IndexType> VeloMirror;
      typedef LAFEM::VectorMirror<MemType, DataType, IndexType> PresMirror;
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

      /// our system gate
      VeloGate gate_velo;
      PresGate gate_pres;
      SystemGate gate_sys;

      /// our coarse-level system muxer
      VeloMuxer coarse_muxer_velo;
      PresMuxer coarse_muxer_pres;
      SystemMuxer coarse_muxer_sys;

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

      /// CTOR
      StokesBlockedSystemLevel() :
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

      virtual ~StokesBlockedSystemLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return (*this->matrix_sys).bytes () + (*this->matrix_s).bytes() + transfer_sys.bytes();
      }

      void compile_system_matrix()
      {
        (*matrix_sys).block_a() = (*matrix_a).clone(LAFEM::CloneMode::Shallow);
        (*matrix_sys).block_b() = (*matrix_b).clone(LAFEM::CloneMode::Shallow);
        (*matrix_sys).block_d() = (*matrix_d).clone(LAFEM::CloneMode::Shallow);
      }

      void compile_system_transfer()
      {
        // clone content into our global transfer matrix
        transfer_sys.get_mat_prol().template at<0,0>() = transfer_velo.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<0,0>() = transfer_velo.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_prol().template at<1,1>() = transfer_pres.get_mat_prol().clone(LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<1,1>() = transfer_pres.get_mat_rest().clone(LAFEM::CloneMode::Shallow);
      }

      template<typename M_, typename D_, typename I_, typename SMA_, typename SMB_, typename SMD_, typename SM_, typename TV_, typename TP_>
      void convert(const StokesBlockedSystemLevel<dim_, M_, D_, I_, SMA_, SMB_, SMD_, SM_, TV_, TP_> & other)
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

        transfer_velo.convert(&coarse_muxer_velo, other.transfer_velo);
        transfer_pres.convert(&coarse_muxer_pres, other.transfer_pres);

        this->compile_system_matrix();
        this->compile_system_transfer();
      }

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

          // create (empty) velocity mirror
          VeloMirror mirror_velo;
          Assembly::MirrorAssembler::assemble_mirror(mirror_velo, space_velo, *halo);

          // create (empty) pressure mirror
          PresMirror mirror_pres;
          Assembly::MirrorAssembler::assemble_mirror(mirror_pres, space_pres, *halo);

          // create a system mirror
          SystemMirror mirror_sys(mirror_velo.clone(), mirror_pres.clone());

          // push mirror into gates
          gate_velo.push(int(rank), std::move(mirror_velo));
          if(!mirror_pres.get_gather().empty())
            gate_pres.push(int(rank), std::move(mirror_pres));
          gate_sys.push(int(rank), std::move(mirror_sys));
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
        LocalVeloTransferMatrix& loc_prol_wrapped = loc_trans.get_mat_prol();
        LocalVeloTransferMatrix& loc_rest_wrapped = loc_trans.get_mat_rest();

        // get the unwrapped types
        typename LocalVeloTransferMatrix::BaseClass& loc_prol = loc_prol_wrapped;
        typename LocalVeloTransferMatrix::BaseClass& loc_rest = loc_rest_wrapped;

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_velo_fine, space_velo_coarse);
        }

        // create a local weight vector
        LocalVeloVector loc_vec_weight = loc_prol_wrapped.create_vector_l();

        // create a scalar weight vector for the assembly
        auto loc_scal_vec_weight = loc_prol.create_vector_l();

        // get the data arrays of the weight vectors
        auto* v_wb = loc_vec_weight.elements();
        auto* v_ws = loc_scal_vec_weight.elements();

        // assemble prolongation matrix
        {
          loc_prol.format();
          loc_scal_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_scal_vec_weight,
            space_velo_fine, space_velo_coarse, cubature);

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
      }

      template<typename SpacePres_, typename Cubature_>
      void assemble_pressure_transfer(const SpacePres_& space_pres_fine, const SpacePres_& space_pres_coarse, const Cubature_& cubature)
      {
        // get local transfer operator
        LocalPresTransfer& loc_trans = this->transfer_pres.local();

        // get local transfer matrices
        LocalPresTransferMatrix& loc_prol = loc_trans.get_mat_prol();
        LocalPresTransferMatrix& loc_rest = loc_trans.get_mat_rest();

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_pres_fine, space_pres_coarse);
        }

        // get local pressure weight vector
        LocalPresVector loc_vec_weight = loc_prol.create_vector_l();

        // assemble prolongation matrix
        {
          loc_prol.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight,
            space_pres_fine, space_pres_coarse, cubature);

          // synchronise weight vector
          this->gate_pres.sync_0(loc_vec_weight);

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol.scale_rows(loc_prol, loc_vec_weight);

          // copy and transpose
          loc_rest = loc_prol.transpose();
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

      template<typename SpaceVelo_, typename SpacePres_, typename Cubature_>
      void assemble_grad_div_matrices(const SpaceVelo_& space_velo, const SpacePres_& space_pres, const Cubature_& cubature)
      {
        Assembly::GradPresDivVeloAssembler::assemble(this->matrix_b.local(), this->matrix_d.local(),
          space_velo, space_pres, cubature);
      }

      template<typename SpaceVelo_>
      void assemble_velo_struct(const SpaceVelo_& space_velo)
      {
        // assemble matrix structure
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_a.local(), space_velo);
      }


      template<typename SpacePres_>
      void assemble_pres_struct(const SpacePres_& space_pres)
      {
        // assemble matrix structure
        Assembly::SymbolicAssembler::assemble_matrix_std1(this->matrix_s.local(), space_pres);
      }
    }; // class StokesBlockedSystemLevel<...>

    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
      typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
      typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrixV_ = LAFEM::SparseMatrixBWrappedCSR<MemType_, DataType_, IndexType_, dim_>,
      typename TransferMatrixP_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
      >
    class StokesBlockedUnitVeloNonePresSystemLevel :
      public StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_>
    {
    public:
      typedef StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_,
        MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_, TransferMatrixV_, TransferMatrixP_> BaseClass;

      // define local filter types
      typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalVeloFilter;
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
        return this->filter_sys.bytes() + BaseClass::bytes();
      }

      void compile_system_filter()
      {
        (*filter_sys).template at<0>() = (*filter_velo).clone(LAFEM::CloneMode::Shallow);
        (*filter_sys).template at<1>() = (*filter_pres).clone(LAFEM::CloneMode::Shallow);
      }
    };

    template<
      typename SpaceVelo_,
      typename SpacePres_>
    class StokesBlockedAssemblerLevel
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
      explicit StokesBlockedAssemblerLevel(DomainLevelType& dom_lvl) :
        domain_level(dom_lvl),
        mesh(domain_level.get_mesh()),
        trafo(mesh),
        space_velo(trafo),
        space_pres(trafo)
      {
      }

      virtual ~StokesBlockedAssemblerLevel()
      {
      }
    };
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STOKES_BLOCKED_HPP
