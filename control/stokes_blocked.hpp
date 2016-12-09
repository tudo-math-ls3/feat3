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
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/gpdv_assembler.hpp>

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
        transfer_sys.get_mat_prol().template at<0,0>().clone(transfer_velo.get_mat_prol(), LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<0,0>().clone(transfer_velo.get_mat_rest(), LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_prol().template at<1,1>().clone(transfer_velo.get_mat_prol(), LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_rest().template at<1,1>().clone(transfer_velo.get_mat_rest(), LAFEM::CloneMode::Shallow);
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
    }; // class StokesBlockedSystemLevel<...>

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
      Cubature::DynamicFactory cubature;

    public:
      explicit StokesBlockedAssemblerLevel(DomainLevelType& dom_lvl) :
        domain_level(dom_lvl),
        mesh(domain_level.get_mesh()),
        trafo(mesh),
        space_velo(trafo),
        space_pres(trafo),
        cubature("auto-degree:" + stringify(Math::sqr(SpaceVeloType::local_degree)+2))
      {
      }

      virtual ~StokesBlockedAssemblerLevel()
      {
      }

      template<typename SystemLevel_>
      void assemble_gates(const DomainLayerType& dom_layer, SystemLevel_& sys_level)
      {
        // create our gates
        typename SystemLevel_::VeloGate& gate_velo = sys_level.gate_velo;
        typename SystemLevel_::PresGate& gate_pres = sys_level.gate_pres;
        typename SystemLevel_::SystemGate& gate_sys  = sys_level.gate_sys;

        // set the gate comm
        gate_velo.set_comm(dom_layer.get_comm());
        gate_pres.set_comm(dom_layer.get_comm());
        gate_sys.set_comm(dom_layer.get_comm());

        // loop over all ranks
        for(Index i(0); i < dom_layer.size(); ++i)
        {
          Index rank = dom_layer.get_rank(i);

          // try to find our halo
          auto* halo = domain_level.find_halo_part(rank);
          if(halo == nullptr)
            throw InternalError("ERROR: Halo not found!");

          // create (empty) velocity mirror
          typename SystemLevel_::VeloMirror mirror_velo;
          Assembly::MirrorAssembler::assemble_mirror(mirror_velo, this->space_velo, *halo);

          // create (empty) pressure mirror
          typename SystemLevel_::PresMirror mirror_pres;
          Assembly::MirrorAssembler::assemble_mirror(mirror_pres, this->space_pres, *halo);

          // create a system mirror
          typename SystemLevel_::SystemMirror mirror_sys(mirror_velo.clone(), mirror_pres.clone());

          // push mirror into gates
          gate_velo.push(int(rank), std::move(mirror_velo));
          gate_pres.push(int(rank), std::move(mirror_pres));
          gate_sys.push(int(rank), std::move(mirror_sys));
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
      void assemble_grad_div_matrices(SystemLevel_& sys_level)
      {
        Assembly::GradPresDivVeloAssembler::assemble(*sys_level.matrix_b, *(sys_level.matrix_d),
          this->space_velo, this->space_pres, this->cubature);
      }

      template<typename SystemLevel_>
      void assemble_pres_laplace(SystemLevel_& sys_level)
      {
        // get the global matrix
        typename SystemLevel_::GlobalSchurMatrix& mat_glob = sys_level.matrix_s;

        // get the local matrix
        typename SystemLevel_::LocalSchurMatrix& mat_loc = *mat_glob;

        // assemble matrix structure?
        if (mat_loc.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_std1(mat_loc, this->space_pres);
        }

        // assemble pressure laplace matrix
        {
          mat_loc.format();
          Assembly::Common::LaplaceOperator laplace_op;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_loc, laplace_op, this->space_pres, this->cubature);
        }
      }

      template<typename SystemLevel_>
      void assemble_velo_struct(SystemLevel_& sys_level)
      {
        // get the global matrix
        typename SystemLevel_::GlobalMatrixBlockA& mat_glob = sys_level.matrix_a;

        // get the local matrix
        typename SystemLevel_::LocalMatrixBlockA& mat_loc = *mat_glob;

        // assemble matrix structure?
        Assembly::SymbolicAssembler::assemble_matrix_std1(mat_loc, this->space_velo);
      }

      template<typename SystemLevel_>
      void assemble_velo_transfer(SystemLevel_& sys_level_fine, StokesBlockedAssemblerLevel& level_coarse)
      {
        // get global transfer operator
        typename SystemLevel_::GlobalVeloTransfer& glob_trans = sys_level_fine.transfer_velo;

        // get local transfer operator
        typename SystemLevel_::LocalVeloTransfer& loc_trans = glob_trans.local();

        // get local transfer matrices
        typename SystemLevel_::LocalVeloTransferMatrix& loc_prol_wrapped = loc_trans.get_mat_prol();
        typename SystemLevel_::LocalVeloTransferMatrix& loc_rest_wrapped = loc_trans.get_mat_rest();

        // get the unwrapped types
        typedef typename SystemLevel_::LocalVeloTransferMatrix WrappedTransfer;
        typename WrappedTransfer::BaseClass& loc_prol = loc_prol_wrapped;
        typename WrappedTransfer::BaseClass& loc_rest = loc_rest_wrapped;

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, this->space_velo, level_coarse.space_velo);
        }

        // create a local weight vector
        typename SystemLevel_::LocalVeloVector loc_vec_weight = loc_prol_wrapped.create_vector_l();

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
            this->space_velo, level_coarse.space_velo, this->cubature);

          // copy weights from scalar to blocked
          for(Index i(0); i < loc_prol.rows(); ++i)
            v_wb[i] = v_ws[i];

          // synchronise blocked weight vector
          sys_level_fine.gate_velo.sync_0(loc_vec_weight);

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

      template<typename SystemLevel_>
      void assemble_pres_transfer(SystemLevel_& sys_level_fine, StokesBlockedAssemblerLevel& level_coarse)
      {
        // get global transfer operator
        typename SystemLevel_::GlobalPresTransfer& glob_trans = sys_level_fine.transfer_pres;

        // get local transfer operator
        typename SystemLevel_::LocalPresTransfer& loc_trans = glob_trans.local();

        // get local transfer matrices
        typename SystemLevel_::LocalPresTransferMatrix& loc_prol = loc_trans.get_mat_prol();
        typename SystemLevel_::LocalPresTransferMatrix& loc_rest = loc_trans.get_mat_rest();

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, this->space_pres, level_coarse.space_pres);
        }

        // get local pressure weight vector
        typename SystemLevel_::LocalPresVector loc_vec_weight = loc_prol.create_vector_l();

        // assemble prolongation matrix
        {
          loc_prol.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight,
            this->space_pres, level_coarse.space_pres, this->cubature);

          // synchronise weight vector
          sys_level_fine.gate_pres.sync_0(loc_vec_weight);

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol.scale_rows(loc_prol, loc_vec_weight);

          // copy and transpose
          loc_rest = loc_prol.transpose();
        }
      }
    };
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STOKES_BLOCKED_HPP
