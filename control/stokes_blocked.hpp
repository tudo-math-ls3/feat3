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
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>>
    class StokesBlockedSystemLevel
    {
    public:
      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;
      static constexpr int dim = dim_;

      /// our compatible velocity sub-system type
      typedef BlockedBasicSystemLevel<dim, MemType_, DataType_, IndexType_, MatrixBlockA_> VeloSystemLevel;
      /// our compatible pressure sub-system type
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> PresSystemLevel;

      /// scalar types
      typedef ScalarMatrix_ LocalScalarMatrix;
      typedef typename LocalScalarMatrix::VectorTypeL LocalScalarVector;

      /// define local blocked matrix type
      typedef MatrixBlockA_ LocalMatrixBlockA;
      typedef MatrixBlockB_ LocalMatrixBlockB;
      typedef MatrixBlockD_ LocalMatrixBlockD;
      typedef LocalScalarMatrix LocalSchurMatrix;
      typedef LAFEM::SaddlePointMatrix<LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockD> LocalSystemMatrix;

      /// define local vector types
      typedef typename LocalMatrixBlockB::VectorTypeL LocalVeloVector;
      typedef typename LocalMatrixBlockD::VectorTypeL LocalPresVector;
      typedef LAFEM::TupleVector<LocalVeloVector, LocalPresVector> LocalSystemVector;

      /// define mirror types
      //typedef LAFEM::VectorMirrorBlocked<MemType_, DataType, IndexType, dim_> VeloMirror;
      typedef LAFEM::VectorMirror<MemType, DataType, IndexType> VeloMirror;
      typedef LAFEM::VectorMirror<MemType, DataType, IndexType> PresMirror;
      typedef LAFEM::TupleMirror<VeloMirror, PresMirror> SystemMirror;

      // define gates
      typedef Global::Gate<LocalVeloVector, VeloMirror> VeloGate;
      typedef Global::Gate<LocalPresVector, PresMirror> PresGate;
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;

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

      /* ***************************************************************************************** */

      /// our system gate
      VeloGate gate_velo;
      PresGate gate_pres;
      SystemGate gate_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;
      GlobalMatrixBlockA matrix_a;
      GlobalMatrixBlockB matrix_b;
      GlobalMatrixBlockD matrix_d;
      GlobalSchurMatrix matrix_s;

      /// CTOR
      StokesBlockedSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys),
        matrix_a(&gate_velo, &gate_velo),
        matrix_b(&gate_velo, &gate_pres),
        matrix_d(&gate_pres, &gate_velo),
        matrix_s(&gate_pres, &gate_pres)
      {
      }

      virtual ~StokesBlockedSystemLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return (*this->matrix_sys).bytes () + (*this->matrix_s).bytes();
      }

      void compile_system_matrix()
      {
        (*matrix_sys).block_a() = (*matrix_a).clone(LAFEM::CloneMode::Shallow);
        (*matrix_sys).block_b() = (*matrix_b).clone(LAFEM::CloneMode::Shallow);
        (*matrix_sys).block_d() = (*matrix_d).clone(LAFEM::CloneMode::Shallow);
      }

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const BlockedBasicSystemLevel<dim_, M_, D_, I_, SM_> & other)
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
    }; // class StokesBlockedSystemLevel<...>


    template<typename SystemLevel_>
    class StokesBlockedTransferLevel
    {
    public:
      typedef SystemLevel_ SystemLevel;
      static constexpr int dim = SystemLevel_::dim;
      typedef typename SystemLevel_::MemType MemType;
      typedef typename SystemLevel_::DataType DataType;
      typedef typename SystemLevel_::IndexType IndexType;

      typedef typename SystemLevel_::LocalScalarVector LocalScalarVector;
      typedef typename SystemLevel_::LocalScalarMatrix LocalScalarMatrix;

      typedef LAFEM::SparseMatrixBWrappedCSR<MemType, DataType, IndexType, dim> LocalVeloTransferMatrix;
      typedef LocalScalarMatrix LocalPresTransferMatrix;
      typedef LAFEM::TupleDiagMatrix<LocalVeloTransferMatrix, LocalPresTransferMatrix> LocalSystemTransferMatrix;

      typedef Global::Matrix<LocalVeloTransferMatrix, typename SystemLevel::VeloMirror, typename SystemLevel::VeloMirror> GlobalVeloTransferMatrix;
      typedef Global::Matrix<LocalPresTransferMatrix, typename SystemLevel::PresMirror, typename SystemLevel::PresMirror> GlobalPresTransferMatrix;
      typedef Global::Matrix<LocalSystemTransferMatrix, typename SystemLevel::SystemMirror, typename SystemLevel::SystemMirror> GlobalSystemTransferMatrix;

      // (global) transfer matrices
      GlobalVeloTransferMatrix prol_velo, rest_velo;
      GlobalPresTransferMatrix prol_pres, rest_pres;
      GlobalSystemTransferMatrix prol_sys, rest_sys;

      explicit StokesBlockedTransferLevel()
      {
      }

      explicit StokesBlockedTransferLevel(SystemLevel_& lvl_coarse, SystemLevel_& lvl_fine) :
        prol_velo(&lvl_fine.gate_velo, &lvl_coarse.gate_velo),
        rest_velo(&lvl_coarse.gate_velo, &lvl_fine.gate_velo),
        prol_pres(&lvl_fine.gate_pres, &lvl_coarse.gate_pres),
        rest_pres(&lvl_coarse.gate_pres, &lvl_fine.gate_pres),
        prol_sys(&lvl_fine.gate_sys, &lvl_coarse.gate_sys),
        rest_sys(&lvl_coarse.gate_sys, &lvl_fine.gate_sys)
      {
      }

      virtual ~StokesBlockedTransferLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return (*prol_sys).bytes() + (*rest_sys).bytes();
      }

      void compile_system_transfer()
      {
        // clone content into our global transfer matrix
        (*prol_sys).template at<0,0>().clone(*prol_velo, LAFEM::CloneMode::Shallow);
        (*prol_sys).template at<1,1>().clone(*prol_pres, LAFEM::CloneMode::Shallow);
        (*rest_sys).template at<0,0>().clone(*rest_velo, LAFEM::CloneMode::Shallow);
        (*rest_sys).template at<1,1>().clone(*rest_pres, LAFEM::CloneMode::Shallow);
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source StokesBasicTransferLevel content as content of current StokesBasicTransferLevel.
       *
       */
      template <typename SL_>
      void convert(SystemLevel_ & lvl_coarse, SystemLevel_ & lvl_fine, const StokesBlockedTransferLevel<SL_> & other)
      {
        prol_velo.convert(&lvl_fine.gate_velo, &lvl_coarse.gate_velo, other.prol_velo);
        rest_velo.convert(&lvl_coarse.gate_velo, &lvl_fine.gate_velo, other.rest_velo);
        prol_pres.convert(&lvl_fine.gate_pres, &lvl_coarse.gate_pres, other.prol_pres);
        rest_pres.convert(&lvl_coarse.gate_pres, &lvl_fine.gate_pres, other.rest_pres);

        (*prol_sys).template at<0,0>().clone(*prol_velo, LAFEM::CloneMode::Shallow);
        (*prol_sys).template at<1,1>().clone(*prol_pres, LAFEM::CloneMode::Shallow);
        (*rest_sys).template at<0,0>().clone(*rest_velo, LAFEM::CloneMode::Shallow);
        (*rest_sys).template at<1,1>().clone(*rest_pres, LAFEM::CloneMode::Shallow);
      }
    }; // struct StokesBlockedTransferLevel<...>


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

      template<typename TransferLevel_>
      void assemble_velo_transfer(TransferLevel_& trans_level, StokesBlockedAssemblerLevel& level_coarse)
      {
        // get global transfer matrices
        typename TransferLevel_::GlobalVeloTransferMatrix& glob_prol = trans_level.prol_velo;
        typename TransferLevel_::GlobalVeloTransferMatrix& glob_rest = trans_level.rest_velo;

        // get local transfer matrices
        typename TransferLevel_::LocalVeloTransferMatrix& loc_prol_wrapped = (*glob_prol);
        typename TransferLevel_::LocalVeloTransferMatrix& loc_rest_wrapped = (*glob_rest);

        // get the unwrapped types
        typedef typename TransferLevel_::LocalVeloTransferMatrix WrappedTransfer;
        typename WrappedTransfer::BaseClass& loc_prol = loc_prol_wrapped;
        typename WrappedTransfer::BaseClass& loc_rest = loc_rest_wrapped;

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, this->space_velo, level_coarse.space_velo);
        }

        // create a global pressure weight vector
        auto glob_vec_weight = glob_prol.create_vector_l();

        // get local pressure weight vector
        auto& loc_vec_weight = (*glob_vec_weight);

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
          glob_vec_weight.sync_0();

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

      template<typename TransferLevel_>
      void assemble_pres_transfer(TransferLevel_& trans_level, StokesBlockedAssemblerLevel& level_coarse)
      {
        // get global transfer matrices
        typename TransferLevel_::GlobalPresTransferMatrix& glob_prol = trans_level.prol_pres;
        typename TransferLevel_::GlobalPresTransferMatrix& glob_rest = trans_level.rest_pres;

        // get local transfer matrices
        typename TransferLevel_::LocalPresTransferMatrix& loc_prol = (*glob_prol);
        typename TransferLevel_::LocalPresTransferMatrix& loc_rest = (*glob_rest);

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, this->space_pres, level_coarse.space_pres);
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
            this->space_pres, level_coarse.space_pres, this->cubature);

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
    };
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STOKES_BLOCKED_HPP
