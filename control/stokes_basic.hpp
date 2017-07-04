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
        transfer_sys.get_mat_trunc().template at<0,0>().clone(transfer_velo.get_mat_trunc(), LAFEM::CloneMode::Shallow);
        transfer_sys.get_mat_trunc().template at<1,1>().clone(transfer_pres.get_mat_trunc(), LAFEM::CloneMode::Shallow);
        transfer_sys.compile();
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

        compile_system_transfer();
      }

      /// \todo find out what to do for disc/cont pressure spaces here...
      template<typename DomainLevel_>
      void assemble_gates(const Domain::VirtualLevel<DomainLevel_>& virt_dom_lvl)
      {
        const auto& dom_level = virt_dom_lvl.level();
        const auto& dom_layer = virt_dom_lvl.layer();
        const auto& space_velo = dom_level.space_velo;
        const auto& space_pres = dom_level.space_pres;

        // set the gate comm
        this->gate_velo.set_comm(dom_layer.comm_ptr());
        this->gate_pres.set_comm(dom_layer.comm_ptr());
        this->gate_sys.set_comm(dom_layer.comm_ptr());

        // loop over all ranks
        for(Index i(0); i < dom_layer.neighbour_count(); ++i)
        {
          int rank = dom_layer.neighbour_rank(i);

          // try to find our halo
          const auto* halo = dom_level.find_halo_part(rank);
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
          this->gate_velo.push(rank, std::move(mirror_velo));
          if(!mirror_pres.empty())
            this->gate_pres.push(rank, std::move(mirror_pres));
          this->gate_sys.push(rank, std::move(mirror_sys));
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

      template<typename DomainLevel_>
      void assemble_coarse_muxers(const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse)
      {
        // assemble muxer parent
        if(virt_lvl_coarse.is_parent())
        {
          XASSERT(virt_lvl_coarse.is_child());

          const auto& layer_c = virt_lvl_coarse.layer_c();
          const DomainLevel_& level_p = virt_lvl_coarse.level_p();

          // loop over all children
          for(Index i(0); i < layer_c.child_count(); ++i)
          {
            const auto* child = level_p.find_patch_part(int(i));
            XASSERT(child != nullptr);
            SystemMirror child_mirror_sys;
            VeloMirror& child_mirror_v = child_mirror_sys.template at<0>();
            PresMirror& child_mirror_p = child_mirror_sys.template at<1>();
            Assembly::MirrorAssembler::assemble_mirror(child_mirror_v._sub_mirror, level_p.space_velo, *child);
            Assembly::MirrorAssembler::assemble_mirror(child_mirror_p, level_p.space_pres, *child);
            this->coarse_muxer_velo.push_child(child_mirror_v.clone(LAFEM::CloneMode::Shallow));
            this->coarse_muxer_pres.push_child(child_mirror_p.clone(LAFEM::CloneMode::Shallow));
            this->coarse_muxer_sys.push_child(std::move(child_mirror_sys));
          }
        }

        // assemble muxer child
        if(virt_lvl_coarse.is_child())
        {
          const auto& layer_c = virt_lvl_coarse.layer_c();
          const DomainLevel_& level_c = virt_lvl_coarse.level_c();

          SystemMirror parent_mirror_sys;
          VeloMirror& parent_mirror_v = parent_mirror_sys.template at<0>();
          PresMirror& parent_mirror_p = parent_mirror_sys.template at<1>();

          // manually set up an identity gather/scatter matrix
          {
            Index n = level_c.space_velo.get_num_dofs();
            parent_mirror_v._sub_mirror = ScalarMirror(n, n);
            auto* idx = parent_mirror_v._sub_mirror.indices();
            for(Index i(0); i < n; ++i)
              idx[i] = i;
          }
          {
            Index n = level_c.space_pres.get_num_dofs();
            parent_mirror_p = ScalarMirror(n, n);
            auto* idx = parent_mirror_p.indices();
            for(Index i(0); i < n; ++i)
              idx[i] = i;
          }

          // set parent and sibling comms
          this->coarse_muxer_velo.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            parent_mirror_v.clone(LAFEM::CloneMode::Shallow)
          );
          this->coarse_muxer_pres.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            parent_mirror_p.clone(LAFEM::CloneMode::Shallow)
          );
          this->coarse_muxer_sys.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            std::move(parent_mirror_sys)
          );

          // compile muxer
          LocalVeloVector tmpl_v(level_c.space_velo.get_num_dofs());
          LocalPresVector tmpl_p(level_c.space_pres.get_num_dofs());
          LocalSystemVector tmpl_s(tmpl_v.clone(), tmpl_p.clone());
          this->coarse_muxer_velo.compile(tmpl_v);
          this->coarse_muxer_pres.compile(tmpl_p);
          this->coarse_muxer_sys.compile(tmpl_s);
        }
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_velocity_transfer(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_velo;
        const auto& space_c = level_c.space_velo;

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
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol_vx, space_f, space_c);

          for(int i(1); i < loc_prol_v.num_row_blocks; ++i)
            loc_prol_v.get(i,i) = loc_prol_vx.clone(LAFEM::CloneMode::Layout);
        }

        // get local velocity weight vector
        LocalVeloVector loc_vec_weight = loc_prol_v.create_vector_l();

        // get local weight vector components
        LocalScalarVector& loc_vec_wx = loc_vec_weight.get(0);

        // assemble prolongation matrix
        {
          loc_prol_vx.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol_vx, loc_vec_wx, space_f, space_c, cubature);

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

        // compile velocity transfer
        this->transfer_velo.compile();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_pressure_transfer(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_pres;
        const auto& space_c = level_c.space_pres;

        // get local transfer operator
        LocalPresTransfer& loc_trans = this->transfer_pres.local();

        // get local transfer matrices
        LocalPresTransferMatrix& loc_prol_p = loc_trans.get_mat_prol();
        LocalPresTransferMatrix& loc_rest_p = loc_trans.get_mat_rest();

        // assemble structure?
        if(loc_prol_p.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol_p, space_f, space_c);
        }

        LocalPresVector loc_vec_weight = loc_prol_p.create_vector_l();

        // assemble prolongation matrix
        {
          loc_prol_p.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol_p, loc_vec_weight, space_f, space_c, cubature);

          // synchronise weight vector
          this->gate_pres.sync_0(loc_vec_weight);

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol_p.scale_rows(loc_prol_p, loc_vec_weight);

          // copy and transpose
          loc_rest_p = loc_prol_p.transpose();
        }

        // compile pressure transfer
        this->transfer_pres.compile();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_transfers(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature)
      {
        this->assemble_velocity_transfer(virt_lvl_fine, virt_lvl_coarse, cubature);
        this->assemble_pressure_transfer(virt_lvl_fine, virt_lvl_coarse, cubature);
        this->compile_system_transfer();
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_velocity_truncation(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const StokesBasicSystemLevel* sys_lvl_coarse = nullptr)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT((sys_lvl_coarse != nullptr) || !virt_lvl_coarse.is_parent());

        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_velo;
        const auto& space_c = level_c.space_velo;

        // get local transfer operator
        LocalVeloTransfer& loc_trans = this->transfer_velo.local();

        // get local transfer matrices
        const LocalVeloTransferMatrix& loc_rest_v = loc_trans.get_mat_rest();
        LocalVeloTransferMatrix& loc_trunc_v = loc_trans.get_mat_trunc();

        // get the matrix blocks
        const LocalScalarTransferMatrix& loc_rest_vx = loc_rest_v.get(0,0);
        LocalScalarTransferMatrix& loc_trunc_vx = loc_trunc_v.get(0,0);

        // restriction matrix must be already assembled
        XASSERTM(loc_rest_vx.size() > Index(0), "you need to call 'assemble_velocity_transfer' first");

        // clone sparsity pattern of restriction matrix
        loc_trunc_vx = loc_rest_vx.clone(LAFEM::CloneMode::Layout);

        // create local weight vector
        LocalVeloVector loc_vec_weight = loc_rest_v.create_vector_l();

        // get local weight vector component
        LocalScalarVector& loc_vec_wx = loc_vec_weight.get(0);

        // format
        loc_trunc_vx.format();
        loc_vec_wx.format();

        // assemble truncation matrix
        Assembly::GridTransfer::assemble_truncation(loc_trunc_vx, loc_vec_wx, space_f, space_c, cubature);

        // expand weight vector
        for(int i(1); i < loc_vec_weight.num_blocks; ++i)
          loc_vec_weight.get(i).copy(loc_vec_wx);

        // We now need to synchronise the weight vector in analogy to the prolongation matrix assembly.
        // Note that the weight vector is now a coarse-level vector, so we need to synchronise over
        // the coarse-level gate. This may be a bit more complicated if the coarse level is a ghost
        // level, as in this case we have to join/split around the synch operation.

        // synchronise weight vector using the muxer/gate
        if(!virt_lvl_coarse.is_child())
        {
          // The coarse level is a simple (non-child) level that exists on all processes,
          // so simply synch over the coarse-level gate:
          sys_lvl_coarse->gate_velo.sync_0(loc_vec_weight);
        }
        else if(virt_lvl_coarse.is_parent())
        {
          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then synch that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a synch over the (non-existent) coarse-level
          // child gate, which is what we actually require here...

          // create temporary vector on parent partitioning
          LocalVeloVector loc_tmp = sys_lvl_coarse->gate_velo._freqs.clone(LAFEM::CloneMode::Allocate);

          // join child weights over muxer
          this->coarse_muxer_velo.join(loc_vec_weight, loc_tmp);

          // sync over coarse gate
          sys_lvl_coarse->gate_velo.sync_0(loc_tmp);

          // split over muxer
          this->coarse_muxer_velo.split(loc_vec_weight, loc_tmp);
        }
        else // ghost
        {
          // The coarse level is a ghost level, i.e. a child but not a parent. In this case, we
          // only have to participate in the join/send operations of the muxer, which are part
          // of the operations executed on the parents handled by the else-if case above.

          this->coarse_muxer_velo.join_send(loc_vec_weight);

          // parent performs sync over its gate here (see above else-if)

          this->coarse_muxer_velo.split_recv(loc_vec_weight);
        }

        // invert components
        loc_vec_weight.component_invert(loc_vec_weight);

        // scale reduction matrix
        loc_trunc_vx.scale_rows(loc_trunc_vx, loc_vec_wx);

        // clone for remaining components
        for(int i(1); i < loc_trunc_v.num_row_blocks; ++i)
          loc_trunc_v.get(i,i) = loc_trunc_vx.clone(LAFEM::CloneMode::Shallow);
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_pressure_truncation(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const StokesBasicSystemLevel* sys_lvl_coarse = nullptr)
      {
        // if the coarse level is a parent, then we need the coarse system level
        XASSERT((sys_lvl_coarse != nullptr) || !virt_lvl_coarse.is_parent());

        // get fine and coarse domain levels
        const DomainLevel_& level_f = *virt_lvl_fine;
        const DomainLevel_& level_c = virt_lvl_coarse.is_child() ? virt_lvl_coarse.level_c() : *virt_lvl_coarse;

        const auto& space_f = level_f.space_pres;
        const auto& space_c = level_c.space_pres;

        // get local transfer operator
        LocalPresTransfer& loc_trans = this->transfer_pres.local();

        // get local transfer matrices
        const LocalPresTransferMatrix& loc_rest = loc_trans.get_mat_rest();
        LocalPresTransferMatrix& loc_trunc = loc_trans.get_mat_trunc();

        // restriction matrix must be already assembled
        XASSERTM(loc_rest.size() > Index(0), "you need to call 'assemble_prescity_transfer' first");

        // clone sparsity pattern of restriction matrix
        loc_trunc = loc_rest.clone(LAFEM::CloneMode::Layout);

        // create local weight vector
        LocalPresVector loc_vec_weight = loc_trunc.create_vector_l();

        // format
        loc_trunc.format();
        loc_vec_weight.format();

        // assemble truncation matrix
        Assembly::GridTransfer::assemble_truncation(loc_trunc, loc_vec_weight, space_f, space_c, cubature);

        // We now need to synchronise the weight vector in analogy to the prolongation matrix assembly.
        // Note that the weight vector is now a coarse-level vector, so we need to synchronise over
        // the coarse-level gate. This may be a bit more complicated if the coarse level is a ghost
        // level, as in this case we have to join/split around the synch operation.

        // synchronise weight vector using the muxer/gate
        if(!virt_lvl_coarse.is_child())
        {
          // The coarse level is a simple (non-child) level that exists on all processes,
          // so simply synch over the coarse-level gate:
          sys_lvl_coarse->gate_pres.sync_0(loc_vec_weight);
        }
        else if(virt_lvl_coarse.is_parent())
        {
          // The coarse level is a child level and this is one of the parent processes which contain
          // the coarse-level gate. So we first need to join the weight vector onto the parent processes
          // first, then synch that joined vector over the parent gate and finally split the result
          // to all child processes -- this "emulates" a synch over the (non-existent) coarse-level
          // child gate, which is what we actually require here...

          // create temporary vector on parent partitioning
          LocalPresVector loc_tmp = sys_lvl_coarse->gate_pres._freqs.clone(LAFEM::CloneMode::Allocate);

          // join child weights over muxer
          this->coarse_muxer_pres.join(loc_vec_weight, loc_tmp);

          // sync over coarse gate
          sys_lvl_coarse->gate_pres.sync_0(loc_tmp);

          // split over muxer
          this->coarse_muxer_pres.split(loc_vec_weight, loc_tmp);
        }
        else // ghost
        {
          // The coarse level is a ghost level, i.e. a child but not a parent. In this case, we
          // only have to participate in the join/send operations of the muxer, which are part
          // of the operations executed on the parents handled by the else-if case above.

          this->coarse_muxer_pres.join_send(loc_vec_weight);

          // parent performs sync over its gate here (see above else-if)

          this->coarse_muxer_pres.split_recv(loc_vec_weight);
        }

        // invert components
        loc_vec_weight.component_invert(loc_vec_weight);

        // scale reduction matrix
        loc_trunc.scale_rows(loc_trunc, loc_vec_weight);
      }

      template<typename DomainLevel_, typename Cubature_>
      void assemble_truncations(
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_fine,
        const Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse,
        const Cubature_& cubature,
        const StokesBasicSystemLevel* sys_lvl_coarse = nullptr)
      {
        this->assemble_velocity_truncation(virt_lvl_fine, virt_lvl_coarse, cubature, sys_lvl_coarse);
        this->assemble_pressure_truncation(virt_lvl_fine, virt_lvl_coarse, cubature, sys_lvl_coarse);
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
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STOKES_BASIC_HPP
