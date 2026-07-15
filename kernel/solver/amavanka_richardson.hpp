// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/solver/iterative.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/solver/amavanka.hpp>
// reuse the recursive Intern::FilterHelper (meta-filter -> flat unit-filter map) from the GMRES variant
#include <kernel/solver/amavanka_gmres.hpp>

#include <kernel/solver/backend_aux/richardson_apply_helper.hpp>

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

#include <vector>

namespace FEAT
{
  namespace Solver
  {
    /// device entry points, implemented in backend_aux/amavanka_richardson.cu
    namespace VankaRichardson
    {
      template<typename DT_, typename IT_, int NC_>
      FEAT::Solver::Status cuda_apply_amavanka_richardson_seq(
        const RichardsonHelper::RichardsonApplyInfo<DT_, IT_, NC_>& apply_info,
        RichardsonHelper::RichardsonStatus<DT_>& status);

      template<typename DT_, typename IT_, int NC_>
      FEAT::Solver::Status cuda_apply_amavanka_richardson_parallel(
        const RichardsonHelper::RichardsonApplyInfo<DT_, IT_, NC_>& apply_info,
        RichardsonHelper::RichardsonStatus<DT_>& status);
    }

    namespace Intern
    {
      // block-size deduction: DenseVectorBlocked exposes ::block_size, plain DenseVector does not (=> 1)
      template<typename, typename = void>
      static constexpr int ricci_dim_vector = 1;
      template<typename T_>
      static constexpr int ricci_dim_vector<T_, std::void_t<decltype(T_::block_size)>> = T_::block_size;
    }

    /**
     * \brief Fused Richardson(Schwarz(AmaVanka)) smoother specialization for a global (meta) system.
     *
     * This solver reproduces the classic \c Richardson(omega) preconditioned by
     * \c Schwarz(AmaVanka) smoother stack, but fuses the whole per-step update
     * \f$ u \gets u + \omega\, f(V\, f(b - A u)) \f$ into a single hand-synchronized device driver.
     *
     * In contrast to \c AmaVankaGMRES (which flattens the meta system into one scalar CSR), this
     * solver walks the \c (k x l) BCSR meta blocks natively via
     * LAFEM::Arch::MatVecMultBCSRDense on device, so no scalar-CSR copy of the system is built.
     *
     * \note Phase 1 limitations (see helper/gotchas): only \c UnitFilter / \c UnitFilterBlocked DOFs
     * are enforced on the device path (Slip/Mean filters are dropped, as in \c AmaVankaGMRES), and
     * only the fixed-step mode (\c min_iter == \c max_iter, no residual norm) is supported on device.
     * The generic (CPU) path applies the full system filter and is used as the reference.
     *
     * \author Maximilian Esser
     */
    template<typename Matrix_, typename Filter_>
    class AmaVankaRichardson;

    template<
      typename LocalMatrix_,
      typename LocalFilter_,
      typename MirrorRow_,
      typename MirrorCol_>
    class AmaVankaRichardson<Global::Matrix<LocalMatrix_, MirrorRow_, MirrorCol_>, Global::Filter<LocalFilter_, MirrorRow_>> :
      public IterativeSolver<typename Global::Matrix<LocalMatrix_, MirrorRow_, MirrorCol_>::VectorTypeR>
    {
    public:
      typedef Global::Matrix<LocalMatrix_, MirrorRow_, MirrorCol_> MatrixType;
      typedef Global::Filter<LocalFilter_, MirrorRow_> FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef typename LocalMatrix_::IndexType IndexType;
      typedef IterativeSolver<VectorType> BaseClass;
      typedef Solver::AmaVanka<LocalMatrix_, LocalFilter_> VankaType;

      /// number of meta components (== number of meta row/col blocks of the square system)
      static constexpr int num_comp = LocalMatrix_::num_row_blocks;
      static_assert(LocalMatrix_::num_row_blocks == LocalMatrix_::num_col_blocks,
        "AmaVankaRichardson requires a square meta system");

    protected:
      /// the system matrix
      const MatrixType& _system_matrix;
      /// the system filter
      const FilterType& _system_filter;
      /// the communicator
      const Dist::Comm* _comm;
      /// the AmaVanka preconditioner (kept alive; its stored matrix is referenced directly)
      std::shared_ptr<VankaType> _vanka;
      /// Richardson damping factor
      DataType _omega;
      /// scratch vectors
      VectorType _vec_def, _vec_cor;

      bool _use_cuda;
      bool _always_use_parallel_kernel;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix   the (global) system matrix
       * \param[in] filter   the (global) system filter
       * \param[in] vanka    the AmaVanka preconditioner (its inverse matrix is used as V)
       * \param[in] omega    the Richardson damping factor
       */
      explicit AmaVankaRichardson(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<VankaType> vanka, DataType omega = DataType(1)) :
        BaseClass("AmaVankaRichardson"),
        _system_matrix(matrix),
        _system_filter(filter),
        _comm(matrix.get_comm()),
        _vanka(std::move(vanka)),
        _omega(omega),
        _vec_def(),
        _vec_cor(),
        _use_cuda(false),
        _always_use_parallel_kernel(false)
      {
        this->_set_comm_by_matrix(matrix);
      }

      virtual ~AmaVankaRichardson() {}

      virtual String name() const override { return "AmaVankaRichardson"; }

      virtual void set_use_cuda(bool use) { this->_use_cuda = use; }
      virtual void set_use_parallel_kernel(bool use) { this->_always_use_parallel_kernel = use; }
      virtual void set_omega(DataType omega) { this->_omega = omega; }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        this->_vanka->init_symbolic();
        this->_vec_def = this->_system_matrix.create_vector_r();
        this->_vec_cor = this->_system_matrix.create_vector_r();
        this->_vec_def.format();
        this->_vec_cor.format();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        Util::ScopedNvtxRange nvtx("Vanka: Init Numeric");
        this->_vanka->init_numeric();
      }

      virtual void done_numeric() override
      {
        this->_vanka->done_numeric();
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        this->_vec_cor.clear();
        this->_vec_def.clear();
        this->_vanka->done_symbolic();
        BaseClass::done_symbolic();
      }

      /// \copydoc SolverBase::apply()
      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        vec_cor.format();
        this->_status = _apply_intern_dispatch(vec_cor, vec_def);
        this->plot_summary();
        return this->_status;
      }

      /// \copydoc IterativeSolver::correct()
      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        this->_status = _apply_intern_dispatch(vec_sol, vec_rhs);
        this->plot_summary();
        return this->_status;
      }

    protected:
      Status _apply_intern_dispatch(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        #ifdef FEAT_HAVE_CUDA
        if(this->_use_cuda)
          return this->_cuda_apply_intern(vec_sol, vec_rhs);
        #endif
        return this->_apply_intern(vec_sol, vec_rhs);
      }

      /**
       * \brief Generic (CPU/reference) Richardson(Schwarz(AmaVanka)) loop.
       *
       * This is the algorithmic ground truth the device driver mirrors; it applies the \b full system
       * filter (including Slip/Mean) and is what the fused CUDA path is validated against on
       * pure-Dirichlet systems.
       */
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        VectorType& vec_def(this->_vec_def);
        VectorType& vec_cor(this->_vec_cor);
        const auto& vanka_mat = this->_vanka->get_vanka_matrix();

        const Index steps = Math::max(this->_max_iter, Index(1));
        for(Index step(0); step < steps; ++step)
        {
          // d = rhs - A*sol   (global apply includes the type-0 synchronization)
          this->_system_matrix.apply(vec_def, vec_sol, vec_rhs, -DataType(1));
          this->_system_filter.filter_def(vec_def);

          // cor = V*d  (local AmaVanka apply) + Schwarz sync_1 + correction filter
          vanka_mat.apply(vec_cor.local(), vec_def.local());
          vec_cor.sync_1();
          this->_system_filter.filter_cor(vec_cor);

          // u += omega * cor
          vec_sol.axpy(vec_cor, this->_omega);
        }

        this->_num_iter = steps;
        this->_def_cur = this->_def_init = DataType(0);
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
        return Status::success;
      }

#ifdef FEAT_HAVE_CUDA
      /// fills a MetaBlock from a concrete SparseMatrixBCSR (any block size)
      template<int bh_, int bw_>
      static void _fill_block(RichardsonHelper::MetaBlock<DataType, IndexType>& mb,
        const LAFEM::SparseMatrixBCSR<DataType, IndexType, bh_, bw_>& m)
      {
        if(m.num_nzes() == Index(0))
          return; // leave the block empty
        // ensure this block owns a persistent cuSPARSE SpMV workspace (built once, reused across
        // applies; required by the cuSPARSE BSR SpMV path on CUDA >= 13) and hand it to the driver
        if(!m.get_workspace(PreferredBackend::cuda))
          m.init_workspace();
        mb.workspace = m.get_workspace(PreferredBackend::cuda).get();
        mb.val = &m.val_arbiter();
        mb.row = &m.row_ptr_arbiter();
        mb.col = &m.col_idx_arbiter();
        mb.block_height = bh_;
        mb.block_width = bw_;
        mb.num_rows = m.num_rows();
        mb.num_cols = m.num_cols();
        mb.num_nnz = m.num_nzes();
      }

      /// extract the system matrix (a SaddlePointMatrix: blocks A=(0,0), B=(0,1), D=(1,0), zero (1,1))
      template<typename MA_, typename MB_, typename MD_>
      static void _extract_matrix(RichardsonHelper::MetaMatrix<DataType, IndexType, num_comp, num_comp>& out,
        const LAFEM::SaddlePointMatrix<MA_, MB_, MD_>& mat)
      {
        static_assert(num_comp == 2, "saddle-point extraction expects two components");
        _fill_block(out.at(0, 0), mat.block_a());
        _fill_block(out.at(0, 1), mat.block_b());
        _fill_block(out.at(1, 0), mat.block_d());
        // (1,1) stays empty
      }

      /// extract the AmaVanka matrix (a full 2x2 TupleMatrix with an empty (1,1) block)
      template<typename Row0_, typename Row1_>
      static void _extract_matrix(RichardsonHelper::MetaMatrix<DataType, IndexType, num_comp, num_comp>& out,
        const LAFEM::TupleMatrix<Row0_, Row1_>& mat)
      {
        static_assert(num_comp == 2, "tuple extraction expects two components");
        _fill_block(out.at(0, 0), mat.template at<0, 0>());
        _fill_block(out.at(0, 1), mat.template at<0, 1>());
        _fill_block(out.at(1, 0), mat.template at<1, 0>());
        _fill_block(out.at(1, 1), mat.template at<1, 1>());
      }

      /// recursively fill per-component vector arbiters, layout, filters and gates
      template<int C_, typename LocGate_>
      void _fill_components(RichardsonHelper::RichardsonApplyInfo<DataType, IndexType, num_comp>& info,
        VectorType& vec_sol, const VectorType& vec_rhs, const LocGate_* gate)
      {
        if constexpr(C_ < num_comp)
        {
          auto& sc = vec_sol.local().template at<C_>();
          const auto& rc = vec_rhs.local().template at<C_>();
          auto& dc = this->_vec_def.local().template at<C_>();
          auto& cc = this->_vec_cor.local().template at<C_>();
          constexpr int bs = Intern::ricci_dim_vector<std::decay_t<decltype(sc)>>;

          info.components[C_] = RichardsonHelper::MetaVectorInfo{bs, sc.size()};
          info.vec_sol[C_] = &sc.elements_arbiter();
          info.vec_rhs[C_] = &rc.elements_arbiter();
          info.vec_def[C_] = &dc.elements_arbiter();
          info.vec_cor[C_] = &cc.elements_arbiter();

          // gate for this component (empty if no row gate)
          auto& gw = info.gates[C_];
          gw.comm = const_cast<Dist::Comm*>(this->_comm);
          if(gate != nullptr)
          {
            const auto& ranks = gate->get_ranks();
            const auto& mirrors = gate->get_mirrors();
            gw.ranks.assign(ranks.begin(), ranks.end());
            gw.freqs = &gate->get_freqs().template at<C_>().elements_arbiter();
            for(std::size_t k(0); k < ranks.size(); ++k)
            {
              const auto& sub = mirrors[k].template at<C_>();
              gw.mirror_idx.push_back(&sub.indices_arbiter());
              gw.mirror_size.push_back(sub.num_indices());
            }
          }
          _fill_components<C_ + 1>(info, vec_sol, vec_rhs, gate);
        }
      }

      virtual Status _cuda_apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        RichardsonHelper::RichardsonApplyInfo<DataType, IndexType, num_comp> info;
        _extract_matrix(info.mat_sys, this->_system_matrix.local());
        _extract_matrix(info.mat_vanka, this->_vanka->get_vanka_matrix());

        // per-component unit-defect filters
        _extract_filters<0>(info);

        // per-component vectors + gates
        const auto* gate = this->_system_matrix.get_row_gate();
        _fill_components<0>(info, vec_sol, vec_rhs, gate);

        info.omega = this->_omega;
        info.min_iter = this->_min_iter;
        info.max_iter = Math::max(this->_max_iter, Index(1));
        info.always_use_parallel_kernel = this->_always_use_parallel_kernel;

        RichardsonHelper::RichardsonStatus<DataType> stat;
        if((this->_comm == nullptr || this->_comm->size() == 1) && !this->_always_use_parallel_kernel)
          this->_status = VankaRichardson::cuda_apply_amavanka_richardson_seq(info, stat);
        else
          this->_status = VankaRichardson::cuda_apply_amavanka_richardson_parallel(info, stat);

        this->_num_iter = stat.num_iter;
        this->_def_cur = stat.def_norm;
        this->_def_init = stat.def_init;
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), this->_status, stat.num_iter));
        return this->_status;
      }

      /// recursively build the flat per-component unit filters from the (meta) system filter
      template<int C_>
      void _extract_filters(RichardsonHelper::RichardsonApplyInfo<DataType, IndexType, num_comp>& info)
      {
        if constexpr(C_ < num_comp)
        {
          auto& vc = this->_vec_def.local().template at<C_>();
          std::map<IndexType, DataType> fmap;
          Intern::FilterHelper<std::decay_t<decltype(this->_system_filter.local().template at<C_>())>>::add_filter(
            fmap, this->_system_filter.local().template at<C_>(), vc);
          // stash the flattened unit filter so its arbiters stay alive for the whole apply
          auto& uf = _unit_filters[C_];
          uf = LAFEM::UnitFilter<DataType, IndexType>(vc.size_raw(), fmap);
          info.def_filters[C_] = RichardsonHelper::UnitFilterInfo<DataType, IndexType>{
            &uf.elements_arbiter(), &uf.indices_arbiter(), uf.num_nzes_raw()};
          _extract_filters<C_ + 1>(info);
        }
      }

      /// flattened per-component unit filters (kept alive across the device apply)
      std::array<LAFEM::UnitFilter<DataType, IndexType>, std::size_t(num_comp)> _unit_filters;
#endif // FEAT_HAVE_CUDA
    }; // class AmaVankaRichardson<...>

    /**
     * \brief Creates a new AmaVankaRichardson smoother object
     *
     * \param[in] matrix   the (global) system matrix
     * \param[in] filter   the (global) system filter
     * \param[in] vanka    the AmaVanka preconditioner
     * \param[in] omega    the Richardson damping factor
     */
    template<typename Matrix_, typename Filter_, typename Vanka_>
    inline std::shared_ptr<AmaVankaRichardson<Matrix_, Filter_>> new_amavanka_richardson(
      const Matrix_& matrix, const Filter_& filter, std::shared_ptr<Vanka_> vanka,
      typename Matrix_::DataType omega = typename Matrix_::DataType(1))
    {
      return std::make_shared<AmaVankaRichardson<Matrix_, Filter_>>(matrix, filter, std::move(vanka), omega);
    }
  } // namespace Solver
} // namespace FEAT
