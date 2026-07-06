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
#include <kernel/lafem/arch/batched_vector_apply.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/mean_filter_blocked.hpp>

#include <kernel/solver/backend_aux/csr_apply_helper.hpp>

#include <kernel/solver/krylov_orthogonalizer.hpp>

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

// includes, FEAT
#include <vector>

namespace FEAT
{
  namespace Solver
  {
    namespace VankaGMRES
    {
      template<typename DT_, typename IT_>
      FEAT::Solver::Status cuda_apply_amavanka_gmres_seq(Memory::Arbiter& vec_sol, const Memory::Arbiter& vec_rhs, const CSRHelper::GMRESApplyInfo<DT_, IT_>& apply_info, CSRHelper::GMRESStatus<DT_>& status);

      template<typename DT_, typename IT_>
      FEAT::Solver::Status cuda_apply_amavanka_gmres_parallel(Memory::Arbiter& vec_sol, const Memory::Arbiter& vec_rhs, const CSRHelper::GMRESApplyInfo<DT_, IT_>& apply_info, CSRHelper::GMRESStatus<DT_>& status);
    }

    namespace Intern
    {
      template<typename, typename = void>
      static constexpr int dim_vector = 1;

      template<typename T>
      static constexpr int dim_vector<T,
      std::void_t<decltype(T::block_size)>> = T::block_size;

      template<typename InMirror_>
      class MirrorHelper;

      template<typename First_, typename... Rest_>
      class MirrorHelper<LAFEM::TupleMirror<First_, Rest_...>>
      {
      public:
        typedef LAFEM::TupleMirror<First_, Rest_...> InMirrorType;

        template<typename OutMirror_, typename FVec_, typename RVec_>
        static std::pair<Index, Index> fill_mirror(OutMirror_& mirror, const InMirrorType& input, const LAFEM::TupleVector<FVec_, RVec_>& tvec, std::pair<Index, Index> offsets = {0u, 0u})
        {
          offsets = MirrorHelper<First_>::fill_mirror(mirror, input.first(), tvec.first(), offsets);
          return MirrorHelper<LAFEM::TupleMirror<Rest_...>>::fill_mirror(mirror, input.rest(), tvec.rest(), offsets);
        }
      };

      template<typename First_>
      class MirrorHelper<LAFEM::TupleMirror<First_>>
      {
      public:
        typedef LAFEM::TupleMirror<First_> InMirrorType;

        template<typename OutMirror_, typename FVec_>
        static std::pair<Index, Index> fill_mirror(OutMirror_& mirror, const InMirrorType& input, const LAFEM::TupleVector<FVec_>& tvec, std::pair<Index, Index> offsets = {0u, 0u})
        {
          offsets = MirrorHelper<First_>::fill_mirror(mirror, input.first(), tvec.first(), offsets);
          return offsets;
        }
      };

      template<typename DT_, typename IT_>
      class MirrorHelper<LAFEM::VectorMirror<DT_, IT_>>
      {
      public:
        typedef LAFEM::VectorMirror<DT_, IT_> InMirrorType;

        template<typename OutMirror_, typename Vec_>
        static std::pair<Index, Index> fill_mirror(OutMirror_& mirror, const InMirrorType& input, const Vec_& vec, std::pair<Index, Index> offsets = {0u, 0u})
        {
          // get indices (host-side views into the mirror index arrays)
          const auto indices_view = input.indices_view_r();
          const auto* indices = indices_view.get_r();
          const Index& dof_offset = offsets.first;
          const Index& ptr_offset = offsets.second;
          auto target_view = mirror.indices_view_w();
          auto* target = target_view.get_w() + ptr_offset;
          constexpr int dim = dim_vector<Vec_>;
          for(Index k = 0; k < input.num_indices(); ++k)
          {
            for(int d = 0; d < dim; ++d)
            {
              target[k*dim + Index(d)] = IT_(indices[k]*dim + Index(d) + dof_offset);
            }
          }
          offsets.first += vec.size_raw();
          offsets.second += input.num_indices()*dim;
          return offsets;
        }
      };

      template<typename InFilterType_>
      class FilterHelper
      {
      public:
        template<typename DT_, typename IT_, typename Vec_>
        static Index add_filter(std::map<IT_, DT_>&, const InFilterType_& DOXY(filter), const Vec_& vec, Index offset = 0u)
        {
          // XABORTM("Not implemented");
          return offset + vec.size_raw();
        }
      };

      template<typename First_, typename... Rest_>
      class FilterHelper<LAFEM::TupleFilter<First_, Rest_...>>
      {
      public:
        typedef LAFEM::TupleFilter<First_, Rest_...> InFilterType;
        /// by default, add nothing from the filter
        template<typename DT_, typename IT_, typename FVec_, typename RVec_>
        static Index add_filter(std::map<IT_, DT_>& filter_map, const InFilterType& input, const LAFEM::TupleVector<FVec_, RVec_>& vec, Index offset = 0u)
        {
          offset = FilterHelper<First_>::add_filter(filter_map, input.first(), vec.first(), offset);
          offset = FilterHelper<LAFEM::TupleFilter<Rest_...>>::add_filter(filter_map, input.rest(), vec.rest(), offset);
          return offset;
        }
      };

      template<typename First_>
      class FilterHelper<LAFEM::TupleFilter<First_>>
      {
      public:
        typedef LAFEM::TupleFilter<First_> InFilterType;
        /// by default, add nothing from the filter
        template<typename DT_, typename IT_, typename FVec_>
        static Index add_filter(std::map<IT_, DT_>& filter_map, const InFilterType& input, const LAFEM::TupleVector<FVec_>& vec, Index offset = 0u)
        {
          offset = FilterHelper<First_>::add_filter(filter_map, input.first(), vec.first(), offset);
          return offset;
        }
      };

      template<typename First_, typename... Rest_>
      class FilterHelper<LAFEM::FilterChain<First_, Rest_...>>
      {
      public:
        typedef LAFEM::FilterChain<First_, Rest_...> InFilterType;
        /// by default, add nothing from the filter
        template<typename DT_, typename IT_, typename Vec_>
        static Index add_filter(std::map<IT_, DT_>& filter_map, const InFilterType& input, const Vec_& vec, Index offset = 0u)
        {
          FilterHelper<First_>::add_filter(filter_map, input.first(), vec, offset);
          return FilterHelper<LAFEM::FilterChain<Rest_...>>::add_filter(filter_map, input.rest(), vec, offset);
        }
      };

      template<typename First_>
      class FilterHelper<LAFEM::FilterChain<First_>>
      {
      public:
        typedef LAFEM::FilterChain<First_> InFilterType;
        /// by default, add nothing from the filter
        template<typename DT_, typename IT_, typename Vec_>
        static Index add_filter(std::map<IT_, DT_>& filter_map, const InFilterType& input, const Vec_& vec, Index offset = 0u)
        {
          return FilterHelper<First_>::add_filter(filter_map, input.first(), vec, offset);
        }
      };

      template<typename FType_>
      class FilterHelper<LAFEM::FilterSequence<FType_>>
      {
      public:
        typedef LAFEM::FilterChain<FType_> InFilterType;
        /// by default, add nothing from the filter
        template<typename DT_, typename IT_, typename Vec_>
        static Index add_filter(std::map<IT_, DT_>& filter_map, const InFilterType& input, const Vec_& vec, Index offset = 0u)
        {
          Index new_offset = offset;
          for(auto& pair : input)
          {
            new_offset = FilterHelper<FType_>::add_filter(filter_map, pair.second, vec, offset);
          }
          return new_offset;
        }
      };

      template<typename DTT_, typename ITT_>
      class FilterHelper<LAFEM::UnitFilter<DTT_, ITT_>>
      {
      public:
        typedef LAFEM::UnitFilter<DTT_, ITT_> InFilterType;
        /// by default, add nothing from the filter
        template<typename DT_, typename IT_>
        static Index add_filter(std::map<IT_, DT_>& filter_map, const InFilterType& input, const LAFEM::DenseVector<DTT_, ITT_>& vec, Index offset = 0u)
        {
          if(input.num_nzes() > 0u)
          {
            const auto& sv = input.get_filter_vector();
            const auto idx_view = sv.indices_view_r();
            const auto val_view = sv.elements_view_r();
            const auto* svi = idx_view.get_r();
            const auto* svd = val_view.get_r();
            for(Index k = 0; k < sv.num_nzes(); ++k)
            {
              filter_map.insert_or_assign(IT_(svi[k]+offset), DT_(svd[k]));
            }
          }
          return offset+vec.size();
        }
      };

      template<typename DTT_, typename ITT_, int dim_>
      class FilterHelper<LAFEM::UnitFilterBlocked<DTT_, ITT_, dim_>>
      {
      public:
        typedef LAFEM::UnitFilterBlocked<DTT_, ITT_, dim_> InFilterType;
        /// by default, add nothing from the filter
        template<typename DT_, typename IT_>
        static Index add_filter(std::map<IT_, DT_>& filter_map, const InFilterType& input, const LAFEM::DenseVectorBlocked<DTT_, ITT_, dim_>& vec, Index offset = 0u)
        {
          if(input.num_nzes() > 0u)
          {
            const auto& sv = input.get_filter_vector();
            const auto idx_view = sv.indices_view_r();
            const auto val_view = sv.elements_view_r();
            const auto* svi = idx_view.get_r();
            const auto* svd = val_view.get_r();
            for(Index k = 0; k < sv.num_nzes(); ++k)
            {
              for(int d = 0; d < dim_; ++d)
                filter_map.insert_or_assign(IT_(dim_)*svi[k]+IT_(d)+IT_(offset), DT_(svd[k][d]));
            }
          }
          return offset+vec.size_raw();
        }
      };
    }
    /**
     * \brief GMRES(k) solver implementation specialization for an Stokes typed GMRES(AmaVanka) implementation
     *
     * This GMRES variant is specialized with a vanka type preconditionier matrix and scalarizes the internal system
     * which allows for very specific optimizations in combination with a kernel fusing of the main loop.
     * Since this approach is memory heavy, it is recommended to only use this solver as a coarse grid solver for
     * the case of many gpus or very small local system sizes.
     * In alle other cases, the standard GMRES with the batched dot variants should be sufficient.
     *
     * \see
     * Y. Saad: A flexible inner-outer preconditioned GMRES algorithm;
     * SIAM Journal on Scientific Computing, Volume 14 Issue 2, pp. 461-469, March 1993
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Maximilian Esser
     */
    template<
      typename Matrix_,
      typename Filter_>
    class AmaVankaGMRES;

    template<
      typename LocalMatrix_,
      typename LocalFilter_,
      typename MirrorRow_,
      typename MirrorCol_>
    class AmaVankaGMRES<Global::Matrix<LocalMatrix_, MirrorRow_, MirrorCol_>, Global::Filter<LocalFilter_, MirrorRow_>>:
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

      // the internal types used
      typedef LAFEM::SparseMatrixCSR<DataType, IndexType> SparseMatType;
      typedef LAFEM::VectorMirror<DataType, IndexType> MirrorType;
      typedef Global::Matrix<SparseMatType, MirrorType, MirrorType> GlobalMatType;
      typedef LAFEM::DenseVector<DataType, IndexType> VecType;
      typedef LAFEM::UnitFilter<DataType, IndexType> FilType;
      typedef Global::Vector<VecType, MirrorType> GlobalVecType;
      typedef Global::Gate<VecType, MirrorType> GateType;


    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// the communicator used
      const Dist::Comm* _comm;
      /// The Amavanka preconditionier we copy our preconditionier matrix from
      std::shared_ptr<VankaType> _vanka;
      /// Our internal gate
      GateType _gate;
      /// Our internal system matrix
      GlobalMatType _internal_system_mat;
      /// Our internal vanka matrix
      SparseMatType _local_vanka_mat;
      /// our local filter object
      FilType _internal_filter;


      /// krylov dimension
      Index _krylov_dim;
      /// inner pseudo-residual scaling factor, see the documentation of this class for details
      DataType _inner_res_scale;
      /// the orthogonalization strategy used to build the Krylov basis; this specialization only
      /// supports the batched classical Gram-Schmidt (\c classical_gram_schmidt) and its (conditional)
      /// reorthogonalization variant (\c cgs2)
      KrylovOrthogonalizerType _ortho_type = KrylovOrthogonalizerType::classical_gram_schmidt;
      /// DGKS norm-ratio threshold below which a corrective reorthogonalization sweep is triggered
      /// (only relevant for cgs2); see CGS2Orthogonalizer
      DataType _reortho_tol = DataType(0.3);
      /// if set, always perform the corrective reorthogonalization sweep in cgs2 mode
      bool _force_reortho = false;
      /// krylov basis vector
      std::vector<GlobalVecType> _vec_v;
      // std::vector<VecType> _vec_v_helper;
      // arbiters of the krylov basis vectors (the kernel acquires its device views from them)
      std::vector<Memory::Arbiter*> _vec_v_elements;
      GlobalVecType _vec_w;
      /// helper vectors
      GlobalVecType _vec_sol, _vec_rhs;
      // VectorType _helper_v;
      /// Givens rotation coefficients
      std::vector<DataType> _c, _s, _q;
      /// Hessenberg matrix
      std::vector<std::vector<DataType>> _h;
      /// scratch for the cgs2 corrective sweep's dot products, kept separate from the accumulated h[k]
      std::vector<DataType> _corr;
      /// device-resident accumulation buffer for the scalar-path batched dots (the kernels acquire
      /// their views from this arbiter)
      Memory::Arbiter _dots_arb;

      Index _global_vec_size;
      Index _global_used_elements;
      Index _internal_flop_counter;

      bool _use_cuda_gmres;
      bool _prefetch_data;
      bool _always_use_parallel_kernel;

      // void _construct_filter(FilType& internal_filter, const LocalFilter_& sys_filter)
      // {
      //   // construct Local vector
      //   VecType loc_tmp;
      //   {
      //     auto tmp_vec = this->_system_matrix.local().create_vector_r();
      //     tmp_vec.format(DataType(1.8));
      //     sys_filter.filter_def(tmp_vec);
      //     loc_tmp.convert(tmp_vec);
      //   }

      //   internal_filter = FilType(loc_tmp.size());

      //   const auto* val = loc_tmp.elements();
      //   for(Index k = 0; k < loc_tmp.size(); ++k)
      //   {
      //     if(Math::abs(loc_tmp[k]) < DataType(0.0001))
      //     {
      //       internal_filter.add(k, DataType(0));
      //     }
      //   }
      // }

    public:

      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] krylov_dim
       * The maximum Krylov subspace dimension. Must be > 0.
       *
       * \param[in] inner_res_scale
       * The scaling factor for the inner GMRES loop residual. See the documentation of this class for details.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       */
      explicit AmaVankaGMRES(const MatrixType& matrix, const FilterType& filter, std::shared_ptr<VankaType> vanka, Index krylov_dim,
        DataType inner_res_scale = DataType(0)) :
        BaseClass("AmaVankaGMRES(" + stringify(krylov_dim) + ")"),
        _system_matrix(matrix),
        _system_filter(filter),
        _comm(_system_matrix.get_comm()),
        _vanka(std::move(vanka)),
        _gate(),
        _internal_system_mat(),
        _local_vanka_mat(),
        _internal_filter(),
        _krylov_dim(krylov_dim),
        _global_vec_size(),
        _global_used_elements(),
        _internal_flop_counter(0u),
        _use_cuda_gmres(false),
        _prefetch_data(false),
        _always_use_parallel_kernel(false)
      {
        if(_comm)
          _gate.set_comm(_comm);
        // we always need to compute defect norms
        this->_force_def_norm_calc = true;

        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
        set_inner_res_scale(inner_res_scale);

        const Index local_mat_rows = _system_matrix.local().num_rows_raw();
        const Index local_mat_cols = _system_matrix.local().num_cols_raw();

        // construct our internal gate
        if(matrix.get_row_gate())
        {
          const auto& row_gate = *matrix.get_row_gate();
          const auto& ranks = row_gate.get_ranks();
          const auto& mirrors = row_gate.get_mirrors();
          const auto& freqs = row_gate.get_freqs();
          for(Index k = 0; k < ranks.size(); ++k)
          {
            //construct mirror
            const Index num_indices = mirrors.at(k).buffer_size(freqs);
            MirrorType mirror(local_mat_rows, num_indices);
            auto offsets = Intern::MirrorHelper<MirrorRow_>::fill_mirror(mirror, mirrors.at(k), freqs);
            XASSERTM(offsets.first == local_mat_rows, "Offset does not match overall mirror size");
            XASSERTM(offsets.second == num_indices, "Offset does not match overall mirror indices");

            this->_gate.push(ranks.at(k), std::move(mirror));
          }
          // convert frequency
          VecType tmp_freq;
          tmp_freq.convert(freqs);
          this->_gate.compile(std::move(tmp_freq));

          // our internal global matrix uses the new gate
          this->_internal_system_mat = GlobalMatType(&this->_gate, &this->_gate, local_mat_rows, local_mat_cols);
          this->_internal_filter = FilType(local_mat_rows, 0u);
        }
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~AmaVankaGMRES()
      {
      }

      /// \copydoc BaseClass::name()
      virtual String name() const override
      {
        return "AmaVankaGMRES";
      }

      virtual void set_use_cuda_gmres(bool use)
      {
        this->_use_cuda_gmres = use;
      }

      virtual void set_prefetch_data(bool use)
      {
        this->_prefetch_data = use;
      }

      virtual void set_use_parallel_kernel(bool use)
      {
        this->_always_use_parallel_kernel = use;
      }
      /**
       * \brief Forces the calculation of defect norms in every iteration (overridden)
       *
       * This solver always requires the calculation of defect norms, so this function has been
       * overridden and calling it has no effect.
       */
      virtual void force_defect_norm_calc(bool DOXY(force)) override
      {
        // force is already applied in constructor
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // initialize internal matrix representation, for this, first initialize vanka preconditioner
        this->_vanka->init_symbolic();

        // convert the global matrix
        this->_internal_system_mat.local().convert(_system_matrix.local());

        {
          // create temp helper vector to assemble filter
          auto tmp_sys_vec = _system_matrix.local().create_vector_r();
          std::map<IndexType, DataType> tmp_map;
          // convert the filter we assume we only require unit defect filter
          const Index in_size = Intern::FilterHelper<std::decay_t<decltype(this->_system_filter.local())>>::add_filter(tmp_map, this->_system_filter.local(), tmp_sys_vec);
          this->_internal_filter = FilType(in_size, tmp_map);
        }


        _c.reserve(_krylov_dim);
        _s.reserve(_krylov_dim);
        _q.reserve(_krylov_dim + 1);
        _h.resize(_krylov_dim);

        for(Index i(0); i < _krylov_dim; ++i)
        {
          _h.at(i).resize(i+1);
        }
        _corr.resize(_krylov_dim);
        _dots_arb = Memory::Arbiter((_krylov_dim + 1) * sizeof(DataType));

        _vec_v.reserve(_krylov_dim+1);
        _vec_v_elements.reserve(_krylov_dim+1);

        _vec_v.emplace_back(this->_internal_system_mat.create_vector_r());
        _vec_v.back().format();
        // _vec_v_helper.emplace_back(_vec_v.back().local().clone(LAFEM::CloneMode::Shallow));
        _vec_v_elements.emplace_back(&_vec_v.back().local().elements_arbiter());
        _vec_w = this->_internal_system_mat.create_vector_r();
        _vec_w.format();
        for(Index i(0); i < _krylov_dim; ++i)
        {
          _vec_v.emplace_back(this->_vec_v.front().clone(LAFEM::CloneMode::Layout));
          _vec_v.back().format();
          // _vec_v_helper.emplace_back(_vec_v.back().local().clone(LAFEM::CloneMode::Shallow));
          _vec_v_elements.emplace_back(&_vec_v.back().local().elements_arbiter());
        }
        _vec_sol = this->_internal_system_mat.create_vector_l();
        _vec_rhs = this->_internal_system_mat.create_vector_r();

        this->_global_used_elements = this->_internal_system_mat.num_nzes_raw();
        this->_global_vec_size = _vec_sol.size();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        this->_vanka->init_numeric();
        this->_local_vanka_mat.convert(this->_vanka->get_vanka_matrix());
        this->_vanka->done_numeric();
        this->_internal_system_mat.local().convert(_system_matrix.local());

        this->_local_vanka_mat.init_workspace();
        this->_internal_system_mat.local().init_workspace();
        this->_dots_arb.format(Memory::Location::main);

      }

      virtual void done_symbolic() override
      {
        this->_global_used_elements = 0u;
        this->_global_vec_size = 0u;
        _vec_sol.clear();
        _vec_rhs.clear();
        _vec_v_elements.clear();
        // _vec_v_helper.clear();
        _vec_v.clear();
        _vec_w.clear();
        this->_internal_filter.clear();
        this->_internal_system_mat.local().clear();
        this->_vanka->done_symbolic();
        BaseClass::done_symbolic();
      }

      virtual void done_numeric() override
      {
        this->_internal_system_mat.local().reset_workspace();
        this->_local_vanka_mat.reset_workspace();
        this->_local_vanka_mat.clear();
      }

      /**
       * \brief Sets the inner Krylov space dimension
       *
       * \param[in] krylov_dim
       * The k in GMRES(k)
       *
       */
      virtual void set_krylov_dim(Index krylov_dim)
      {
        XASSERT(krylov_dim > Index(0));
        _krylov_dim = krylov_dim;
      }

      /**
       * \brief Sets the inner residual scale
       *
       * \param[in] inner_res_scale
       * Scaling parameter for convergence of the inner iteration
       */
      virtual void set_inner_res_scale(DataType inner_res_scale)
      {
        XASSERT(inner_res_scale >= DataType(0));
        _inner_res_scale = inner_res_scale;
      }

      /**
       * \brief Sets the Krylov subspace orthogonalization strategy
       *
       * \param[in] ortho_type
       * The orthogonalization strategy to use. This specialization only supports the batched
       * classical Gram-Schmidt (\c classical_gram_schmidt) and its (conditional) reorthogonalization
       * variant (\c cgs2); passing \c modified_gram_schmidt is not supported.
       *
       * \see KrylovOrthogonalizerType
       */
      virtual void set_orthogonalizer(KrylovOrthogonalizerType ortho_type)
      {
        XASSERTM(ortho_type == KrylovOrthogonalizerType::classical_gram_schmidt
          || ortho_type == KrylovOrthogonalizerType::cgs2,
          "AmaVankaGMRES only supports classical_gram_schmidt or cgs2");
        _ortho_type = ortho_type;
      }

      /**
       * \brief Sets the DGKS norm-ratio threshold that triggers a corrective reorthogonalization sweep
       *
       * Only relevant in \c cgs2 mode. See CGS2Orthogonalizer.
       */
      virtual void set_reortho_tol(DataType tol)
      {
        XASSERT(tol >= DataType(0));
        _reortho_tol = tol;
      }

      /**
       * \brief Forces (or unforces) a corrective reorthogonalization sweep in every step
       *
       * Only relevant in \c cgs2 mode. See CGS2Orthogonalizer.
       */
      virtual void set_force_reortho(bool force)
      {
        _force_reortho = force;
      }

      /// returns the overall distributed flops performed by the last apply call
      /// do not comm -> reduce this number
      Index get_global_flops() const
      {
        return _internal_flop_counter;
      }

      /// \copydoc IterativeSolver::apply()
      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // marshal the external rhs into the internal vector (replaces the removed set_vec)
        this->_vec_rhs.local().copy(vec_rhs.local());
        this->_vec_v.at(0).local().copy(vec_rhs.local());

        // clear solution vector
        _vec_sol.format();

        // apply
        this->_status = _apply_intern_dispatch(this->_vec_sol, this->_vec_rhs);
        // marshal the internal solution back out (replaces the removed set_vec_inv)
        this->_vec_sol.local().copy_to(vec_sol.local());
        this->plot_summary();
        return this->_status;
      }

      /// \copydoc SolverBase::correct()
      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // marshal the external vectors into the internal ones (replaces the removed set_vec)
        this->_vec_sol.local().copy(vec_sol.local());
        this->_vec_rhs.local().copy(vec_rhs.local());

        this->_internal_system_mat.apply(this->_vec_v.at(0), this->_vec_sol, this->_vec_rhs, -DataType(1));
        this->_internal_filter.filter_def(this->_vec_v.at(0).local());

        // apply
        this->_status = _apply_intern_dispatch(this->_vec_sol, this->_vec_rhs);
        // marshal the internal solution back out (replaces the removed set_vec_inv)
        this->_vec_sol.local().copy_to(vec_sol.local());
        this->plot_summary();
        return this->_status;
      }

    protected:
      static std::uint64_t _start_cuda_range([[maybe_unused]] const char* name)
      {
        std::uint64_t id = 0;
        #ifdef FEAT_HAVE_CUDA
        id = Util::cuda_start_range(name);
        #endif
        return id;
      }

      static void _stop_cuda_range([[maybe_unused]] std::uint64_t id)
      {
        #ifdef FEAT_HAVE_CUDA
        Util::cuda_stop_range(id);
        #endif
      }

      template<typename VecT_>
      static void _calc_dots(std::vector<DataType>& loc_dots, const VecT_& target, const std::vector<Memory::Arbiter>& vec_arb, const Index num_dots)
      {
        const auto id_dot = _start_cuda_range("Calc Dots");
        for(Index i = 0; i < num_dots; ++i)
        {
          loc_dots[i] = DataType(0);
        }
        Intern::template BatchedHelper<VecT_>::batched_dot(loc_dots.data(), target, vec_arb.data(), num_dots, num_dots);
        _stop_cuda_range(id_dot);
      }

      /// \brief Batched dot: computes \p num_dots projection dots into \p loc_dots and returns the
      /// appended self-dot \f$\langle w,w\rangle\f$ in the same batched reduction.
      ///
      /// \p vec_arb holds all \c krylov_dim+1 basis vectors, so slot \p num_dots of the reduction is
      /// the target vector itself and thus yields the self-dot without a separate norm reduction.
      DataType _calc_dots_s(std::vector<DataType>& loc_dots, const VecType& vec_target, const std::vector<Memory::Arbiter*>& vec_arb, const Index num_dots)
      {
        const auto id_dot = _start_cuda_range("Calc Dots");
        // zero the accumulation arbiter (host), run the batched dot (device-accumulates for cuda), read back
        _dots_arb.format();
        LAFEM::Arch::BatchedDotProduct::exec<DataType, DataType>(this->_dots_arb, vec_target.elements_arbiter(), vec_target.size(), vec_arb.data(), num_dots + 1u);
        DataType self_dot(0);
        {
          Memory::TypedView<DataType> dv(this->_dots_arb.view(Memory::Location::main, Memory::Access::read));
          const DataType* dp = dv.get_r();
          for(Index i = 0; i < num_dots; ++i)
            loc_dots[i] = dp[i];
          self_dot = dp[num_dots];
        }
        _stop_cuda_range(id_dot);
        return self_dot;
      }

      /// \brief Distributed batched (frequency-weighted) dot with appended self-dot; the single
      /// all-reduce covers the \p num_dots projections \b and the self-dot together.
      DataType _calc_triple_dots_s(std::vector<DataType>& loc_dots, const VecType& vec_target, const std::vector<Memory::Arbiter*>& vec_arb, const VecType& freqs, const Index num_dots)
      {
        const auto id_dot = _start_cuda_range("Calc Dots");
        _dots_arb.format();
        LAFEM::Arch::BatchedTripleDotProduct::exec<DataType, DataType>(this->_dots_arb, vec_target.elements_arbiter(), vec_target.size(), vec_arb.data(), freqs.elements_arbiter(), num_dots + 1u);
        DataType self_dot(0);
        {
          // acquire the host view (device->host copy if needed) and all-reduce projections + self-dot in one go
          Memory::TypedView<DataType> dv(this->_dots_arb.view(Memory::Location::main, Memory::Access::read_write));
          DataType* dp = dv.get_w();
          this->_comm->allreduce(dp, dp, num_dots + 1u, Dist::op_sum);
          for(Index i = 0; i < num_dots; ++i)
            loc_dots[i] = dp[i];
          self_dot = dp[num_dots];
        }
        _stop_cuda_range(id_dot);
        return self_dot;
      }

      /**
       * \brief Orthogonalizes v[i+1] against v[0..i] and returns its post-orthogonalization norm.
       *
       * Recreates the classical Gram-Schmidt (\c classical_gram_schmidt) / conditional CGS2 (\c cgs2)
       * strategies of KrylovSubspaceOrthogonalizer: a single batched reduction supplies both the
       * projection coefficients (written into \c _h.at(i)) and the appended self-dot, from which the
       * norm is recovered algebraically as \f$\sqrt{\langle w,w\rangle - \sum_k h_k^2}\f$ (with a
       * norm2 fallback on heavy cancellation). In cgs2 mode a corrective sweep is added when forced
       * or when the DGKS criterion \f$\text{nrm\_sqr} < \tau^2\langle w,w\rangle\f$ triggers.
       */
      DataType _orthogonalize_cpu(const Index i)
      {
        VecType& target = this->_vec_v.at(i+1).local();
        const Index num_basis = i + 1u;
        const bool distributed = this->_comm && (this->_comm->size() > 1);

        // first classical Gram-Schmidt sweep: batched projection dots + appended self-dot
        DataType self_dot = distributed
          ? _calc_triple_dots_s(this->_h.at(i), target, this->_vec_v_elements, this->_gate.get_freqs(), num_basis)
          : _calc_dots_s(this->_h.at(i), target, this->_vec_v_elements, num_basis);
        _internal_flop_counter += (num_basis + 1u) * 2u * _global_vec_size;

        DataType sum_sqr(0);
        for(Index k(0); k < num_basis; ++k)
          sum_sqr += Math::sqr(this->_h.at(i).at(k));

        // subtract the projection in a single batched cumulative axpy: target -= sum_k h[k] * v[k]
        _add_vecs_s(target, this->_vec_v_elements, this->_h.at(i), num_basis, DataType(-1));
        _internal_flop_counter += num_basis * 2u * _global_vec_size;

        DataType nrm_sqr = self_dot - sum_sqr;

        // conditional corrective reorthogonalization sweep (cgs2); the second batched reduction
        // supplies both the corrections and the updated self-dot, so no extra reduction is needed
        if((this->_ortho_type == KrylovOrthogonalizerType::cgs2)
          && (this->_force_reortho || (nrm_sqr < Math::sqr(this->_reortho_tol) * self_dot)))
        {
          DataType self_dot2 = distributed
            ? _calc_triple_dots_s(this->_corr, target, this->_vec_v_elements, this->_gate.get_freqs(), num_basis)
            : _calc_dots_s(this->_corr, target, this->_vec_v_elements, num_basis);
          _internal_flop_counter += (num_basis + 1u) * 2u * _global_vec_size;
          DataType sum_sqr2(0);
          for(Index k(0); k < num_basis; ++k)
          {
            this->_h.at(i).at(k) += this->_corr.at(k);
            sum_sqr2 += Math::sqr(this->_corr.at(k));
          }
          _add_vecs_s(target, this->_vec_v_elements, this->_corr, num_basis, DataType(-1));
          _internal_flop_counter += num_basis * 2u * _global_vec_size;
          nrm_sqr = self_dot2 - sum_sqr2;
        }

        // recover the norm algebraically; fall back to an explicit norm2 on heavy cancellation
        // (guarantees a positive result so the subsequent normalization cannot produce a NaN)
        if(nrm_sqr <= Math::Limits<DataType>::epsilon() * DataType(100))
        {
          _internal_flop_counter += 2u * _global_vec_size;
          return this->_vec_v.at(i+1).norm2();
        }
        return Math::sqrt(nrm_sqr);
      }

      template<typename VecT_>
      static void _add_vecs(VecT_& target_vec, const std::vector<Memory::Arbiter*>& vec_arb, const std::vector<DataType>& loc_scales, const Index num_axpy, DataType alpha)
      {
        const auto id_vec_add = _start_cuda_range("Add vecs");
        Intern::template BatchedHelper<VecT_>::commulative_axpy(target_vec, vec_arb.data(), loc_scales.data(), num_axpy, alpha, num_axpy);
        _stop_cuda_range(id_vec_add);
      }

      static void _add_vecs_s(VecType& target_vec, const std::vector<Memory::Arbiter*>& vec_arb, const std::vector<DataType>& loc_scales, const Index num_axpy, DataType alpha)
      {
        const auto id_vec_add = _start_cuda_range("Add vecs");
        LAFEM::Arch::CumulativeAxpy::exec<DataType>(target_vec.elements_arbiter(), target_vec.size(), vec_arb.data(), loc_scales.data(), num_axpy, alpha);
        _stop_cuda_range(id_vec_add);
      }

      Status _set_initial_defectp(const GlobalVecType& vec_def, const GlobalVecType&)
      {
        // reset iteration variables
        this->_def_init = this->_def_prev = this->_def_cur = DataType(0);
        this->_num_iter = Index(0);
        this->_num_stag_iter = Index(0);

        // first, let's see if we have to compute the defect at all
        bool calc_def = this->_force_def_norm_calc
          || (this->_min_iter < this->_max_iter)
          || (this->_plot_iter() || this->_plot_summary())
          || (this->_min_stag_iter > Index(0));

        // if we don't have to compute the initial defect, we can exit here
        if(!calc_def)
          return Status::progress;

        // compute new initial defect
        this->_def_init = this->_def_cur = this->_def_prev = vec_def.norm2();
        Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_init, this->get_num_iter()));
        _internal_flop_counter += 2u * _global_vec_size;

        // plot iteration line?
        if(this->_plot_iter())
          this->_plot_iter_line(Index(0), this->_def_init, this->_def_init);

        // ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(this->_def_init))
          return Status::aborted;

        // check against low absolute tolerance; this may be fulfilled if we have only a handful of dofs
        if(this->_def_init < this->_tol_abs_low)
          return Status::success;

        // check if the initial defect is zero; we test against eps^2 here
        if(this->_def_init <= Math::sqr(Math::eps<DataType>()))
          return Status::success;

        // continue iterating
        return Status::progress;
      }

      Status _set_new_defectp(const GlobalVecType& vec_def, const GlobalVecType&)
      {
        // increase iteration count
        ++this->_num_iter;

        // store previous defect
        this->_def_prev = this->_def_cur;

        // first, let's see if we have to compute the defect at all
        bool calc_def = this->_force_def_norm_calc
          || (this->_min_iter < this->_max_iter)
          || (this->_plot_iter() || this->_plot_summary())
          || (this->_min_stag_iter > Index(0));

        // compute new defect
        if(calc_def)
        {
          this->_def_cur = vec_def.norm2();
          Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_cur, this->get_num_iter()));
          _internal_flop_counter += 2u * _global_vec_size;
        }

        // analyze defect
        Status status = this->_analyze_defect(this->_num_iter, this->_def_cur, this->_def_prev, true);

        // plot defect?
        if(this->_plot_iter(status))
          this->_plot_iter_line(this->_num_iter, this->_def_cur, this->_def_prev);

        // return status
        return status;
      }

      void _apply_vanka(GlobalVecType& v, const GlobalVecType& w)
      {
        this->_local_vanka_mat.apply(v.local(), w.local());
        v.sync_1();
        this->_internal_filter.filter_cor(v.local());
        _internal_flop_counter += 2u * _global_used_elements;
      }

      template<typename DT_, typename IT_>
      CSRHelper::CRSMat<DT_, IT_> _extract_mat(const LAFEM::SparseMatrixCSR<DT_, IT_>& mat)
      {
        // hand over the data arbiters; the kernel acquires its device views from them
        return CSRHelper::CRSMat<DT_, IT_>{&mat.val_arbiter(), &mat.row_ptr_arbiter(), &mat.col_idx_arbiter(), mat.num_cols(), mat.num_cols(), mat.num_nzes()};
      }

      template<typename DT_, typename IT_>
      CSRHelper::UnitFilterInfo<DT_, IT_> _extract_filter(const LAFEM::UnitFilter<DT_, IT_>& filter)
      {
        return CSRHelper::UnitFilterInfo<DT_, IT_>{&filter.elements_arbiter(), &filter.indices_arbiter(), filter.num_nzes_raw()};
      }

      Status _apply_intern_dispatch(GlobalVecType& vec_sol, const GlobalVecType& vec_rhs)
      {
        _internal_flop_counter = 0u;
        #ifdef FEAT_HAVE_CUDA
        if(this->_use_cuda_gmres)
        {
          return this->_cuda_apply_intern(vec_sol, vec_rhs);
        }
        #endif
        return this->_apply_intern(vec_sol, vec_rhs);
      }

      #ifdef FEAT_HAVE_CUDA
      virtual Status _cuda_apply_intern(GlobalVecType& vec_sol, const GlobalVecType& vec_rhs)
      {
        IterationStats pre_iter(*this);
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
        pre_iter.destroy();
        CSRHelper::GMRESApplyInfo<DataType, IndexType> apply_info;
        {
          apply_info.mat_sys = _extract_mat(this->_internal_system_mat.local());
          apply_info.mat_vanka = _extract_mat(this->_local_vanka_mat);
          apply_info.def_filter = _extract_filter(this->_internal_filter);
          apply_info.vec_v = &this->_vec_v_elements;
          apply_info.vec_w = &this->_vec_w.local().elements_arbiter();
          apply_info.gate = &this->_gate;
          apply_info.size = this->_internal_system_mat.local().num_rows();

          apply_info.res_scale = this->_inner_res_scale;
          apply_info.abs_tol = this->_tol_abs;
          apply_info.rel_tol = this->_tol_rel;
          apply_info.krylov_dim = this->_krylov_dim;
          apply_info.max_iter = this->_max_iter;
          apply_info.min_iter = this->_min_iter;
          apply_info.plot_interval = this->_plot_interval;
          apply_info.plot_mode = this->_plot_mode;
          apply_info.prefetch_data = this->_prefetch_data;
          apply_info.ortho_type = this->_ortho_type;
          apply_info.reortho_tol = this->_reortho_tol;
          apply_info.force_reortho = this->_force_reortho;
        }

        CSRHelper::GMRESStatus<DataType> stat;
        {
          stat.def_norm = DataType(0);
          stat.def_init = DataType(0);
          stat.num_iter = Index(0);
        }

        if((this->_comm->size() == 1) && !this->_always_use_parallel_kernel)
          this->_status = VankaGMRES::cuda_apply_amavanka_gmres_seq(vec_sol.local().elements_arbiter(), vec_rhs.local().elements_arbiter(), apply_info, stat);
        else
          this->_status = VankaGMRES::cuda_apply_amavanka_gmres_parallel(vec_sol.local().elements_arbiter(), vec_rhs.local().elements_arbiter(), apply_info, stat);
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), this->_status, stat.num_iter));
        this->_def_cur = stat.def_norm;
        this->_def_init = stat.def_init;
        this->_num_iter = stat.num_iter;
        return this->_status;
      }
      #endif

      virtual Status _apply_intern(GlobalVecType& vec_sol, const GlobalVecType& vec_rhs)
      {
        const auto id_total = _start_cuda_range("Internal Apply");
        IterationStats pre_iter(*this);
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
        const GlobalMatType& matrix(this->_internal_system_mat);
        const FilType& filter(this->_internal_filter);

        // compute initial defect
        Status status = this->_set_initial_defectp(this->_vec_v.at(0), vec_sol);

        pre_iter.destroy();



        // outer GMRES loop
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          _q.clear();
          _s.clear();
          _c.clear();
          _q.push_back(this->_def_cur);

          // normalize v[0]
          this->_vec_v.at(0).scale(this->_vec_v.at(0), DataType(1) / _q.back());
          _internal_flop_counter += 2u*this->_vec_v.at(0).size_raw();

          const auto id_inner = _start_cuda_range("Inner GMRES Loop");
          // inner GMRES loop
          Index i(0);
          while(i < this->_krylov_dim)
          {
            // if(i == 0)
            // {
            //   std::cout << "Rank " << _comm->rank() << " :" << this->_vec_v[0].local() << "\n-------------\n";
            // }
            // apply vanka matrix
            {
              const auto id_vanka = _start_cuda_range("Solve Vanka");
              _apply_vanka(this->_vec_w, this->_vec_v.at(i));
              _stop_cuda_range(id_vanka);
            }
            // // sync by calculating mean
            // if(i == 0)
            // {
            //   std::cout << "Rank " << _comm->rank() << " :" << this->_vec_w.local() << "\n-------------\n";
            // }

            // v[i+1] := A*w
            {
              const auto id_matvec = _start_cuda_range("A*w");
              matrix.apply(this->_vec_v.at(i+1), this->_vec_w);
              filter.filter_def(this->_vec_v.at(i+1).local());
              _stop_cuda_range(id_matvec);
              _internal_flop_counter += 2u*_global_used_elements;
            }

            // Gram-Schmidt process (classical / conditional CGS2), with the self-dot fused into the
            // batched reduction so the orthogonalized norm is recovered without a separate norm2
            const DataType alpha = _orthogonalize_cpu(i);

            // normalize v[i+1]
            const auto id_scale= _start_cuda_range("Norm scale");
            this->_vec_v.at(i+1).scale(this->_vec_v.at(i+1), DataType(1) / alpha);
            _internal_flop_counter += 2u * _global_vec_size;
            _stop_cuda_range(id_scale);

            {
              const auto id_givens = _start_cuda_range("Givens Rotation");
              // apply Givens rotations
              for(Index k(0); k < i; ++k)
              {
                DataType t(this->_h.at(i).at(k));
                this->_h.at(i).at(k  ) = this->_c.at(k) * t + this->_s.at(k) * this->_h.at(i).at(k+1);
                this->_h.at(i).at(k+1) = this->_s.at(k) * t - this->_c.at(k) * this->_h.at(i).at(k+1);
              }
              _stop_cuda_range(id_givens);
            }


            const auto id_def_push = _start_cuda_range("Push defect");
            // compute beta
            const DataType beta = Math::sqrt(Math::sqr(this->_h.at(i).at(i)) + Math::sqr(alpha));

            // compute next plane rotation
            _s.push_back(alpha / beta);
            _c.push_back(this->_h.at(i).at(i) / beta);

            this->_h.at(i).at(i) = beta;
            this->_q.push_back(this->_s.back() * this->_q.at(i));
            this->_q.at(i) *= this->_c.back();

            _stop_cuda_range(id_def_push);
            // push our new defect
            if(++i < this->_krylov_dim)
            {
              // get the absolute defect
              DataType def_cur = Math::abs(this->_q.back());

              // did we diverge?
              if((def_cur > this->_div_abs) || (def_cur > (this->_div_rel * this->_def_init)))
                break;

              // minimum number of iterations performed?
              if(!(this->_num_iter < this->_min_iter))
              {
                // did we converge?
                if((def_cur <= _inner_res_scale * this->_tol_abs) &&
                  ((def_cur <= _inner_res_scale * (this->_tol_rel * this->_def_init)) || (def_cur <= _inner_res_scale * this->_tol_abs_low) ))
                  break;

                // maximum number of iterations performed?
                if(this->_num_iter >= this->_max_iter)
                  break;
              }

              // set our pseudo defect
              // increase iteration count
              ++this->_num_iter;

              // compute new defect
              this->_def_cur = def_cur;

              // plot?
              if(this->_plot_iter())
              {
                String msg = this->_plot_name
                  +  "* " + stringify(this->_num_iter).pad_front(this->_iter_digits)
                  + " : " + stringify_fp_sci(this->_def_cur);
                this->_print_line(msg);
              }
            }
          }

          _stop_cuda_range(id_inner);


          {
            const Index n = Math::min(i, this->_krylov_dim);
            const auto id_h_solve = _start_cuda_range("Solve Hessenberg");
            // solve H*q = q
            for(Index k(n); k > 0;)
            {
              --k;
              this->_q.at(k) /= this->_h.at(k).at(k);
              for(Index j(k); j > 0;)
              {
                --j;
                this->_q.at(j) -= this->_h.at(k).at(j) * this->_q.at(k);
              }
            }
            _stop_cuda_range(id_h_solve);

            // Runtime::SyncGuard sync_guard;
            // update solution
            this->_vec_w.format();
            // _add_vecs(this->_vec_w, this->_vec_v, this->_q, n, DataType(1));
            _add_vecs_s(this->_vec_w.local(), this->_vec_v_elements, this->_q, n, DataType(1));
            _internal_flop_counter += n * 2u * _global_vec_size;
          }

          {
            const auto id_vanka = _start_cuda_range("Solve Vanka");
            _apply_vanka(this->_vec_v.front(), this->_vec_w);
            vec_sol.axpy(this->_vec_v.front());
            _internal_flop_counter += _global_vec_size;
            _stop_cuda_range(id_vanka);
          }

          // compute "real" residual
          {
            const auto id_residual = _start_cuda_range("Residual Calc");
            matrix.apply(this->_vec_v.at(0), vec_sol, vec_rhs, -DataType(1));
            filter.filter_def(this->_vec_v.at(0).local());
            _stop_cuda_range(id_residual);
            _internal_flop_counter += 2u*_global_used_elements;
          }

          // set the current defect
          status = this->_set_new_defectp(this->_vec_v.at(0), vec_sol);
        }

        _stop_cuda_range(id_total);
        // finished
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
        return status;
      }
    }; // class AmaVankaGMRES<...>

    /**
     * \brief Creates a new GMRES solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] krylov_dim
     * The maximum Krylov subspace dimension. Must be > 0.
     *
     * \param[in] inner_res_scale
     * The scaling factor for the inner GMRES loop residual.
     * See the details of the GMRES class documentation for details
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new GMRES object.
     */
    template<typename Matrix_, typename Filter_, typename VankaType_>
    inline std::shared_ptr<AmaVankaGMRES<Matrix_, Filter_>> new_amavanka_gmres(
      const Matrix_& matrix, const Filter_& filter, std::shared_ptr<VankaType_> vanka, Index krylov_dim,
      typename Matrix_::DataType inner_res_scale = typename Matrix_::DataType(0))
    {
      return std::make_shared<AmaVankaGMRES<Matrix_, Filter_>>(matrix, filter, std::move(vanka), krylov_dim, inner_res_scale);
    }

  } // namespace Solver
} // namespace FEAT
