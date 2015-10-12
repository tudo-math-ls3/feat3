#pragma once
#ifndef KERNEL_SOLVER_FGMRES_HPP
#define KERNEL_SOLVER_FGMRES_HPP 1

// includes, FEAST
#include <kernel/solver/iterative.hpp>

// includes, FEAST
#include <vector>

namespace FEAST
{
  namespace Solver
  {
    /**
     * \brief FGMRES(k) solver implementation
     *
     * This class implements the Restarted Flexible Generalised Minimal Residual solver (FGMRES(k)).
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
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_>
    class FGMRES :
      public PreconditionedIterativeSolver<typename Matrix_::VectorTypeR>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<VectorType> BaseClass;

      typedef SolverBase<VectorType> PrecondType;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// krylov dimension
      Index _krylov_dim;
      /// inner pseudo-residual scaling factor
      DataType _inner_res_scale;
      /// krylov basis vectors
      std::vector<VectorType> _vec_v, _vec_z;
      /// Givens rotation coefficients
      std::vector<DataType> _c, _s, _q;
      /// Hessenberg matrix
      std::vector<std::vector<DataType>> _h;

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
       * The scaling factor for the inner GMRES loop residual.
       * Set this to zero unless you know what you are doing.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       */
      explicit FGMRES(const MatrixType& matrix, const FilterType& filter, Index krylov_dim,
        DataType inner_res_scale = DataType(0), std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("FGMRES(" + stringify(krylov_dim) + ")", precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _krylov_dim(krylov_dim),
        _inner_res_scale(inner_res_scale)
      {
        _c.reserve(krylov_dim);
        _s.reserve(krylov_dim);
        _q.reserve(krylov_dim);
        _h.resize(krylov_dim);
        for(Index i(0); i < krylov_dim; ++i)
          _h.at(i).resize(i+1);
      }

      virtual String name() const override
      {
        return "FGMRES";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        _vec_v.push_back(this->_system_matrix.create_vector_r());
        for(Index i(0); i < _krylov_dim; ++i)
        {
          _vec_v.push_back(this->_vec_v.front().clone(LAFEM::CloneMode::Layout));
          _vec_z.push_back(this->_vec_v.front().clone(LAFEM::CloneMode::Layout));
        }
      }

      virtual void done_symbolic() override
      {
        _vec_v.clear();
        _vec_z.clear();
        BaseClass::done_symbolic();
      }

      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // save input rhs vector as initial defect
        this->_vec_v.at(0).copy(vec_rhs);
        //this->_system_filter.filter_def(this->_vec_v.at(0));

        // clear solution vector
        vec_sol.format();

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute initial defect
        this->_system_matrix.apply(this->_vec_v.at(0), vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_v.at(0));

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);

        // compute initial defect
        Status status = this->_set_initial_defect(this->_vec_v.at(0), vec_sol);

        // outer GMRES loop
        while(status == Status::progress)
        {
          _q.clear();
          _s.clear();
          _c.clear();
          _q.push_back(this->_def_cur);

          // normalise v[0]
          this->_vec_v.at(0).scale(this->_vec_v.at(0), DataType(1) / _q.back());

          // inner GMRES loop
          Index i(0);
          while(i < this->_krylov_dim)
          {
            // apply preconditioner
            if(!this->_apply_precond(this->_vec_z.at(i), this->_vec_v.at(i), filter))
              return Status::aborted;
            //filter.filter_cor(this->_vec_z.at(i));

            // v[i+1] := A*z[i]
            matrix.apply(this->_vec_v.at(i+1), this->_vec_z.at(i));
            filter.filter_def(this->_vec_v.at(i+1));

            // Gram-Schmidt process
            for(Index k(0); k <= i; ++k)
            {
              this->_h.at(i).at(k) = this->_vec_v.at(i+1).dot(this->_vec_v.at(k));
              this->_vec_v.at(i+1).axpy(this->_vec_v.at(k), this->_vec_v.at(i+1), -this->_h.at(i).at(k));
            }

            // normalise v[i+1]
            DataType alpha = this->_vec_v.at(i+1).norm2();
            this->_vec_v.at(i+1).scale(this->_vec_v.at(i+1), DataType(1) / alpha);

            // apply Givens rotations
            for(Index k(0); k < i; ++k)
            {
              DataType t(this->_h.at(i).at(k));
              this->_h.at(i).at(k  ) = this->_c.at(k) * t + this->_s.at(k) * this->_h.at(i).at(k+1);
              this->_h.at(i).at(k+1) = this->_s.at(k) * t - this->_c.at(k) * this->_h.at(i).at(k+1);
            }

            // compute beta
            DataType beta = Math::sqrt(Math::sqr(this->_h.at(i).at(i)) + Math::sqr(alpha));

            // compute next plane rotation
            _s.push_back(alpha / beta);
            _c.push_back(this->_h.at(i).at(i) / beta);

            this->_h.at(i).at(i) = beta;
            this->_q.push_back(this->_s.back() * this->_q.at(i));
            this->_q.at(i) *= this->_c.back();

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
                if((def_cur <= _inner_res_scale * this->_tol_abs) && (def_cur <= _inner_res_scale * (this->_tol_rel * this->_def_init)))
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
              if(this->_plot)
              {
                std::cout << this->_plot_name
                  <<  "* " << stringify(this->_num_iter).pad_front(this->_iter_digits)
                  << " : " << scientify(this->_def_cur) << std::endl;
              }
            }
          }

          Index n = Math::min(i, this->_krylov_dim);

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

          // update solution
          for(Index k(0); k < n; ++k)
            vec_sol.axpy(this->_vec_z.at(k), vec_sol, this->_q.at(k));

          // compute "real" residual
          matrix.apply(this->_vec_v.at(0), vec_sol, vec_rhs, -DataType(1));
          filter.filter_def(this->_vec_v.at(0));

          // set the current defect
          status = this->_set_new_defect(this->_vec_v.at(0), vec_sol);
        }

        // finished
        return status;
      }
    }; // class FGMRES<...>

    /**
     * \brief Creates a new FGMRES solver object
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
     * Set this to zero unless you know what you are doing.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new FGMRES object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAST_COMPILER_GNU) && (FEAST_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<FGMRES<Matrix_, Filter_>> new_fgmres(
      const Matrix_& matrix, const Filter_& filter, Index krylov_dim,
      typename Matrix_::DataType inner_res_scale = typename Matrix_::DataType(0))
    {
      return std::make_shared<FGMRES<Matrix_, Filter_>>(matrix, filter, krylov_dim, inner_res_scale);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<FGMRES<Matrix_, Filter_>> new_fgmres(
      const Matrix_& matrix, const Filter_& filter, Index krylov_dim,
      typename Matrix_::DataType inner_res_scale,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<FGMRES<Matrix_, Filter_>>(matrix, filter, krylov_dim, inner_res_scale, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<FGMRES<Matrix_, Filter_>> new_fgmres(
      const Matrix_& matrix, const Filter_& filter, Index krylov_dim,
      typename Matrix_::DataType inner_res_scale = typename Matrix_::DataType(0),
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<FGMRES<Matrix_, Filter_>>(matrix, filter, krylov_dim, inner_res_scale, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_FGMRES_HPP
