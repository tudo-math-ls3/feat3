// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_FGMRES_HPP
#define KERNEL_SOLVER_FGMRES_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

// includes, FEAT
#include <vector>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief FGMRES(k) solver implementation
     *
     * This class implements the Restarted Flexible Generalized Minimal Residual solver (FGMRES(k)), which is
     * a variation of the "normal" preconditioned Restarted GMRES (GMRES(k)) solver that supports varying
     * preconditioners, i.e. it allows the use of preconditioner \f$M^{-1}\f$ which does not necessarily have to
     * fulfill the linearity condition \f$M^{-1}(\alpha x + y) = \alpha M^{-1}x + M^{-1}y\f$.
     * If used with a preconditioner that does fulfill the linearity condition or if used without a preconditioner,
     * the FGMRES(k) and GMRES(k) solvers are equivalent with the exception of rounding errors, of course.
     *
     * <u><b>Details about the inner resolution scaling factor:</b></u>\n
     * In contrast to other popular iterative Krylov subspace methods like PCG or BiCGStab, the GMRES method and all
     * of its variants -- including the FMGRES(k) solver implemented by this class -- do \b not compute or update the
     * residual/defect vector after each iteration! This is not an oversight, but it is a fundamental property of the
     * GMRES methods that both the solution vector as well as the residual/defect vectors are only computed after the
     * GMRES method has terminated or, in the case of restarted [F]GMRES(k) variants, after each restart, which happens
     * after at most k iterations. Saad describes this problem and proposes a solution in his book \cite Saad03 in
     * Proposition 6.9, which also states that the euclidean norm of the residual \f$\|b-Ax_m\|_2\f$ is equal to
     * a scalar quantity \f$|\gamma_{m+1}|\f$, which is easily computed from the other quantities that have to be
     * updated in each GMRES iteration. Therefore, one might simply use this scalar to implement a absolute or relative
     * defect tolerance stopping criterion without the need to explicitly compute the residual/defect vector itself.
     *
     * However, in the real world case of finite precision arithmetic, the identity is merely an approximation, i.e.
     * we actually have \f$|\gamma_{m+1}|\approx\|b-Ax_m\|_2\f$ and therefore checking \f$|\gamma_{m+1}|\f$ against
     * the stopping criterion may lead to false positives and therefore premature termination of the inner FGMRES(k)
     * loop. At the end of each inner FGMRES(k) loop, the solution and residual/defect vectors are updated and now
     * the real residual/defect norm is computed and checked against the tolerances and it may now happen that
     * \f$|\gamma_{m+1}| < \textnormal{tol}_{\textnormal{abs}} < \|b-Ax_m\|_2\f\f$ (and analogously for the relative
     * tolerance, of course), therefore leading to another restart of the inner FGMRES(k) loop. In the worst case,
     * this restarted inner loop may again produce a false positive after 1 iteration, therefore triggering another
     * restart, thus effectively leading to an extremely slowly converging FGMRES(1) iteration.
     *
     * To deal with this problem, this implementation includes an inner residual scaling factor \f$\delta\geq 0\f$,
     * which is used to check \f$\delta|\gamma_{m+1}| < \textnormal{tol}_{\textnormal{abs}}\f$ for the absolute and
     * \f$\delta|\gamma_{m+1}| < \textnormal{tol}_{\textnormal{rel}}\|r_0\|_2\f$ for the relative stopping criterion,
     * respectively. By choosing a \f$0 < \delta < 1\f$, one can ensure that the inner GMRES loop has to fulfill a
     * tighter tolerance than the outer loop, therefore reducing the risk of false positives. It is also possible to
     * set \f$\delta = 0\f$, in which case the inner GMRES loop will only terminate prematurely if it detects that it
     * has arrived at the exact solution and no more inner iterations must be performed to avoid division by zero.
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
      /// inner pseudo-residual scaling factor, see the documentation of this class for details
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
       * The scaling factor for the inner GMRES loop residual. See the documentation of this class for details.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       */
      explicit FGMRES(const MatrixType& matrix, const FilterType& filter, Index krylov_dim,
        DataType inner_res_scale = DataType(0), std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("FGMRES(" + stringify(krylov_dim) + ")", precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _krylov_dim(krylov_dim)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
        set_inner_res_scale(inner_res_scale);
      }

      explicit FGMRES(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter, std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("FGMRES", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _inner_res_scale(DataType(0))
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);

        // Set _inner_res_scale by parameter or use default value
        auto inner_res_scale_p = section->query("inner_res_scale");
        if(inner_res_scale_p.second && (!inner_res_scale_p.first.parse(this->_inner_res_scale) || (this->_inner_res_scale < DataType(0))))
          throw ParseError(section_name + ".inner_res_scale", inner_res_scale_p.first, "a non-negative float");

        // Check if we have set _krylov_vim
        auto krylov_dim_p = section->query("krylov_dim");
        if(!krylov_dim_p.second)
          throw ParseError("FGMRES config section is missing the mandatory krylov_dim!");

        if(!krylov_dim_p.first.parse(this->_krylov_dim) || (this->_krylov_dim <= Index(0)))
          throw ParseError(section_name + ".krylov_dim", krylov_dim_p.first, "a positive integer");

        this->set_plot_name("FGMRES("+stringify(_krylov_dim)+")");
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~FGMRES()
      {
      }

      /// \copydoc BaseClass::name()
      virtual String name() const override
      {
        return "FGMRES";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        _c.reserve(_krylov_dim);
        _s.reserve(_krylov_dim);
        _q.reserve(_krylov_dim + 1);
        _h.resize(_krylov_dim);

        for(Index i(0); i < _krylov_dim; ++i)
        {
          _h.at(i).resize(i+1);
        }

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

      /**
       * \brief Sets the inner Krylov space dimension
       *
       * \param[in] krylov_dim
       * The k in FGMRES(k)
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

      /// \copydoc IterativeSolver::apply()
      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // save input rhs vector as initial defect
        this->_vec_v.at(0).copy(vec_rhs);

        // clear solution vector
        vec_sol.format();

        // apply
        this->_status = _apply_intern(vec_sol, vec_rhs);
        this->plot_summary();
        return this->_status;
      }

      /// \copydoc SolverBase::correct()
      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute initial defect
        this->_system_matrix.apply(this->_vec_v.at(0), vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_v.at(0));

        // apply
        this->_status = _apply_intern(vec_sol, vec_rhs);
        this->plot_summary();
        return this->_status;
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        IterationStats pre_iter(*this);
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);

        // compute initial defect
        Status status = this->_set_initial_defect(this->_vec_v.at(0), vec_sol);

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

          // inner GMRES loop
          Index i(0);
          while(i < this->_krylov_dim)
          {
            // apply preconditioner
            if(!this->_apply_precond(this->_vec_z.at(i), this->_vec_v.at(i), filter))
            {
              stat.destroy();
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }
            //filter.filter_cor(this->_vec_z.at(i));

            // v[i+1] := A*z[i]
            matrix.apply(this->_vec_v.at(i+1), this->_vec_z.at(i));
            filter.filter_def(this->_vec_v.at(i+1));

            // Gram-Schmidt process
            for(Index k(0); k <= i; ++k)
            {
              this->_h.at(i).at(k) = this->_vec_v.at(i+1).dot(this->_vec_v.at(k));
              this->_vec_v.at(i+1).axpy(this->_vec_v.at(k), -this->_h.at(i).at(k));
            }

            // normalize v[i+1]
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
            vec_sol.axpy(this->_vec_z.at(k), this->_q.at(k));

          // compute "real" residual
          matrix.apply(this->_vec_v.at(0), vec_sol, vec_rhs, -DataType(1));
          filter.filter_def(this->_vec_v.at(0));

          // set the current defect
          status = this->_set_new_defect(this->_vec_v.at(0), vec_sol);
        }

        // finished
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
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
     * See the details of the FGMRES class documentation for details
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new FGMRES object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<FGMRES<Matrix_, Filter_>> new_fgmres(
      const Matrix_& matrix, const Filter_& filter, Index krylov_dim,
      typename Matrix_::DataType inner_res_scale = typename Matrix_::DataType(0),
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<FGMRES<Matrix_, Filter_>>(matrix, filter, krylov_dim, inner_res_scale, precond);
    }

    /**
     * \brief Creates a new FGMRES solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new FGMRES object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<FGMRES<Matrix_, Filter_>> new_fgmres(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<FGMRES<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_FGMRES_HPP
