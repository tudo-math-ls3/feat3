#pragma once
#ifndef KERNEL_SOLVER_BICGSTAB_HPP
#define KERNEL_SOLVER_BICGSTAB_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {

    /**
     * \brief Enum for diffent preconditioning variants for BiCGStab
     */
    enum class BiCGStabPreconVariant
    {
      undefined = 0,
      left,
      right
    };

    /// \cond internal
    /**
     * \brief Streaming operator for BiCGStabPreconVariants
     */
    inline std::ostream& operator<<(std::ostream& os, BiCGStabPreconVariant precon_variant)
    {
      switch(precon_variant)
      {
        case BiCGStabPreconVariant::undefined:
          return os << "undefined";
        case BiCGStabPreconVariant::left:
          return os << "left";
        case BiCGStabPreconVariant::right:
          return os << "right";
        default:
          return os << "-unknown-";
      }
    }

    inline void operator<<(BiCGStabPreconVariant& precon_variant, const String& precon_variant_name)
    {
        if(precon_variant_name == "undefined")
          precon_variant = BiCGStabPreconVariant::undefined;
        else if(precon_variant_name == "left")
          precon_variant = BiCGStabPreconVariant::left;
        else if(precon_variant_name == "right")
          precon_variant = BiCGStabPreconVariant::right;
        else
          throw InternalError(__func__, __FILE__, __LINE__, "Unknown BiCGStabPreconVariant identifier string "
              +precon_variant_name);
    }

    /// \endcond

    /**
     * \brief (Preconditioned) Bi-Conjugate Gradient Stabilized solver implementation
     *
     * This class implements the BiCGStab solver from \cite Vor92 in a modified form to make it more readable,
     * similar in PCG in structure and variables and to save one application of the preconditioner in each iteration.
     * The solver has a left-only and a right-only preconditioned variant, which can be switched at runtime. They
     * differ only in the choice of the scalar product to compute various quantities.
     *
     * \note Which variant gives the better convergence is highly problem dependent.
     *
     * \note Even if no preconditioning is used, the variant has to be set to left or right.
     *
     * \see BiCGStabPreconVariant
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * This is the algorithm:
     * \f{align*}{
     *   r_0 & := b - A x_0 \\
     *   \hat{r}_0 & := r_0 \\
     *   \tilde{r}_0 & := M^{-1} r \\
     *   \tilde{p}_0 & := \tilde{r}_0 \\
     *   \rho_0 & :=
     *   \begin{cases}
     *     \left< \hat{r}_0, \tilde{r}_0 \right>, & \mathrm{~left~preconditioned} \\
     *     \left< \hat{r}_0, r_0 \right>, & \mathrm{~right~preconditioned}
     *   \end{cases}\\
     *   \intertext{For $k = 0, \dots, \mathrm{max\_iter}$}
     *   q_k & = A \tilde{p}_k & \mathrm{(~first~application~of~the~matrix)}\\
     *   \tilde{q}_k & = M^{-1} q_k &\mathrm{(~first~application~of~the~preconditioner)}\\
     *   \alpha_k & =
     *   \begin{cases}
     *     \frac{\rho_k}{\left< \hat{r}_0, q_k \right>}, & \mathrm{~left~preconditioned} \\
     *     \frac{\rho_k}{\left< \hat{r}_0, \tilde{q}_k \right>}, & \mathrm{~right~preconditioned} \\
     *   \end{cases}\\
     *   x_{k+1/2} & = x_k + \alpha_k \tilde{p}\\
     *   r_{k+1/2} & = r_k - \alpha_k q_k \\
     *   \mathrm{If~} x_{k+1/2}, r_{k+1/2}& \mathrm{~are~good~enough:~break}\\
     *   \tilde{r}_{k+1/2} & = \tilde{r}_k - \alpha_k \tilde{q}_k \\
     *   t_k & = A \tilde{r}_{k+1/2} & \mathrm{(~second~application~of~the~matrix)}\\
     *   \tilde{t}_k & = M^{-1} t_k & \mathrm{(~second~application~of~the~preconditioner)}\\
     *   \omega_k& =
     *   \begin{cases}
     *     \frac{\left< \tilde{t}_k, \tilde{r}_{k+1/2} \right>}{\left< \tilde{t}_k, \tilde{t}_k \right>},
     *     & \mathrm{~left~preconditioned} \\
     *     \frac{\left< t_k, r_{k+1/2} \right>}{\left< t_k, t_k \right>}, & \mathrm{~right~preconditioned} \\
     *   \end{cases}\\
     *   x_{k+1} & = x_{k+1/2} + \omega_k \tilde{r}_{k+1/2} \\
     *   r_{k+1} & = r_{k+1/2} - \omega_k t_{k} \\
     *   \mathrm{If~} x_{k+1}, r_{k+1} & \mathrm{~are~good~enough:~break}\\
     *   \tilde{r}_{k+1} & = \tilde{r}_{k+1/2} - \omega_k \tilde{t}_{k} \\
     *   \rho_{k+1} & :=
     *   \begin{cases}
     *     \left< \hat{r}_0, \tilde{r}_k \right>, & \mathrm{~left~preconditioned} \\
     *     \left< \hat{r}_0, r_k \right>, & \mathrm{~right~preconditioned}
     *   \end{cases}\\
     *   \beta_k & = \frac{\rho_k \omega_k}{\rho_{k+1} \alpha_k} \\
     *   \tilde{p}_{k+1} & = \tilde{r}_{k+1} + \beta_k (\tilde{p}_k - \omega \tilde{q}_k)
     * \f}
     *
     * \author Jordi Paul and Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class BiCGStab :
      public PreconditionedIterativeSolver<typename Matrix_::VectorTypeR>
    {
      public:
        /// The matrix type this solver can be applied to
        typedef Matrix_ MatrixType;
        /// The filter to apply
        typedef Filter_ FilterType;
        /// The vector type this solver can be applied to
        typedef typename MatrixType::VectorTypeR VectorType;
        /// The floating point precision
        typedef typename MatrixType::DataType DataType;
        /// Our base class
        typedef PreconditionedIterativeSolver<VectorType> BaseClass;
        /// The type of the preconditioner
        typedef SolverBase<VectorType> PrecondType;

      protected:
        /// the matrix for the solver
        const MatrixType& _system_matrix;
        /// the filter for the solver
        const FilterType& _system_filter;
        /// The preconditioned primal search direction
        VectorType _vec_p_tilde;
        /// q~[k] = M^{-1} A q[k]
        VectorType _vec_q_tilde;
        /// The defect vector
        VectorType _vec_r;
        /// The preconditioned defect vector
        VectorType _vec_r_tilde;
        /// The arbitrary starting vector
        VectorType _vec_r_hat_0;
        /// t~[k] = M^{-1} A r~[k+1/2]
        VectorType _vec_t_tilde;
        /// Temporary vector
        VectorType _vec_tmp;
        /// Precondition from left or right?
        BiCGStabPreconVariant _precon_variant;

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
         * \param[in] precond
         * A pointer to the preconditioner. May be \c nullptr.
         *
         * \param[in] precon_variant
         * Which preconditioning variant to use, defaults to left
         */
        explicit BiCGStab(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr,
        BiCGStabPreconVariant precon_variant = BiCGStabPreconVariant::left) :
          BaseClass("BiCGStab", precond),
          _system_matrix(matrix),
          _system_filter(filter),
          _precon_variant(precon_variant)
          {
          }

        /**
         * \brief Constructor using a PropertyMap
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
         */
        explicit BiCGStab(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter, std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass("BiCGStab", section_name, section, precond),
          _system_matrix(matrix),
          _system_filter(filter),
          _precon_variant(BiCGStabPreconVariant::left)
          {
            // Check if we have set _p_variant
            auto p_variant_p = section->query("precon_variant");
            if(p_variant_p.second)
            {
              BiCGStabPreconVariant precon_variant;
              precon_variant << p_variant_p.first;
              set_precon_variant(precon_variant);
            }
          }

        /// \copydoc SolverBase::name()
        virtual String name() const override
        {
          return "BiCGStab";
        }

        /// \copydoc SolverBase::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();
          // create all temporary vectors
          _vec_p_tilde = this->_system_matrix.create_vector_r();
          _vec_q_tilde = this->_system_matrix.create_vector_r();
          _vec_r = this->_system_matrix.create_vector_r();
          _vec_r_hat_0 = this->_system_matrix.create_vector_r();
          _vec_r_tilde = this->_system_matrix.create_vector_r();
          _vec_t_tilde = this->_system_matrix.create_vector_r();
          _vec_tmp = this->_system_matrix.create_vector_r();
        }

        /// \copydoc SolverBase::done_symbolic()
        virtual void done_symbolic() override
        {
          _vec_p_tilde.clear();
          _vec_q_tilde.clear();
          _vec_r.clear();
          _vec_r_hat_0.clear();
          _vec_r_tilde.clear();
          _vec_t_tilde.clear();
          _vec_tmp.clear();

          BaseClass::done_symbolic();
        }

        /// \copydoc IterativeSolver::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
        {
          // compute initial defect
          this->_system_matrix.apply(this->_vec_r, vec_sol, vec_rhs, -DataType(1));
          this->_system_filter.filter_def(this->_vec_r);

          // apply
          return _apply_intern(vec_sol, vec_rhs);
        }

        /// \copydoc SolverBase::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // save rhs vector as initial defect
          this->_vec_r.copy(vec_def);
          //this->_system_filter.filter_def(this->_vec_r);

          // format solution vector
          vec_cor.format();

          return _apply_intern(vec_cor, vec_def);
        }


      /**
       * \brief Sets the preconditioning variant (left or right)
       *
       * \param[in] precon_variant
       * Which preconditioning variant to use.
       */
      virtual void set_precon_variant(BiCGStabPreconVariant precon_variant)
      {
        XASSERT(precon_variant == BiCGStabPreconVariant::left || precon_variant == BiCGStabPreconVariant::right);
        _precon_variant = precon_variant;
      }

      /// \copydoc SolverBase::write_config()
      virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name) const override
      {
        XASSERT(parent != nullptr);

        PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

        my_section->add_entry("precon_variant", stringify(_precon_variant));

        return my_section;
      }

      protected:
      /**
       * \brief Internal function, applies the solver
       *
       * \param[in] vec_sol
       * The current solution vector, gets overwritten
       *
       * \param[in] vec_rhs
       * The right hand side vector. This is unused in this function, as the intial defect was alredy computed and
       * stored in _vec_r.
       *
       * \returns A status code.
       */
        Status _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
          VectorType& vec_p_tilde  (_vec_p_tilde);
          VectorType& vec_r        (_vec_r);
          VectorType& vec_r_tilde  (_vec_r_tilde);
          VectorType& vec_t        (_vec_tmp);
          VectorType& vec_t_tilde  (_vec_t_tilde);
          VectorType& vec_q        (_vec_tmp);
          VectorType& vec_q_tilde  (_vec_q_tilde);

          const MatrixType& mat_sys(this->_system_matrix);
          const FilterType& fil_sys(this->_system_filter);
          Status status(Status::progress);

          // Apply preconditioner to initial defect and save it to p_0
          if(!this->_apply_precond(vec_p_tilde, _vec_r, fil_sys))
          {
            Statistics::add_solver_expression(
              std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          // We set the arbitrary vector r^[0] = r[0]
          _vec_r_hat_0.copy(vec_r);
          // r~[0] = M^{-1} r[0]
          vec_r_tilde.copy(vec_p_tilde);

          const VectorType& vec_r_hat_0(_vec_r_hat_0);

          DataType rho(0);

          // Left preconditioned: rho[k] = <r^[0], r~[k]>
          if(_precon_variant == BiCGStabPreconVariant::left)
          {
            rho = vec_r_hat_0.dot(vec_r_tilde);
          }
          // Right preconditioned: rho[k] = <r^[0], r[k]>
          else
          {
            rho = vec_r_hat_0.dot(vec_r);
          }

          // Compute initial defect
          status = this->_set_initial_defect(vec_r, vec_sol);

          while(status == Status::progress)
          {
            IterationStats stat(*this);

            // q[k] = A p~[k]
            mat_sys.apply(vec_q, vec_p_tilde);
            fil_sys.filter_def(vec_q);

            // Apply preconditioner to primal search direction
            // q~[k] = M^{-1} q[k]
            if(!this->_apply_precond(vec_q_tilde, vec_q, fil_sys))
            {
              Statistics::add_solver_expression(
                std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            DataType alpha(0);
            // Left preconditioned: alpha[k] = rho[k] / <r^[0], q~[k]>
            if(_precon_variant == BiCGStabPreconVariant::left)
            {
              alpha = rho / vec_r_hat_0.dot(vec_q_tilde);
            }
            // Right preconditioned: alpha[k] = rho[k] / <r^[0], q[k]>
            else
            {
              alpha = rho / vec_r_hat_0.dot(vec_q);
            }

            // First "half" update
            // x[k+1/2] = x[k] + alpha[k] p~[k]
            vec_sol.axpy(vec_p_tilde, vec_sol, alpha);

            // r[k+1/2] = r[k] - alpha[k] q[k]
            vec_r.axpy(vec_q, vec_r, -alpha);

            // Check if we are already converged or failed after the "half" update
            {
              Status status_half(Status::progress);

              DataType def_old(this->_def_cur);
              DataType def_half(vec_r.norm2());

              // ensure that the defect is neither NaN nor infinity
              if(!Math::isfinite(def_half))
              {
                status_half = Status::aborted;
              }

              // is diverged?
              if(this->is_diverged(def_half))
              {
                this->_def_cur = def_half;
                status_half = Status::diverged;
              }

              // is converged?
              if(this->is_converged(def_half))
              {
                this->_def_cur = def_half;
                status_half = Status::success;
              }

              // If we break early, we still count this iteration and plot it if necessary
              if(status_half != Status::progress)
              {
                ++this->_num_iter;
                // plot?
                if(this->_plot)
                {
                  std::cout << this->_plot_name
                  <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
                  << " : " << stringify_fp_sci(this->_def_cur)
                  << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
                  << " / " << stringify_fp_fix(this->_def_cur / def_old)
                  << std::endl;
                }
                Statistics::add_solver_expression(
                  std::make_shared<ExpressionEndSolve>(this->name(), status_half, this->get_num_iter()));

                return status_half;
              }
            }

            // Update preconditioned defect: r~[k+1/2] = r~[k] - alpha q[k]
            vec_r_tilde.axpy(vec_q_tilde, vec_r_tilde, -alpha);

            // t = A r~[k+1/2]
            mat_sys.apply(vec_t, vec_r_tilde);
            fil_sys.filter_def(vec_t);

            // Apply preconditioner for correction direction
            // t~[k] = M^{-1} t
            if(!this->_apply_precond(vec_t_tilde, vec_t, fil_sys))
            {
              Statistics::add_solver_expression(
                std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            DataType omega(0);
            // Left preconditioned: omega[k] = <t~[k], r~[k+1/2] / <t~[k], t~[k]>
            if(_precon_variant == BiCGStabPreconVariant::left)
            {
              omega = vec_t_tilde.dot(vec_r_tilde) / vec_t_tilde.dot(vec_t_tilde);
            }
            // Right preconditioned: omega[k] = <t[k], r[k+1/2] / <t[k], t[k]>
            else
            {
              omega = vec_t.dot(vec_r) / vec_t.dot(vec_t);
            }

            XASSERTM(Math::isfinite(omega),"BiCGStab breakdown!");

            // Second "half" update
            // x[k+1] = x[k+1/2] + omega r~[k+1/2]
            vec_sol.axpy(vec_r_tilde, vec_sol, omega);

            // Upate defect
            // r[k+1] = r[k] - omega t[k]
            vec_r.axpy(vec_t, vec_r, -omega);

            // Compute defect norm
            status = this->_set_new_defect(vec_r, vec_sol);

            if(status != Status::progress)
            {
              Statistics::add_solver_expression(
                std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
              return status;
            }

            // Update preconditioned defect
            // r~[k+1] = r~[k] - omega t~[k]
            vec_r_tilde.axpy(vec_t_tilde, vec_r_tilde, -omega);

            // Save old rho
            DataType rho2(rho);
            // Left preconditioned: rho[k+1] = <r^[0], r~[k+1]>
            if(_precon_variant == BiCGStabPreconVariant::left)
            {
              rho = vec_r_hat_0.dot(vec_r_tilde);
            }
            // Right preconditioned: rho[k+1] = <r^[0], r[k+1]>
            else
            {
              rho = vec_r_hat_0.dot(vec_r);
            }

            // beta[k] = (rho[k+1]*alpha[k]) / (rho[k]*omega[k])
            DataType beta((rho/rho2)*(alpha/omega));
            XASSERTM(Math::isfinite(beta),"BiCGStab  breakdown!");

            // p~[k+1] = r~[k+1] + beta(p~[k] - omega[k] q~[k])
            vec_p_tilde.axpy(vec_q_tilde, vec_p_tilde, -omega);
            vec_p_tilde.axpy(vec_p_tilde, vec_r_tilde, beta);

          }

          // we should never reach this point...
          Statistics::add_solver_expression(
            std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return Status::undefined;
        }
    }; // class BiCGStab<...>

    /**
     * \brief Creates a new BiCGStab solver object
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
     * \param[in] precon_variant
     * Which preconditioning variant to use, defaults to left
     *
     * \returns
     * A shared pointer to a new BiCGStab object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const Matrix_& matrix, const Filter_& filter)
      {
        return std::make_shared<BiCGStab<Matrix_, Filter_>>(matrix, filter, nullptr);
      }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond, BiCGStabPreconVariant precon_variant = BiCGStabPreconVariant::left)
      {
        return std::make_shared<BiCGStab<Matrix_, Filter_>>(matrix, filter, precond, precon_variant);
      }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr,
      BiCGStabPreconVariant precon_variant = BiCGStabPreconVariant::left)
      {
        return std::make_shared<BiCGStab<Matrix_, Filter_>>(matrix, filter, precond, precon_variant);
      }
#endif

    /**
     * \brief Creates a new BiCGStab solver object using a PropertyMap
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
     * A shared pointer to a new BiCGStab object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
      {
        return std::make_shared<BiCGStab<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
      }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
      {
        return std::make_shared<BiCGStab<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
      }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
      {
        return std::make_shared<BiCGStab<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
      }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_BICGSTAB_HPP
