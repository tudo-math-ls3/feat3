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
     * This class implements the BiCGStab solver from \cite Vor92.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * Note that the preconditioner can be applied left-only or right-only - the algorithm differs only in the
     * computation of the parameter \f$ \omega \f$. In the right-preconditioned variant, it is computed to minimise
     * the residual of the unpreconditioned system. This saves one call to the preconditioner and the storage of one
     * internal vector. However, even if both matrix and preconditioner are near-symmetric, this variant destroys
     * the near-symmetry of the preconditioned operator.
     *
     * So if you suspect your system to be near symmetric and you have a symmetric preconditioner (such as Jacobi),
     * it might be beneficial to set_precon_variant(left). The right-variant is cheaper, though.
     *
     * \author Jordi Paul
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
        /// The unpreconditioned primal search direction
        VectorType _vec_p;
        /// The defect vector
        VectorType _vec_r;
        /// The "dual" intial defect vector
        VectorType _vec_r_tilde_0;
        /// The defect vector after the "half" update
        VectorType _vec_s;
        /// t = A z
        VectorType _vec_t;
        /// v = A y
        VectorType _vec_v;
        /// The precontitioned _vec_t
        VectorType _vec_w;
        /// The preconditioned primal search direction
        VectorType _vec_y;
        /// The preconditioned defect vector after the "half" update
        VectorType _vec_z;
        /// Use the p_variant?
        BiCGStabPreconVariant _precon_variant;
        /// We need to know if init_symbolic was called because the _precon_variant cannot be switched then
        bool _have_init_symbolic;

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
         * Which preconditioning variant to use, defaults to right
         */
        explicit BiCGStab(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr, BiCGStabPreconVariant precon_variant = BiCGStabPreconVariant::right) :
          BaseClass("BiCGStab", precond),
          _system_matrix(matrix),
          _system_filter(filter),
          _precon_variant(precon_variant),
          _have_init_symbolic(false)
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
          _have_init_symbolic(false)
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
          _vec_p = this->_system_matrix.create_vector_r();
          _vec_r = this->_system_matrix.create_vector_r();
          _vec_r_tilde_0 = this->_system_matrix.create_vector_r();
          _vec_s = this->_system_matrix.create_vector_r();
          _vec_t = this->_system_matrix.create_vector_r();
          _vec_v = this->_system_matrix.create_vector_r();
          _vec_y = this->_system_matrix.create_vector_r();
          _vec_z = this->_system_matrix.create_vector_r();

          // If the p-variant is used, we skip one application of the preconditioner so this can be a shallow copy
          if(_precon_variant == BiCGStabPreconVariant::right)
          {
            _vec_w.clone(_vec_t, LAFEM::CloneMode::Shallow);
          }
          else
          {
            _vec_w = this->_system_matrix.create_vector_r();
          }

          _have_init_symbolic = true;
        }

        /// \copydoc SolverBase::done_symbolic()
        virtual void done_symbolic() override
        {
          _vec_p.clear();
          _vec_r.clear();
          _vec_r_tilde_0.clear();
          _vec_s.clear();
          _vec_t.clear();
          _vec_v.clear();
          _vec_y.clear();
          _vec_z.clear();
          // If the right-variant is used, we skip one application of the preconditioner so this is a shallow copy
          // and must not be cleared
          if(_precon_variant != BiCGStabPreconVariant::right)
          {
            _vec_w.clear();
          }
          BaseClass::done_symbolic();
          _have_init_symbolic = false;
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
        XASSERTM(!_have_init_symbolic,"precon_variant can only be switched before init_symbolic() is called!");
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
        virtual Status _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
          VectorType& vec_p        (_vec_p);
          VectorType& vec_r        (_vec_r);
          _vec_r_tilde_0.copy(_vec_r);
          const VectorType& vec_r_tilde_0(_vec_r_tilde_0);
          VectorType& vec_s        (_vec_s);
          VectorType& vec_t        (_vec_t);
          VectorType& vec_v        (_vec_v);
          VectorType& vec_w        (_vec_w);
          VectorType& vec_y        (_vec_y);
          VectorType& vec_z        (_vec_z);
          const MatrixType& mat_sys(this->_system_matrix);
          const FilterType& fil_sys(this->_system_filter);
          Status status(Status::progress);

          DataType alpha(1);
          DataType rho(1);
          DataType omega(1);

          DataType rho2(rho);

          vec_v.format();
          vec_p.format();

          // Compute initial defect
          status = this->_set_initial_defect(vec_r, vec_sol);

          while(status == Status::progress)
          {
            IterationStats stat(*this);

            rho2 = rho;
            rho = vec_r_tilde_0.dot(vec_r);

            DataType beta((rho/rho2)*(alpha/omega));
            XASSERTM(Math::isfinite(beta),"BiCGStab  breakdown!");

            // p[k] = r[k] + beta(p[k-1] - omega[k-1]v[k-1])
            vec_p.axpy(vec_v, vec_p, -omega);
            vec_p.axpy(vec_p, vec_r, beta);

            // Apply preconditioner to primal search direction
            if(!this->_apply_precond(vec_y, vec_p, fil_sys))
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            // v[k] = Ay[k]
            mat_sys.apply(vec_v, vec_y);
            fil_sys.filter_def(vec_v);

            alpha = rho / vec_r_tilde_0.dot(vec_v);

            // s = r[k] - alpha v[k]
            vec_s.axpy(vec_v, vec_r, -alpha);

            // First "half" update
            // x[k+1/2] = x[k] + alpha y
            vec_sol.axpy(vec_y, vec_sol, alpha);

            // Check if we are already converged or failed after the "half" update
            {
              Status status_half(Status::progress);

              DataType def_old(this->_def_cur);
              DataType def_half(vec_s.norm2());

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

            // Apply preconditioner for correction direction
            if(!this->_apply_precond(vec_z, vec_s, fil_sys))
            {
              Statistics::add_solver_expression(
                std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            // t = Az
            mat_sys.apply(vec_t, vec_z);
            fil_sys.filter_def(vec_t);

            // If the right-variant is used, we set omega[k] = < t, s > / < t, t >
            // Otherwise (meaning the left-variant), compute w = M^{-1} t and set omega[k] = < w, s> / < w, w >
            if(_precon_variant == BiCGStabPreconVariant::left)
            {
              // Apply preconditioner for the computation of omega
              if(!this->_apply_precond(vec_w, vec_t, fil_sys))
              {
                Statistics::add_solver_expression(
                  std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
                return Status::aborted;
              }
            }

            DataType normsqr_vec_w(vec_w.norm2sqr());
            omega = vec_w.dot(vec_z) / normsqr_vec_w;

            XASSERTM(Math::isfinite(omega),"BiCGStab breakdown!");

            // Second "half" update
            vec_sol.axpy(vec_z, vec_sol, omega);

            // Upate defect
            vec_r.axpy(vec_t, vec_s, -omega);

            // compute defect norm
            status = this->_set_new_defect(vec_r, vec_sol);

            if(status != Status::progress)
            {
              Statistics::add_solver_expression(
                std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
              return status;
            }

          }

          // we should never reach this point...
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
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
     * Which preconditioning variant to use, defaults to right
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
      std::shared_ptr<Precond_> precond, BiCGStabPreconVariant precon_variant = BiCGStabPreconVariant::right)
      {
        return std::make_shared<BiCGStab<Matrix_, Filter_>>(matrix, filter, precond, precon_variant);
      }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr,
      BiCGStabPreconVariant precon_variant = BiCGStabPreconVariant::right)
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
