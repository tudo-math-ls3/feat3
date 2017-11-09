#pragma once
#ifndef KERNEL_SOLVER_BICGSTABL_HPP
#define KERNEL_SOLVER_BICGSTABL_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {

    /**
     * \brief Enum for diffent preconditioning variants for BiCGStabL
     */

    enum class BiCGStabLPreconVariant
    {
      left,
      right
    };

    /// \cond internal
    /**
     * \brief Streaming operator for BiCGStabLPreconVariants
     */
    inline std::ostream& operator<<(std::ostream& os, BiCGStabLPreconVariant precon_variant)
    {
      switch(precon_variant)
      {
        case BiCGStabLPreconVariant::left:
          return os << "left";
        case BiCGStabLPreconVariant::right:
          return os << "right";
        default:
          return os << "-unknown-";
      }
    }

    inline void operator<<(BiCGStabLPreconVariant& precon_variant, const String& precon_variant_name)
    {
      if(precon_variant_name == "left")
        precon_variant = BiCGStabLPreconVariant::left;
      else if(precon_variant_name == "right")
        precon_variant = BiCGStabLPreconVariant::right;
      else
        throw InternalError(__func__, __FILE__, __LINE__, "Unknown BiCGStabLPreconVariant identifier string "
            +precon_variant_name);
    }

    /// \endcond

    /**
     * \brief (Preconditioned) BiCGStab(l) solver implementation
     *
     * This class implements the BiCGStab(l) solver from \cite Sleijpen1993 with left or right preconditioning.
     * The algorithm is identical to the algorithm given by Sleijpen and Fokkema, but implements preconditioning
     * in an intuitive way :
     *
     * Right precondtioning : solve AMy=b for y and set x=My
     *
     * Left preconditioning : solve MAx=Mb (the abort criterion controls the real residuals b-Ax)
     *
     *
     * \note Which variant gives the better convergence is highly problem dependent.
     *
     * \note Even if no preconditioning is used, the variant has to be set to left or right.
     *
     * \note In problems with few unknowns and a good preconditioner (MultiGrid) where the first outer iteration would converge
     *       a high polynomial degree (l) can lead to stability problems, because the convergence criterion is not calculated
     *       until l inner iterations. Reducing the polynomial degree should help in those cases.
     *
     * \see BiCGStabLPreconVariant
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     *
     * \author Jonas Duennebacke
     */
    template<typename Matrix_, typename Filter_>
    class BiCGStabL :
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

        /// arbitrary vector for the algorithm
        VectorType _vec_r0_tilde;
        /// vector for the algorithm
        VectorType _vec_u0;
        /// auxiliary vector for preconditioning
        VectorType _vec_pc;

        /// vector 'list' for the algorithm (size = l+1);
        std::vector<VectorType> _vec_uj_hat;

        /// vector 'list' for the algorithm (size = l+1); r[0] is the residual vector of the preconitioned system
        std::vector<VectorType> _vec_rj_hat;

        /// Precondition from left or right?
        BiCGStabLPreconVariant _precon_variant;

        /// Parameter l configuring the solver
        int _l;

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
         * \param[in] l
         * A parameter for the solver configuration, defaults to 2.
         *
         * \param[in] precond
         * A pointer to the preconditioner. May be \c nullptr.
         *
         * \param[in] precon_variant
         * Which preconditioning variant to use, defaults to left.
         */
        explicit BiCGStabL(const MatrixType& matrix, const FilterType& filter,
          int l = 2,                                                             //default parameter 2 (1 would lead to standard BiCGStab algorithm)
          std::shared_ptr<PrecondType> precond = nullptr,
          BiCGStabLPreconVariant precon_variant = BiCGStabLPreconVariant::left)
           :
          BaseClass("BiCGStab(" + stringify(l) +")", precond),
          _system_matrix(matrix),
          _system_filter(filter),
          _precon_variant(precon_variant),
          _l(l)
        {
          XASSERT(precon_variant == BiCGStabLPreconVariant::left || precon_variant == BiCGStabLPreconVariant::right);
          XASSERTM(_l > 0, "bicgstabl polynomial degree must be greather than zero!");
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

        explicit BiCGStabL(const String& section_name, PropertyMap* section,
          const MatrixType& matrix, const FilterType& filter,
          std::shared_ptr<PrecondType> precond = nullptr)
           :
          BaseClass("BiCGStabL", section_name, section, precond),
          _system_matrix(matrix),
          _system_filter(filter),
          _precon_variant(BiCGStabLPreconVariant::left)
        {
          // Check if we have set _p_variant
          auto p_variant_p = section->query("precon_variant");
          if(p_variant_p.second)
          {
            BiCGStabLPreconVariant precon_variant;
            precon_variant << p_variant_p.first;
            set_precon_variant(precon_variant);
          }

          //Check if we have set the solver parameter _l
          std::pair<String , bool> l_pair = section->query("polynomial_degree");

          int l_pm = 2;                                                 //use default l, if the parameter is not given
          if (l_pair.second)
          {
            l_pair.first.trim_me();

            l_pm = std::stoi(l_pair.first);

            XASSERTM(l_pm > 0, "bicgstabl polynomial degree must be greather than zero!");
          }
          this->_l = l_pm;
          if(this->_plot_name == "BiCGStabL")
            this->set_plot_name("BiCGStab(" + stringify(_l) + ")");
        }

        /// \copydoc SolverBase::name()
        virtual String name() const override
        {
          return "BiCGStabL";
        }

        /// \copydoc SolverBase::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();

          // create all temporary vectors
          _vec_r0_tilde = this->_system_matrix.create_vector_r();
          _vec_u0       = this->_system_matrix.create_vector_r();
          _vec_pc       = this->_system_matrix.create_vector_r();

          _vec_uj_hat.resize(Index(_l+1));
          _vec_rj_hat.resize(Index(_l+1));
          for (int i=0; i <= _l; i++)
          {
            _vec_rj_hat.at(Index(i)) = this->_system_matrix.create_vector_r();
            _vec_uj_hat.at(Index(i)) = this->_system_matrix.create_vector_r();
          }
        }

        /// \copydoc SolverBase::done_symbolic()
        virtual void done_symbolic() override
        {
          _vec_r0_tilde.clear();
          _vec_u0.clear();
          _vec_pc.clear();
          for (int i=0; i <= _l; i++)
          {
            _vec_rj_hat.at(Index(i)).clear();
            _vec_uj_hat.at(Index(i)).clear();
          }
        BaseClass::done_symbolic();
        }

        /// \copydoc IterativeSolver::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
        {
          // compute initial defect
          this->_system_matrix.apply(this->_vec_rj_hat.at(0), vec_sol, vec_rhs, -DataType(1));
          this->_system_filter.filter_def(this->_vec_rj_hat.at(0));

          // apply
          Status st(_apply_intern(vec_sol, vec_rhs));

          this->plot_summary(st);

          return st;
        }

        /// \copydoc SolverBase::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // save rhs vector as initial defect
          this->_vec_rj_hat.at(0).copy(vec_def);
          //this->_system_filter.filter_def(this->_vec_rj_hat[0]);

          // format solution vector
          vec_cor.format();

          Status st(_apply_intern(vec_cor, vec_def));

          this->plot_summary(st);

          return st;
        }


        /**
         * \brief Sets the preconditioning variant (left or right)
         *
         * \param[in] precon_variant
         * Which preconditioning variant to use.
         */
        virtual void set_precon_variant(BiCGStabLPreconVariant precon_variant)
        {
          XASSERT(precon_variant == BiCGStabLPreconVariant::left || precon_variant == BiCGStabLPreconVariant::right);
          _precon_variant = precon_variant;
        }

        /// \copydoc SolverBase::write_config()
        virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name) const override
        {
          XASSERT(parent != nullptr);

          PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

          my_section->add_entry("precon_variant", stringify(_precon_variant));
          my_section->add_entry("polynomial_degree", stringify(_l));
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
         * The right hand side vector. This is used for the calculation of the real residual in the case of right preconditioning
         *
         * \returns A status code.
         */
        Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
        {
          IterationStats pre_iter(*this);
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));


          int l (_l);

          //choose  vec_r0_tilde as vec_r0_hat (vec_rj_hat[0])
          _vec_r0_tilde.copy(_vec_rj_hat.at(0));

          std::vector<DataType>  tau               ( Index( l*(l-1) ) , DataType(0));             //vector for tau_ij (j from 0 to l-1 , i from 0 to j-1, so LESS THAN HALF of the vector gets used)
                                                                                                  //tau_ij = tau[j*(l-1) + i] (column major ordering)

          std::vector<DataType>  sigma             ( Index(l) , DataType(0) );                    //0-offset in implementation; 1-offset in pseudo-code
          std::vector<DataType>  gamma_j           ( Index(l) , DataType(0) );                    //0-offset in implementation; 1-offset in pseudo-code
          std::vector<DataType>  gamma_prime       ( Index(l) , DataType(0) );                    //0-offset in implementation; 1-offset in pseudo-code
          std::vector<DataType>  gamma_prime_prime ( Index(l) , DataType(0) );                    //0-offset in implementation; 1-offset in pseudo-code

          DataType rho_0(1), alpha(0), omega(1), rho_1(0), beta(0), gamma(0);


          const MatrixType& mat_sys(this->_system_matrix);
          const FilterType& fil_sys(this->_system_filter);
          Status status(Status::progress);

          // Compute initial defect
          status = this->_set_initial_defect(_vec_rj_hat.at(0), vec_sol);

          if (_precon_variant ==  BiCGStabLPreconVariant::left)
          {
            _vec_pc.copy(_vec_rj_hat.at(0));
            if(!this->_apply_precond(_vec_rj_hat.at(0), _vec_pc, fil_sys))
            {
              Statistics::add_solver_expression(
                std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }
          }


          pre_iter.destroy();
          while(status == Status::progress)
          {
            IterationStats stat(*this);
            rho_0 *= (-omega);


            for( int j(0); j < l ; j++)
            {
              rho_1 = _vec_rj_hat.at( Index(j) ).dot(_vec_r0_tilde);
              beta  = alpha * rho_1 / rho_0;
              rho_0 = rho_1;
              for ( int i(0) ; i <= j ; i++)
              {
                _vec_uj_hat.at( Index(i) ).axpy(_vec_uj_hat.at( Index(i) ), _vec_rj_hat.at( Index(i) ), -beta );
              }



              if (_precon_variant ==  BiCGStabLPreconVariant::left)
              {
                mat_sys.apply(_vec_pc, _vec_uj_hat.at(Index(j)) );
                fil_sys.filter_def(_vec_pc);
                if(!this->_apply_precond(_vec_uj_hat.at( Index(j+1) ), _vec_pc, fil_sys))
                {
                  Statistics::add_solver_expression(
                    std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
                  return Status::aborted;
                }
              }
              else
              {
                if(!this->_apply_precond(_vec_pc, _vec_uj_hat.at( Index(j) ), fil_sys))
                {
                  Statistics::add_solver_expression(
                    std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
                  return Status::aborted;
                }
                mat_sys.apply(_vec_uj_hat.at( Index(j+1) ), _vec_pc);
                fil_sys.filter_def(_vec_uj_hat.at( Index(j+1) ));
              }

              gamma = _vec_uj_hat.at( Index(j+1) ).dot(_vec_r0_tilde);
              alpha = rho_0 / gamma;

              for ( int i(0) ; i <= j ; i++)
              {
                _vec_rj_hat.at(Index(i)).axpy( _vec_uj_hat.at( Index(i+1) ) , _vec_rj_hat.at( Index(i) ), -alpha );
              }

              if (_precon_variant ==  BiCGStabLPreconVariant::left)
              {

                mat_sys.apply( _vec_pc, _vec_rj_hat.at(Index(j)) );
                fil_sys.filter_def(_vec_pc);
                if(!this->_apply_precond(_vec_rj_hat.at( Index(j+1) ) , _vec_pc, fil_sys))
                {
                  Statistics::add_solver_expression(
                    std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
                  return Status::aborted;
                }
              }
              else
              {
                if(!this->_apply_precond(_vec_pc, _vec_rj_hat.at(Index(j)), fil_sys))
                {
                  Statistics::add_solver_expression(
                    std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
                  return Status::aborted;
                }
                mat_sys.apply(_vec_rj_hat.at(Index(j+1)), _vec_pc);
                fil_sys.filter_def( _vec_rj_hat.at(Index(j+1)) );
              }

              vec_sol.axpy(_vec_uj_hat.at(0), vec_sol, alpha);
            }
            for ( int j(0) ; j < l ; j++)       //changed counter in comparison to pseudo code
            {

              for ( int i(0); i<j ; i++)        //changed counter
              {
                tau.at( Index(j*(l-1) + i) ) = ( _vec_rj_hat.at( Index(i+1) ).dot( _vec_rj_hat.at( Index(j+1) ) ) ) / sigma.at(Index(i));
                _vec_rj_hat.at(Index(j+1)).axpy(_vec_rj_hat.at(Index(i+1)), _vec_rj_hat.at(Index(j+1)), -tau.at( Index(j*(l-1) + i) ) );
              }
              sigma.at(Index(j)) = _vec_rj_hat.at( Index(j+1) ).dot( _vec_rj_hat.at(Index(j+1)) );
              gamma_prime.at(Index(j)) = ( _vec_rj_hat.at(0).dot( _vec_rj_hat.at(Index(j+1)) ) ) / sigma.at(Index(j));
            }
            gamma_j.at(Index(l-1)) = gamma_prime.at(Index(l-1));
            omega = gamma_j.at(Index(l-1));

            for ( int j(l-2) ; j >= 0 ; j--)  //changed counter
            {
              gamma_j.at(Index(j)) = gamma_prime.at(Index(j));
              for ( int i(j+1) ; i < l ; i++)
              {
                gamma_j.at(Index(j)) -= tau.at( Index(i*(l-1) + j) ) * gamma_j.at(Index(i));
              }
            }

            for ( int j(0); j < l-1 ; j++)                                      //changed counter
            {
              gamma_prime_prime.at( Index(j) ) = gamma_j.at( Index(j+1) );
              for ( int i(j+1) ; i < l-1 ; i++)                                  //changed counter
              {
                gamma_prime_prime.at(Index(j)) += tau.at(Index(i*(l-1) + j)) * gamma_j.at(Index(i+1));
              }
            }
                      vec_sol.axpy(_vec_rj_hat.at(0),           vec_sol,      gamma_j.at(0));
            _vec_rj_hat.at(0).axpy(_vec_rj_hat.at(Index(l)), _vec_rj_hat.at(0), -gamma_prime.at(Index(l-1)));
            _vec_uj_hat.at(0).axpy(_vec_uj_hat.at(Index(l)), _vec_uj_hat.at(0),     -gamma_j.at(Index(l-1)));

            for ( int j(1) ; j<l; j++)
            {
              _vec_uj_hat.at(0).axpy(_vec_uj_hat.at(Index(j)), _vec_uj_hat.at(0),           -gamma_j.at(Index(j-1)));
                        vec_sol.axpy(_vec_rj_hat.at(Index(j)),           vec_sol,  gamma_prime_prime.at(Index(j-1)));
              _vec_rj_hat.at(0).axpy(_vec_rj_hat.at(Index(j)), _vec_rj_hat.at(0),       -gamma_prime.at(Index(j-1)));
            }
            // Compute defect norm
            if (_precon_variant ==  BiCGStabLPreconVariant::left)
            {
              mat_sys.apply(_vec_pc, vec_sol);
              fil_sys.filter_def(_vec_pc);
              _vec_pc.axpy( vec_rhs, _vec_pc , -1);
              fil_sys.filter_def(_vec_pc);
              status = this->_set_new_defect(_vec_pc, vec_sol);
            }
            else
            {
              status = this->_set_new_defect(_vec_rj_hat.at(0), vec_sol);
            }

            if(status != Status::progress)
            {
              stat.destroy();
              Statistics::add_solver_expression(
                std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
              if (_precon_variant == BiCGStabLPreconVariant::right)
              {
                if(!this->_apply_precond(_vec_pc, vec_sol, fil_sys))
                {
                  Statistics::add_solver_expression(
                    std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
                  return Status::aborted;
                }
                vec_sol.copy(_vec_pc);
              }
              return status;
            }

          }
          // we should never reach this point...
          Statistics::add_solver_expression(
            std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return Status::undefined;
        }
    };

    /**
     * \brief Creates a new BiCGStabL solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] l
     * A parameter for the solver configuration
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \param[in] precon_variant
     * Which preconditioning variant to use, defaults to left
     *
     * \returns
     * A shared pointer to a new BiCGStabL object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStabL<Matrix_, Filter_>> new_bicgstabl(
      const Matrix_& matrix, const Filter_& filter, int l = 2 )
    {
      return std::make_shared<BiCGStabL<Matrix_, Filter_>>(matrix, filter, l, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<BiCGStabL<Matrix_, Filter_>> new_bicgstabl(
      const Matrix_& matrix, const Filter_& filter, int l,
      std::shared_ptr<Precond_> precond, BiCGStabLPreconVariant precon_variant = BiCGStabLPreconVariant::left)
    {
      return std::make_shared<BiCGStabL<Matrix_, Filter_>>(matrix, filter, l, precond, precon_variant);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStabL<Matrix_, Filter_>> new_bicgstabl(
      const Matrix_& matrix, const Filter_& filter, int l = 2,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr,
      BiCGStabLPreconVariant precon_variant = BiCGStabLPreconVariant::left)
    {
      return std::make_shared<BiCGStabL<Matrix_, Filter_>>(matrix, filter, l, precond, precon_variant);
    }
#endif

    /**
     * \brief Creates a new BiCGStabL solver object using a PropertyMap
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
     * A shared pointer to a new BiCGStabL object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr


#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStabL<Matrix_, Filter_>> new_bicgstabl(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<BiCGStabL<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<BiCGStabL<Matrix_, Filter_>> new_bicgstabl(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<BiCGStabL<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStabL<Matrix_, Filter_>> new_bicgstabl(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<BiCGStabL<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif

  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_BICGSTABL_HPP
