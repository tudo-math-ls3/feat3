#pragma once
#ifndef KERNEL_SOLVER_POLYNOMIAL_PRECOND_HPP
#define KERNEL_SOLVER_POLYNOMIAL_PRECOND_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Polynomial preconditioner implementation
     *
     * This class represents the Neumann-Polynomial-Preconditioner \f$M^{-1} = \sum_{k=0}^m (I - \tilde M^{-1}A)^k \tilde M^{-1}\f$
     *
     * \note As of now, M is a hardcoded Jacobi (main diagonal inverse)
     *
     */
    template<typename Matrix_, typename Filter_>
    class PolynomialPrecond :
      public SolverBase<typename Matrix_::VectorTypeL>
    {
    public:
      /// The matrix type
      typedef Matrix_ MatrixType;
      /// The filter type
      typedef Filter_ FilterType;
      /// The type of vector this solver can be applied to
      typedef typename MatrixType::VectorTypeL VectorType;
      /// The floating point precision
      typedef typename MatrixType::DataType DataType;
      /// Our base class
      typedef SolverBase<VectorType> BaseClass;

    protected:
      /// The system matrix
      const MatrixType& _matrix;
      /// The filter for projecting solution and defect to subspaces
      const FilterType& _filter;
      /// order m of preconditioner
      Index _m;
      /// The damping parameter for the internal jacobi preconditioner.
      DataType _omega;
      /// The component-wise inverted diagonal of _matrix
      VectorType _inv_diag;
      /// auxilary vectors
      VectorType _aux1, _aux2, _aux3;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The source matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] m
       * The order of the polynom
       *
       * \param[in] omega
       * The damping parameter for the internal jacobi preconditioner.
       */
      explicit PolynomialPrecond(const MatrixType& matrix, const FilterType& filter, Index m, DataType omega = DataType(1)) :
        _matrix(matrix),
        _filter(filter),
        _m(m),
        _omega(omega)
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
       * \returns
       * A shared pointer to a new PolynomialPrecond object.
       */
      explicit PolynomialPrecond(const String& section_name, PropertyMap* section,
      const MatrixType& matrix, const FilterType& filter) :
        BaseClass(section_name, section),
        _matrix(matrix),
        _filter(filter),
        _m(1),
        _omega(1)
        {
          auto omega_p = section->query("omega");
          if(omega_p.second)
          {
            set_omega(DataType(std::stod(omega_p.first)));
          }
          else
          {
            throw InternalError(__func__,__FILE__,__LINE__,
            name()+" config section is missing the mandatory omega key!");
          }

          auto m_p = section->query("m");
          if(m_p.second)
          {
            set_m(std::stoul(m_p.first));
          }
          else
          {
            throw InternalError(__func__,__FILE__,__LINE__,
            name()+" config section is missing the mandatory m key!");
          }
        }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~PolynomialPrecond()
      {
      }

      /// \copydoc SolverBase::name()
      virtual String name() const override
      {
        return "Polynomial";
      }

      /// \copydoc SolverBase::init_symbolic()
      virtual void init_symbolic() override
      {
        _inv_diag = _matrix.create_vector_r();
        _aux1 = _matrix.create_vector_r();
        _aux2 = _matrix.create_vector_r();
        _aux3 = _matrix.create_vector_r();
      }

      /// \copydoc SolverBase::done_symbolic()
      virtual void done_symbolic() override
      {
        _inv_diag.clear();
        _aux1.clear();
        _aux2.clear();
        _aux3.clear();
      }

      /// \copydoc SolverBase::init_numeric()
      virtual void init_numeric() override
      {
        // extract matrix diagonal
        _matrix.extract_diag(_inv_diag);

        // invert diagonal elements
        _inv_diag.component_invert(_inv_diag, _omega);
      }

      /**
       * \brief Sets the damping parameter
       *
       * \param[in] omega
       * The new damping parameter.
       *
       */
      void set_omega(DataType omega)
      {
        XASSERT(omega > DataType(0));
        _omega = omega;
      }

      /**
       * \brief Sets the polynomial order parameter
       *
       * \param[in] m
       * The new polynomial order parameter.
       *
       */
      void set_m(Index m)
      {
        XASSERT(m > DataType(0));
        _m = m;
      }

      /// \copydoc SolverBase::write_config()
      virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name) const override
      {

        PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

        my_section->add_entry("omega", stringify_fp_sci(_omega));
        my_section->add_entry("m", stringify(_m));

        return my_section;
      }

      /// \copydoc SolverBase::apply()
      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        /*
         * preconditioner is given by
         *   \f$ M^-1 = \left(I + (I - \tilde M^{-1}A) + ... + (I - \tilde M^{-1 }A)^m\right) \tilde M^{-1} \f$
         *
         * the preconditioner only works, if
         *   ||I - \tilde M^{-1} A||_2 < 1.
         */

        vec_cor.component_product(_inv_diag, vec_def);
        _aux3.copy(vec_cor);

        for (Index i = 1; i <= _m; ++i)
        {
          _matrix.apply(_aux1, vec_cor);
          _filter.filter_def(_aux1);
          _aux2.component_product(_inv_diag, _aux1);
          vec_cor.axpy(vec_cor, _aux3);
          vec_cor.axpy(_aux2, vec_cor, DataType(-1.0));
        }

        this->_filter.filter_cor(vec_cor);

        return Status::success;
      }
    }; // class PolynomialPrecond<...>

    /**
     * \brief Creates a new PoynomialPrecond solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] m
     * The order of the polynom
     *
     * \param[in] omega
     * The damping parameter for the internal jacobi preconditioner.
     *
     * \returns
     * A shared pointer to a new PolynomialPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PolynomialPrecond<Matrix_, Filter_>> new_polynomial_precond(
      const Matrix_& matrix, const Filter_& filter,
      const Index m, const typename Matrix_::DataType omega = typename Matrix_::DataType(1))
    {
      return std::make_shared<PolynomialPrecond<Matrix_, Filter_>>(matrix, filter, m, omega);
    }

    /**
     * \brief Creates a new PolynomialPrecond solver object using a PropertyMap
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
     * \returns
     * A shared pointer to a new PolynomialPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PolynomialPrecond<Matrix_, Filter_>> new_polynomial_precond(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PolynomialPrecond<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_POLYNOMIAL_PRECOND_HPP
