// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/solver/base.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Jacobi preconditioner implementation
     *
     * This class implements a simple damped Jacobi preconditioner.
     *
     * This implementation works for the following matrix types and combinations thereof:
     * - all LAFEM::SparseMatrix types
     * - LAFEM::DenseMatrix
     * - LAFEM::PowerDiagMatrix
     * - LAFEM::PowerFullMatrix
     * - LAFEM::TupleDiagMatrix
     * - Global::Matrix
     *
     * Moreover, this implementation supports all data and index types.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class JacobiPrecond :
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
      /// The damping parameter
      DataType _omega;
      /// The component-wise inverted diagonal of _matrix
      VectorType _inv_diag;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix whose main diagonal is to be used.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] omega
       * The damping parameter for the preconditioner.
       */
      explicit JacobiPrecond(const MatrixType& matrix, const FilterType& filter, DataType omega = DataType(1)) :
        _matrix(matrix),
        _filter(filter),
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
       * A shared pointer to a new JacobiPrecond object.
       */
      explicit JacobiPrecond(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter) :
        BaseClass(section_name, section),
        _matrix(matrix),
        _filter(filter),
        _omega(1)
      {
        auto omega_p = section->query("omega");
        if(omega_p.second && (!omega_p.first.parse(this->_omega) || (this->_omega <= DataType(0))))
          throw ParseError(section_name + ".omega", omega_p.first, "a positive float");
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~JacobiPrecond()
      {
      }

      /// \copydoc SolverBase::name()
      virtual String name() const override
      {
        return "Jacobi";
      }

      /// \copydoc SolverBase::init_symbolic()
      virtual void init_symbolic() override
      {
        _inv_diag = _matrix.create_vector_r();
      }

      /// \copydoc SolverBase::done_symbolic()
      virtual void done_symbolic() override
      {
        _inv_diag.clear();
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

      /// \copydoc SolverBase::apply()
      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        vec_cor.component_product(_inv_diag, vec_def);
        this->_filter.filter_cor(vec_cor);
        return Status::success;
      }
    }; // class JacobiPrecond<...>

    /**
     * \brief Creates a new JacobiPrecond solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] omega
     * The damping parameter for Jacobi.
     *
     * \returns
     * A shared pointer to a new JacobiPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<JacobiPrecond<Matrix_, Filter_>> new_jacobi_precond(
      const Matrix_& matrix, const Filter_& filter,
      const typename Matrix_::DataType omega = typename Matrix_::DataType(1))
    {
      return std::make_shared<JacobiPrecond<Matrix_, Filter_>>(matrix, filter, omega);
    }

    /**
     * \brief Creates a new JacobiPrecond solver object using a PropertyMap
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
     * \returns
     * A shared pointer to a new JacobiPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<JacobiPrecond<Matrix_, Filter_>> new_jacobi_precond(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<JacobiPrecond<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }

  } // namespace Solver
} // namespace FEAT
