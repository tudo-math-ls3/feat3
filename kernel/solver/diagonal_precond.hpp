#pragma once
#ifndef KERNEL_SOLVER_DIAGONAL_PRECOND_HPP
#define KERNEL_SOLVER_DIAGONAL_PRECOND_HPP 1

// includes, FEAST
#include <kernel/solver/base.hpp>

namespace FEAST
{
  namespace Solver
  {
    /**
     * \brief Diagonal preconditioner implementation
     *
     * This class implements a simple Diagonal preconditioner.
     *
     * \author Peter Zajac
     */
    template<typename Vector_, typename Filter_>
    class DiagonalPrecond :
      public SolverBase<Vector_>
    {
    public:
      typedef Vector_ VectorType;
      typedef Filter_ FilterType;

    protected:
      const VectorType& _diag;
      const FilterType& _filter;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] diag
       * The vector representing a diagonal matrix.
       *
       * \param[in] filter
       * The filter.
       */
      explicit DiagonalPrecond(const VectorType& diag, const FilterType& filter) :
        _diag(diag),
        _filter(filter)
      {
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        vec_cor.component_product(_diag, vec_def);
        this->_filter.filter_cor(vec_cor);
        return Status::success;
      }
    }; // class DiagonalPrecond<...>

    /**
     * \brief Creates a new DiagonalPrecond solver object
     *
     * \param[in] diag
     * The vector representing the diagonal matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] omega
     * The damping parameter for Jacobi.
     *
     * \returns
     * A shared pointer to a new DiagonalPrecond object.
     */
    template<typename Vector_, typename Filter_>
    inline std::shared_ptr<DiagonalPrecond<Vector_, Filter_>> new_diagonal_precond(
      const Vector_& diag, const Filter_& filter)
    {
      return std::make_shared<DiagonalPrecond<Vector_, Filter_>>(diag, filter);
    }
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_DIAGONAL_PRECOND_HPP