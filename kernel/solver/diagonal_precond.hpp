// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_DIAGONAL_PRECOND_HPP
#define KERNEL_SOLVER_DIAGONAL_PRECOND_HPP 1

// includes, FEAT
#include <kernel/solver/base.hpp>

namespace FEAT
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
      typedef SolverBase<Vector_> BaseClass;
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

      explicit DiagonalPrecond(const String& section_name, PropertyMap* section,
      const VectorType& diag, const FilterType& filter) :
        BaseClass(section_name, section),
        _diag(diag),
        _filter(filter)
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Diagonal";
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
} // namespace FEAT

#endif // KERNEL_SOLVER_DIAGONAL_PRECOND_HPP
