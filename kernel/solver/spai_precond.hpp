#pragma once
#ifndef KERNEL_SOLVER_SPAI_PRECOND_HPP
#define KERNEL_SOLVER_SPAI_PRECOND_HPP 1

// includes, FEAT
#include <kernel/solver/precon_wrapper.hpp>

namespace FEAT
{
  namespace Solver
  {
    /// \todo reimplement this
    template<typename Matrix_, typename Filter_>
    using SPAIPrecond = Solver::PreconWrapper<Matrix_, Filter_, Solver::SPAIPreconditioner>;

    /**
     * \brief Creates a new SPAIPrecond solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \todo adjust arguments to reimplemented class
     *
     * \returns
     * A shared pointer to a new SPAIPrecond object.
     */
    template<typename Matrix_, typename Filter_, typename... Args_>
    inline std::shared_ptr<SPAIPrecond<Matrix_, Filter_>> new_spai_precond(
      const Matrix_& matrix, const Filter_& filter, Args_&&... args)
    {
      return std::make_shared<SPAIPrecond<Matrix_, Filter_>>
        (filter, matrix, std::forward<Args_>(args)...);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_SPAI_PRECOND_HPP
