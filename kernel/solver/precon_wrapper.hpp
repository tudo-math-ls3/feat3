#pragma once
#ifndef KERNEL_SOLVER_PRECON_WRAPPER_HPP
#define KERNEL_SOLVER_PRECON_WRAPPER_HPP 1

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/lafem/preconditioner.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Provisional LAFEM Preconditioner Wrapper class template
     *
     * This class template acts as a wrapper around the preconditioners implemented in the
     * <c>lafem/preconditioner.hpp</c> header file.
     *
     * \tparam Matrix_
     * The matrix class; is passed as the first parameter to the preconditioner class template.
     *
     * \tparam Precond_
     * The preconditioner class template.
     *
     * <b>Example</b>:\n
     * To use the ILUPreconditioner class for CSR-matrices, one would have to use the following class template
     * combination:
     * <c>PreconWrapper<SparseMatrixCSR<Mem::Main,double>, ILUPreconditioner></c>.
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_,
      template<typename,typename> class Precon_>
    class PreconWrapper :
      public SolverBase<typename Matrix_::VectorTypeR>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef typename MatrixType::VectorTypeR VectorType;

    protected:
      /// the filter object
      const Filter_& _filter;
      /// the actual preconditioner object
      Precon_<MatrixType, VectorType> _precond;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] args
       * The arguments which are passed to the preconditioner object constructor.
       * For the required set of arguments, see the documentation of the corresponding
       * preconditioner class template.
       */
      template<typename... Args_>
      explicit PreconWrapper(const Filter_& filter, Args_&&... args) :
        _filter(filter),
        _precond(std::forward<Args_>(args)...)
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return _precond.name();
      }

      /// Applies the preconditioner.
      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        _precond.apply(vec_cor, vec_def);
        this->_filter.filter_cor(vec_cor);
        return Status::success;
      }
    }; // class PreconWrapper<...>
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PRECON_WRAPPER_HPP
