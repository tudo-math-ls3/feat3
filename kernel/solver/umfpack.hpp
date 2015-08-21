#pragma once
#ifndef KERNEL_SOLVER_UMFPACK_HPP
#define KERNEL_SOLVER_UMFPACK_HPP 1

// includes, FEAST
#include <kernel/solver/base.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace FEAST
{
  namespace Solver
  {
#if defined(FEAST_HAVE_UMFPACK) || defined(DOXYGEN)
    /**
     * \brief UMFPACK solver class
     *
     * This class provides an implementation of the SolverBase interface using the
     * direct solver UMFPACK for doing the actual dirty work.
     *
     * \attention
     * This class is only declared if FEAST was configured to build and link against
     * the \c UMFPACK third-party library.
     *
     * \author Peter Zajac
     */
    class Umfpack :
      public SolverBase<LAFEM::DenseVector<Mem::Main, double, Index>>
    {
    public:
      /// compatible matrix type
      typedef LAFEM::SparseMatrixCSR<Mem::Main, double, Index> MatrixType;
      /// compatible vector type
      typedef LAFEM::DenseVector<Mem::Main, double, Index> VectorType;

      /// our base class
      typedef SolverBase<VectorType> BaseClass;

    private:
      /// system matrix
      const MatrixType& _system_matrix;
      /// umfpack control array
      double* _umf_control;
      /// umfpack symbolic factorisation pointer
      void* _umf_symbolic;
      /// umfpack numeric factorisation pointer
      void* _umf_numeric;

      /// symbolic peak memory size
      std::size_t _sym_peak_size;
      /// symbolic factorisation memory size
      std::size_t _sym_mem_size;
      /// numeric factorisation memory size
      std::size_t _num_mem_size;
      /// total peak memory size
      std::size_t _umf_peak_size;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] system_matrix
       * A reference to the system matrix to be factorised.
       */
      explicit Umfpack(const MatrixType& system_matrix);

      /// virtual destructor
      virtual ~Umfpack();

      virtual void init_symbolic() override;
      virtual void done_symbolic() override;
      virtual void init_numeric() override;
      virtual void done_numeric() override;

      /**
       * \brief Solves a linear system with the factorised system matrix.
       *
       * \param[in,out] x
       * A reference to the solution vector. The vector must be allocated to the correct length, but its
       * initial contents are ignored.
       *
       * \param[in] b
       * A reference to the right-hand-side of the linear system.
       */
      virtual Status apply(VectorType& x, const VectorType& b) override;
    }; // class Umfpack
#endif // defined(FEAST_HAVE_UMFPACK) || defined(DOXYGEN)
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_UMFPACK_HPP
