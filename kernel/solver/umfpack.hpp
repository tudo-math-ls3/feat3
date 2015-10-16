#pragma once
#ifndef KERNEL_SOLVER_UMFPACK_HPP
#define KERNEL_SOLVER_UMFPACK_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/mean_filter.hpp>

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

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Umfpack";
      }

      virtual void init_symbolic() override;
      virtual void done_symbolic() override;
      virtual void init_numeric() override;
      virtual void done_numeric() override;

      /**
       * \brief Solves a linear system with the factorised system matrix.
       *
       * \param[in,out] vec_sol
       * A reference to the solution vector. The vector must be allocated to the correct length, but its
       * initial contents are ignored.
       *
       * \param[in] vec_rhs
       * A reference to the right-hand-side of the linear system.
       */
      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override;
    }; // class Umfpack

    /**
     * \brief UMFPACK Mean solver class
     *
     * This class implements a variant of the Umfpack solver, which is capable of solving linear
     * systems including the integral mean contraint or any other constraint that can be expressed
     * as a scalar Lagrange multiplier to the original system matrix.
     *
     * As Umfpack is a direct solver, it cannot directly utilise the MeanFilter class, which is used
     * by iterative solvers to enforce the integral mean constraint. Therefore, this class implements
     * an algorithm, which extends the original system matrix by a Lagrange multiplier vector and
     * applied the Umfpack solver onto this extended linear system.
     *
     * This class offers two constructors:
     * - A CTOR which takes a SparseMatrixCSR and a MeanFilter as input
     * - A CTOR which takes a SparseMatrixCSR and a DenseVector representing the Lagrange multiplier
     *   as input.
     *
     * \author Peter Zajac
     */
    class UmfpackMean :
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
      /// weight vector
      const VectorType& _weight_vector;
      /// our extended system matrix
      MatrixType _solver_matrix;
      /// two temporary extended vectors
      VectorType _vec_x, _vec_b;
      /// our internal Umfpack solver object
      Umfpack _umfpack;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] system_matrix
       * A reference to the system matrix to be factorised.
       *
       * \param[in] weight_vector
       * The weight vector to be used as a Lagrange multiplier.
       */
      explicit UmfpackMean(const MatrixType& system_matrix, const VectorType& weight_vector);

      /**
       * \brief Constructor
       *
       * \param[in] system_matrix
       * A reference to the system matrix to be factorised.
       *
       * \param[in] mean_filter
       * A reference to the mean filter containing the weight vector.
       */
      explicit UmfpackMean(
        const MatrixType& system_matrix,
        const LAFEM::MeanFilter<Mem::Main, double, Index>& mean_filter) :
        UmfpackMean(system_matrix, mean_filter.get_vec_dual())
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "UmfpackMean";
      }

      virtual void init_symbolic() override;
      virtual void done_symbolic() override;
      virtual void init_numeric() override;
      virtual void done_numeric() override;

      /**
       * \brief Solves a linear system with the factorised system matrix.
       *
       * \param[in,out] vec_sol
       * A reference to the solution vector. The vector must be allocated to the correct length, but its
       * initial contents are ignored.
       *
       * \param[in] vec_rhs
       * A reference to the right-hand-side of the linear system.
       */
      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override;
    }; // class UmfpackMean
#endif // defined(FEAST_HAVE_UMFPACK) || defined(DOXYGEN)
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_UMFPACK_HPP
