#pragma once
#ifndef KERNEL_SOLVER_UMFPACK_HPP
#define KERNEL_SOLVER_UMFPACK_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/mean_filter.hpp>

namespace FEAT
{
  namespace Solver
  {
#if defined(FEAT_HAVE_UMFPACK) || defined(DOXYGEN)
    /**
     * \brief UMFPACK solver class
     *
     * This class provides an implementation of the SolverBase interface using the
     * direct solver UMFPACK for doing the actual dirty work.
     *
     * \note
     * This solver can only be applied onto SparseMatrixCSR<Mem::Main,double,Index> matrices.
     * If you want to apply UMFPACK on other matrix types, use the GenericUmfpack solver instead.
     *
     * \attention
     * This class is only declared if FEAT was configured to build and link against
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

    /**
     * \brief Generic UMFPACK solver class
     *
     * This class effectively wraps around an Umfpack solver object
     *
     * \attention
     * This class is only declared if FEAT was configured to build and link against
     * the \c UMFPACK third-party library.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_>
    class GenericUmfpack :
      public SolverBase<typename Matrix_::VectorTypeL>
    {
    public:
      /// our base class
      typedef SolverBase<typename Matrix_::VectorTypeL> BaseClass;
      /// our matrix type
      typedef Matrix_ MatrixType;
      /// our vector type
      typedef typename MatrixType::VectorTypeL VectorType;

    protected:
      /// our matrix
      const MatrixType& _matrix;
      /// the matrix for our Umfpack solver (SparseMatrixCSR<Mem::Main, double, Index>)
      typename Umfpack::MatrixType _umf_matrix;
      /// the vectors for our Umfpack solver (DenseVector<Mem::Main, double, Index>)
      typename Umfpack::VectorType _umf_vsol, _umf_vrhs;
      /// the actual Umpfack solver object
      Umfpack _umfpack;

    public:
      explicit GenericUmfpack(const MatrixType& matrix) :
        _matrix(matrix),
        _umfpack(_umf_matrix)
      {
      }

      /// virtual destructor
      virtual ~GenericUmfpack()
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "GenericUmfpack";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // convert matrix to obtain the structure
        _umf_matrix.convert(_matrix);

        // create vectors
        _umf_vsol = _umf_matrix.create_vector_l();
        _umf_vrhs = _umf_matrix.create_vector_l();

        // factorise symbolic
        _umfpack.init_symbolic();
      }

      virtual void done_symbolic() override
      {
        _umfpack.done_symbolic();

        _umf_vrhs.clear();
        _umf_vsol.clear();
        _umf_matrix.clear();

        BaseClass::done_symbolic();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // convert our system matrix to <Mem::Main, double, Index> (if necessary)
        typename MatrixType::template ContainerType<Mem::Main, double, Index> mat_main;
        mat_main.convert(_matrix);

        // get the array of our CSR matrix
        const Index* row_ptr = _umf_matrix.row_ptr();
        /*const*/ Index* col_idx = _umf_matrix.col_ind();
        double* val = _umf_matrix.val();

        // copy entries into our CSR matrix
        for(Index i(0); i < _umf_matrix.rows(); ++i)
          mat_main.set_line(i, val + row_ptr[i], col_idx + row_ptr[i], 0);

        // factorise
        _umfpack.init_numeric();
      }

      virtual void done_numeric() override
      {
        _umfpack.done_numeric();
      }

      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // convert RHS vector
        _umf_vrhs.copy(vec_rhs);

        // solve
        Status status = _umfpack.apply(_umf_vsol, _umf_vrhs);

        // convert sol vector
        _umf_vsol.copy_inv(vec_sol);

        return status;
      }
    }; // class GenericUmfpack<...>

    /**
     * \brief Creates a new GenericUmfpack solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \returns
     * A shared pointer to a new GenericUmfpack object.
     */
    template<typename Matrix_>
    inline std::shared_ptr<GenericUmfpack<Matrix_>> new_generic_umfpack(const Matrix_& matrix)
    {
      return std::make_shared<GenericUmfpack<Matrix_>> (matrix);
    }

#endif // defined(FEAT_HAVE_UMFPACK) || defined(DOXYGEN)
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_UMFPACK_HPP
