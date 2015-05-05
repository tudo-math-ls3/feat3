#pragma once
#ifndef KERNEL_LAFEM_UMFPACK_HPP
#define KERNEL_LAFEM_UMFPACK_HPP 1
#include <kernel/base_header.hpp>

#ifdef FEAST_HAVE_UMFPACK
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/exception.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief UMFPACK solver wrapper class
     *
     * \author Peter Zajac
     */
    class Umfpack
    {
    public:
      /// compatible matrix type
      typedef SparseMatrixCSR<Mem::Main, double> MatrixType;
      /// compatible vector type
      typedef DenseVector<Mem::Main, double> VectorType;

      /**
       * \brief Base-class for UMFPACK errors
       *
       * An instance of this class may be thrown if there is no derived class
       * which suits the error that occurred.
       */
      class Error : public Exception
      {
      public:
        explicit Error(const String& msg) : Exception(msg) {}
      };

      /**
       * \brief Out-Of-Memory error
       *
       * This exception is self-explanatory: UMFPACK went out of memory.
       */
      class OutOfMemoryError : public Error
      {
      public:
        OutOfMemoryError() : Error("out of memory") {}
      };

      /**
       * \brief Invalid matrix structure error
       *
       * This exception is thrown if the matrix structure is invalid.\n
       * According to UMFPACK's documentation, this may have one of the following reasons:
       * - Number of non-zero entries is negative
       * - row_ptr[0] is nonzero
       * - row_ptr[i+1] < row_ptr[i] for some row i
       * - a row_ptr entry is out-of-bounds
       * - duplicate entries in a matrix row
       * - matrix entries are not sorted by ascending column indices
       */
      class InvalidMatrixError : public Error
      {
      public:
        InvalidMatrixError() : Error("invalid matrix structure") {}
      };

      /**
       * \brief Singular matrix error
       *
       * This error is self-explanatory: The matrix is singular.
       */
      class SingularMatrixError : public Error
      {
      public:
        SingularMatrixError() : Error("singular matrix") {}
      };

    private:
      /// system matrix pointer
      const MatrixType* _system_matrix;
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
      /// constructor
      explicit Umfpack();

      /**
       * \brief Constructor
       *
       * This constructor automatically calls init() to factorise the input matrix.
       *
       * \param[in] system_matrix
       * A pointer to the system matrix to be factorised.
       */
      explicit Umfpack(const MatrixType* system_matrix);

      /// virtual destructor
      virtual ~Umfpack();

      /**
       * \brief Performs symbolic factorisation.
       *
       * \param[in] system_matrix
       * A pointer to the system matrix to be factorised. Must be valid until free_symbolic() is called.
       *
       * \param[in] ignore_data
       * Specifies whether the symbolic factorisation shall ignore the value data of the matrix.
       * Should be set to \c true if the value array of the matrix is not yet initialised.
       */
      void init_symbolic(const MatrixType* system_matrix, bool ignore_data = false);

      /**
       * \brief Performs numeric factorisation.
       *
       * This function performs the numeric factorisation of the system matrix that has been assigned
       * by a prior init_symbolic call.
       */
      void init_numeric();

      /**
       * \brief Performs both symbolic and numeric factorisation.
       *
       * \param[in] system_matrix
       * A pointer to the system matrix to be factorised. Must be valid until free_symbolic() is called.
       */
      void init(const MatrixType* system_matrix);

      /**
       * \brief Frees the numeric factorisation.
       */
      void free_numeric();

      /**
       * \brief Frees the symbolic factorisation.
       */
      void free_symbolic();

      /**
       * \brief Frees both the numeric and symbolic factorisation.
       */
      void free();

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
      void solve(VectorType& x, const VectorType& b);
    }; // class Umfpack
  } // namespace LAFEM
} // namespace FEAST

#endif // FEAST_HAVE_UMFPACK
#endif // KERNEL_LAFEM_UMFPACK_HPP
