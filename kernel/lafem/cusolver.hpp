#pragma once
#ifndef KERNEL_LAFEM_CUSOLVER_HPP
#define KERNEL_LAFEM_CUSOLVER_HPP 1
#include <kernel/base_header.hpp>

#ifdef FEAST_HAVE_CUSOLVER
#ifndef  __CUDACC__
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#endif

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Cuda LU sparse solver wrapper class
     *
     * \author Dirk Ribbrock
     */
    class CuSolverLU
    {
      public:
#ifndef  __CUDACC__
      static void solve(DenseVector<Mem::Main, double, unsigned int> & x, SparseMatrixCSR<Mem::Main, double, unsigned int> & A,
        DenseVector<Mem::Main, double, unsigned int> & b)
      {
        solve_intern((int)x.size(), (int)A.used_elements(), A.val(), (const int*) A.row_ptr(), (const int*) A.col_ind(), b.elements(), x.elements());
      }
#endif

      static void solve_intern(int n, int nnzA, const double * csrValA, const int * csrRowPtrA, const int * csrColIndA,
        const double * b, double * x);
    };

    /**
     * \brief Cuda QR sparse solver wrapper class
     *
     * \author Dirk Ribbrock
     */
    class CuSolverQR
    {
      public:
#ifndef  __CUDACC__
      static void solve(DenseVector<Mem::CUDA, double, unsigned int> & x, SparseMatrixCSR<Mem::CUDA, double, unsigned int> & A,
        DenseVector<Mem::CUDA, double, unsigned int> & b)
      {
        solve_intern((int)x.size(), (int)A.used_elements(), A.val(), (const int*) A.row_ptr(), (const int*) A.col_ind(), b.elements(), x.elements());
      }
#endif

      static void solve_intern(int m, int nnz, const double * csrValA, const int * csrRowPtrA, const int * csrColIndA,
        const double * b, double * x);
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // FEAST_HAVE_CUSOLVER
#endif // KERNEL_LAFEM_CUSOLVER_HPP
