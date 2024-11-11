// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include <kernel/base_header.hpp>

#include <kernel/solver/base.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
        int cuda_lu(int n, int nnzA, const double * csrValA, const int * csrRowPtrA, const int * csrColIndA,
        const double * b, double * x);

        int cuda_qr(int m, int nnz, const double * csrValA, const int * csrRowPtrA, const int * csrColIndA,
        const double * b, double * x);
    }

    /**
     * \brief CuSolverLU solver class
     *
     * This class provides an implementation of the SolverBase interface using the
     * direct cuda solver cusolverSpDcsrlsvluHost for doing the actual dirty work.
     *
     * \attention
     * This class is only declared if FEAT was configured to build and link against
     * cuda and the cuda sdk ships with the cusolver library (cuda version >= 7).
     */
    class CuSolverLU :
      public SolverBase<LAFEM::DenseVector<double, unsigned int>>
    {
      private:
        /// system matrix
        const LAFEM::SparseMatrixCSR<double, unsigned int> & _system_matrix;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] system_matrix
         * A reference to the system matrix to be factorized.
         */
        explicit CuSolverLU(const LAFEM::SparseMatrixCSR<double, unsigned int> & system_matrix) :
          _system_matrix(system_matrix)
          {
          }

        /// Returns the name of the solver.
        virtual String name() const override
        {
          return "CuSolverLU";
        }

        /**
         * \brief Solves a linear system with the factorized system matrix.
         *
         * \param[in,out] x
         * A reference to the solution vector. The vector must be allocated to the correct length, but its
         * initial contents are ignored.
         *
         * \param[in] b
         * A reference to the right-hand-side of the linear system.
         */
        virtual Status apply(LAFEM::DenseVector<double, unsigned int> & x, const LAFEM::DenseVector<double, unsigned int> & b) override
          {
            int status = Intern::cuda_lu((int)x.size(), (int)_system_matrix.used_elements(), _system_matrix.val(), (const int*) _system_matrix.row_ptr(), (const int*) _system_matrix.col_ind(), b.elements(), x.elements());

            return (status == 0) ? Status::success :  Status::aborted;
          }
    };

    /**
     * \brief CuSolverQR solver class
     *
     * This class provides an implementation of the SolverBase interface using the
     * direct cuda solver cusolverSpDcsrlsvqr for doing the actual dirty work.
     *
     * \attention
     * This class is only declared if FEAT was configured to build and link against
     * cuda and the cuda sdk ships with the cusolver library (cuda version >= 7).
     */
    class CuSolverQR :
      public SolverBase<LAFEM::DenseVector<double, unsigned int>>
    {
      private:
        /// system matrix
        const LAFEM::SparseMatrixCSR<double, unsigned int> & _system_matrix;
      public:
        /**
         * \brief Constructor
         *
         * \param[in] system_matrix
         * A reference to the system matrix to be factorized.
         */
        explicit CuSolverQR(const LAFEM::SparseMatrixCSR<double, unsigned int> & system_matrix) :
          _system_matrix(system_matrix)
        {
        }

        /// Returns the name of the solver.
        virtual String name() const override
        {
          return "CuSolverQR";
        }

        virtual Status apply(LAFEM::DenseVector<double, unsigned int> & x, const LAFEM::DenseVector<double, unsigned int> & b) override
        {
          int status = Intern::cuda_qr((int)x.size(), (int)_system_matrix.used_elements(), _system_matrix.val(), (const int*) _system_matrix.row_ptr(), (const int*) _system_matrix.col_ind(), b.elements(), x.elements());

          return (status == 0) ? Status::success :  Status::aborted;
        }
    };
  } // namespace LAFEM
} // namespace FEAT
