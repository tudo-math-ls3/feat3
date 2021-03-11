// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/solver/spai_precond.hpp>

#include <mkl.h>

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      void spai_mkl_solve_minimization_problem(float *Ahat, unsigned long lenJk_in, unsigned long rowlen_in, unsigned long mrl_in, float *rhs, float *sol)
      {
        MKL_INT lenJk = (MKL_INT)lenJk_in;
        MKL_INT rowlen = (MKL_INT)rowlen_in;
        MKL_INT mrl = (MKL_INT)mrl_in;

        // leading dimension of A (max allocated column size)
        const MKL_INT lda = rowlen*mrl;

        // storage for the qr deflections
        LAFEM::DenseVector<Mem::Main, float, unsigned long> tau(rowlen_in);

        // QR decomposition
        MKL_INT ret_error = LAPACKE_sgeqrf(LAPACK_COL_MAJOR, lenJk, rowlen, Ahat, lda, tau.elements());
        if(ret_error)
          XABORTM("SPAI: sgeqrf: Invalid argument no. " + stringify(-ret_error));

        // rhs <- Qt * rhs
        ret_error = LAPACKE_sormqr (LAPACK_COL_MAJOR, 'L', 'T', lenJk, 1, rowlen, Ahat, lda, tau.elements(), rhs, lenJk);
        if(ret_error)
          XABORTM("SPAI: sormqr: Invalid argument no. " + stringify(-ret_error));

        // rhs <- solution of R*x = rhs
        cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, rowlen, 1, 1.0f, Ahat, lda, rhs, rowlen);

        // copy result to SPAI
        cblas_scopy(rowlen, rhs, 1, sol, 1);
      }

      void spai_mkl_solve_minimization_problem(double *Ahat, unsigned long lenJk_in, unsigned long rowlen_in, unsigned long mrl_in, double *rhs, double *sol)
      {
        MKL_INT lenJk = (MKL_INT)lenJk_in;
        MKL_INT rowlen = (MKL_INT)rowlen_in;
        MKL_INT mrl = (MKL_INT)mrl_in;

        // leading dimension of A (max allocated column size)
        const MKL_INT lda = rowlen*mrl;

        // storage for the qr deflections
        LAFEM::DenseVector<Mem::Main, double, unsigned long> tau(rowlen_in);

        // QR decomposition
        MKL_INT ret_error = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, lenJk, rowlen, Ahat, lda, tau.elements());
        if(ret_error)
          XABORTM("SPAI: dgeqrf: Invalid argument no. " + stringify(-ret_error));

        // rhs <- Qt * rhs
        ret_error = LAPACKE_dormqr (LAPACK_COL_MAJOR, 'L', 'T', lenJk, 1, rowlen, Ahat, lda, tau.elements(), rhs, lenJk);
        if(ret_error)
          XABORTM("SPAI: dormqr: Invalid argument no. " + stringify(-ret_error));

        // rhs <- solution of R*x = rhs
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, rowlen, 1, 1.0f, Ahat, lda, rhs, rowlen);

        // copy result to SPAI
        cblas_dcopy(rowlen, rhs, 1, sol, 1);
      }
    } // namespace Intern
  } // namespace Solver
} // namespace FEAT
