// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_csr_dense.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/exception.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
#include <mkl_spblas.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  // the code below casts Index* to MKL_INT*
  static_assert(sizeof(Index) == sizeof(MKL_INT), "size mismatch between MKL_INT and FEAT::Index");

  void MatVecMultCSRDense::exec_mkl_impl(float * r, const float alpha, const float * const x, const float * const y, const float * const val, const Index * const col_idx, const Index * const row_ptr, const Index rows, const Index cols, const Index, const bool transposed)
  {
    MKL_INT mrows = (MKL_INT)rows;
    MKL_INT mcolumns = (MKL_INT)cols;

    if(y == nullptr)
      Memory::memset_main(r, 0, sizeof(float) * (transposed ? cols : rows));
    else if(r != y)
    {
      MKL_INT one = 1;
      MKL_INT n = transposed ? mcolumns : mrows;
      scopy(&n, y, &one, r, &one);
    }

    float beta = 1.0f;

  #ifndef FEAT_USE_MKL_LEGACY_SPMV

    sparse_operation_t opt = (transposed ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE);

    sparse_matrix_t A;
    FEAT_DISABLE_WARNINGS
    sparse_status_t status = mkl_sparse_s_create_csr(&A, SPARSE_INDEX_BASE_ZERO, mrows, mcolumns, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_idx, (float*) val);
    FEAT_RESTORE_WARNINGS
    if (status != SPARSE_STATUS_SUCCESS)
      XABORTM("MKL Sparse Error occurred in execution!\n");
    matrix_descr md;
    md.type = SPARSE_MATRIX_TYPE_GENERAL;
    status = mkl_sparse_s_mv(opt, alpha, A, md, x, beta, r);
    if (status != SPARSE_STATUS_SUCCESS)
      XABORTM("MKL Sparse Error occurred in execution!\n");
    mkl_sparse_destroy(A);

  #else

    char matdescra[6];
    matdescra[0] = 'G';
    matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
    matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
    matdescra[3] = 'C';

    char trans = (transposed ? 'T' : 'N');

    FEAT_DISABLE_WARNINGS
    mkl_scsrmv(&trans, &mrows, &mcolumns, (const float*)&alpha, matdescra, (const float*) val, (const MKL_INT*)col_idx, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const float*)x, (const float*)&beta, r);
    FEAT_RESTORE_WARNINGS

  #endif
  }

  void MatVecMultCSRDense::exec_mkl_impl(double * r, const double alpha, const double * const x, const double * const y, const double * const val, const Index * const col_idx, const Index * const row_ptr, const Index rows, const Index cols, const Index, const bool transposed)
  {
    MKL_INT mrows = (MKL_INT)rows;
    MKL_INT mcolumns = (MKL_INT)cols;

    if(y == nullptr)
      Memory::memset_main(r, 0, sizeof(double) * (transposed ? cols : rows));
    else if(r != y)
    {
      MKL_INT one = 1;
      MKL_INT n = transposed ? mcolumns : mrows;
      dcopy(&n, y, &one, r, &one);
    }

    double beta = 1.0;

  #ifndef FEAT_USE_MKL_LEGACY_SPMV

    sparse_operation_t opt = (transposed ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE);

    sparse_matrix_t A;
    FEAT_DISABLE_WARNINGS
    sparse_status_t status = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, mrows, mcolumns, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_idx, (double*) val);
    FEAT_RESTORE_WARNINGS
    if (status != SPARSE_STATUS_SUCCESS)
      XABORTM("MKL Sparse Error occurred in execution!\n");
    matrix_descr md;
    md.type = SPARSE_MATRIX_TYPE_GENERAL;
    status = mkl_sparse_d_mv(opt, alpha, A, md, x, beta, r);
    if (status != SPARSE_STATUS_SUCCESS)
      XABORTM("MKL Sparse Error occurred in execution!\n");
    mkl_sparse_destroy(A);

  #else

    char matdescra[6];
    matdescra[0] = 'G';
    matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
    matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
    matdescra[3] = 'C';

    char trans = (transposed ? 'T' : 'N');

    FEAT_DISABLE_WARNINGS
    mkl_dcsrmv(&trans, &mrows, &mcolumns, (const double*)&alpha, matdescra, (const double*) val, (const MKL_INT*)col_idx, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const double*)x, (const double*)&beta, r);
    FEAT_RESTORE_WARNINGS

  #endif
  }
} // namespace FEAT::LAFEM::Arch
