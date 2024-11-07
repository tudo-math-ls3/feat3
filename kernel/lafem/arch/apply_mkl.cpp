// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/apply.hpp>
#include <kernel/util/exception.hpp>

#include <cstring>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
#include <mkl_spblas.h>
FEAT_RESTORE_WARNINGS

/// \todo remove this pragma once MKL is switched to sparse executor interface
// disable 'deprecated' warnings for MKL functions, because this file won't compile otherwise
#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(disable: 4996)
#endif

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Apply::csr_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;

  if (r != y)
  {
    MKL_INT one = 1;
    if (transposed)
      scopy(&mcolumns, (const float*)y, &one, r, &one);
    else
      scopy(&mrows, (const float*)y, &one, r, &one);
  }

#ifdef FEAT_USE_MKL_SPARSE_EXECUTOR

  sparse_operation_t opt;
  if (transposed)
    opt = SPARSE_OPERATION_TRANSPOSE;
  else
    opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS
  sparse_status_t status = mkl_sparse_s_create_csr(&A, SPARSE_INDEX_BASE_ZERO, mrows, mcolumns, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (float*) val);
  FEAT_RESTORE_WARNINGS
  if (status != SPARSE_STATUS_SUCCESS)
    XABORTM("MKL Sparse Error occurred in execution!\n");
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  status = mkl_sparse_s_mv(opt, a, A, md, x, b, r);
  if (status != SPARSE_STATUS_SUCCESS)
    XABORTM("MKL Sparse Error occurred in execution!\n");
  mkl_sparse_destroy(A);

#else

  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';

  char trans;
  if (transposed)
    trans = 'T';
  else
    trans = 'N';

  FEAT_DISABLE_WARNINGS
  mkl_scsrmv(&trans, &mrows, &mcolumns, (const float*)&a, matdescra, (const float*) val, (const MKL_INT*)col_ind, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const float*)x, (const float*)&b, r);
  FEAT_RESTORE_WARNINGS

#endif
}

void Apply::csr_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;

  if (r != y)
  {
    MKL_INT one = 1;
    if (transposed)
      dcopy(&mcolumns, (const double*)y, &one, r, &one);
    else
      dcopy(&mrows, (const double*)y, &one, r, &one);
  }

#ifdef FEAT_USE_MKL_SPARSE_EXECUTOR

  sparse_operation_t opt;
  if (transposed)
    opt = SPARSE_OPERATION_TRANSPOSE;
  else
    opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS
  sparse_status_t status = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, mrows, mcolumns, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (double*) val);
  FEAT_RESTORE_WARNINGS
  if (status != SPARSE_STATUS_SUCCESS)
    XABORTM("MKL Sparse Error occurred in execution!\n");
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  status = mkl_sparse_d_mv(opt, a, A, md, x, b, r);
  if (status != SPARSE_STATUS_SUCCESS)
    XABORTM("MKL Sparse Error occurred in execution!\n");
  mkl_sparse_destroy(A);

#else

  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';

  char trans;
  if (transposed)
    trans = 'T';
  else
    trans = 'N';

  FEAT_DISABLE_WARNINGS
  mkl_dcsrmv(&trans, &mrows, &mcolumns, (const double*)&a, matdescra, (const double*) val, (const MKL_INT*)col_ind, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const double*)x, (const double*)&b, r);
  FEAT_RESTORE_WARNINGS

#endif
}

void Apply::bcsr_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mrows * mblocksize;

  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mcopysize, (const float*)y, &one, r, &one);
  }


#ifdef FEAT_USE_MKL_SPARSE_EXECUTOR

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS
  mkl_sparse_s_create_bsr(&A, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, mrows, mcolumns, mblocksize, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (float*) val);
  FEAT_RESTORE_WARNINGS
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_s_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);

#else

  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';

  FEAT_DISABLE_WARNINGS
  mkl_sbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (const float*)&a, matdescra, (const float*) val, (const MKL_INT*)col_ind, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const float*)x, (const float*)&b, r);
  FEAT_RESTORE_WARNINGS

#endif
}

void Apply::bcsr_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mrows * mblocksize;

  if (r != y)
  {
    MKL_INT one = 1;
    dcopy(&mcopysize, (const double*)y, &one, r, &one);
  }


#ifdef FEAT_USE_MKL_SPARSE_EXECUTOR

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS
  mkl_sparse_d_create_bsr(&A, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, mrows, mcolumns, mblocksize, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (double*) val);
  FEAT_RESTORE_WARNINGS
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_d_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);

#else

  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';

  FEAT_DISABLE_WARNINGS
  mkl_dbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (const double*)&a, matdescra, (const double*) val, (const MKL_INT*)col_ind, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const double*)x, (const double*)&b, r);
  FEAT_RESTORE_WARNINGS

#endif
}

void Apply::bcsr_transposed_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mcolumns * mblocksize;

  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mcopysize, (const float*)y, &one, r, &one);
  }


#ifdef FEAT_USE_MKL_SPARSE_EXECUTOR

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS
    mkl_sparse_s_create_bsr(&A, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, mrows, mcolumns, mblocksize, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (float*) val);
  FEAT_RESTORE_WARNINGS
    matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_s_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);

#else

  char trans = 'T';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';

  FEAT_DISABLE_WARNINGS
    mkl_sbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (const float*)&a, matdescra, (const float*) val, (const MKL_INT*)col_ind, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const float*)x, (const float*)&b, r);
  FEAT_RESTORE_WARNINGS

#endif
}

void Apply::bcsr_transposed_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mcolumns * mblocksize;

  if (r != y)
  {
    MKL_INT one = 1;
    dcopy(&mcopysize, (const double*)y, &one, r, &one);
  }


#ifdef FEAT_USE_MKL_SPARSE_EXECUTOR

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS
    mkl_sparse_d_create_bsr(&A, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, mrows, mcolumns, mblocksize, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (double*) val);
  FEAT_RESTORE_WARNINGS
    matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_d_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);

#else

  char trans = 'T';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';

  FEAT_DISABLE_WARNINGS
    mkl_dbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (const double*)&a, matdescra, (const double*) val, (const MKL_INT*)col_ind, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const double*)x, (const double*)&b, r);
  FEAT_RESTORE_WARNINGS

#endif
}

void Apply::dense_mkl(float * r, const float alpha, const float beta, const float * const y, const float * const val, const float * const x, const Index rows, const Index columns)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;

  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mrows, (const float*)y, &one, r, &one);
  }
  cblas_sgemv(CblasRowMajor, CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
}

void Apply::dense_mkl(double * r, const double alpha, const double beta, const double * const y, const double * const val, const double * const x, const Index rows, const Index columns)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;

  if (r != y)
  {
    MKL_INT one = 1;
    dcopy(&mrows, (const double*)y, &one, r, &one);
  }
  cblas_dgemv(CblasRowMajor, CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
}
