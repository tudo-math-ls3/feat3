// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/apply.hpp>

#include <cstring>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Apply<Mem::Main>::csr_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;

  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mrows, (const float*)y, &one, r, &one);
  }

  sparse_operation_t opt;
  if (transposed)
    opt = SPARSE_OPERATION_TRANSPOSE;
  else
    opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS;
  mkl_sparse_s_create_csr(&A, SPARSE_INDEX_BASE_ZERO, mrows, mcolumns, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (float*) val);
  FEAT_RESTORE_WARNINGS;
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_s_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);
}

void Apply<Mem::Main>::csr_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;

  if (r != y)
  {
    MKL_INT one = 1;
    dcopy(&mrows, (const double*)y, &one, r, &one);
  }

  sparse_operation_t opt;
  if (transposed)
    opt = SPARSE_OPERATION_TRANSPOSE;
  else
    opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS;
  mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, mrows, mcolumns, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (double*) val);
  FEAT_RESTORE_WARNINGS;
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_d_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);
}

void Apply<Mem::Main>::csrb_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
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

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS;
  mkl_sparse_s_create_bsr(&A, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, mrows, mcolumns, mblocksize, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (float*) val);
  FEAT_RESTORE_WARNINGS;
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_s_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);
}

void Apply<Mem::Main>::csrb_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
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

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS;
  mkl_sparse_d_create_bsr(&A, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, mrows, mcolumns, mblocksize, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_ind, (double*) val);
  FEAT_RESTORE_WARNINGS;
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_d_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);
}

void Apply<Mem::Main>::coo_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows,
  const Index columns, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT ue = (MKL_INT)used_elements;

  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mrows, (const float*)y, &one, r, &one);
  }

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS;
  mkl_sparse_s_create_coo(&A, SPARSE_INDEX_BASE_ZERO, mrows, mcolumns, ue, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, (float*) val);
  FEAT_RESTORE_WARNINGS;
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_s_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);
}

void Apply<Mem::Main>::coo_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const row_ptr, const Index * const col_ptr,
  const Index rows, const Index columns, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT ue = (MKL_INT)used_elements;

  if (r != y)
  {
    MKL_INT one = 1;
    dcopy(&mrows, (const double*)y, &one, r, &one);
  }

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;

  sparse_matrix_t A;
  FEAT_DISABLE_WARNINGS;
  mkl_sparse_d_create_coo(&A, SPARSE_INDEX_BASE_ZERO, mrows, mcolumns, ue, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, (double*) val);
  FEAT_RESTORE_WARNINGS;
  matrix_descr md;
  md.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_sparse_d_mv(opt, a, A, md, x, b, r);
  mkl_sparse_destroy(A);
}

void Apply<Mem::Main>::dense_mkl(float * r, const float alpha, const float beta, const float * const y, const float * const val, const float * const x, const Index rows, const Index columns)
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

void Apply<Mem::Main>::dense_mkl(double * r, const double alpha, const double beta, const double * const y, const double * const val, const double * const x, const Index rows, const Index columns)
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
