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
  char trans;
  if (transposed)
    trans = 'T';
  else
    trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mrows, (float*)y, &one, r, &one);
  }
  mkl_scsrmv(&trans, &mrows, &mcolumns, (float*)&a, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, (float*)&b, r);
}

void Apply<Mem::Main>::csr_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  char trans;
  if (transposed)
    trans = 'T';
  else
    trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r != y)
  {
    MKL_INT one = 1;
    dcopy(&mrows, (double*)y, &one, r, &one);
  }
  mkl_dcsrmv(&trans, &mrows, &mcolumns, (double*)&a, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, (double*)&b, r);
}

void Apply<Mem::Main>::csrb_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mrows * mblocksize;
  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mcopysize, (float*)y, &one, r, &one);
  }
  mkl_sbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (float*)&a, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, (float*)&b, r);
}

void Apply<Mem::Main>::csrb_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mrows * mblocksize;
  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r != y)
  {
    MKL_INT one = 1;
    dcopy(&mcopysize, (double*)y, &one, r, &one);
  }
  mkl_dbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (double*)&a, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, (double*)&b, r);
}

void Apply<Mem::Main>::coo_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows,
  const Index columns, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mrows, (float*)y, &one, r, &one);
  }
  mkl_scoomv(&trans, &mrows, &mcolumns, (float*)&a, matdescra, (float*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float*)x, (float*)&b, r);
}

void Apply<Mem::Main>::coo_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const row_ptr, const Index * const col_ptr,
  const Index rows, const Index columns, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r != y)
  {
    MKL_INT one = 1;
    dcopy(&mrows, (double*)y, &one, r, &one);
  }
  mkl_dcoomv(&trans, &mrows, &mcolumns, (double*)&a, matdescra, (double*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double*)x, (double*)&b, r);
}

void Apply<Mem::Main>::dense_mkl(float * r, const float alpha, const float beta, const float * const y, const float * const val, const float * const x, const Index rows, const Index columns)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;

  if (r != y)
  {
    MKL_INT one = 1;
    scopy(&mrows, (float*)y, &one, r, &one);
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
    dcopy(&mrows, (double*)y, &one, r, &one);
  }
  cblas_dgemv(CblasRowMajor, CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
}