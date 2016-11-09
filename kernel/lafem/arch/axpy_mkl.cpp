// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/axpy.hpp>

#include <cstring>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Axpy<Mem::Main>::dv_mkl(float * r, const float a, const float * const x, const float * const y, const Index size)
{
  if (r == y)
  {
    cblas_saxpy((MKL_INT)size, a, x, 1, r, 1);
  }
  else if (r == x)
  {
    float * t = new float[size];
    memcpy(t, x, sizeof(float) * size);
    memcpy(r, y, sizeof(float) * size);
    cblas_saxpy((MKL_INT)size, a, t, 1, r, 1);
    delete[] t;
  }
  else
  {
    memcpy(r, y, sizeof(float) * size);
    cblas_saxpy((MKL_INT)size, a, x, 1, r, 1);
  }
}

void Axpy<Mem::Main>::dv_mkl(double * r, const double a, const double * const x, const double * const y, const Index size)
{
  if (r == y)
  {
    cblas_daxpy((MKL_INT)size, a, x, 1, r, 1);
  }
  else if (r == x)
  {
    double * t = new double[size];
    memcpy(t, x, sizeof(double) * size);
    memcpy(r, y, sizeof(double) * size);
    cblas_daxpy((MKL_INT)size, a, t, 1, r, 1);
    delete[] t;
  }
  else
  {
    memcpy(r, y, sizeof(double) * size);
    cblas_daxpy((MKL_INT)size, a, x, 1, r, 1);
  }
}

void Axpy<Mem::Main>::csr_mkl(float * r, const float a, const float * const x, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  float beta = 1.f;
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
  if (r == y)
  {
    mkl_scsrmv(&trans, &mrows, &mcolumns, (float*)&a, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    scopy(&mrows, (float*)y, &one, r, &one);
    mkl_scsrmv(&trans, &mrows, &mcolumns, (float*)&a, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, &beta, r);
  }
}

void Axpy<Mem::Main>::csr_mkl(double * r, const double a, const double * const x, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  double beta = 1;
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
  if (r == y)
  {
    mkl_dcsrmv(&trans, &mrows, &mcolumns, (double*)&a, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    dcopy(&mrows, (double*)y, &one, r, &one);
    mkl_dcsrmv(&trans, &mrows, &mcolumns, (double*)&a, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, &beta, r);
  }
}

void Axpy<Mem::Main>::csrb_mkl(float * r, const float a, const float * const x, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mrows * mblocksize;
  float beta = 1.f;
  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == y)
  {
    mkl_sbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (float*)&a, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    scopy(&mcopysize, (float*)y, &one, r, &one);
    mkl_sbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (float*)&a, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, &beta, r);
  }
}

void Axpy<Mem::Main>::csrb_mkl(double * r, const double a, const double * const x, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mrows * mblocksize;
  double beta = 1.;
  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == y)
  {
    mkl_dbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (double*)&a, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    dcopy(&mcopysize, (double*)y, &one, r, &one);
    mkl_dbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (double*)&a, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, &beta, r);
  }
}

void Axpy<Mem::Main>::coo_mkl(float * r, const float a, const float * const x, const float * const y, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows,
  const Index columns, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  float beta = 1.f;
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == y)
  {
    mkl_scoomv(&trans, &mrows, &mcolumns, (float*)&a, matdescra, (float*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    scopy(&mrows, (float*)y, &one, r, &one);
    mkl_scoomv(&trans, &mrows, &mcolumns, (float*)&a, matdescra, (float*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float*)x, &beta, r);
  }
}

void Axpy<Mem::Main>::coo_mkl(double * r, const double a, const double * const x, const double * const y, const double * const val, const Index * const row_ptr, const Index * const col_ptr,
  const Index rows, const Index columns, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  double beta = 1.;
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == y)
  {
    mkl_dcoomv(&trans, &mrows, &mcolumns, (double*)&a, matdescra, (double*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    dcopy(&mrows, (double*)y, &one, r, &one);
    mkl_dcoomv(&trans, &mrows, &mcolumns, (double*)&a, matdescra, (double*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double*)x, &beta, r);
  }
}

void Axpy<Mem::Main>::dense_mkl(float * r, const float alpha, const float * const y, const float * const val, const float * const x, const Index rows, const Index columns)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  float beta = 1.;

  if (r == y)
  {
    cblas_sgemv(CblasRowMajor, CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
  }
  else
  {
    MKL_INT one = 1;
    scopy(&mrows, (float*)y, &one, r, &one);
    cblas_sgemv(CblasRowMajor, CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
  }
}

void Axpy<Mem::Main>::dense_mkl(double * r, const double alpha, const double * const y, const double * const val, const double * const x, const Index rows, const Index columns)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  double beta = 1.;

  if (r == y)
  {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
  }
  else
  {
    MKL_INT one = 1;
    dcopy(&mrows, (double*)y, &one, r, &one);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
  }
}
