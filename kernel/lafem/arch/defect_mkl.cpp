// includes, FEAT
#include <kernel/lafem/arch/defect.hpp>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Defect<Mem::Main>::csr_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const float * const x, const Index rows, const Index columns, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  float alpha = -1.f;
  float beta = 1.f;
  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == rhs)
  {
    mkl_scsrmv(&trans, &mrows, &mcolumns, &alpha, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    scopy(&mrows, (float*)rhs, &one, r, &one);
    mkl_scsrmv(&trans, &mrows, &mcolumns, &alpha, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, &beta, r);
  }
}

void Defect<Mem::Main>::csr_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const double * const x, const Index rows, const Index columns, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  double alpha = -1.;
  double beta = 1.;
  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == rhs)
  {
    mkl_dcsrmv(&trans, &mrows, &mcolumns, &alpha, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    dcopy(&mrows, (double*)rhs, &one, r, &one);
    mkl_dcsrmv(&trans, &mrows, &mcolumns, &alpha, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, &beta, r);
  }
}

void Defect<Mem::Main>::csrb_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const float * const x, const Index rows, const Index columns, const Index, const int blocksize)
{

  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mrows * mblocksize;
  float alpha = -1.f;
  float beta = 1.f;
  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == rhs)
  {
    mkl_sbsrmv(&trans, &mrows, &mcolumns, &mblocksize, &alpha, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    scopy(&mcopysize, (float*)rhs, &one, r, &one);
    mkl_sbsrmv(&trans, &mrows, &mcolumns, &mblocksize, &alpha, matdescra, (float*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (float*)x, &beta, r);
  }
}

void Defect<Mem::Main>::csrb_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const double * const x, const Index rows, const Index columns, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  MKL_INT mcopysize = mrows * mblocksize;
  double alpha = -1.;
  double beta = 1.;
  char trans = 'N';
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == rhs)
  {
    mkl_dbsrmv(&trans, &mrows, &mcolumns, &mblocksize, &alpha, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    dcopy(&mcopysize, (double*)rhs, &one, r, &one);
    mkl_dbsrmv(&trans, &mrows, &mcolumns, &mblocksize, &alpha, matdescra, (double*) val, (MKL_INT*)col_ind, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr) + 1, (double*)x, &beta, r);
  }
}

void Defect<Mem::Main>::coo_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const float * const x, const Index rows, const Index columns, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  float alpha = -1.f;
  float beta = 1.f;
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == rhs)
  {
    mkl_scoomv(&trans, &mrows, &mcolumns, &alpha, matdescra, (float*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    scopy(&mrows, (float*)rhs, &one, r, &one);
    mkl_scoomv(&trans, &mrows, &mcolumns, &alpha, matdescra, (float*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float*)x, &beta, r);
  }
}

void Defect<Mem::Main>::coo_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const double * const x, const Index rows, const Index columns, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  double alpha = -1.;
  double beta = 1.;
  char matdescra[6];
  matdescra[0] = 'G';
  matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
  matdescra[3] = 'C';
  if (r == rhs)
  {
    mkl_dcoomv(&trans, &mrows, &mcolumns, &alpha, matdescra, (double*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double*)x, &beta, r);
  }
  else
  {
    MKL_INT one = 1;
    dcopy(&mrows, (double*)rhs, &one, r, &one);
    mkl_dcoomv(&trans, &mrows, &mcolumns, &alpha, matdescra, (double*)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double*)x, &beta, r);
  }
}
