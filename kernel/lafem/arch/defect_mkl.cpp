// includes, FEAST
#include <kernel/lafem/arch/defect.hpp>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

void Defect<Mem::Main>::csr_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const float * const x, const Index rows, const Index /*columns*/, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  float * t = new float[rows];
  mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, t);
  vsSub(mrows, rhs, t, r);
  delete[] t;
}

void Defect<Mem::Main>::csr_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const double * const x, const Index rows, const Index /*columns*/, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  double * t = new double[rows];
  mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, t);
  vdSub(mrows, rhs, t, r);
  delete[] t;
}

void Defect<Mem::Main>::csrb_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const float * const x, const Index rows, const Index /*columns*/, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  char trans = 'N';
  float * t = new float[mrows * mblocksize];
  mkl_cspblas_sbsrgemv(&trans, &mrows, &mblocksize, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, t);
  vsSub(mrows * mblocksize, rhs, t, r);
  delete[] t;
}

void Defect<Mem::Main>::csrb_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const double * const x, const Index rows, const Index /*columns*/, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  char trans = 'N';
  double * t = new double[mrows * mblocksize];
  mkl_cspblas_dbsrgemv(&trans, &mrows, &mblocksize, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, t);
  vdSub(mrows * mblocksize, rhs, t, r);
  delete[] t;
}

void Defect<Mem::Main>::coo_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const float * const x, const Index rows, const Index /*columns*/, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  float * t = new float[rows];
  mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float *)x, t);
  vsSub(mrows, rhs, t, r);
  delete[] t;
}

void Defect<Mem::Main>::coo_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const double * const x, const Index rows, const Index /*columns*/, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  double * t = new double[rows];
  mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double *)x, t);
  vdSub(mrows, rhs, t, r);
  delete[] t;
}
