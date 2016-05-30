// includes, FEAT
#include <kernel/lafem/arch/product_matvec.hpp>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void ProductMatVec<Mem::Main>::csr_mkl(float * r, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const float * const x, const Index rows, const Index, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, r);
}

void ProductMatVec<Mem::Main>::csr_mkl(double * r, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const double * const x, const Index rows, const Index, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, r);
}

void ProductMatVec<Mem::Main>::csrb_mkl(float * r, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const float * const x, const Index rows, const Index, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  char trans = 'N';
  mkl_cspblas_sbsrgemv(&trans, &mrows, &mblocksize, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, r);
}

void ProductMatVec<Mem::Main>::csrb_mkl(double * r, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const double * const x, const Index rows, const Index, const Index, const int blocksize)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mblocksize = (MKL_INT)blocksize;
  char trans = 'N';
  mkl_cspblas_dbsrgemv(&trans, &mrows, &mblocksize, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, r);
}

void ProductMatVec<Mem::Main>::coo_mkl(float * r, const float * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const float * const x, const Index rows, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float *)x, r);
}

void ProductMatVec<Mem::Main>::coo_mkl(double * r, const double * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const double * const x, const Index rows, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double *)x, r);
}
