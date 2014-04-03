// includes, FEAST
#include <kernel/lafem/arch/defect.hpp>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

void Defect<Mem::Main, Algo::MKL>::csr(float * r, const float * const rhs, const float * const val, const Index * const col_ind, const Index * const row_ptr, const float * const x, const Index rows, const Index, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, r);
  vsSub(mrows, rhs, r, r);
}

void Defect<Mem::Main, Algo::MKL>::csr(double * r, const double * const rhs, const double * const val, const Index * const col_ind, const Index * const row_ptr, const double * const x, const Index rows, const Index, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, r);
  vdSub(mrows, rhs, r, r);
}

void Defect<Mem::Main, Algo::MKL>::coo(float * r, const float * const rhs, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const float * const x, const Index rows, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float *)x, r);
  vsSub(mrows, rhs, r, r);
}

void Defect<Mem::Main, Algo::MKL>::coo(double * r, const double * const rhs, const double * const val, const Index * const row_ptr, const Index * const col_ptr, const double * const x, const Index rows, const Index used_elements)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';
  MKL_INT ue = (MKL_INT)used_elements;
  mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double *)x, r);
  vdSub(mrows, rhs, r, r);
}
