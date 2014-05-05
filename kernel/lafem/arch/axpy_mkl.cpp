// includes, FEAST
#include <kernel/lafem/arch/axpy.hpp>

#include <cstring>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

void Axpy<Mem::Main, Algo::MKL>::dv(float * r, const float a, const float * const x, const float * const y, const Index size)
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

void Axpy<Mem::Main, Algo::MKL>::dv(double * r, const double a, const double * const x, const double * const y, const Index size)
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

void Axpy<Mem::Main, Algo::MKL>::dv(float * r, const float * const a, const float * const x, const float * const y, const Index size)
{
  if (r == y)
  {
    float * t = new float[size];
    vsMul((MKL_INT)size, x, a, t);
    vsAdd((MKL_INT)size, t, y, r);
    delete[] t;
  }
  else if (r == x)
  {
    vsMul((MKL_INT)size, x, a, r);
    vsAdd((MKL_INT)size, x, y, r);
  }
  else
  {
    vsMul((MKL_INT)size, x, a, r);
    vsAdd((MKL_INT)size, r, y, r);
  }
}

void Axpy<Mem::Main, Algo::MKL>::dv(double * r, const double * const a, const double * const x, const double * const y, const Index size)
{
  if (r == y)
  {
    double * t = new double[size];
    vdMul((MKL_INT)size, x, a, t);
    vdAdd((MKL_INT)size, t, y, r);
    delete[] t;
  }
  else if (r == x)
  {
    vdMul((MKL_INT)size, x, a, r);
    vdAdd((MKL_INT)size, x, y, r);
  }
  else
  {
    vdMul((MKL_INT)size, x, a, r);
    vdAdd((MKL_INT)size, r, y, r);
  }
}

void Axpy<Mem::Main, Algo::MKL>::csr(float * r, const float a, const float * const x, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';

  if (r == y)
  {
    float * t = new float[rows];
    mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, t);
    cblas_saxpy((MKL_INT)mrows, a, t, 1, r, 1);
    delete[] t;
  }
  else if (r == x)
  {
    float * t = new float[rows];
    memcpy(t, x, sizeof(float) * rows);
    mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, t, r);
    memcpy(t, y, sizeof(float) * rows);
    cblas_saxpy(mrows, a, r, 1, t, 1);
    memcpy(r, t, sizeof(float) * rows);
    delete[] t;
  }
  else
  {
    mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, r);
    float * t = new float[rows];
    memcpy(t, y, sizeof(float) * rows);
    cblas_saxpy(mrows, a, r, 1, t, 1);
    memcpy(r, t, sizeof(float) * rows);
    delete[] t;
  }
}

void Axpy<Mem::Main, Algo::MKL>::csr(double * r, const double a, const double * const x, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';

  if (r == y)
  {
    double * t = new double[rows];
    mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, t);
    cblas_daxpy((MKL_INT)mrows, a, t, 1, r, 1);
    delete[] t;
  }
  else if (r == x)
  {
    double * t = new double[rows];
    memcpy(t, x, sizeof(double) * rows);
    mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, t, r);
    memcpy(t, y, sizeof(double) * rows);
    cblas_daxpy(mrows, a, r, 1, t, 1);
    memcpy(r, t, sizeof(double) * rows);
    delete[] t;
  }
  else
  {
    mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, r);
    double * t = new double[rows];
    memcpy(t, y, sizeof(double) * rows);
    cblas_daxpy(mrows, a, r, 1, t, 1);
    memcpy(r, t, sizeof(double) * rows);
    delete[] t;
  }
}

void Axpy<Mem::Main, Algo::MKL>::csr(float * r, const float * const a, const float * const x, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';

  if (r == y)
  {
    float * t = new float[rows];
    mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, t);
    vsMul(mrows, t, a, t);
    vsAdd(mrows, x, y, r);
    delete[] t;
  }
  else if (r == x)
  {
    float * t = new float[rows];
    memcpy(t, x, sizeof(float) * rows);
    mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, t, r);
    vsMul(mrows, r, a, r);
    vsAdd(mrows, r, y, r);
    delete[] t;
  }
  else
  {
    mkl_cspblas_scsrgemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (float *)x, r);
    vsMul(mrows, r, a, r);
    vsAdd(mrows, r, y, r);
  }
}

void Axpy<Mem::Main, Algo::MKL>::csr(double * r, const double * const a, const double * const x, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index, const Index)
{
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';

  if (r == y)
  {
    double * t = new double[rows];
    mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, t);
    vdMul(mrows, t, a, t);
    vdAdd(mrows, x, y, r);
    delete[] t;
  }
  else if (r == x)
  {
    double * t = new double[rows];
    memcpy(t, x, sizeof(double) * rows);
    mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, t, r);
    vdMul(mrows, r, a, r);
    vdAdd(mrows, r, y, r);
    delete[] t;
  }
  else
  {
    mkl_cspblas_dcsrgemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ind, (double *)x, r);
    vdMul(mrows, r, a, r);
    vdAdd(mrows, r, y, r);
  }
}

void Axpy<Mem::Main, Algo::MKL>::coo(float * r, const float a, const float * const x, const float * const y, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements)
{
  MKL_INT ue = (MKL_INT)used_elements;
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';

  if (r == y)
  {
    float * t = new float[rows];
    mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float *)x, t);
    cblas_saxpy((MKL_INT)mrows, a, t, 1, r, 1);
    delete[] t;
  }
  else if (r == x)
  {
    float * t = new float[rows];
    memcpy(t, x, sizeof(float) * rows);
    mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, t, r);
    memcpy(t, y, sizeof(float) * rows);
    cblas_saxpy(mrows, a, r, 1, t, 1);
    memcpy(r, t, sizeof(float) * rows);
    delete[] t;
  }
  else
  {
    mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float *)x, r);
    float * t = new float[rows];
    memcpy(t, y, sizeof(float) * rows);
    cblas_saxpy(mrows, a, r, 1, t, 1);
    memcpy(r, t, sizeof(float) * rows);
    delete[] t;
  }
}

void Axpy<Mem::Main, Algo::MKL>::coo(double * r, const double a, const double * const x, const double * const y, const double * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements)
{
  MKL_INT ue = (MKL_INT)used_elements;
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';

  if (r == y)
  {
    double * t = new double[rows];
    mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double *)x, t);
    cblas_daxpy((MKL_INT)mrows, a, t, 1, r, 1);
    delete[] t;
  }
  else if (r == x)
  {
    double * t = new double[rows];
    memcpy(t, x, sizeof(double) * rows);
    mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, t, r);
    memcpy(t, y, sizeof(double) * rows);
    cblas_daxpy(mrows, a, r, 1, t, 1);
    memcpy(r, t, sizeof(double) * rows);
    delete[] t;
  }
  else
  {
    mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double *)x, r);
    double * t = new double[rows];
    memcpy(t, y, sizeof(double) * rows);
    cblas_daxpy(mrows, a, r, 1, t, 1);
    memcpy(r, t, sizeof(double) * rows);
    delete[] t;
  }
}

void Axpy<Mem::Main, Algo::MKL>::coo(float * r, const float * const a, const float * const x, const float * const y, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements)
{
  MKL_INT ue = (MKL_INT)used_elements;
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';

  if (r == y)
  {
    float * t = new float[rows];
    mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float *)x, t);
    vsMul(mrows, t, a, t);
    vsAdd(mrows, x, y, r);
    delete[] t;
  }
  else if (r == x)
  {
    float * t = new float[rows];
    memcpy(t, x, sizeof(float) * rows);
    mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, t, r);
    vsMul(mrows, r, a, r);
    vsAdd(mrows, r, y, r);
    delete[] t;
  }
  else
  {
    mkl_cspblas_scoogemv(&trans, &mrows, (float *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (float *)x, r);
    vsMul(mrows, r, a, r);
    vsAdd(mrows, r, y, r);
  }
}

void Axpy<Mem::Main, Algo::MKL>::coo(double * r, const double * const a, const double * const x, const double * const y, const double * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements)
{
  MKL_INT ue = (MKL_INT)used_elements;
  MKL_INT mrows = (MKL_INT)rows;
  char trans = 'N';

  if (r == y)
  {
    double * t = new double[rows];
    mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double *)x, t);
    vdMul(mrows, t, a, t);
    vdAdd(mrows, x, y, r);
    delete[] t;
  }
  else if (r == x)
  {
    double * t = new double[rows];
    memcpy(t, x, sizeof(double) * rows);
    mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, t, r);
    vdMul(mrows, r, a, r);
    vdAdd(mrows, r, y, r);
    delete[] t;
  }
  else
  {
    mkl_cspblas_dcoogemv(&trans, &mrows, (double *)val, (MKL_INT*)row_ptr, (MKL_INT*)col_ptr, &ue, (double *)x, r);
    vdMul(mrows, r, a, r);
    vdAdd(mrows, r, y, r);
  }
}
