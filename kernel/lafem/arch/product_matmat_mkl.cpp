// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/product_matmat.hpp>

#include <mkl_blas.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void ProductMatMat<Mem::Main>::dense_mkl(float * r, const float * const x, const float * const y, const Index rows, const Index columns, const Index inner)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT minner = (MKL_INT)inner;
  float one = 1.f;
  float zero = 0.f;
  char trans = 'N';
  sgemm(&trans, &trans, &mcolumns, &mrows, &minner, &one, y, &mcolumns, x, &minner, &zero, r, &mcolumns);
}

void ProductMatMat<Mem::Main>::dense_mkl(double * r, const double * const x, const double * const y, const Index rows, const Index columns, const Index inner)
{
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT minner = (MKL_INT)inner;
  double one = 1.;
  double zero = 0.;
  char trans = 'N';
  dgemm(&trans, &trans, &mcolumns, &mrows, &minner, &one, y, &mcolumns, x, &minner, &zero, r, &mcolumns);
}
