// includes, FEAST
#include <kernel/lafem/product_matvec.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void ProductMatVec<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const SparseMatrixCSR<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b)
{
  MKL_INT rows = a.rows();
  char trans = 'N';
  mkl_cspblas_scsrgemv(&trans, &rows, (float *)a.val(), (MKL_INT*)a.row_ptr(), (MKL_INT*)a.col_ind(), (float *)b.elements(), r.elements());
}

void ProductMatVec<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const SparseMatrixCSR<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b)
{
  MKL_INT rows = a.rows();
  char trans = 'N';
  mkl_cspblas_dcsrgemv(&trans, &rows, (double *)a.val(), (MKL_INT*)a.row_ptr(), (MKL_INT*)a.col_ind(), (double *)b.elements(), r.elements());
}
