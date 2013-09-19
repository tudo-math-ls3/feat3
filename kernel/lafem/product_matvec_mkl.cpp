// includes, FEAST
#include <kernel/lafem/product_matvec.hpp>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void ProductMatVec<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const SparseMatrixCSR<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b)
{
  MKL_INT rows = (MKL_INT)a.rows();
  char trans = 'N';
  mkl_cspblas_scsrgemv(&trans, &rows, (float *)a.val(), (MKL_INT*)a.row_ptr(), (MKL_INT*)a.col_ind(), (float *)b.elements(), r.elements());
}

void ProductMatVec<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const SparseMatrixCSR<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b)
{
  MKL_INT rows = (MKL_INT)a.rows();
  char trans = 'N';
  mkl_cspblas_dcsrgemv(&trans, &rows, (double *)a.val(), (MKL_INT*)a.row_ptr(), (MKL_INT*)a.col_ind(), (double *)b.elements(), r.elements());
}

void ProductMatVec<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const SparseMatrixCOO<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b)
{
  MKL_INT rows = (MKL_INT)a.rows();
  char trans = 'N';
  MKL_INT ue = (MKL_INT)a.used_elements();
  mkl_cspblas_scoogemv(&trans, &rows, (float *)a.val(), (MKL_INT*)a.row(), (MKL_INT*)a.column(), &ue, (float *)b.elements(), r.elements());
}

void ProductMatVec<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const SparseMatrixCOO<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b)
{
  MKL_INT rows = (MKL_INT)a.rows();
  char trans = 'N';
  MKL_INT ue = (MKL_INT)a.used_elements();
  mkl_cspblas_dcoogemv(&trans, &rows, (double *)a.val(), (MKL_INT*)a.row(), (MKL_INT*)a.column(), &ue, (double *)b.elements(), r.elements());
}
