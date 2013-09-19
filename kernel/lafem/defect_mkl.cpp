// includes, FEAST
#include <kernel/lafem/defect.hpp>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void Defect<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & rhs, const SparseMatrixCSR<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b)
{
  MKL_INT rows = (MKL_INT)a.rows();
  char trans = 'N';
  mkl_cspblas_scsrgemv(&trans, &rows, (float *)a.val(), (MKL_INT*)a.row_ptr(), (MKL_INT*)a.col_ind(), (float *)b.elements(), r.elements());
  vsSub((MKL_INT)r.size(), rhs.elements(), r.elements(), r.elements());
}

void Defect<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & rhs, const SparseMatrixCSR<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b)
{
  MKL_INT rows = (MKL_INT)a.rows();
  char trans = 'N';
  mkl_cspblas_dcsrgemv(&trans, &rows, (double *)a.val(), (MKL_INT*)a.row_ptr(), (MKL_INT*)a.col_ind(), (double *)b.elements(), r.elements());
  vdSub((MKL_INT)r.size(), rhs.elements(), r.elements(), r.elements());
}

void Defect<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & rhs, const SparseMatrixCOO<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b)
{
  MKL_INT rows = (MKL_INT)a.rows();
  char trans = 'N';
  MKL_INT ue = (MKL_INT)a.used_elements();
  mkl_cspblas_scoogemv(&trans, &rows, (float *)a.val(), (MKL_INT*)a.row(), (MKL_INT*)a.column(), &ue, (float *)b.elements(), r.elements());
  vsSub((MKL_INT)r.size(), rhs.elements(), r.elements(), r.elements());
}

void Defect<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & rhs, const SparseMatrixCOO<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b)
{
  MKL_INT rows = (MKL_INT)a.rows();
  char trans = 'N';
  MKL_INT ue = (MKL_INT)a.used_elements();
  mkl_cspblas_dcoogemv(&trans, &rows, (double *)a.val(), (MKL_INT*)a.row(), (MKL_INT*)a.column(), &ue, (double *)b.elements(), r.elements());
  vdSub((MKL_INT)r.size(), rhs.elements(), r.elements(), r.elements());
}
