// includes, FEAST
#include <kernel/lafem/scale.hpp>
#include <kernel/lafem/algorithm.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void Scale<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const float s)
{
  if (r.elements() == x.elements())
  {
    cblas_sscal(r.size(), s, r.elements(), 1);
  }
  else
  {
    copy(r, x);
    cblas_sscal(r.size(), s, r.elements(), 1);
  }
}

void Scale<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const double s)
{
  if (r.elements() == x.elements())
  {
    cblas_dscal(r.size(), s, r.elements(), 1);
  }
  else
  {
    copy(r, x);
    cblas_dscal(r.size(), s, r.elements(), 1);
  }
}

void Scale<Algo::MKL>::value(SparseMatrixCOO<Mem::Main, float> & r, const SparseMatrixCOO<Mem::Main, float> & x, const float s)
{
  if (r.val() == x.val())
  {
    cblas_sscal(r.used_elements(), s, r.val(), 1);
  }
  else
  {
    MemoryPool<Mem::Main>::copy(r.val(), x.val(), r.used_elements() * sizeof(float));
    cblas_sscal(r.used_elements(), s, r.val(), 1);
  }
}

void Scale<Algo::MKL>::value(SparseMatrixCOO<Mem::Main, double> & r, const SparseMatrixCOO<Mem::Main, double> & x, const double s)
{
  if (r.val() == x.val())
  {
    cblas_dscal(r.used_elements(), s, r.val(), 1);
  }
  else
  {
    MemoryPool<Mem::Main>::copy(r.val(), x.val(), r.used_elements() * sizeof(double));
    cblas_dscal(r.used_elements(), s, r.val(), 1);
  }
}

void Scale<Algo::MKL>::value(SparseMatrixCSR<Mem::Main, float> & r, const SparseMatrixCSR<Mem::Main, float> & x, const float s)
{
  if (r.val() == x.val())
  {
    cblas_sscal(r.used_elements(), s, r.val(), 1);
  }
  else
  {
    MemoryPool<Mem::Main>::copy(r.val(), x.val(), r.used_elements() * sizeof(float));
    cblas_sscal(r.used_elements(), s, r.val(), 1);
  }
}

void Scale<Algo::MKL>::value(SparseMatrixCSR<Mem::Main, double> & r, const SparseMatrixCSR<Mem::Main, double> & x, const double s)
{
  if (r.val() == x.val())
  {
    cblas_dscal(r.used_elements(), s, r.val(), 1);
  }
  else
  {
    MemoryPool<Mem::Main>::copy(r.val(), x.val(), r.used_elements() * sizeof(double));
    cblas_dscal(r.used_elements(), s, r.val(), 1);
  }
}

void Scale<Algo::MKL>::value(SparseMatrixELL<Mem::Main, float> & r, const SparseMatrixELL<Mem::Main, float> & x, const float s)
{
  if (r.Ax() == x.Ax())
  {
    cblas_sscal(r.stride() * r.num_cols_per_row(), s, r.Ax(), 1);
  }
  else
  {
    MemoryPool<Mem::Main>::copy(r.Ax(), x.Ax(), r.stride() * r.num_cols_per_row() * sizeof(float));
    cblas_sscal(r.stride() * r.num_cols_per_row(), s, r.Ax(), 1);
  }
}

void Scale<Algo::MKL>::value(SparseMatrixELL<Mem::Main, double> & r, const SparseMatrixELL<Mem::Main, double> & x, const double s)
{
  if (r.Ax() == x.Ax())
  {
    cblas_dscal(r.stride() * r.num_cols_per_row(), s, r.Ax(), 1);
  }
  else
  {
    MemoryPool<Mem::Main>::copy(r.Ax(), x.Ax(), r.stride() * r.num_cols_per_row() * sizeof(double));
    cblas_dscal(r.stride() * r.num_cols_per_row(), s, r.Ax(), 1);
  }
}
