// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>

#include <cstring>

#ifdef FEAST_GMP
#include <gmpxx.h>
#include <mpfr.h>
#endif

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void ProductMatVec<Mem::Main, Algo::Generic>::csr(DT_ * r, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const DT_ * const x, const Index rows)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    DT_ sum(0);
    const Index end(row_ptr[row + 1]);
    for (Index i(row_ptr[row]) ; i < end ; ++i)
    {
      sum += val[i] * x[col_ind[i]];
    }
    r[row] = sum;
  }
}

template void ProductMatVec<Mem::Main, Algo::Generic>::csr(float *, const float * const, const Index * const, const Index * const, const float * const, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::csr(double *, const double * const, const Index * const, const Index * const, const double * const, const Index);
#ifdef FEAST_GMP
template void ProductMatVec<Mem::Main, Algo::Generic>::csr(mpf_class *, const mpf_class * const, const Index * const, const Index * const, const mpf_class * const, const Index);
#endif

template <typename DT_>
void ProductMatVec<Mem::Main, Algo::Generic>::ell(DT_ * r, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    const Index * tAj(Aj);
    const DT_ * tAx(Ax);
    DT_ sum(0);
    tAj += row;
    tAx += row;

    const Index max(Arl[row]);
    for(Index n(0); n < max ; n++)
    {
      const DT_ A_ij = *tAx;

      const Index col = *tAj;
      sum += A_ij * x[col];

      tAj += stride;
      tAx += stride;
    }
    r[row] = sum;
  }
}

template void ProductMatVec<Mem::Main, Algo::Generic>::ell(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::ell(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
#ifdef FEAST_GMP
template void ProductMatVec<Mem::Main, Algo::Generic>::ell(mpf_class *, const mpf_class * const, const Index * const, const Index * const, const mpf_class * const, const Index, const Index);
#endif

template <typename DT_>
void ProductMatVec<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ * const val, const Index * const row_ptr, const Index * const col_ptr, const DT_ * const x, const Index rows, const Index used_elements)
{
  memset(r, 0, sizeof(DT_) * rows);
  for (Index i(0) ; i < used_elements ; ++i)
  {
    r[row_ptr[i]] += val[i] * x[col_ptr[i]];
  }
}

template void ProductMatVec<Mem::Main, Algo::Generic>::coo(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::coo(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
#ifdef FEAST_GMP
template void ProductMatVec<Mem::Main, Algo::Generic>::coo(mpf_class *, const mpf_class * const, const Index * const, const Index * const, const mpf_class * const, const Index, const Index);
#endif
