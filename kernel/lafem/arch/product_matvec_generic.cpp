// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>

#include <cstring>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_, typename IT_>
void ProductMatVec<Mem::Main, Algo::Generic>::csr(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    DT_ sum(0);
    const IT_ end(row_ptr[row + 1]);
    for (IT_ i(row_ptr[row]) ; i < end ; ++i)
    {
      sum += val[i] * x[col_ind[i]];
    }
    r[row] = sum;
  }
}

template void ProductMatVec<Mem::Main, Algo::Generic>::csr(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::csr(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void ProductMatVec<Mem::Main, Algo::Generic>::ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    const IT_ * tAj(Aj);
    const DT_ * tAx(Ax);
    DT_ sum(0);
    tAj += row;
    tAx += row;

    const IT_ max(Arl[row]);
    for(IT_ n(0); n < max ; n++)
    {
      const DT_ A_ij = *tAx;

      const IT_ col = *tAj;
      sum += A_ij * x[col];

      tAj += stride;
      tAx += stride;
    }
    r[row] = sum;
  }
}

template void ProductMatVec<Mem::Main, Algo::Generic>::ell(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::ell(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

template <typename DT_>
void ProductMatVec<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ * const val, const Index * const row_ptr, const Index * const col_ptr, const DT_ * const x, const Index rows, const Index used_elements)
{
  Index iter(0);
  for (Index row(0); row < rows; ++row)
  {
    DT_ sum(DT_(0));
    while (iter < used_elements && row_ptr[iter] == row)
    {
      sum += val[iter] * x[col_ptr[iter]];
      ++iter;
    }
    r[row] = sum;
  }
}

template void ProductMatVec<Mem::Main, Algo::Generic>::coo(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::coo(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
