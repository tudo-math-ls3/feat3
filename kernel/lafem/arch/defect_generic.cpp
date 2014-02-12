// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/difference.hpp>

#include <cstring>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Defect<Mem::Main, Algo::Generic>::csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const DT_ * const x, const Index rows)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    DT_ sum(0);
    const Index end(row_ptr[row + 1]);
    for (Index i(row_ptr[row]) ; i < end ; ++i)
    {
      sum += val[i] * x[col_ind[i]];
    }
    r[row] = rhs[row] - sum;
  }
}

template void Defect<Mem::Main, Algo::Generic>::csr(float *, const float * const, const float * const, const Index * const, const Index * const, const float * const, const Index);
template void Defect<Mem::Main, Algo::Generic>::csr(double *, const double * const, const double * const, const Index * const, const Index * const, const double * const, const Index);

template <typename DT_>
void Defect<Mem::Main, Algo::Generic>::ell(DT_ * r, const DT_ * const rhs, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const DT_ * const x, const Index stride, const Index rows)
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
    r[row] = rhs[row] - sum;
  }
}

template void Defect<Mem::Main, Algo::Generic>::ell(float *, const float * const, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void Defect<Mem::Main, Algo::Generic>::ell(double *,const double * const,  const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);

template <typename DT_>
void Defect<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ * const rhs, const DT_ * const val, const Index * const row_ptr, const Index * const col_ptr, const DT_ * const x, const Index rows, const Index used_elements)
{
  Index iter(0);
  for (Index row(0); row < rows; ++row)
  {
    DT_ sum(DT_(0));
    while (row_ptr[iter] == row && iter < used_elements)
    {
      sum += val[iter] * x[col_ptr[iter]];
      ++iter;
    }
    r[row] = rhs[row] - sum;
  }
}

template void Defect<Mem::Main, Algo::Generic>::coo(float *, const float * const, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void Defect<Mem::Main, Algo::Generic>::coo(double *, const double * const, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
