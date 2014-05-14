#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
  #error "Do not include this implementation-only header file directly!"
#endif

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_, typename IT_>
void ScaleRows<Mem::Main, Algo::Generic>::csr(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    const IT_ end(row_ptr[row + 1]);
    for (IT_ i(row_ptr[row]) ; i < end ; ++i)
    {
      r[i] = a[i] * x[row];
    }
  }
}

template <typename DT_, typename IT_>
void ScaleRows<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ * const a, const IT_ * const row_idx, const IT_ * const /*col_idx*/, const DT_ * const x, const Index, const Index, const Index used_elements)
{
  for (Index i(0) ; i < used_elements ; ++i)
  {
      r[i] = a[i] * x[row_idx[i]];
  }
}

template <typename DT_, typename IT_>
void ScaleRows<Mem::Main, Algo::Generic>::ell(DT_ * r, const DT_ * const Ax, const IT_ * const /*Aj*/, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    const DT_ * tAx(Ax);
    DT_ * tr(r);
    tAx += row;
    tr += row;

    const IT_ max(Arl[row]);
    for(IT_ n(0); n < max ; n++)
    {
      *tr = *tAx * x[row];

      tAx += stride;
      tr += stride;
    }
  }
}

// ***********************************************

template <typename DT_, typename IT_>
void ScaleCols<Mem::Main, Algo::Generic>::csr(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    const IT_ end(row_ptr[row + 1]);
    for (IT_ i(row_ptr[row]) ; i < end ; ++i)
    {
      r[i] = a[i] * x[col_ind[i]];
    }
  }
}

template <typename DT_, typename IT_>
void ScaleCols<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ * const a, const IT_ * const /*row_idx*/, const IT_ * const col_idx, const DT_ * const x, const Index, const Index, const Index used_elements)
{
  for (Index i(0) ; i < used_elements ; ++i)
  {
    r[i] = a[i] * x[col_idx[i]];
  }
}

template <typename DT_, typename IT_>
void ScaleCols<Mem::Main, Algo::Generic>::ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    const IT_ * tAj(Aj);
    const DT_ * tAx(Ax);
    DT_ * tr(r);
    tAj += row;
    tAx += row;
    tr += row;

    const IT_ max(Arl[row]);
    for(IT_ n(0); n < max ; n++)
    {
      *tr = *tAx * x[*tAj];

      tAj += stride;
      tAx += stride;
      tr += stride;
    }
  }
}

#endif // KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP
