#pragma once
#ifndef KERNEL_LAFEM_ARCH_DEFECT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_DEFECT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_DEFECT_HPP
  #error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_, typename IT_>
void Defect<Mem::Main, Algo::Generic>::csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
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

template <typename DT_, typename IT_, Index BlockHeight_, Index BlockWidth_>
void Defect<Mem::Main, Algo::Generic>::csrb(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
{
  Tiny::Vector<DT_, BlockHeight_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockHeight_> *>(r));
  const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const bval(reinterpret_cast<const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const>(val));
  const Tiny::Vector<DT_, BlockWidth_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockWidth_> * const>(x));
  const Tiny::Vector<DT_, BlockHeight_> * const brhs(reinterpret_cast<const Tiny::Vector<DT_, BlockHeight_> * const>(rhs));

  for (Index row(0) ; row < rows ; ++row)
  {
    Tiny::Vector<DT_, BlockHeight_> bsum(0);
    const IT_ end(row_ptr[row + 1]);
    for (IT_ i(row_ptr[row]) ; i < end ; ++i)
    {
      for (Index h(0) ; h < BlockHeight_ ; ++h)
      {
        for (Index w(0) ; w < BlockWidth_ ; ++w)
        {
          bsum[h] += bval[i][h][w] * bx[col_ind[i]][w];
        }
      }
    }
    br[row] = brhs[row] - bsum;
  }
}

template <typename DT_, typename IT_>
void Defect<Mem::Main, Algo::Generic>::ell(DT_ * r, const DT_ * const rhs, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    const IT_ * tAj(Aj);
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

template <typename DT_, typename IT_>
void Defect<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const row_ptr, const IT_ * const col_ptr, const DT_ * const x, const Index rows, const Index used_elements)
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
    r[row] = rhs[row] - sum;
  }
}

template <typename DT_, typename IT_>
void Defect<Mem::Main, Algo::Generic>::banded(DT_ * r, const DT_ * const y, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns)
{
#ifdef START_OFFSET
#warning Overwriting definition of START_OFFSET
#undef START_OFFSET
#endif

#ifdef END_OFFSET
#error Overwriting definition of END_OFFSET
#undef END_OFFSET
#endif

#define START_OFFSET(j) ((j == Index(-1)) ? rows : ((j == k) ? 0 : rows - offsets[j] - 1))
#define END_OFFSET(j) ((j == Index(-1)) ? rows : ((j == num_of_offsets) ? 0 : columns + rows - offsets[j] - 1))

  // Search first offset of the upper triangular matrix
  Index k(0);
  while (k < num_of_offsets && offsets[k] + 1 < rows)
  {
    ++k;
  }

  // iteration over all offsets of the lower triangular matrix
  for (Index i(k + 1); i > 0;)
  {
    --i;

    // iteration over all offsets of the upper triangular matrix
    for (Index j(num_of_offsets + 1); j > 0;)
    {
      --j;

      // iteration over all rows which contain the offsets between offset i and offset j
      for (Index l(Math::max(START_OFFSET(i), END_OFFSET(j))); l < Math::min(START_OFFSET(i-1), END_OFFSET(j-1)); ++l)
      {
        DT_ s(0);
        for (Index a(i); a < j; ++a)
        {
          s += val[a * rows + l] * x[l + offsets[a] + 1 - rows];
        }
        r[l] = y[l] - s;
      }
    }
  }

#undef START_OFFSET
#undef END_OFFSET
}

#endif // KERNEL_LAFEM_ARCH_DEFECT_GENERIC_HPP
