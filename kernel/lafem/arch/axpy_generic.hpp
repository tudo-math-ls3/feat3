#pragma once
#ifndef KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_AXPY_HPP
  #error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Axpy<Mem::Main, Algo::Generic>::dv(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
{
  if (r == y)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] += a * x[i];
    }
  }
  else if (r == x)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] *= a;
      r[i]+= y[i];
    }
  }
  else
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] = a * x[i] + y[i];
    }
  }
}


template <typename DT_>
void Axpy<Mem::Main, Algo::Generic>::dv(DT_ * r, const DT_  * const a, const DT_ * const x, const DT_ * const y, const Index size)
{
  if (r == y)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] += a[i] * x[i];
    }
  }
  else if(r == x)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] *= a[i];
      r[i]+= y[i];
    }
  }
  else if(r == a)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] *= x[i];
      r[i] += y[i];
    }
  }
  else
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] = a[i] * x[i] + y[i];
    }
  }
}


template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    DT_ sum(0);
    const Index end(row_ptr[row + 1]);
    for (Index i(row_ptr[row]) ; i < end ; ++i)
    {
      sum += val[i] * x[col_ind[i]];
    }
    r[row] = (sum * a)+ y[row];
  }
}

template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::csr(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index)
{
  for (Index row(0) ; row < rows ; ++row)
  {
    DT_ sum(0);
    const Index end(row_ptr[row + 1]);
    for (Index i(row_ptr[row]) ; i < end ; ++i)
    {
      sum += val[i] * x[col_ind[i]];
    }
    r[row] = (sum * a[row])+ y[row];
  }
}

template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::ell(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const Index stride, const Index rows)
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
    r[row] = (sum * a) + y[row];
  }
}

template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::ell(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const Index stride, const Index rows)
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
    r[row] = (sum * a[row]) + y[row];
  }
}

template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const row_ptr, const IT_ * const col_ptr, const Index rows, const Index used_elements)
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
    r[row] = a * sum + y[row];
  }
}

template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const row_ptr, const IT_ * const col_ptr, const Index rows, const Index used_elements)
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
    r[row] = a[row] * sum + y[row];
  }
}

template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::banded(DT_ * r, const DT_ * const y, const DT_ alpha, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns)
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
        r[l] = y[l] + alpha * s;
      }
    }
  }

#undef START_OFFSET
#undef END_OFFSET
}


template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::banded(DT_ * r, const DT_ * const y, const DT_ * const alpha, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns)
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
          s += val[a * rows + l] * x[l - columns + offsets[a]];
        }
        r[l] = y[l] + alpha[l] * s;
      }
    }
  }

#undef START_OFFSET
#undef END_OFFSET
}

#endif // KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
