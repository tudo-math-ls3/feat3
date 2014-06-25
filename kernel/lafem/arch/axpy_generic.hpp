#pragma once
#ifndef KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_AXPY_HPP
  #error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>

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
    const IT_ end(row_ptr[row + 1]);
    for (IT_ i(row_ptr[row]) ; i < end ; ++i)
    {
      sum += val[i] * x[col_ind[i]];
    }
    r[row] = (sum * a) + y[row];
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
    r[row] = (sum * a[row]) + y[row];
  }
}

template <typename DT_, typename IT_, Index BlockHeight_, Index BlockWidth_>
void Axpy<Mem::Main, Algo::Generic>::csrb(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index)
{
  Tiny::Vector<DT_, BlockHeight_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockHeight_> *>(r));
  const Tiny::Vector<DT_, BlockHeight_> * const by(reinterpret_cast<const Tiny::Vector<DT_, BlockHeight_> * const>(y));
  const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const bval(reinterpret_cast<const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const>(val));
  const Tiny::Vector<DT_, BlockWidth_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockWidth_> * const>(x));

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
    br[row] = (bsum * a) + by[row];
  }
}

template <typename DT_, typename IT_, Index BlockHeight_, Index BlockWidth_>
void Axpy<Mem::Main, Algo::Generic>::csrb(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index)
{
  Tiny::Vector<DT_, BlockHeight_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockHeight_> *>(r));
  const Tiny::Vector<DT_, BlockHeight_> * const ba(reinterpret_cast<const Tiny::Vector<DT_, BlockHeight_> * const>(a));
  const Tiny::Vector<DT_, BlockHeight_> * const by(reinterpret_cast<const Tiny::Vector<DT_, BlockHeight_> * const>(y));
  const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const bval(reinterpret_cast<const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const>(val));
  const Tiny::Vector<DT_, BlockWidth_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockWidth_> * const>(x));

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
    br[row] = (bsum * ba[row]) + by[row];
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

    const IT_ max(Arl[row]);
    for(IT_ n(0); n < max ; n++)
    {
      const DT_ A_ij = *tAx;

      const IT_ col = *tAj;
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
  for (IT_ row(0); row < IT_(rows); ++row)
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


namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      namespace Intern
      {
        namespace AxpyBanded
        {
          template <typename DT_, typename IT_, Index noo, Index i, Index j>
          struct Iteration_Left
          {
            static void f(DT_ * const r, const DT_ * y, const DT_ alpha,
                          const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Index start(Math::max(Intern::ProductMatVecBanded::start_offset(j-1, offsets, rows, columns, noo),
                                    Intern::ProductMatVecBanded::end_offset(i-1, offsets, rows, columns, noo) + 1));
              Index end  (Math::min(Intern::ProductMatVecBanded::start_offset(j-2, offsets, rows, columns, noo),
                                    Intern::ProductMatVecBanded::end_offset(i-2, offsets, rows, columns, noo) + 1));

              FEAST_IVDEP
                for (Index l(start); l < end; ++l)
                {
                  r[l] = y[l] + alpha * Intern::ProductMatVecBanded::Single_Entry<DT_, IT_, i-j>::f(val + (j-1) * rows + l,
                                                                                                    offsets + (j-1), x + l + 1 - rows, rows);
                }

              Iteration_Left<DT_, IT_, noo, i, j-1>::f(r, y, alpha, val, offsets, x, rows, columns);
            }

            static void f(DT_ * const r, const DT_ * y, const DT_ * const alpha,
                          const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Index start(Math::max(Intern::ProductMatVecBanded::start_offset(j-1, offsets, rows, columns, noo),
                                    Intern::ProductMatVecBanded::end_offset(i-1, offsets, rows, columns, noo) + 1));
              Index end  (Math::min(Intern::ProductMatVecBanded::start_offset(j-2, offsets, rows, columns, noo),
                                    Intern::ProductMatVecBanded::end_offset(i-2, offsets, rows, columns, noo) + 1));

              FEAST_IVDEP
                for (Index l(start); l < end; ++l)
                {
                  r[l] = y[l] + alpha[l] * Intern::ProductMatVecBanded::Single_Entry<DT_, IT_, i-j>::f(val + (j-1) * rows + l,
                                                                                                       offsets + (j-1), x + l + 1 - rows, rows);
                }

              Iteration_Left<DT_, IT_, noo, i, j-1>::f(r, y, alpha, val, offsets, x, rows, columns);
            }
};

          template <typename DT_, typename IT_, Index noo, Index i>
          struct Iteration_Left<DT_, IT_, noo, i, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ * const /*y*/, const DT_ /*alpha*/,
                          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }

            static void f(DT_ * const /*r*/, const DT_ * const /*y*/, const DT_ * const /*alpha*/,
                          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }
};

          /********************************************************************/

          template <typename DT_, typename IT_, Index noo, Index i>
          struct Iteration_Right
          {
            static void f(DT_ * const r, const DT_ * const y, const DT_ alpha,
                          const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Iteration_Left<DT_, IT_, noo, i, i-1>::f(r, y, alpha, val, offsets, x, rows, columns);
              Iteration_Right<DT_, IT_, noo, i-1>::f(r, y, alpha, val, offsets, x, rows, columns);
            }

            static void f(DT_ * const r, const DT_ * const y, const DT_ * const alpha,
                          const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Iteration_Left<DT_, IT_, noo, i, i-1>::f(r, y, alpha, val, offsets, x, rows, columns);
              Iteration_Right<DT_, IT_, noo, i-1>::f(r, y, alpha, val, offsets, x, rows, columns);
            }
};

          template <typename DT_, typename IT_, Index noo>
          struct Iteration_Right<DT_, IT_, noo, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ * const /*y*/, const DT_ /*alpha*/,
                          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }

            static void f(DT_ * const /*r*/, const DT_ * const /*y*/, const DT_ * const /*alpha*/,
                          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }
};

          /********************************************************************/

          template <typename DT_, typename IT_>
          void axpy_banded_generic(DT_ * r, const DT_ * const y, const DT_ alpha, const DT_ * const val,
                                   const IT_ * const offsets, const DT_ * const x,
                                   const Index num_of_offsets, const Index rows, const Index columns)
          {
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
                const Index start(Math::max(Intern::ProductMatVecBanded::start_offset(  i, offsets, rows, columns, num_of_offsets),
                                            Intern::ProductMatVecBanded::end_offset  (  j, offsets, rows, columns, num_of_offsets) + 1));
                const Index stop (Math::min(Intern::ProductMatVecBanded::start_offset(i-1, offsets, rows, columns, num_of_offsets),
                                            Intern::ProductMatVecBanded::end_offset  (j-1, offsets, rows, columns, num_of_offsets) + 1));
                for (Index l(start); l < stop; ++l)
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
          }

          template <typename DT_, typename IT_>
          void axpy_banded_generic(DT_ * r, const DT_ * const y, const DT_ * const alpha, const DT_ * const val,
                                   const IT_ * const offsets, const DT_ * const x,
                                   const Index num_of_offsets, const Index rows, const Index columns)
          {
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
                const Index start(Math::max(Intern::ProductMatVecBanded::start_offset(  i, offsets, rows, columns, num_of_offsets),
                                            Intern::ProductMatVecBanded::end_offset  (  j, offsets, rows, columns, num_of_offsets) + 1));
                const Index stop (Math::min(Intern::ProductMatVecBanded::start_offset(i-1, offsets, rows, columns, num_of_offsets),
                                            Intern::ProductMatVecBanded::end_offset  (j-1, offsets, rows, columns, num_of_offsets) + 1));
                for (Index l(start); l < stop; ++l)
                {
                  DT_ s(0);
                  for (Index a(i); a < j; ++a)
                  {
                    s += val[a * rows + l] * x[l + offsets[a] + 1 - rows];
                  }
                  r[l] = y[l] + alpha[l] * s;
                }
              }
            }
          }
        } // namespace AxpyBanded
      } // namespace Intern
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::banded(DT_ * r, const DT_ * const y, const DT_ alpha,
                                            const DT_ * const val, const IT_ * const offsets, const DT_ * const x,
                                            const Index num_of_offsets, const Index rows, const Index columns)
{
  switch (num_of_offsets)
  {
  case 3:
    Intern::AxpyBanded::Iteration_Right<DT_, IT_, 3, 4>::f(r, y, alpha, val, offsets, x, rows, columns);
    break;
  case 5:
    Intern::AxpyBanded::Iteration_Right<DT_, IT_, 5, 6>::f(r, y, alpha, val, offsets, x, rows, columns);
    break;
  case 9:
    Intern::AxpyBanded::Iteration_Right<DT_, IT_, 9, 10>::f(r, y, alpha, val, offsets, x, rows, columns);
    break;
  case 25:
    Intern::AxpyBanded::Iteration_Right<DT_, IT_, 25, 26>::f(r, y, alpha, val, offsets, x, rows, columns);
    break;
  default:
#if DEBUG
    std::cout << "Axpy not optimized for " << num_of_offsets << " offsets!" << std::endl;
#endif
    Intern::AxpyBanded::axpy_banded_generic(r, y, alpha, val, offsets, x, num_of_offsets, rows, columns);
  }
}

template <typename DT_, typename IT_>
void Axpy<Mem::Main, Algo::Generic>::banded(DT_ * r, const DT_ * const y, const DT_ * const alpha,
                                            const DT_ * const val, const IT_ * const offsets, const DT_ * const x,
                                            const Index num_of_offsets, const Index rows, const Index columns)
{
  switch (num_of_offsets)
  {
  case 3:
    Intern::AxpyBanded::Iteration_Right<DT_, IT_, 3, 4>::f(r, y, alpha, val, offsets, x, rows, columns);
    break;
  case 5:
    Intern::AxpyBanded::Iteration_Right<DT_, IT_, 5, 6>::f(r, y, alpha, val, offsets, x, rows, columns);
    break;
  case 9:
    Intern::AxpyBanded::Iteration_Right<DT_, IT_, 9, 10>::f(r, y, alpha, val, offsets, x, rows, columns);
    break;
  case 25:
    Intern::AxpyBanded::Iteration_Right<DT_, IT_, 25, 26>::f(r, y, alpha, val, offsets, x, rows, columns);
    break;
  default:
#if DEBUG
    std::cout << "Axpy not optimized for " << num_of_offsets << " offsets!" << std::endl;
#endif
    Intern::AxpyBanded::axpy_banded_generic(r, y, alpha, val, offsets, x, num_of_offsets, rows, columns);
  }
}

#endif // KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
