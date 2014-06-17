#pragma once
#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATVEC_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_PRODUCT_MATVEC_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATVEC_HPP
  #error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>

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

template <typename DT_, typename IT_, Index BlockHeight_, Index BlockWidth_>
void ProductMatVec<Mem::Main, Algo::Generic>::csrb(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
{
  Tiny::Vector<DT_, BlockHeight_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockHeight_> *>(r));
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
    br[row] = bsum;
  }
}

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

template <typename DT_, typename IT_>
void ProductMatVec<Mem::Main, Algo::Generic>::coo(DT_ * r, const DT_ * const val, const IT_ * const row_ptr, const IT_ * const col_ptr, const DT_ * const x, const Index rows, const Index used_elements)
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
    r[row] = sum;
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
        namespace ProductMatVecBanded
        {
          template <typename DT_, typename IT_, Index i>
          struct Single_Entry
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return Single_Entry<DT_, IT_, i-9>::f(val + 9*rows, offsets + 9, x , rows) +
                *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[2 * rows] * x[offsets[2]] +
                val[3 * rows] * x[offsets[3]] +
                val[4 * rows] * x[offsets[4]] +
                val[5 * rows] * x[offsets[5]] +
                val[6 * rows] * x[offsets[6]] +
                val[7 * rows] * x[offsets[7]] +
                val[8 * rows] * x[offsets[8]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 25>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[ 2 * rows] * x[offsets[ 2]] +
                val[ 3 * rows] * x[offsets[ 3]] +
                val[ 4 * rows] * x[offsets[ 4]] +
                val[ 5 * rows] * x[offsets[ 5]] +
                val[ 6 * rows] * x[offsets[ 6]] +
                val[ 7 * rows] * x[offsets[ 7]] +
                val[ 8 * rows] * x[offsets[ 8]] +
                val[ 9 * rows] * x[offsets[ 9]] +
                val[10 * rows] * x[offsets[10]] +
                val[11 * rows] * x[offsets[11]] +
                val[12 * rows] * x[offsets[12]] +
                val[13 * rows] * x[offsets[13]] +
                val[14 * rows] * x[offsets[14]] +
                val[15 * rows] * x[offsets[15]] +
                val[16 * rows] * x[offsets[16]] +
                val[17 * rows] * x[offsets[17]] +
                val[18 * rows] * x[offsets[18]] +
                val[19 * rows] * x[offsets[19]] +
                val[20 * rows] * x[offsets[20]] +
                val[21 * rows] * x[offsets[21]] +
                val[22 * rows] * x[offsets[22]] +
                val[23 * rows] * x[offsets[23]] +
                val[24 * rows] * x[offsets[24]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 1>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index /*rows*/)
            {
              return *val * x[*offsets];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 2>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 3>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[2 * rows] * x[offsets[2]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 4>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[2 * rows] * x[offsets[2]] +
                val[3 * rows] * x[offsets[3]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 5>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[2 * rows] * x[offsets[2]] +
                val[3 * rows] * x[offsets[3]] +
                val[4 * rows] * x[offsets[4]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 6>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[2 * rows] * x[offsets[2]] +
                val[3 * rows] * x[offsets[3]] +
                val[4 * rows] * x[offsets[4]] +
                val[5 * rows] * x[offsets[5]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 7>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[2 * rows] * x[offsets[2]] +
                val[3 * rows] * x[offsets[3]] +
                val[4 * rows] * x[offsets[4]] +
                val[5 * rows] * x[offsets[5]] +
                val[6 * rows] * x[offsets[6]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 8>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[2 * rows] * x[offsets[2]] +
                val[3 * rows] * x[offsets[3]] +
                val[4 * rows] * x[offsets[4]] +
                val[5 * rows] * x[offsets[5]] +
                val[6 * rows] * x[offsets[6]] +
                val[7 * rows] * x[offsets[7]];
            }
          };

          template <typename DT_, typename IT_>
          struct Single_Entry<DT_, IT_, 9>
          {
            static inline DT_ f(const DT_ * const val, const IT_ * const offsets,
                                const DT_ * const x, const Index rows)
            {
              return *val * x[*offsets] +
                val[rows] * x[offsets[1]] +
                val[2 * rows] * x[offsets[2]] +
                val[3 * rows] * x[offsets[3]] +
                val[4 * rows] * x[offsets[4]] +
                val[5 * rows] * x[offsets[5]] +
                val[6 * rows] * x[offsets[6]] +
                val[7 * rows] * x[offsets[7]] +
                val[8 * rows] * x[offsets[8]];
            }
          };

          /********************************************************************/

          template <typename IT_>
          inline Index start_offset(const Index i, const IT_ * const offsets,
                                    const Index rows, const Index columns, const Index noo)
          {
            if (i == Index(-1))
            {
              return rows;
            }
            else if (i == noo)
            {
              return Index(0);
            }
            else
            {
              return std::max(columns + Index(1), rows + columns - offsets[i]) - columns - Index(1);
            }
          }

          template <typename IT_>
          inline Index end_offset(const Index i, const IT_ * const offsets,
                                  const Index rows, const Index columns, const Index noo)
          {
            if (i == Index (-1))
            {
              return rows;
            }
            else if (i == noo)
            {
              return Index(0);
            }
            else
            {
              return std::min(rows, columns + rows - offsets[i] - Index(1));
            }
          }

          /********************************************************************/

          template <typename DT_, typename IT_, Index noo, Index i, Index j>
          struct Iteration_Left
          {
            static void f(DT_ * const r, const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Index start(std::max(start_offset(j-1, offsets, rows, columns, noo),
                                   end_offset(i-1, offsets, rows, columns, noo)));
              Index end  (std::min(start_offset(j-2, offsets, rows, columns, noo),
                                   end_offset(i-2, offsets, rows, columns, noo)));

              FEAST_IVDEP
              for (Index l(start); l < end; ++l)
              {
                r[l] = Single_Entry<DT_, IT_, i-j>::f(val + (j-1) * rows + l,
                                                      offsets + (j-1), x + l + 1 - rows, rows);
              }

              Iteration_Left<DT_, IT_, noo, i, j-1>::f(r, val, offsets, x, rows, columns);
            }
          };

          template <typename DT_, typename IT_, Index noo, Index i>
          struct Iteration_Left<DT_, IT_, noo, i, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }
          };

          /********************************************************************/

          template <typename DT_, typename IT_, Index noo, Index i>
          struct Iteration_Right
          {
            static void f(DT_ * const r, const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Iteration_Left<DT_, IT_, noo, i, i-1>::f(r, val, offsets, x, rows, columns);
              Iteration_Right<DT_, IT_, noo, i-1>::f(r, val, offsets, x, rows, columns);
            }
          };

          template <typename DT_, typename IT_, Index noo>
          struct Iteration_Right<DT_, IT_, noo, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }
          };

          /********************************************************************/

          template <typename DT_, typename IT_>
          void product_matvec_banded_generic(DT_ * r, const DT_ * const val,
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
                const Index start(Math::max(start_offset(  i, offsets, rows, columns, num_of_offsets),
                                            end_offset  (  j, offsets, rows, columns, num_of_offsets)));
                const Index stop (Math::min(start_offset(i-1, offsets, rows, columns, num_of_offsets),
                                            end_offset  (j-1, offsets, rows, columns, num_of_offsets)));
                for (Index l(start); l < stop; ++l)
                {
                  DT_ s(0);
                  for (Index a(i); a < j; ++a)
                  {
                    s += val[a * rows + l] * x[l + offsets[a] + 1 - rows];
                  }
                  r[l] = s;
                }
              }
            }
          }
        } // namespace ProductMatVecBanded
      } // namespace Intern
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

template <typename DT_, typename IT_>
void ProductMatVec<Mem::Main, Algo::Generic>::banded(DT_ * r, const DT_ * const val,
                                                     const IT_ * const offsets, const DT_ * const x,
                                                     const Index num_of_offsets, const Index rows, const Index columns)
{
  switch (num_of_offsets)
  {
  case 3:
    Intern::ProductMatVecBanded::Iteration_Right<DT_, IT_, 3, 4>::f(r, val, offsets, x, rows, columns);
    break;
  case 5:
    Intern::ProductMatVecBanded::Iteration_Right<DT_, IT_, 5, 6>::f(r, val, offsets, x, rows, columns);
    break;
  case 9:
    Intern::ProductMatVecBanded::Iteration_Right<DT_, IT_, 9,10>::f(r, val, offsets, x, rows, columns);
    break;
  case 25:
    Intern::ProductMatVecBanded::Iteration_Right<DT_, IT_,25,26>::f(r, val, offsets, x, rows, columns);
    break;
  default:
#if DEBUG
    std::cout << "ProductMatVec not optimized for " << num_of_offsets << " offsets!" << std::endl;
#endif
    Intern::ProductMatVecBanded::product_matvec_banded_generic(r, val, offsets, x, num_of_offsets, rows, columns);
  }
}

#endif // KERNEL_LAFEM_ARCH_PRODUCT_MATVEC_GENERIC_HPP
