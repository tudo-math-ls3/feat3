#pragma once
#ifndef KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_AXPY_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {

      template <typename DT_>
      void Axpy<Mem::Main>::dv_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
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

      template <typename DT_, typename IT_>
      void Axpy<Mem::Main>::csr_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
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

      template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void Axpy<Mem::Main>::csrb_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
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
            for (int h(0) ; h < BlockHeight_ ; ++h)
            {
              for (int w(0) ; w < BlockWidth_ ; ++w)
              {
                bsum[h] += bval[i][h][w] * bx[col_ind[i]][w];
              }
            }
          }
          br[row] = (bsum * a) + by[row];
        }
      }

      namespace Intern
      {
        namespace AxpyELL
        {
          template <typename DT_>
          FORCE_INLINE void single_entry_axpy(Index k, DT_ * const r, const DT_ * const b, const DT_ a, const DT_ * const y)
          {
            r[k] = b[k] * a + y[k];
          }

          template <typename DT_, typename IT_, Index C_>
          struct AxpySpecialisation
          {
            static void f(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y,
                          const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs,
                          const IT_ * const cl, const Index /*C*/, const Index rows)
            {
              DT_ tmp[C_];
              const DT_ * const ctmp(static_cast<const DT_ * const>(tmp));

              for (Index i(0) ; i < rows/C_ ; ++i)
              {
                Intern::LoopUnroller<0, C_>::step(Intern::ProductMatVecELL::zero_entry, tmp);

                for (Index j(0) ; j < cl[i] ; ++j)
                {
                  Intern::LoopUnroller<0, C_>::step(Intern::ProductMatVecELL::single_matrix_entry, tmp, val + cs[i] + j*C_, x, col_ind + cs[i] + j*C_);
                }

                Intern::LoopUnroller<0, C_>::step(single_entry_axpy, r + i*C_, ctmp, a, y + i*C_);
              }

              Index i(rows/C_);
              {
                for (Index k(0) ; k < rows%C_ ; ++k)
                {
                  tmp[k] = DT_(0);

                  for (Index j(0) ; j < cl[i] ; ++j)
                  {
                    tmp[k] += val[cs[i]+j*C_+k] * x[col_ind[cs[i]+j*C_+k]];
                  }

                  r[i*C_+k] = tmp[k] * a + y[i*C_+k];
                }
              }
            }
          };

          template <typename DT_, typename IT_>
          struct AxpyGeneric
          {
            static void f(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y,
                          const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs,
                          const IT_ * const cl, const Index C, const Index rows)
            {
              for (Index i(0) ; i < rows/C ; ++i)
              {
                for (Index k(0); k < C; ++k)
                {
                  r[i*C + k] = DT_(0);
                }

                for (Index j(0) ; j < cl[i] ; ++j)
                {
                  for (Index k(0); k < C; ++k)
                  {
                    r[i*C+k] += val[cs[i]+j*C+k] * x[col_ind[cs[i]+j*C+k]];
                  }
                }

                for (Index k(0); k < C; ++k)
                {
                  r[i*C+k] = r[i*C+k] * a + y[i*C+k];
                }
              }

              Index i(rows/C);
              {
                for (Index k(0) ; k < rows%C ; ++k)
                {
                  r[i*C+k] = DT_(0);

                  for (Index j(0) ; j < cl[i] ; ++j)
                  {
                    r[i*C+k] += val[cs[i]+j*C+k] * x[col_ind[cs[i]+j*C+k]];
                  }

                  r[i*C+k] = r[i*C+k] * a + y[i*C+k];
                }
              }
            }
          };
        } // namespace AxpyELL
      } // namespace Intern

      template <typename DT_, typename IT_>
      void Axpy<Mem::Main>::ell_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
                                               const IT_ * const col_ind, const IT_ * const cs,
                                               const IT_ * const cl,
                                               const Index C, const Index rows)
      {
        switch (C)
        {
        case 2:
          Intern::AxpyELL::AxpySpecialisation<DT_, IT_, 2>::f(r, a, x, y, val, col_ind, cs, cl, C, rows);
          break;
        case 4:
          Intern::AxpyELL::AxpySpecialisation<DT_, IT_, 4>::f(r, a, x, y, val, col_ind, cs, cl, C, rows);
          break;
        case 8:
          Intern::AxpyELL::AxpySpecialisation<DT_, IT_, 8>::f(r, a, x, y, val, col_ind, cs, cl, C, rows);
          break;
        case 16:
          Intern::AxpyELL::AxpySpecialisation<DT_, IT_,16>::f(r, a, x, y, val, col_ind, cs, cl, C, rows);
          break;
        case 32:
          Intern::AxpyELL::AxpySpecialisation<DT_, IT_,32>::f(r, a, x, y, val, col_ind, cs, cl, C, rows);
          break;
        default:
#ifdef DEBUG
          /// \todo print warning in feast log file
          std::cout << "Warning: Axpy not optimized for chunk size = " << C << "!" << std::endl;
#endif
          Intern::AxpyELL::AxpyGeneric<DT_, IT_>::f(r, a, x, y, val, col_ind, cs, cl, C, rows);
        }
      }

      template <typename DT_, typename IT_>
      void Axpy<Mem::Main>::coo_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
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
                  DT_ tmp(0);
                  Intern::LoopUnroller<0, i-j>::step(Intern::ProductMatVecBanded::single_matrix_entry, &tmp, val + (j-1) * rows + l,
                                                     offsets + (j-1), x + l + 1 - rows, rows);
                  r[l] = y[l] + alpha * tmp;
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
          };

          template <typename DT_, typename IT_, Index noo>
          struct Iteration_Right<DT_, IT_, noo, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ * const /*y*/, const DT_ /*alpha*/,
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
        } // namespace AxpyBanded
      } // namespace Intern

      template <typename DT_, typename IT_>
      void Axpy<Mem::Main>::banded_generic(DT_ * r, const DT_ * const y, const DT_ alpha,
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
#ifdef DEBUG
          /// \todo print warning in feast log file
          std::cout << "Warning: Axpy not optimized for " << num_of_offsets << " offsets!" << std::endl;
#endif
          Intern::AxpyBanded::axpy_banded_generic(r, y, alpha, val, offsets, x, num_of_offsets, rows, columns);
        }
      }

      template <typename DT_>
      void Axpy<Mem::Main>::dense_generic(DT_ * r, const DT_ alpha, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
      {
        for (Index row(0) ; row < rows ; ++row)
        {
          DT_ sum(0);
          for (Index col(0); col < columns; ++col)
          {
            sum += val[row * columns + col] * x[col];
          }
          r[row] = y[row] + alpha * sum;
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
