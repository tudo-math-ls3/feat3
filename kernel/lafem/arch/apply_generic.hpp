#pragma once
#ifndef KERNEL_LAFEM_ARCH_APPLY_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_APPLY_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_APPLY_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_>
      void Apply<Mem::Main>::csr_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                                               const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
      {
        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool<Mem::Main>::set_memory(r, DT_(0), (transposed?columns:rows));
        }
        else if (r != y)
        {
          MemoryPool<Mem::Main>::copy(r, y, (transposed?columns:rows));
        }

        if (transposed)
        {
          for (Index col(0) ; col < columns ; ++col)
          {
            r[col] = b * r[col];
          }
          for (Index row(0) ; row < rows ; ++row)
          {
            for (Index i(row_ptr[row]) ; i < row_ptr[row+1] ; ++i)
            {
              r[col_ind[i]] += val[i] * x[row] * a;
            }
          }
        }
        else
        {
          for (Index row(0) ; row < rows ; ++row)
          {
            DT_ sum(0);
            const IT_ end(row_ptr[row + 1]);
            for (IT_ i(row_ptr[row]) ; i < end ; ++i)
            {
              sum += val[i] * x[col_ind[i]];
            }
            r[row] = (sum * a) + (b * r[row]);
          }
        }
      }

      template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void Apply<Mem::Main>::csrb_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                                                const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index)
      {
        Tiny::Vector<DT_, BlockHeight_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockHeight_> *>(r));
        const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const bval(reinterpret_cast<const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *>(val));
        const Tiny::Vector<DT_, BlockWidth_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockWidth_> *>(x));

        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool<Mem::Main>::set_memory(r, DT_(0), /*(transposed?columns:rows)*/ rows * BlockHeight_);
        }
        else if (r != y)
        {
          MemoryPool<Mem::Main>::copy(r, y, /*(transposed?columns:rows)*/ rows * BlockHeight_);
        }

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
          br[row] = (bsum * a) + (b * br[row]);
        }
      }

      template <typename DT_, typename IT_, int BlockSize_>
      void Apply<Mem::Main>::csrsb_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index)
      {
        Tiny::Vector<DT_, BlockSize_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockSize_> *>(r));
        const Tiny::Vector<DT_, BlockSize_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockSize_> *>(x));

        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool<Mem::Main>::set_memory(r, DT_(0), /*(transposed?columns:rows)*/ rows * BlockSize_);
        }
        else if (r != y)
        {
          MemoryPool<Mem::Main>::copy(r, y, /*(transposed?columns:rows)*/ rows * BlockSize_);
        }

        for (Index row(0) ; row < rows ; ++row)
        {
          Tiny::Vector<DT_, BlockSize_> bsum(0);
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i(row_ptr[row]) ; i < end ; ++i)
          {
            bsum += val[i] * bx[col_ind[i]];
          }
          br[row] = (bsum * a) + (b * br[row]);
        }
      }

      namespace Intern
      {
        template <Index Start, Index End, Index Step = 1>
        struct LoopUnroller
        {
          template <typename ... Params>
          static FORCE_INLINE void step(void f(Index, Params ...), Params... parameters)
          {
            f(Start, parameters ...);
            LoopUnroller<Start+Step, End, Step>::step(f, parameters ...);
          }
        };

        template <Index Start, Index Step>
        struct LoopUnroller<Start, Start, Step>
        {
          template <typename ... Params>
          static FORCE_INLINE void step(void (Index, Params ...), Params ...)
          {
          }
        };

        namespace ApplyELL
        {
          template <typename DT_>
          FORCE_INLINE void single_entry_apply(Index k, DT_ * const r, const DT_ * const b, const DT_ a, const DT_ * const y, const DT_ beta)
          {
            r[k] = b[k] * a + y[k] * beta;
          }

          template <typename DT_>
          FORCE_INLINE void zero_entry(Index k, DT_ * r)
          {
            r[k] = DT_(0);
          }

          template <typename DT_, typename IT_>
          FORCE_INLINE void single_matrix_entry(Index k, DT_ * const r, const DT_ * const val,
                                                const DT_ * const x, const IT_ * const col_ind)
          {
            r[k] += val[k] * x[col_ind[k]];
          }

          template <typename DT_, typename IT_, Index C_>
          struct ProductMatVecSpecialisation
          {
            static void f(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs,
                          const IT_ * const cl, const DT_ * const x, const Index /*C*/, const Index rows)
            {
              for (Index i(0) ; i < rows/C_ ; ++i)
              {
                Intern::LoopUnroller<0, C_>::step(zero_entry, r + i*C_);

                for (Index j(0) ; j < cl[i] ; ++j)
                {
                  Intern::LoopUnroller<0, C_>::step(single_matrix_entry, r + i*C_, val + cs[i] + j*C_, x, col_ind + cs[i] + j*C_);
                }
              }

              Index i(rows/C_);
              {
                for (Index k(0) ; k < rows%C_ ; ++k)
                {
                  r[i*C_+k] = DT_(0);

                  for (Index j(0) ; j < cl[i] ; ++j)
                  {
                    r[i*C_+k] += val[cs[i]+j*C_+k] * x[col_ind[cs[i]+j*C_+k]];
                  }
                }
              }
            }
          };

          template <typename DT_, typename IT_, Index C_>
          struct ApplySpecialisation
          {
            static void f(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b,
                          const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs,
                          const IT_ * const cl, const Index /*C*/, const Index rows)
            {
              DT_ tmp[C_];
              const DT_ * const ctmp(static_cast<DT_ *>(tmp));

              for (Index i(0) ; i < rows/C_ ; ++i)
              {
                Intern::LoopUnroller<0, C_>::step(Intern::ApplyELL::zero_entry, tmp);

                for (Index j(0) ; j < cl[i] ; ++j)
                {
                  Intern::LoopUnroller<0, C_>::step(Intern::ApplyELL::single_matrix_entry, tmp, val + cs[i] + j*C_, x, col_ind + cs[i] + j*C_);
                }

                Intern::LoopUnroller<0, C_>::step(single_entry_apply, r + i*C_, ctmp, a, ((const DT_ *)r) + i*C_, b);
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

                  r[i*C_+k] = tmp[k] * a + b * r[i*C_+k];
                }
              }
            }
          };

          template <typename DT_, typename IT_>
          struct ApplyGeneric
          {
            static void f(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b,
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
                  r[i*C+k] = r[i*C+k] * a + b * r[i*C+k];
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

                  r[i*C+k] = r[i*C+k] * a + b * r[i*C+k];
                }
              }
            }
          };
        } // namespace ApplyELL
      } // namespace Intern

      template <typename DT_, typename IT_>
      void Apply<Mem::Main>::ell_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                                               const IT_ * const col_ind, const IT_ * const cs,
                                               const IT_ * const cl,
                                               const Index C, const Index rows)
      {
        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool<Mem::Main>::set_memory(r, DT_(0), rows);
        }
        else if (r != y)
        {
          MemoryPool<Mem::Main>::copy(r, y, rows);
        }

        switch (C)
        {
        case 2:
          Intern::ApplyELL::ApplySpecialisation<DT_, IT_, 2>::f(r, a, x, b, val, col_ind, cs, cl, C, rows);
          break;
        case 4:
          Intern::ApplyELL::ApplySpecialisation<DT_, IT_, 4>::f(r, a, x, b, val, col_ind, cs, cl, C, rows);
          break;
        case 8:
          Intern::ApplyELL::ApplySpecialisation<DT_, IT_, 8>::f(r, a, x, b, val, col_ind, cs, cl, C, rows);
          break;
        case 16:
          Intern::ApplyELL::ApplySpecialisation<DT_, IT_,16>::f(r, a, x, b, val, col_ind, cs, cl, C, rows);
          break;
        case 32:
          Intern::ApplyELL::ApplySpecialisation<DT_, IT_,32>::f(r, a, x, b, val, col_ind, cs, cl, C, rows);
          break;
        default:
#ifdef DEBUG
          /// \todo print warning in feat log file
          std::cout << "Warning: Apply not optimized for chunk size = " << C << "!" << std::endl;
#endif
          Intern::ApplyELL::ApplyGeneric<DT_, IT_>::f(r, a, x, b, val, col_ind, cs, cl, C, rows);
        }
      }

      template <typename DT_, typename IT_>
      void Apply<Mem::Main>::coo_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                                               const IT_ * const row_ptr, const IT_ * const col_ptr, const Index rows, const Index, const Index used_elements)
      {
        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool<Mem::Main>::set_memory(r, DT_(0), rows);
        }
        else if (r != y)
        {
          MemoryPool<Mem::Main>::copy(r, y, rows);
        }

        Index iter(0);
        for (IT_ row(0); row < IT_(rows); ++row)
        {
          DT_ sum(DT_(0));
          while (iter < used_elements && row_ptr[iter] == row)
          {
            sum += val[iter] * x[col_ptr[iter]];
            ++iter;
          }
          r[row] = a * sum + b * r[row];
        }
      }

      namespace Intern
      {
        namespace ApplyBanded
        {
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
              return Math::max(columns + Index(1), rows + columns - Index(offsets[i])) - columns - Index(1);
            }
          }

          template <typename IT_>
          inline Index end_offset(const Index i, const IT_ * const offsets,
                                  const Index rows, const Index columns, const Index noo)
          {
            if (i == Index (-1))
            {
              return rows - 1;
            }
            else if (i == noo)
            {
              return Index(-1);
            }
            else
            {
              return Math::min(rows, columns + rows - Index(offsets[i]) - Index(1)) - Index(1);
            }
          }

          template <typename DT_, typename IT_>
          FORCE_INLINE void single_matrix_entry(Index k, DT_ * const res, const DT_ * const val,
                                                const IT_ * const offsets, const DT_ * const x, Index rows)
          {
            *res += val[k * rows] * x[offsets[k]];
          }

          template <typename DT_, typename IT_, Index noo, Index i, Index j>
          struct Iteration_Left
          {
            static void f(DT_ * const r, const DT_ alpha, const DT_ beta,
                          const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Index start(Math::max(Intern::ApplyBanded::start_offset(j-1, offsets, rows, columns, noo),
                                    Intern::ApplyBanded::end_offset(i-1, offsets, rows, columns, noo) + 1));
              Index end  (Math::min(Intern::ApplyBanded::start_offset(j-2, offsets, rows, columns, noo),
                                    Intern::ApplyBanded::end_offset(i-2, offsets, rows, columns, noo) + 1));

              FEAT_IVDEP
                for (Index l(start); l < end; ++l)
                {
                  DT_ tmp(0);
                  Intern::LoopUnroller<0, i-j>::step(Intern::ApplyBanded::single_matrix_entry, &tmp, val + (j-1) * rows + l,
                                                     offsets + (j-1), x + l + 1 - rows, rows);
                  r[l] = beta * r[l] + alpha * tmp;
                }

              Iteration_Left<DT_, IT_, noo, i, j-1>::f(r, alpha, beta, val, offsets, x, rows, columns);
            }
          };

          template <typename DT_, typename IT_, Index noo, Index i>
          struct Iteration_Left<DT_, IT_, noo, i, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ /*alpha*/, const DT_ /*beta*/,
                          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }
          };

          /********************************************************************/

          template <typename DT_, typename IT_, Index noo, Index i>
          struct Iteration_Right
          {
            static void f(DT_ * const r, const DT_ alpha, const DT_ beta,
                          const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Iteration_Left<DT_, IT_, noo, i, i-1>::f(r, alpha, beta, val, offsets, x, rows, columns);
              Iteration_Right<DT_, IT_, noo, i-1>::f(r, alpha, beta, val, offsets, x, rows, columns);
            }
          };

          template <typename DT_, typename IT_, Index noo>
          struct Iteration_Right<DT_, IT_, noo, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ /*alpha*/, const DT_ /*beta*/,
                          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }
          };

          /********************************************************************/

          template <typename DT_, typename IT_>
          void apply_banded_generic(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta,
                                   const DT_ * const val, const IT_ * const offsets,
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
                const Index start(Math::max(Intern::ApplyBanded::start_offset(  i, offsets, rows, columns, num_of_offsets),
                                            Intern::ApplyBanded::end_offset  (  j, offsets, rows, columns, num_of_offsets) + 1));
                const Index stop (Math::min(Intern::ApplyBanded::start_offset(i-1, offsets, rows, columns, num_of_offsets),
                                            Intern::ApplyBanded::end_offset  (j-1, offsets, rows, columns, num_of_offsets) + 1));
                for (Index l(start); l < stop; ++l)
                {
                  DT_ s(0);
                  for (Index a(i); a < j; ++a)
                  {
                    s += val[a * rows + l] * x[l + offsets[a] + 1 - rows];
                  }
                  r[l] = beta * r[l] + alpha * s;
                }
              }
            }
          }
        } // namespace ApplyBanded
      } // namespace Intern

      template <typename DT_, typename IT_>
      void Apply<Mem::Main>::banded_generic(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y,
                                                  const DT_ * const val, const IT_ * const offsets,
                                                  const Index num_of_offsets, const Index rows, const Index columns)
      {
        if (Math::abs(beta) < Math::eps<DT_>())
        {
          MemoryPool<Mem::Main>::set_memory(r, DT_(0), rows);
        }
        else if (r != y)
        {
          MemoryPool<Mem::Main>::copy(r, y, rows);
        }
#ifdef FEAT_UNROLL_BANDED
        switch (num_of_offsets)
        {
        case 3:
          Intern::ApplyBanded::Iteration_Right<DT_, IT_, 3, 4>::f(r, alpha, beta, val, offsets, x, rows, columns);
          break;
        case 5:
          Intern::ApplyBanded::Iteration_Right<DT_, IT_, 5, 6>::f(r, alpha, beta, val, offsets, x, rows, columns);
          break;
        case 9:
          Intern::ApplyBanded::Iteration_Right<DT_, IT_, 9, 10>::f(r, alpha, beta, val, offsets, x, rows, columns);
          break;
        case 25:
          Intern::ApplyBanded::Iteration_Right<DT_, IT_, 25, 26>::f(r, alpha, beta, val, offsets, x, rows, columns);
          break;
        default:
#ifdef DEBUG
          /// \todo print warning in feat log file
          std::cout << "Warning: Apply not optimized for " << num_of_offsets << " offsets!" << std::endl;
#endif
          Intern::ApplyBanded::apply_banded_generic(r, alpha, x, beta, val, offsets, num_of_offsets, rows, columns);
        }
#else
        Intern::ApplyBanded::apply_banded_generic(r, alpha, x, beta, val, offsets, num_of_offsets, rows, columns);
#endif //FEAT_UNROLL_BANDED
      }

      template <typename DT_>
      void Apply<Mem::Main>::dense_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
      {
        if (Math::abs(beta) < Math::eps<DT_>())
        {
          MemoryPool<Mem::Main>::set_memory(r, DT_(0), rows);
        }
        else if (r != y)
        {
          MemoryPool<Mem::Main>::copy(r, y, rows);
        }

        for (Index row(0) ; row < rows ; ++row)
        {
          DT_ sum(0);
          for (Index col(0); col < columns; ++col)
          {
            sum += val[row * columns + col] * x[col];
          }
          r[row] = beta * r[row] + alpha * sum;
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_APPLY_GENERIC_HPP
