// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
      void Apply::csr_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                                               const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index, const bool transposed)
      {
        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool::set_memory(r, DT_(0), (transposed?columns:rows));
        }
        else if (r != y)
        {
          MemoryPool::copy(r, y, (transposed?columns:rows));
        }

        if (transposed)
        {
          DT_ ba = b/a;
          FEAT_PRAGMA_OMP(parallel for)
          for (Index col = 0 ; col < columns ; ++col)
          {
            r[col] = ba * r[col];
          }

          for (Index row = 0 ; row < rows ; ++row)
          {
            for (Index i = row_ptr[row] ; i < row_ptr[row+1] ; ++i)
            {
              r[col_ind[i]] += val[i] * x[row];
            }
          }
          FEAT_PRAGMA_OMP(parallel for)
          for (Index col = 0 ; col < columns ; ++col)
          {
            r[col] = a * r[col];
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index row = 0 ; row < rows ; ++row)
          {
            DT_ sum(0);
            const IT_ end(row_ptr[row + 1]);
            for (IT_ i = row_ptr[row] ; i < end ; ++i)
            {
              sum += val[i] * x[col_ind[i]];
            }
            r[row] = (sum * a) + (b * r[row]);
          }
        }
      }

      template <typename DT_, typename IT_>
      void Apply::cscr_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                                               const IT_ * const col_ind, const IT_ * const row_ptr, const IT_ * const row_numbers,
                                               const Index used_rows, const Index rows, const Index columns, const Index, const bool transposed)
      {
        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool::set_memory(r, DT_(0), (transposed?columns:rows));
        }
        else if (r != y)
        {
          MemoryPool::copy(r, y, (transposed?columns:rows));
        }

        if (transposed)
        {
          DT_ ba = b/a;
          FEAT_PRAGMA_OMP(parallel for)
          for (Index col = 0 ; col < columns ; ++col)
          {
            r[col] = ba * r[col];
          }

          for (Index nzrow(0) ; nzrow < used_rows ; ++nzrow)
          {
            const Index row(row_numbers[nzrow]);
            for (Index i(row_ptr[nzrow]) ; i < row_ptr[nzrow+1] ; ++i)
            {
              r[col_ind[i]] += val[i] * x[row];
            }
          }
          FEAT_PRAGMA_OMP(parallel for)
          for (Index col = 0 ; col < columns ; ++col)
          {
            r[col] = a * r[col];
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index nzrow = 0 ; nzrow < used_rows ; ++nzrow)
          {
            const Index row(row_numbers[nzrow]);
            DT_ sum(0);
            const IT_ end(row_ptr[nzrow + 1]);
            for (IT_ i = row_ptr[nzrow] ; i < end ; ++i)
            {
              sum += val[i] * x[col_ind[i]];
            }
            r[row] = (sum * a) + (b * r[row]);
          }
        }
      }

      template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
      void Apply::bcsr_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                                                const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index)
      {
        Tiny::Vector<DT_, BlockHeight_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockHeight_> *>(r));
        const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const bval(reinterpret_cast<const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *>(val));
        const Tiny::Vector<DT_, BlockWidth_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockWidth_> *>(x));

        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool::set_memory(r, DT_(0), /*(transposed?columns:rows)*/ rows * BlockHeight_);
        }
        else if (r != y)
        {
          MemoryPool::copy(r, y, /*(transposed?columns:rows)*/ rows * BlockHeight_);
        }

        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          Tiny::Vector<DT_, BlockHeight_> bsum(0);
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i = row_ptr[row] ; i < end ; ++i)
          {
            bsum.add_mat_vec_mult(bval[i], bx[col_ind[i]]);
          }
          br[row] = (bsum * a) + (b * br[row]);
        }
      }

      template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
      void Apply::bcsr_transposed_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
        const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index)
      {
        Tiny::Vector<DT_, BlockWidth_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockWidth_> *>(r));
        const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const bval(reinterpret_cast<const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *>(val));
        const Tiny::Vector<DT_, BlockHeight_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockHeight_> *>(x));

        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool::set_memory(r, DT_(0), columns * BlockWidth_);
        }
        else if (r != y)
        {
          MemoryPool::copy(r, y, columns * BlockWidth_);
        }

        DT_ ba = b/a;
        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < columns ; ++col)
        {
          br[col] = ba * br[col];
        }
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          for (Index i(row_ptr[row]) ; i < row_ptr[row+1] ; ++i)
          {
            br[col_ind[i]].add_vec_mat_mult(bx[row], bval[i]);
          }
        }
        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < columns ; ++col)
        {
          br[col] = a * br[col];
        }
      }

      template <int BlockSize_, typename DT_, typename IT_>
      void Apply::csrsb_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index)
      {
        Tiny::Vector<DT_, BlockSize_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockSize_> *>(r));
        const Tiny::Vector<DT_, BlockSize_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockSize_> *>(x));

        if (Math::abs(b) < Math::eps<DT_>())
        {
          MemoryPool::set_memory(r, DT_(0), /*(transposed?columns:rows)*/ rows * BlockSize_);
        }
        else if (r != y)
        {
          MemoryPool::copy(r, y, /*(transposed?columns:rows)*/ rows * BlockSize_);
        }

        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
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

      } // namespace Intern

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

              FEAT_PRAGMA_IVDEP
              for (Index l = start; l < end; ++l)
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
      void Apply::banded_generic(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y,
                                                  const DT_ * const val, const IT_ * const offsets,
                                                  const Index num_of_offsets, const Index rows, const Index columns)
      {
        if (Math::abs(beta) < Math::eps<DT_>())
        {
          MemoryPool::set_memory(r, DT_(0), rows);
        }
        else if (r != y)
        {
          MemoryPool::copy(r, y, rows);
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
          std::cout << "Warning: Apply not optimized for " << num_of_offsets << " offsets!" << "\n";
#endif
          Intern::ApplyBanded::apply_banded_generic(r, alpha, x, beta, val, offsets, num_of_offsets, rows, columns);
        }
#else
        Intern::ApplyBanded::apply_banded_generic(r, alpha, x, beta, val, offsets, num_of_offsets, rows, columns);
#endif //FEAT_UNROLL_BANDED
      }

      template <typename DT_, typename IT_>
      void Apply::banded_transposed_generic(DT_ * DOXY(r), const DT_ DOXY(alpha), const DT_ * const DOXY(x), const DT_ DOXY(beta), const DT_ * const DOXY(y),
        const DT_ * const DOXY(val), const IT_ * const DOXY(offsets),
        const Index DOXY(num_of_offsets), const Index DOXY(rows), const Index DOXY(columns))
      {
        XABORTM("not implemented");
      }

      template <typename DT_>
      void Apply::dense_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
      {
        if (Math::abs(beta) < Math::eps<DT_>())
        {
          MemoryPool::set_memory(r, DT_(0), rows);
        }
        else if (r != y)
        {
          MemoryPool::copy(r, y, rows);
        }

        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          DT_ sum(0);
          for (Index col(0); col < columns; ++col)
          {
            sum += val[row * columns + col] * x[col];
          }
          r[row] = beta * r[row] + alpha * sum;
        }
      }

      template <typename DT_>
      void Apply::dense_transposed_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const y,
        const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
      {
        if (Math::abs(beta) < Math::eps<DT_>())
        {
          MemoryPool::set_memory(r, DT_(0), columns);
        }
        else if (r != y)
        {
          MemoryPool::copy(r, y, columns);
        }

        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < columns ; ++col)
        {
          DT_ sum(0);
          for (Index row(0); row < rows; ++row)
          {
            sum += val[row * columns + col] * x[row];
          }
          r[col] = beta * r[col] + alpha * sum;
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_APPLY_GENERIC_HPP
