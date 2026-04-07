// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MATVECMULT_BANDED_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
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

    namespace ApplyBanded
    {
      template <typename IT_>
      inline Index start_offset(const Index i, const IT_ * const offsets,
        const Index rows, const Index cols, const Index noo)
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
          return Math::max(cols + Index(1), rows + cols - Index(offsets[i])) - cols - Index(1);
        }
      }

      template <typename IT_>
      inline Index end_offset(const Index i, const IT_ * const offsets,
        const Index rows, const Index cols, const Index noo)
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
          return Math::min(rows, cols + rows - Index(offsets[i]) - Index(1)) - Index(1);
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
        void f(DT_ * const r, const DT_ alpha,
          const DT_ * const val, const IT_ * const offsets,
          const DT_ * const x, const Index rows, const Index cols)
        {
          Index start(Math::max(Intern::ApplyBanded::start_offset(j-1, offsets, rows, cols, noo),
            Intern::ApplyBanded::end_offset(i-1, offsets, rows, cols, noo) + 1));
          Index end  (Math::min(Intern::ApplyBanded::start_offset(j-2, offsets, rows, cols, noo),
            Intern::ApplyBanded::end_offset(i-2, offsets, rows, cols, noo) + 1));

          FEAT_PRAGMA_IVDEP
          for (Index l = start; l < end; ++l)
          {
            DT_ tmp(0);
            Intern::LoopUnroller<0, i-j>::step(Intern::ApplyBanded::single_matrix_entry, &tmp, val + (j-1) * rows + l,
              offsets + (j-1), x + l + 1 - rows, rows);
            r[l] += alpha * tmp;
          }

          Iteration_Left<DT_, IT_, noo, i, j-1>::f(r, alpha, val, offsets, x, rows, cols);
        }
      };

      template <typename DT_, typename IT_, Index noo, Index i>
      struct Iteration_Left<DT_, IT_, noo, i, 0>
      {
        void f(DT_ * const /*r*/, const DT_ /*alpha*/,
          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
          const DT_ * const /*x*/, const Index /*rows*/, const Index /*cols*/)
        {
        }
      };

      /********************************************************************/

      template <typename DT_, typename IT_, Index noo, Index i>
      struct Iteration_Right
      {
        void f(DT_ * const r, const DT_ alpha,
          const DT_ * const val, const IT_ * const offsets,
          const DT_ * const x, const Index rows, const Index cols)
        {
          Iteration_Left<DT_, IT_, noo, i, i-1>::f(r, alpha, val, offsets, x, rows, cols);
          Iteration_Right<DT_, IT_, noo, i-1>::f(r, alpha, val, offsets, x, rows, cols);
        }
      };

      template <typename DT_, typename IT_, Index noo>
      struct Iteration_Right<DT_, IT_, noo, 0>
      {
        void f(DT_ * const /*r*/, const DT_ /*alpha*/,
          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
          const DT_ * const /*x*/, const Index /*rows*/, const Index /*cols*/)
        {
        }
      };

      /********************************************************************/

      template <typename DT_, typename IT_>
      void apply_banded_generic(DT_ * r, const DT_ alpha, const DT_ * const x,
        const DT_ * const val, const IT_ * const offsets,
        const Index bands, const Index rows, const Index cols)
      {
        // Search first offset of the upper triangular matrix
        Index k(0);
        while (k < bands && offsets[k] + 1 < rows)
        {
          ++k;
        }

        // iteration over all offsets of the lower triangular matrix
        for (Index i(k + 1); i > 0;)
        {
          --i;

          // iteration over all offsets of the upper triangular matrix
          for (Index j(bands + 1); j > 0;)
          {
            --j;

            // iteration over all rows which contain the offsets between offset i and offset j
            const Index start(Math::max(Intern::ApplyBanded::start_offset(  i, offsets, rows, cols, bands),
              Intern::ApplyBanded::end_offset  (  j, offsets, rows, cols, bands) + 1));
            const Index stop (Math::min(Intern::ApplyBanded::start_offset(i-1, offsets, rows, cols, bands),
              Intern::ApplyBanded::end_offset  (j-1, offsets, rows, cols, bands) + 1));
            for (Index l(start); l < stop; ++l)
            {
              DT_ s(0);
              for (Index a(i); a < j; ++a)
              {
                s += val[a * rows + l] * x[l + offsets[a] + 1 - rows];
              }
              r[l] += alpha * s;
            }
          }
        }
      }
    } // namespace ApplyBanded
  } // namespace Intern

  template <typename DT_, typename IT_>
  void MatVecMultBandedDense::exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y,
    const DT_ * const val, const IT_ * const offsets, const Index bands, const Index rows, const Index cols)
  {
    if(y == nullptr)
      Memory::memset_main(r, 0, sizeof(DT_) * rows);
    else if(r != y)
      Memory::memcopy_main(r, y, sizeof(DT_) * rows);

#ifdef FEAT_UNROLL_BANDED
    switch (bands)
    {
    case 3:
      Intern::ApplyBanded::Iteration_Right<DT_, IT_, 3, 4>::f(r, alpha, val, offsets, x, rows, cols);
      break;
    case 5:
      Intern::ApplyBanded::Iteration_Right<DT_, IT_, 5, 6>::f(r, alpha, val, offsets, x, rows, cols);
      break;
    case 9:
      Intern::ApplyBanded::Iteration_Right<DT_, IT_, 9, 10>::f(r, alpha, val, offsets, x, rows, cols);
      break;
    case 25:
      Intern::ApplyBanded::Iteration_Right<DT_, IT_, 25, 26>::f(r, alpha, val, offsets, x, rows, cols);
      break;
    default:
#ifdef DEBUG
      /// \todo print warning in feat log file
      std::cout << "Warning: Apply not optimized for " << bands << " offsets!\n";
#endif
      Intern::ApplyBanded::apply_banded_generic(r, alpha, x, val, offsets, bands, rows, cols);
    }
#else
    Intern::ApplyBanded::apply_banded_generic(r, alpha, x, val, offsets, bands, rows, cols);
#endif //FEAT_UNROLL_BANDED
  }
} // namespace FEAT::LAFEM::Arch
