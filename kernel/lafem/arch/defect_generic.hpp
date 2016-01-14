#pragma once
#ifndef KERNEL_LAFEM_ARCH_DEFECT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_DEFECT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_DEFECT_HPP
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

      template <typename DT_, typename IT_>
      void Defect<Mem::Main>::csr_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
      {
        for (Index row(0) ; row < rows ; ++row)
        {
          DT_ sum(0);
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i(row_ptr[row]) ; i < end ; ++i)
          {
            sum += val[i] * x[col_ind[i]];
          }
          r[row] = rhs[row] - sum;
        }
      }

      template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void Defect<Mem::Main>::csrb_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
      {
        Tiny::Vector<DT_, BlockHeight_> * br(reinterpret_cast<Tiny::Vector<DT_, BlockHeight_> *>(r));
        const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * const bval(reinterpret_cast<const Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *>(val));
        const Tiny::Vector<DT_, BlockWidth_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, BlockWidth_> *>(x));
        const Tiny::Vector<DT_, BlockHeight_> * const brhs(reinterpret_cast<const Tiny::Vector<DT_, BlockHeight_> *>(rhs));

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
          br[row] = brhs[row] - bsum;
        }
      }

      namespace Intern
      {
        namespace DefectELL
        {
          template <typename DT_>
          FORCE_INLINE void single_entry_defect(Index k, DT_ * const r, const DT_ * const b, const DT_ * const rhs)
          {
            r[k] = rhs[k] - b[k];
          }

          template <typename DT_, typename IT_, Index C_>
          struct DefectSpecialisation
          {
            static void f(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs,
                          const IT_ * const cl, const DT_ * const x, const Index /*C*/, const Index rows)
            {
              DT_ tmp[C_];
              const DT_ * const ctmp(static_cast<DT_ *>(tmp));

              for (Index i(0) ; i < rows/C_ ; ++i)
              {
                Intern::LoopUnroller<0, C_>::step(Intern::ProductMatVecELL::zero_entry, tmp);

                for (Index j(0) ; j < cl[i] ; ++j)
                {
                  Intern::LoopUnroller<0, C_>::step(Intern::ProductMatVecELL::single_matrix_entry, tmp, val + cs[i] + j*C_, x, col_ind + cs[i] + j*C_);
                }

                Intern::LoopUnroller<0, C_>::step(single_entry_defect, r + i*C_, ctmp, rhs + i*C_);
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

                  r[i*C_+k] = rhs[i*C_+k] - tmp[k];
                }
              }
            }
          };

          template <typename DT_, typename IT_>
          struct DefectGeneric
          {
            static void f(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs,
                          const IT_ * const cl, const DT_ * const x, const Index C, const Index rows)
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
                  r[i*C + k] = rhs[i*C + k] - r[i*C + k];
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

                  r[i*C+k] = rhs[i*C+k] - r[i*C+k];
                }
              }
            }
          };
        } // namespace DefectELL
      } // namespace Intern

      template <typename DT_, typename IT_>
      void Defect<Mem::Main>::ell_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val,
                                                 const IT_ * const col_ind, const IT_ * const cs,
                                                 const IT_ * const cl, const DT_ * const x,
                                                 const Index C, const Index rows)
      {
        switch (C)
        {
        case 2:
          Intern::DefectELL::DefectSpecialisation<DT_, IT_, 2>::f(r, rhs, val, col_ind, cs, cl, x, C, rows);
          break;
        case 4:
          Intern::DefectELL::DefectSpecialisation<DT_, IT_, 4>::f(r, rhs, val, col_ind, cs, cl, x, C, rows);
          break;
        case 8:
          Intern::DefectELL::DefectSpecialisation<DT_, IT_, 8>::f(r, rhs, val, col_ind, cs, cl, x, C, rows);
          break;
        case 16:
          Intern::DefectELL::DefectSpecialisation<DT_, IT_,16>::f(r, rhs, val, col_ind, cs, cl, x, C, rows);
          break;
        case 32:
          Intern::DefectELL::DefectSpecialisation<DT_, IT_,32>::f(r, rhs, val, col_ind, cs, cl, x, C, rows);
          break;
        default:
#ifdef DEBUG
          /// \todo print warning in feast log file
          std::cout << "Warning: Defect not optimized for chunk size = " << C << "!" << std::endl;
#endif
          Intern::DefectELL::DefectGeneric<DT_, IT_>::f(r, rhs, val, col_ind, cs, cl, x, C, rows);
        }
      }

      template <typename DT_, typename IT_>
      void Defect<Mem::Main>::coo_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const row_ptr, const IT_ * const col_ptr, const DT_ * const x, const Index rows, const Index, const Index used_elements)
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
          r[row] = rhs[row] - sum;
        }
      }


      namespace Intern
      {
        namespace DefectBanded
        {
          template <typename DT_, typename IT_, Index noo, Index i, Index j>
          struct Iteration_Left
          {
            static void f(DT_ * const r, const DT_ * y,
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
                  r[l] = y[l] - tmp;
                }

              Iteration_Left<DT_, IT_, noo, i, j-1>::f(r, y, val, offsets, x, rows, columns);
            }
          };

          template <typename DT_, typename IT_, Index noo, Index i>
          struct Iteration_Left<DT_, IT_, noo, i, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ * const /*y*/,
                          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }
          };

          /********************************************************************/

          template <typename DT_, typename IT_, Index noo, Index i>
          struct Iteration_Right
          {
            static void f(DT_ * const r, const DT_ * const y,
                          const DT_ * const val, const IT_ * const offsets,
                          const DT_ * const x, const Index rows, const Index columns)
            {
              Iteration_Left<DT_, IT_, noo, i, i-1>::f(r, y, val, offsets, x, rows, columns);
              Iteration_Right<DT_, IT_, noo, i-1>::f(r, y, val, offsets, x, rows, columns);
            }
          };

          template <typename DT_, typename IT_, Index noo>
          struct Iteration_Right<DT_, IT_, noo, 0>
          {
            static void f(DT_ * const /*r*/, const DT_ * const /*y*/,
                          const DT_ * const /*val*/, const IT_ * const /*offsets*/,
                          const DT_ * const /*x*/, const Index /*rows*/, const Index /*columns*/)
            {
            }
          };

          /********************************************************************/

          template <typename DT_, typename IT_>
          void defect_banded_generic(DT_ * r, const DT_ * const y, const DT_ * const val,
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
                  r[l] = y[l] - s;
                }
              }
            }
          }
        } // namespace DefectBanded
      } // namespace Intern



      template <typename DT_, typename IT_>
      void Defect<Mem::Main>::banded_generic(DT_ * r, const DT_ * const rhs,
                                                    const DT_ * const val, const IT_ * const offsets, const DT_ * const x,
                                                    const Index num_of_offsets, const Index rows, const Index columns)
      {
#ifdef FEAST_UNROLL_BANDED
        switch (num_of_offsets)
        {
        case 3:
          Intern::DefectBanded::Iteration_Right<DT_, IT_, 3, 4>::f(r, rhs, val, offsets, x, rows, columns);
          break;
        case 5:
          Intern::DefectBanded::Iteration_Right<DT_, IT_, 5, 6>::f(r, rhs, val, offsets, x, rows, columns);
          break;
        case 9:
          Intern::DefectBanded::Iteration_Right<DT_, IT_, 9, 10>::f(r, rhs, val, offsets, x, rows, columns);
          break;
        case 25:
          Intern::DefectBanded::Iteration_Right<DT_, IT_, 25, 26>::f(r, rhs, val, offsets, x, rows, columns);
          break;
        default:
#ifdef DEBUG
          /// \todo print warning in feast log file
          std::cout << "Warning: Defect not optimized for " << num_of_offsets << " offsets!" << std::endl;
#endif
          Intern::DefectBanded::defect_banded_generic(r, rhs, val, offsets, x, num_of_offsets, rows, columns);
        }
#else
        Intern::DefectBanded::defect_banded_generic(r, rhs, val, offsets, x, num_of_offsets, rows, columns);
#endif //FEAST_UNROLL_BANDED
      }

      template <typename DT_>
      void Defect<Mem::Main>::dense_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
      {
        for (Index row(0) ; row < rows ; ++row)
        {
          DT_ sum(0);
          for (Index col(0); col < columns; ++col)
          {
            sum += val[row * columns + col] * x[col];
          }
          r[row] = rhs[row] - sum;
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_DEFECT_GENERIC_HPP
