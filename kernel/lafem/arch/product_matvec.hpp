#pragma once
#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATVEC_HPP
#define KERNEL_LAFEM_ARCH_PRODUCT_MATVEC_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <iostream> // TODO: Muss noch entfernt werden

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_, typename Algo_>
      struct ProductMatVec;

      template <>
      struct ProductMatVec<Mem::Main, Algo::Generic>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
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

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
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
        static void coo(DT_ * r, const DT_ * const val, const IT_ * const row_ptr, const IT_ * const col_ptr, const DT_ * const x, const Index rows, const Index used_elements)
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
            r[row] = sum;
          }
        }

        template <typename DT_>
        static void band_q2_d2(DT_ * r, const DT_ * const val, const DT_ * const x, const Index nodes_per_row, const Index nodes_per_column)
        {
          const Index rows(nodes_per_row * nodes_per_column);
          const Index m(nodes_per_row);

          r[0] = val[4 * rows] * x[0]
            + val[5 * rows] * x[1]
            + val[6 * rows] * x[m - 1]
            + val[7 * rows] * x[m]
            + val[8 * rows] * x[m + 1];

          for (Index i(1); i < m - 1; ++i)
          {
            r[i] = val[3 * rows + i] * x[i - 1]
              + val[4 * rows + i] * x[i]
              + val[5 * rows + i] * x[i + 1]
              + val[6 * rows + i] * x[m + i - 1]
              + val[7 * rows + i] * x[m + i]
              + val[8 * rows + i] * x[m + i + 1];
          }

          r[m - 1] = val[2 * rows + m - 1] * x[0]
            + val[3 * rows + m - 1] * x[m - 2]
            + val[4 * rows + m - 1] * x[m - 1]
            + val[5 * rows + m - 1] * x[m]
            + val[6 * rows + m - 1] * x[2 * m - 2]
            + val[7 * rows + m - 1] * x[2 * m - 1]
            + val[8 * rows + m - 1] * x[2 * m];

          r[m] = val[1 * rows + m] * x[0]
            + val[2 * rows + m] * x[1]
            + val[3 * rows + m] * x[m - 1]
            + val[4 * rows + m] * x[m]
            + val[5 * rows + m] * x[m + 1]
            + val[6 * rows + m] * x[2 * m - 1]
            + val[7 * rows + m] * x[2 * m]
            + val[8 * rows + m] * x[2 * m + 1];

          for (Index i(1); i < rows - 2 * m - 1; ++i)
          {
          r[m + i] = val[m + i] * x[i - 1]
            + val[1 * rows + m + i] * x[i]
            + val[2 * rows + m + i] * x[i + 1]
            + val[3 * rows + m + i] * x[m + i - 1]
            + val[4 * rows + m + i] * x[m + i]
            + val[5 * rows + m + i] * x[m + i + 1]
            + val[6 * rows + m + i] * x[2 * m + i - 1]
            + val[7 * rows + m + i] * x[2 * m + i]
            + val[8 * rows + m + i] * x[2 * m + i + 1];
          }

          r[rows - m - 1] = val[rows - m - 1] * x[rows - 2 * m - 2]
            + val[2 * rows - m - 1] * x[rows - 2 * m - 1]
            + val[3 * rows - m - 1] * x[rows - 2 * m]
            + val[4 * rows - m - 1] * x[rows - m - 2]
            + val[5 * rows - m - 1] * x[rows - m - 1]
            + val[6 * rows - m - 1] * x[rows - m]
            + val[7 * rows - m - 1] * x[rows - 2]
            + val[8 * rows - m - 1] * x[rows - 1];

          r[rows - m] = val[rows - m] * x[rows - 2 * m - 1]
            + val[2 * rows - m] * x[rows - 2 * m]
            + val[3 * rows - m] * x[rows - 2 * m + 1]
            + val[4 * rows - m] * x[rows - m - 1]
            + val[5 * rows - m] * x[rows - m]
            + val[6 * rows - m] * x[rows - m + 1]
            + val[7 * rows - m] * x[rows - 1];

          for (Index i(1); i < m - 1; ++i)
          {
            r[rows - m + i] = val[rows - m + i] * x[rows - 2 * m + i - 1]
              + val[2 * rows - m + i] * x[rows - 2 * m + i]
              + val[3 * rows - m + i] * x[rows - 2 * m + i + 1]
              + val[4 * rows - m + i] * x[rows - m + i - 1]
              + val[5 * rows - m + i] * x[rows - m + i]
              + val[6 * rows - m + i] * x[rows - m + i + 1];
          }

          r[rows - 1] = val[rows - 1] * x[rows - m - 2]
            + val[2 * rows - 1] * x[rows - m - 1]
            + val[3 * rows - 1] * x[rows - m]
            + val[4 * rows - 1] * x[rows - 2]
            + val[5 * rows - 1] * x[rows - 1];
        }
      };

      extern template void ProductMatVec<Mem::Main, Algo::Generic>::csr(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::csr(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ProductMatVec<Mem::Main, Algo::Generic>::ell(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::ell(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      extern template void ProductMatVec<Mem::Main, Algo::Generic>::coo(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::coo(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::coo(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::coo(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);

      extern template void ProductMatVec<Mem::Main, Algo::Generic>::band_q2_d2(float *, const float * const, const float * const, const Index, const Index);
      extern template void ProductMatVec<Mem::Main, Algo::Generic>::band_q2_d2(double *, const double * const, const double * const, const Index, const Index);

      template <>
      struct ProductMatVec<Mem::Main, Algo::MKL>
      {
        static void csr(float * r, const float * const val, const Index * const col_ind, const Index * const row_ptr, const float * const x, const Index rows, const Index columns, const Index used_elements);
        static void csr(double * r, const double * const val, const Index * const col_ind, const Index * const row_ptr, const double * const x, const Index rows, const Index columns, const Index used_elements);

        static void coo(float * r, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const float * const x, const Index rows, const Index used_elements);
        static void coo(double * r, const double * const val, const Index * const row_ptr, const Index * const col_ptr, const double * const x, const Index rows, const Index used_elements);
      };

      template <>
      struct ProductMatVec<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements);

        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_PRODUCT_MATVEC_HPP
