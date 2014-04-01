#pragma once
#ifndef KERNEL_LAFEM_ARCH_AXPY_HPP
#define KERNEL_LAFEM_ARCH_AXPY_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_, typename Algo_>
      struct Axpy;

      template <>
      struct Axpy<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static void dv(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
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
        static void dv(DT_ * r, const DT_  * const a, const DT_ * const x, const DT_ * const y, const Index size)
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
        static void csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
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
        static void csr(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
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
        static void ell(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const Index stride, const Index rows)
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
        static void ell(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const Index stride, const Index rows)
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
        static void coo(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
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
        static void coo(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val,
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

        template <typename DT_>
        static void banded_q1_d2(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const Index nodes_per_row, const Index nodes_per_column)
        {
          const Index rows(nodes_per_row * nodes_per_column);
          const Index m(nodes_per_row);

          r[0] = y[0] + a * (val[4 * rows] * x[0]
                             + val[5 * rows] * x[1]
                             + val[6 * rows] * x[m - 1]
                             + val[7 * rows] * x[m]
                             + val[8 * rows] * x[m + 1]);

          for (Index i(1); i < m - 1; ++i)
          {
            r[i] = y[i] + a * (val[3 * rows + i] * x[i - 1]
                               + val[4 * rows + i] * x[i]
                               + val[5 * rows + i] * x[i + 1]
                               + val[6 * rows + i] * x[m + i - 1]
                               + val[7 * rows + i] * x[m + i]
                               + val[8 * rows + i] * x[m + i + 1]);
          }

          r[m - 1] = y[m - 1] + a * (val[2 * rows + m - 1] * x[0]
                                     + val[3 * rows + m - 1] * x[m - 2]
                                     + val[4 * rows + m - 1] * x[m - 1]
                                     + val[5 * rows + m - 1] * x[m]
                                     + val[6 * rows + m - 1] * x[2 * m - 2]
                                     + val[7 * rows + m - 1] * x[2 * m - 1]
                                     + val[8 * rows + m - 1] * x[2 * m]);

          r[m] = y[m] + a * (val[1 * rows + m] * x[0]
                             + val[2 * rows + m] * x[1]
                             + val[3 * rows + m] * x[m - 1]
                             + val[4 * rows + m] * x[m]
                             + val[5 * rows + m] * x[m + 1]
                             + val[6 * rows + m] * x[2 * m - 1]
                             + val[7 * rows + m] * x[2 * m]
                             + val[8 * rows + m] * x[2 * m + 1]);

          for (Index i(1); i < rows - 2 * m - 1; ++i)
          {
            r[m + i] = y[m + i] + a * (val[m + i] * x[i - 1]
                                       + val[1 * rows + m + i] * x[i]
                                       + val[2 * rows + m + i] * x[i + 1]
                                       + val[3 * rows + m + i] * x[m + i - 1]
                                       + val[4 * rows + m + i] * x[m + i]
                                       + val[5 * rows + m + i] * x[m + i + 1]
                                       + val[6 * rows + m + i] * x[2 * m + i - 1]
                                       + val[7 * rows + m + i] * x[2 * m + i]
                                       + val[8 * rows + m + i] * x[2 * m + i + 1]);
          }

          r[rows - m - 1] = y[rows - m - 1] + a * (val[rows - m - 1] * x[rows - 2 * m - 2]
                                                   + val[2 * rows - m - 1] * x[rows - 2 * m - 1]
                                                   + val[3 * rows - m - 1] * x[rows - 2 * m]
                                                   + val[4 * rows - m - 1] * x[rows - m - 2]
                                                   + val[5 * rows - m - 1] * x[rows - m - 1]
                                                   + val[6 * rows - m - 1] * x[rows - m]
                                                   + val[7 * rows - m - 1] * x[rows - 2]
                                                   + val[8 * rows - m - 1] * x[rows - 1]);

          r[rows - m] = y[rows - m] + a * (val[rows - m] * x[rows - 2 * m - 1]
                                           + val[2 * rows - m] * x[rows - 2 * m]
                                           + val[3 * rows - m] * x[rows - 2 * m + 1]
                                           + val[4 * rows - m] * x[rows - m - 1]
                                           + val[5 * rows - m] * x[rows - m]
                                           + val[6 * rows - m] * x[rows - m + 1]
                                           + val[7 * rows - m] * x[rows - 1]);

          for (Index i(1); i < m - 1; ++i)
          {
            r[rows - m + i] = y[rows - m + i] + a * (val[rows - m + i] * x[rows - 2 * m + i - 1]
                                                     + val[2 * rows - m + i] * x[rows - 2 * m + i]
                                                     + val[3 * rows - m + i] * x[rows - 2 * m + i + 1]
                                                     + val[4 * rows - m + i] * x[rows - m + i - 1]
                                                     + val[5 * rows - m + i] * x[rows - m + i]
                                                     + val[6 * rows - m + i] * x[rows - m + i + 1]);
          }

          r[rows - 1] = y[rows - 1] + a * (val[rows - 1] * x[rows - m - 2]
                                           + val[2 * rows - 1] * x[rows - m - 1]
                                           + val[3 * rows - 1] * x[rows - m]
                                           + val[4 * rows - 1] * x[rows - 2]
                                           + val[5 * rows - 1] * x[rows - 1]);
        }

        template <typename DT_>
        static void banded_q1_d2(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const Index nodes_per_row, const Index nodes_per_column)
        {
          const Index rows(nodes_per_row * nodes_per_column);
          const Index m(nodes_per_row);

          r[0] = y[0] + a[0] * (val[4 * rows] * x[0]
                                + val[5 * rows] * x[1]
                                + val[6 * rows] * x[m - 1]
                                + val[7 * rows] * x[m]
                                + val[8 * rows] * x[m + 1]);

          for (Index i(1); i < m - 1; ++i)
          {
            r[i] = y[i] + a[i] * (val[3 * rows + i] * x[i - 1]
                                  + val[4 * rows + i] * x[i]
                                  + val[5 * rows + i] * x[i + 1]
                                  + val[6 * rows + i] * x[m + i - 1]
                                  + val[7 * rows + i] * x[m + i]
                                  + val[8 * rows + i] * x[m + i + 1]);
          }

          r[m - 1] = y[m - 1] + a[m - 1] * (val[2 * rows + m - 1] * x[0]
                                            + val[3 * rows + m - 1] * x[m - 2]
                                            + val[4 * rows + m - 1] * x[m - 1]
                                            + val[5 * rows + m - 1] * x[m]
                                            + val[6 * rows + m - 1] * x[2 * m - 2]
                                            + val[7 * rows + m - 1] * x[2 * m - 1]
                                            + val[8 * rows + m - 1] * x[2 * m]);

          r[m] = y[m] + a[m] * (val[1 * rows + m] * x[0]
                                + val[2 * rows + m] * x[1]
                                + val[3 * rows + m] * x[m - 1]
                                + val[4 * rows + m] * x[m]
                                + val[5 * rows + m] * x[m + 1]
                                + val[6 * rows + m] * x[2 * m - 1]
                                + val[7 * rows + m] * x[2 * m]
                                + val[8 * rows + m] * x[2 * m + 1]);

          for (Index i(1); i < rows - 2 * m - 1; ++i)
          {
            r[m + i] = y[m + i] + a[m + i] * (val[m + i] * x[i - 1]
                                              + val[1 * rows + m + i] * x[i]
                                              + val[2 * rows + m + i] * x[i + 1]
                                              + val[3 * rows + m + i] * x[m + i - 1]
                                              + val[4 * rows + m + i] * x[m + i]
                                              + val[5 * rows + m + i] * x[m + i + 1]
                                              + val[6 * rows + m + i] * x[2 * m + i - 1]
                                              + val[7 * rows + m + i] * x[2 * m + i]
                                              + val[8 * rows + m + i] * x[2 * m + i + 1]);
          }

          r[rows - m - 1] = y[rows - m - 1] + a[rows - m - 1] * (val[rows - m - 1] * x[rows - 2 * m - 2]
                                                                 + val[2 * rows - m - 1] * x[rows - 2 * m - 1]
                                                                 + val[3 * rows - m - 1] * x[rows - 2 * m]
                                                                 + val[4 * rows - m - 1] * x[rows - m - 2]
                                                                 + val[5 * rows - m - 1] * x[rows - m - 1]
                                                                 + val[6 * rows - m - 1] * x[rows - m]
                                                                 + val[7 * rows - m - 1] * x[rows - 2]
                                                                 + val[8 * rows - m - 1] * x[rows - 1]);

          r[rows - m] = y[rows - m] + a[rows - m] * (val[rows - m] * x[rows - 2 * m - 1]
                                                     + val[2 * rows - m] * x[rows - 2 * m]
                                                     + val[3 * rows - m] * x[rows - 2 * m + 1]
                                                     + val[4 * rows - m] * x[rows - m - 1]
                                                     + val[5 * rows - m] * x[rows - m]
                                                     + val[6 * rows - m] * x[rows - m + 1]
                                                     + val[7 * rows - m] * x[rows - 1]);

          for (Index i(1); i < m - 1; ++i)
          {
            r[rows - m + i] = y[rows - m + i] + a[rows - m + i] * (val[rows - m + i] * x[rows - 2 * m + i - 1]
                                                                   + val[2 * rows - m + i] * x[rows - 2 * m + i]
                                                                   + val[3 * rows - m + i] * x[rows - 2 * m + i + 1]
                                                                   + val[4 * rows - m + i] * x[rows - m + i - 1]
                                                                   + val[5 * rows - m + i] * x[rows - m + i]
                                                                   + val[6 * rows - m + i] * x[rows - m + i + 1]);
          }

          r[rows - 1] = y[rows - 1] + a[rows - 1] * (val[rows - 1] * x[rows - m - 2]
                                                     + val[2 * rows - 1] * x[rows - m - 1]
                                                     + val[3 * rows - 1] * x[rows - m]
                                                     + val[4 * rows - 1] * x[rows - 2]
                                                     + val[5 * rows - 1] * x[rows - 1]);
        }
      };

      extern template void Axpy<Mem::Main, Algo::Generic>::dv(float *, const float, const float * const, const float * const, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::dv(double *, const double, const double * const, const double * const, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::dv(float *, const float * const, const float * const, const float * const, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::dv(double *, const double * const, const double * const, const double * const, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::csr(float *, const float, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::csr(float *, const float, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::csr(double *, const double, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::csr(double *, const double, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::csr(float *, const float * const, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::csr(float *, const float * const, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::csr(double *, const double * const, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::csr(double *, const double * const, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::ell(float *, const float, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::ell(float *, const float, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::ell(double *, const double, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::ell(double *, const double, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::ell(float *, const float * const, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::ell(float *, const float * const, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::ell(double *, const double * const, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::ell(double *, const double * const, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::coo(float *, const float, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::coo(float *, const float, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::coo(double *, const double, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::coo(double *, const double, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::coo(float *, const float * const, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::coo(float *, const float * const, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::coo(double *, const double * const, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::coo(double *, const double * const, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::banded_q1_d2(float *, const float, const float * const, const float * const, const float * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::banded_q1_d2(double *, const double, const double * const, const double * const, const double * const, const Index, const Index);

      extern template void Axpy<Mem::Main, Algo::Generic>::banded_q1_d2(float *, const float * const, const float * const, const float * const, const float * const, const Index, const Index);
      extern template void Axpy<Mem::Main, Algo::Generic>::banded_q1_d2(double *, const double * const, const double * const, const double * const, const double * const, const Index, const Index);

      template <>
      struct Axpy<Mem::Main, Algo::MKL>
      {
        static void dv(float * r, const float a, const float * const x, const float * const y, const Index size);
        static void dv(double * r, const double a, const double * const x, const double * const y, const Index size);

        static void dv(float * r, const float * const a, const float * const x, const float * const y, const Index size);
        static void dv(double * r, const double * const a, const double * const x, const double * const y, const Index size);

        static void csr(float * r, const float a, const float * const x, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index, const Index);
        static void csr(double * r, const double a, const double * const x, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index, const Index);

        static void csr(float * r, const float * const a, const float * const x, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index, const Index);
        static void csr(double * r, const double * const a, const double * const x, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index, const Index);

        static void coo(float * r, const float a, const float * const x, const float * const y, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements);
        static void coo(double * r, const double a, const double * const x, const double * const y, const double * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements);

        static void coo(float * r, const float * const a, const float * const x, const float * const y, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements);
        static void coo(double * r, const double * const a, const double * const x, const double * const y, const double * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements);
      };

      template <>
      struct Axpy<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static void dv(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static void dv(DT_ * r, const DT_ * a, const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static void csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows, const Index columns, const Index used_elements);
        template <typename DT_>
        static void csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows, const Index columns, const Index used_elements);

        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows, const Index columns, const Index used_elements);
        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows, const Index columns, const Index used_elements);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const Index stride, const Index rows);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const Index stride, const Index rows);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_AXPY_HPP
