#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
#define KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP 1

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
      struct ScaleRows;

      template <>
      struct ScaleRows<Mem::Main, Algo::Generic>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
        {
          for (Index row(0) ; row < rows ; ++row)
          {
            const IT_ end(row_ptr[row + 1]);
            for (IT_ i(row_ptr[row]) ; i < end ; ++i)
            {
              r[i] = a[i] * x[row];
            }
          }
        }

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const Ax, const IT_ * const /*Aj*/, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
        {
          for (Index row(0) ; row < rows ; ++row)
          {
            const DT_ * tAx(Ax);
            DT_ * tr(r);
            tAx += row;
            tr += row;

            const IT_ max(Arl[row]);
            for(IT_ n(0); n < max ; n++)
            {
              *tr = *Ax * x[row];

              tAx += stride;
              tr += stride;
            }
          }
        }
      };

      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleRows<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      template <>
      struct ScaleRows<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows);
      };


      // ***********************************************

      template <typename Mem_, typename Algo_>
      struct ScaleCols;

      template <>
      struct ScaleCols<Mem::Main, Algo::Generic>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
        {
          for (Index row(0) ; row < rows ; ++row)
          {
            const IT_ end(row_ptr[row + 1]);
            for (IT_ i(row_ptr[row]) ; i < end ; ++i)
            {
              r[i] = a[i] * x[col_ind[i]];
            }
          }
        }

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
        {
          for (Index row(0) ; row < rows ; ++row)
          {
            const IT_ * tAj(Aj);
            const DT_ * tAx(Ax);
            DT_ * tr(r);
            tAj += row;
            tAx += row;
            tr += row;

            const IT_ max(Arl[row]);
            for(IT_ n(0); n < max ; n++)
            {
              *tr = *Ax * x[*tAj];

              tAj += stride;
              tAx += stride;
              tr += stride;
            }
          }
        }
      };

      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleCols<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      template <>
      struct ScaleCols<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST


#endif // KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
