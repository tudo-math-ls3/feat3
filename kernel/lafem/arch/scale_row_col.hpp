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

      template <typename Mem_, typename Algo_>
      struct ScaleCols;

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
      };

      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index, const Index);

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
      };

      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index, const Index);

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST


#endif // KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
