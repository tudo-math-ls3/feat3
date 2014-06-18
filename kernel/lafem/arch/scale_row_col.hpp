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
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void coo(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const /*rl*/, const DT_ * const x, const Index C, const Index rows);
};

      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleRows<Mem::Main, Algo::Generic>::coo(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::coo(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::coo(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::coo(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleRows<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      template <>
      struct ScaleRows<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const /*rl*/, const DT_ * const x, const Index C, const Index rows);
      };


      // ***********************************************

      template <typename Mem_, typename Algo_>
      struct ScaleCols;

      template <>
      struct ScaleCols<Mem::Main, Algo::Generic>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void coo(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const /*rl*/, const DT_ * const x, const Index C, const Index rows);
};

      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleCols<Mem::Main, Algo::Generic>::coo(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::coo(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::coo(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::coo(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleCols<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main, Algo::Generic>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      template <>
      struct ScaleCols<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const /*rl*/, const DT_ * const x, const Index C, const Index rows);
};

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/scale_row_col_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
