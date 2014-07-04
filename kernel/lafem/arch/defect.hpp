#pragma once
#ifndef KERNEL_LAFEM_ARCH_DEFECT_HPP
#define KERNEL_LAFEM_ARCH_DEFECT_HPP 1

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
      struct Defect;

      template <>
      struct Defect<Mem::Main, Algo::Generic>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_, Index BlockHeight_, Index BlockWidth_>
        static void csrb(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const rhs, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows);

        template <typename DT_, typename IT_>
        static void coo(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const row_ptr, const IT_ * const col_ptr, const DT_ * const x, const Index rows, const Index used_elements);

        template <typename DT_, typename IT_>
        static void banded(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns);
      };

      extern template void Defect<Mem::Main, Algo::Generic>::csr(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::csr(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::csr(double *, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::csr(double *, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void Defect<Mem::Main, Algo::Generic>::ell(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::ell(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::ell(double *,const double * const,  const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::ell(double *,const double * const,  const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      extern template void Defect<Mem::Main, Algo::Generic>::coo(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::coo(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::coo(double *, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::coo(double *, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      extern template void Defect<Mem::Main, Algo::Generic>::banded(float *, const float * const, const float * const, const Index * const, const float * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main, Algo::Generic>::banded(double *, const double * const, const double * const, const Index * const, const double * const, const Index, const Index, const Index);

      template <>
        struct Defect<Mem::Main, Algo::MKL>
      {
        static void csr(float * r, const float * const rhs, const float * const val, const Index * const col_ind, const Index * const row_ptr, const float * const x, const Index rows, const Index, const Index);
        static void csr(double * r, const double * const rhs, const double * const val, const Index * const col_ind, const Index * const row_ptr, const double * const x, const Index rows, const Index, const Index);

        static void coo(float * r, const float * const rhs, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const float * const x, const Index rows, const Index used_elements);
        static void coo(double * r, const double * const rhs, const double * const val, const Index * const row_ptr, const Index * const col_ptr, const double * const x, const Index rows, const Index used_elements);
      };

      template <>
      struct Defect<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const rhs, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows);

        template <typename DT_, typename IT_>
        static void banded(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/defect_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_DEFECT_HPP
