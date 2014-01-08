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
      struct Defect
      {
      };

      template <>
      struct Defect<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const DT_ * const x, const Index rows);

        template <typename DT_>
        static void ell(DT_ * r, const DT_ * const rhs, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const DT_ * const x, const Index stride, const Index rows);

        template <typename DT_>
        static void coo(DT_ * r, const DT_ * const rhs, const DT_ * const val, const Index * const row_ptr, const Index * const col_ptr, const DT_ * const x, const Index rows, const Index used_elements);
      };

      template <>
      struct Defect<Mem::Main, Algo::MKL>
      {
        static void csr(float * r, const float * const rhs, const float * const val, const Index * const col_ind, const Index * const row_ptr, const float * const x, const Index rows);
        static void csr(double * r, const double * const rhs, const double * const val, const Index * const col_ind, const Index * const row_ptr, const double * const x, const Index rows);

        static void coo(float * r, const float * const rhs, const float * const val, const Index * const row_ptr, const Index * const col_ptr, const float * const x, const Index rows, const Index used_elements);
        static void coo(double * r, const double * const rhs, const double * const val, const Index * const row_ptr, const Index * const col_ptr, const double * const x, const Index rows, const Index used_elements);
      };

      template <>
      struct Defect<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const DT_ * const x, const Index rows);

        template <typename DT_>
        static void ell(DT_ * r, const DT_ * const rhs, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const DT_ * const x, const Index stride, const Index rows);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_DEFECT_HPP
