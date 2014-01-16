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
      struct Axpy
      {
      };

      template <>
      struct Axpy<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static void dv(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static void dv(DT_ * r, const DT_ * a, const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static void csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows);

        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows);

        template <typename DT_>
        static void ell(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const Index stride, const Index rows);

        template <typename DT_>
        static void ell(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const Index stride, const Index rows);

        template <typename DT_>
        static void coo(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements);

        template <typename DT_>
        static void coo(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const Index * const row_ptr, const Index * const col_ptr, const Index rows, const Index used_elements);
      };

      template <>
      struct Axpy<Mem::Main, Algo::MKL>
      {
        static void dv(float * r, const float a, const float * const x, const float * const y, const Index size);
        static void dv(double * r, const double a, const double * const x, const double * const y, const Index size);

        static void dv(float * r, const float * const a, const float * const x, const float * const y, const Index size);
        static void dv(double * r, const double * const a, const double * const x, const double * const y, const Index size);

        static void csr(float * r, const float a, const float * const x, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows);
        static void csr(double * r, const double a, const double * const x, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows);

        static void csr(float * r, const float * const a, const float * const x, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows);
        static void csr(double * r, const double * const a, const double * const x, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows);

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
        static void csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows);

        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows);

        template <typename DT_>
        static void ell(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const Index stride, const Index rows);

        template <typename DT_>
        static void ell(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const Index stride, const Index rows);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_AXPY_HPP
