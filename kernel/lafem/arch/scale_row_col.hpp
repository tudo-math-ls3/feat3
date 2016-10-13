#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
#define KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct ScaleRows;

      template <>
      struct ScaleRows<Mem::Main>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
          csr_generic(r, a, col_ind, row_ptr, x, rows, columns, used_elements);
        }

        template <typename DT_, typename IT_>
        static void coo(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
          coo_generic(r, a, col_ind, row_ptr, x, rows, columns, used_elements);
        }

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const rl, const DT_ * const x, const Index C, const Index rows)
        {
          ell_generic(r, a, col_ind, cs, cl, rl, x, C, rows);
        }

        template <typename DT_, typename IT_>
        static void csr_generic(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void coo_generic(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell_generic(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const /*rl*/, const DT_ * const x, const Index C, const Index rows);
      };

#ifdef FEAT_EICKT
      extern template void ScaleRows<Mem::Main>::csr_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::csr_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::csr_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::csr_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleRows<Mem::Main>::coo_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::coo_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::coo_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::coo_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleRows<Mem::Main>::ell_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::ell_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::ell_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void ScaleRows<Mem::Main>::ell_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);
#endif

      template <>
      struct ScaleRows<Mem::CUDA>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const /*rl*/, const DT_ * const x, const Index C, const Index rows);
      };


      // ***********************************************

      template <typename Mem_>
      struct ScaleCols;

      template <>
      struct ScaleCols<Mem::Main>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
          csr_generic(r, a, col_ind, row_ptr, x, rows, columns, used_elements);
        }

        template <typename DT_, typename IT_>
        static void coo(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
          coo_generic(r, a, col_ind, row_ptr, x, rows, columns, used_elements);
        }

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const rl, const DT_ * const x, const Index C, const Index rows)
        {
          ell_generic(r, a, col_ind, cs, cl, rl, x, C, rows);
        }

        template <typename DT_, typename IT_>
        static void csr_generic(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void coo_generic(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell_generic(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const /*rl*/, const DT_ * const x, const Index C, const Index rows);
      };

#ifdef FEAT_EICKT
      extern template void ScaleCols<Mem::Main>::csr_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::csr_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::csr_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::csr_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleCols<Mem::Main>::coo_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::coo_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::coo_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::coo_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void ScaleCols<Mem::Main>::ell_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::ell_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::ell_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void ScaleCols<Mem::Main>::ell_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);
#endif

      template <>
      struct ScaleCols<Mem::CUDA>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const a, const IT_ * const col_ind, const IT_ * const cs,
                        const IT_ * const cl, const IT_ * const /*rl*/, const DT_ * const x, const Index C, const Index rows);
      };
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/scale_row_col_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
