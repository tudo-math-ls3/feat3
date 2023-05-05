// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_HPP
#define KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct ProductMatMat
      {
        template <typename DT_>
        static void dense(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index columns, const Index inner)
        {
          dense_generic(r, alpha, beta, x, y, z, rows, columns, inner);
        }

#ifdef FEAT_HAVE_HALFMATH
        static void dense(Half * r, const Half alpha, const Half beta, const Half * const x, const Half * const y, const Half * const z, const Index rows, const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dense_cuda, dense_generic, dense_generic, r, alpha, beta, x, y, z, rows, columns, inner)
        }
#endif

        static void dense(float * r, const float alpha, const float beta, const float * const x, const float * const y, const float * const z, const Index rows, const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dense_cuda, dense_mkl, dense_generic, r, alpha, beta, x, y, z, rows, columns, inner)
        }

        static void dense(double * r, const double alpha, const double beta, const double * const x, const double * const y, const double * const z, const Index rows, const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dense_cuda, dense_mkl, dense_generic, r, alpha, beta, x, y, z, rows, columns, inner)
        }

        template <typename DT_, typename IT_>
        static void dsd(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index used_elements,
                                         const DT_ * y, const Index rows, const Index columns, const Index inner)
        {
          dsd_generic(r, alpha, beta, val, col_ind, row_ptr, used_elements, y, rows, columns, inner);
        }

#ifdef FEAT_HAVE_HALFMATH
        static void dsd(Half * r, const Half alpha, const Half beta, const Half * const val, const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index used_elements,
                                         const Half * y, const Index rows,  const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dsd_cuda, dsd_generic, dsd_generic, r, alpha, beta, val, col_ind, row_ptr, used_elements, y, rows, columns, inner)
        }
#endif

        static void dsd(float * r, const float alpha, const float beta, const float * const val, const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index used_elements,
                                         const float * y, const Index rows,  const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dsd_cuda, dsd_generic, dsd_generic, r, alpha, beta, val, col_ind, row_ptr, used_elements, y, rows, columns, inner)
        }

        static void dsd(double * r, const double alpha, const double beta, const double * const val, const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index used_elements,
                                         const double * y, const Index rows,  const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dsd_cuda, dsd_generic, dsd_generic, r, alpha, beta, val, col_ind, row_ptr, used_elements, y, rows, columns, inner)
        }

#ifdef FEAT_HAVE_HALFMATH
        static void dsd(Half * r, const Half alpha, const Half beta, const Half * const val, const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index used_elements,
                                         const Half * y, const Index rows,  const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dsd_cuda, dsd_generic, dsd_generic, r, alpha, beta, val, col_ind, row_ptr, used_elements, y, rows, columns, inner)
        }
#endif

        static void dsd(float * r, const float alpha, const float beta, const float * const val, const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index used_elements,
                                         const float * y, const Index rows,  const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dsd_cuda, dsd_generic, dsd_generic, r, alpha, beta, val, col_ind, row_ptr, used_elements, y, rows, columns, inner)
        }

        static void dsd(double * r, const double alpha, const double beta, const double * const val, const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index used_elements,
                                         const double * y, const Index rows,  const Index columns, const Index inner)
        {
          BACKEND_SKELETON_VOID(dsd_cuda, dsd_generic, dsd_generic, r, alpha, beta, val, col_ind, row_ptr, used_elements, y, rows, columns, inner)
        }

        template <typename DT_>
        static void dense_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index columns, const Index inner);

        static void dense_mkl(float * r, const float alpha, const float beta, const float * const x, const float * const y, const float * const z, const Index rows, const Index columns, const Index inner);
        static void dense_mkl(double * r, const double alpha, const double beta, const double * const x, const double * const y, const double * const z, const Index rows, const Index columns, const Index inner);

        template <typename DT_>
        static void dense_cuda(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index columns, const Index inner);

        template <typename DT_, typename IT_>
        static void dsd_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index used_elements,
                                         const DT_ * y, const Index rows,  const Index columns, const Index inner);

        template <typename DT_, typename IT_>
        static void dsd_cuda(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index used_elements,
                                         const DT_ * y, const Index rows,  const Index columns, const Index inner);

      };

#ifdef FEAT_EICKT
      extern template void ProductMatMat::dense_generic(float *, const float, const float, const float * const, const float * const, const float * const, const Index, const Index, const Index);
      extern template void ProductMatMat::dense_generic(double *, const double, const double, const double * const, const double * const, const double * const, const Index, const Index, const Index);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/product_matmat_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_HPP
