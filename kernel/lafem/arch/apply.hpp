// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_APPLY_HPP
#define KERNEL_LAFEM_ARCH_APPLY_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/backend.hpp>
#include <kernel/lafem/arch/product_matmat.hpp>
#include <kernel/util/half.hpp>

#include <typeinfo>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct Apply
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                        const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns,
                        const Index used_elements, const bool transposed)
        {
          csr_generic(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, transposed);
        }

#ifdef FEAT_HAVE_HALFMATH
        static void csr(Half * r, const Half a, const Half * const x, const Half b, const Half * const y, const Half * const val,
                        const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index rows, const Index columns,
                        const Index used_elements, const bool transposed)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, transposed)
        }
#endif

        static void csr(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val,
                        const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index rows, const Index columns,
                        const Index used_elements, const bool transposed)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_mkl, csr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, transposed)
        }

        static void csr(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val,
                        const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index rows, const Index columns,
                        const Index used_elements, const bool transposed)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_mkl, csr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, transposed)
        }

#ifdef FEAT_HAVE_HALFMATH
        static void csr(Half * r, const Half a, const Half * const x, const Half b, const Half * const y, const Half * const val,
                        const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index rows, const Index columns,
                        const Index used_elements, const bool transposed)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, transposed)
        }
#endif

        static void csr(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val,
                        const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index rows, const Index columns,
                        const Index used_elements, const bool transposed)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, transposed)
        }

        static void csr(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val,
                        const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index rows, const Index columns,
                        const Index used_elements, const bool transposed)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, transposed)
        }

        template <typename DT_, typename IT_>
        static void cscr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                        const IT_ * const col_ind, const IT_ * const row_ptr, const IT_ * const row_numbers, const Index used_rows, const Index rows, const Index columns,
                        const Index used_elements, const bool transposed)
        {
          cscr_generic(r, a, x, b, y, val, col_ind, row_ptr, row_numbers, used_rows, rows, columns, used_elements, transposed);
        }

        template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
        static void bcsr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                         const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns,
                         const Index used_elements)
        {
          bcsr_generic<BlockHeight_, BlockWidth_, DT_, IT_>(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements);
        }

        template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
        static void bcsr_transposed(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
          const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns,
          const Index used_elements)
        {
          bcsr_transposed_generic<BlockHeight_, BlockWidth_, DT_, IT_>(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements);
        }

        template <int BlockHeight_, int BlockWidth_>
        static void bcsr(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val,
                         const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index rows, const Index columns,
                         const Index used_elements)
        {
          if (BlockHeight_ == BlockWidth_)
            BACKEND_SKELETON_VOID_T2(BlockHeight_, BlockWidth_, bcsr_cuda, bcsr_mkl, bcsr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements)
          else
            BACKEND_SKELETON_VOID_T2(BlockHeight_, BlockWidth_, bcsr_cuda, bcsr_generic, bcsr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements)
        }

        template <int BlockHeight_, int BlockWidth_>
        static void bcsr(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val,
                         const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index rows, const Index columns,
                         const Index used_elements)
        {
          if (BlockHeight_ == BlockWidth_)
            BACKEND_SKELETON_VOID_T2(BlockHeight_, BlockWidth_, bcsr_cuda, bcsr_mkl, bcsr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements)
          else
            BACKEND_SKELETON_VOID_T2(BlockHeight_, BlockWidth_, bcsr_cuda, bcsr_generic, bcsr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements)
        }

        template <int BlockHeight_, int BlockWidth_>
        static void bcsr(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val,
                         const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index rows, const Index columns,
                         const Index used_elements)
        {
          BACKEND_SKELETON_VOID_T2(BlockHeight_, BlockWidth_, bcsr_cuda, bcsr_generic, bcsr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements)
        }

        template <int BlockHeight_, int BlockWidth_>
        static void bcsr(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val,
                         const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index rows, const Index columns,
                         const Index used_elements)
        {
          BACKEND_SKELETON_VOID_T2(BlockHeight_, BlockWidth_, bcsr_cuda, bcsr_generic, bcsr_generic, r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements)
        }

        template <int BlockSize_, typename DT_, typename IT_>
        static void csrsb(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements)
        {
          csrsb_generic<BlockSize_, DT_, IT_>(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements);
        }

        template <typename DT_, typename IT_>
        static void banded(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y, const DT_ * const val, const IT_ * const offsets,  const Index num_of_offsets, const Index rows, const Index columns)
        {
          banded_generic(r, alpha, x, beta, y, val, offsets, num_of_offsets, rows, columns);
        }

        template <typename DT_, typename IT_>
        static void banded_transposed(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y, const DT_ * const val, const IT_ * const offsets,  const Index num_of_offsets, const Index rows, const Index columns)
        {
          banded_transposed_generic(r, alpha, x, beta, y, val, offsets, num_of_offsets, rows, columns);
        }

        static void banded(float * r, const float alpha, const float * const x, const float beta, const float * const y, const float * const val, const std::uint64_t * const offsets,  const Index num_of_offsets, const Index rows, const Index columns)
        {
          BACKEND_SKELETON_VOID(banded_cuda, banded_generic, banded_generic, r, alpha, x, beta, y, val, offsets, num_of_offsets, rows, columns)
        }

        static void banded(double * r, const double alpha, const double * const x, const double beta, const double * const y, const double * const val, const std::uint64_t * const offsets,  const Index num_of_offsets, const Index rows, const Index columns)
        {
          BACKEND_SKELETON_VOID(banded_cuda, banded_generic, banded_generic, r, alpha, x, beta, y, val, offsets, num_of_offsets, rows, columns)
        }

        static void banded(float * r, const float alpha, const float * const x, const float beta, const float * const y, const float * const val, const std::uint32_t * const offsets,  const Index num_of_offsets, const Index rows, const Index columns)
        {
          BACKEND_SKELETON_VOID(banded_cuda, banded_generic, banded_generic, r, alpha, x, beta, y, val, offsets, num_of_offsets, rows, columns)
        }

        static void banded(double * r, const double alpha, const double * const x, const double beta, const double * const y, const double * const val, const std::uint32_t * const offsets,  const Index num_of_offsets, const Index rows, const Index columns)
        {
          BACKEND_SKELETON_VOID(banded_cuda, banded_generic, banded_generic, r, alpha, x, beta, y, val, offsets, num_of_offsets, rows, columns)
        }

        template <typename DT_>
        static void dense(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
        {
          dense_generic(r, alpha, beta, y, val, x, rows, columns);
        }

        template <typename DT_>
        static void dense_transposed(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
        {
          dense_transposed_generic(r, alpha, beta, y, val, x, rows, columns);
        }

#ifdef FEAT_HAVE_HALFMATH
        static void dense(Half * r, const Half alpha, const Half beta, const Half * const y, const Half * const val, const Half * const x, const Index rows, const Index columns)
        {
          switch(Backend::get_preferred_backend())
          {
            //no cuda half implementation exists, thus we use the gemm version
            case PreferredBackend::cuda:
              //ProductMatMat::dense_cuda(r, alpha, beta, x, val, y, rows, 1, columns);
              ProductMatMat::dense_cuda(r, alpha, beta, val, x, y, rows, 1, columns);
              break;
            case PreferredBackend::generic:
            default:
              dense_generic(r, alpha, beta, y, val, x, rows, columns);
          }
        }
#endif

        static void dense(float * r, const float alpha, const float beta, const float * const y, const float * const val, const float * const x, const Index rows, const Index columns)
        {
          BACKEND_SKELETON_VOID(dense_cuda, dense_mkl, dense_generic, r, alpha, beta, y, val, x, rows, columns)
        }

        static void dense(double * r, const double alpha, const double beta, const double * const y, const double * const val, const double * const x, const Index rows, const Index columns)
        {
          BACKEND_SKELETON_VOID(dense_cuda, dense_mkl, dense_generic, r, alpha, beta, y, val, x, rows, columns)
        }


        template <typename DT_, typename IT_>
        static void csr_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                        const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index, const bool);

        template <typename DT_, typename IT_>
        static void cscr_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                        const IT_ * const col_ind, const IT_ * const row_ptr, const IT_ * const row_numbers, const Index used_rows,
                        const Index rows, const Index, const Index, const bool);

        template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
        static void bcsr_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
                         const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index);

        template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
        static void bcsr_transposed_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val,
          const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index);

        template <int BlockSize_, typename DT_, typename IT_>
        static void csrsb_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void banded_generic(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y, const DT_ * const val, const IT_ * const offsets,  const Index num_of_offsets, const Index rows, const Index columns);

        template <typename DT_, typename IT_>
        static void banded_transposed_generic(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y, const DT_ * const val, const IT_ * const offsets,  const Index num_of_offsets, const Index rows, const Index columns);

        template <typename DT_>
        static void dense_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const rhs, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns);

        template <typename DT_>
        static void dense_transposed_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const rhs, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns);

        static void csr_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool);
        static void csr_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const bool);

        template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
        static void bcsr_mkl(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements)
        {
          XASSERTM(BlockHeight_ == BlockWidth_, "MKL only supports square blocks!");
          bcsr_mkl(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, BlockHeight_);
        }

        static void bcsr_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize);
        static void bcsr_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize);

        template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
        static void bcsr_transposed_mkl(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements)
        {
          XASSERTM(BlockHeight_ == BlockWidth_, "MKL only supports square blocks!");
          bcsr_mkl(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, BlockHeight_);
        }

        static void bcsr_transposed_mkl(float * r, const float a, const float * const x, const float b, const float * const y, const float * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize);
        static void bcsr_transposed_mkl(double * r, const double a, const double * const x, const double b, const double * const y, const double * const val, const Index * const col_ind, const Index * const row_ptr, const Index rows, const Index columns, const Index, const int blocksize);

        static void dense_mkl(float * r, const float alpha, const float beta, const float * const y, const float * const val, const float * const x, const Index rows, const Index columns);
        static void dense_mkl(double * r, const double alpha, const double beta, const double * const y, const double * const val, const double * const x, const Index rows, const Index columns);

        template <typename DT_, typename IT_>
        static void csr_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const bool transposed);

        template <int BlockHeight_, int BlockWidth_, typename DT_, typename IT_>
        static void bcsr_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements)
        {
          XASSERTM(BlockHeight_ < 10, "The generic cuda bcsr kernel does not support BlockHeight greather than 9!");
          bcsr_wrapper_cuda(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, BlockHeight_, BlockWidth_);
        }

        template <typename DT_, typename IT_>
        static void bcsr_wrapper_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const int BlockHeight, const int BlockWidth);

        template <typename DT_, typename IT_>
        static void bcsr_intern_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const int BlockSize);

        template <typename DT_, typename IT_>
        static void bcsr_intern_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const int BlockHeight, const int BlockWidth);

        template <int BlockSize_, typename DT_, typename IT_>
        static void csrsb_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements);

        template <typename DT_, typename IT_>
        static void banded_cuda(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y, const DT_ * const val, const IT_ * const offsets, const Index num_of_offsets, const Index rows, const Index columns);

        template <typename DT_>
        static void dense_cuda(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns);
      };

#ifdef FEAT_EICKT
      extern template void Apply::csr_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
      extern template void Apply::csr_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
      extern template void Apply::csr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
      extern template void Apply::csr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);

      extern template void Apply::cscr_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const Index, const bool);
      extern template void Apply::cscr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const Index, const bool);
      extern template void Apply::cscr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const Index, const bool);
      extern template void Apply::cscr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const Index, const bool);

      extern template void Apply::banded_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint64_t * const, const Index, const Index, const Index);
      extern template void Apply::banded_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint32_t * const, const Index, const Index, const Index);
      extern template void Apply::banded_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint64_t * const, const Index, const Index, const Index);
      extern template void Apply::banded_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const Index, const Index, const Index);

      extern template void Apply::dense_generic(float *, const float, const float, const float * const, const float * const, const float * const, const Index, const Index);
      extern template void Apply::dense_generic(double *, const double, const double, const double * const, const double * const, const double * const, const Index, const Index);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/apply_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_APPLY_HPP
