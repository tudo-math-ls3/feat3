// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_LUMPING_HPP
#define KERNEL_LAFEM_ARCH_LUMPING_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct Lumping
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * lump, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
        {
          csr_generic(lump, val, col_ind, row_ptr, rows);
        }

        static void csr(float * lump, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, lump, val, col_ind, row_ptr, rows)
        }

        static void csr(double * lump, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, lump, val, col_ind, row_ptr, rows)
        }

        static void csr(float * lump, const float * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, lump, val, col_ind, row_ptr, rows)
        }

        static void csr(double * lump, const double * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, lump, val, col_ind, row_ptr, rows)
        }

        template <typename DT_, typename IT_>
        static void csr_generic(DT_ * lump, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);

        template <typename DT_, typename IT_>
        static void bcsr(DT_ * lump, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          bcsr_generic(lump, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);
        }

        static void bcsr(float * lump, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda, bcsr_generic, bcsr_generic, lump, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr(double * lump, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda, bcsr_generic, bcsr_generic, lump, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr(float * lump, const float * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda, bcsr_generic, bcsr_generic, lump, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr(double * lump, const double * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda, bcsr_generic, bcsr_generic, lump, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        template <typename DT_, typename IT_>
        static void bcsr_generic(DT_* lump, const DT_* const val, const IT_* const col_ind, const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);

        template <typename DT_, typename IT_>
        static void csr_cuda(DT_ * lump, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);

        template <typename DT_, typename IT_>
        static void bcsr_cuda(DT_ * lump, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);
      };

#ifdef FEAT_EICKT
      extern template void Lumping::csr_generic(float *, const float * const, const Index * const, const Index * const, const Index);
      extern template void Lumping::csr_generic(double *, const double * const, const Index * const, const Index * const, const Index);

      extern template void Lumping::bcsr_generic(float *, const float * const, const Index * const, const Index * const, const Index, const int, const int);
      extern template void Lumping::bcsr_generic(double *, const double * const, const Index * const, const Index * const, const Index, const int, const int);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/lumping_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_LUMPING_HPP
