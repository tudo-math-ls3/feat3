// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_DIAGONAL_HPP
#define KERNEL_LAFEM_ARCH_DIAGONAL_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct Diagonal
      {
        template <typename IT_>
        static void csr(IT_ * diag, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
        {
          csr_generic(diag, col_ind, row_ptr, rows);
        }

        static void csr(unsigned long * diag, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, diag, col_ind, row_ptr, rows)
        }

        static void csr(unsigned int * diag, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, diag, col_ind, row_ptr, rows)
        }

        template <typename IT_>
        static void csr_generic(IT_ * diag, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);

        template <typename IT_>
        static void csr_cuda(IT_ * diag, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);

      };

#ifdef FEAT_EICKT
      extern template void Diagonal::csr_generic(unsigned long *, const unsigned long * const, const unsigned long * const, const Index);
      extern template void Diagonal::csr_generic(unsigned int *, const unsigned int * const, const unsigned int * const, const Index);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/diagonal_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_DIAGONAL_HPP
