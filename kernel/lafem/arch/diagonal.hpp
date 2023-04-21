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

        static void csr(std::uint64_t * diag, const std::uint64_t * const col_ind, const std::uint64_t * const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, diag, col_ind, row_ptr, rows)
        }

        static void csr(std::uint32_t * diag, const std::uint32_t * const col_ind, const std::uint32_t * const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda, csr_generic, csr_generic, diag, col_ind, row_ptr, rows)
        }

        template <typename IT_>
        static void csr_generic(IT_ * diag, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);

        template <typename IT_>
        static void csr_cuda(IT_ * diag, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);

      };

#ifdef FEAT_EICKT
      extern template void Diagonal::csr_generic(std::uint64_t *, const std::uint64_t * const, const std::uint64_t * const, const Index);
      extern template void Diagonal::csr_generic(std::uint32_t *, const std::uint32_t * const, const std::uint32_t * const, const Index);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/diagonal_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_DIAGONAL_HPP
