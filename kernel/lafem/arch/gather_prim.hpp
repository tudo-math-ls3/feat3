#pragma once
#ifndef KERNEL_LAFEM_ARCH_GATHER_PRIM_HPP
#define KERNEL_LAFEM_ARCH_GATHER_PRIM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct GatherPrim;

      template <>
      struct GatherPrim<Mem::Main>
      {
        template <typename DT_>
        static void dv_csr(DT_* b,
                           const DT_* v,
                           const Index* col_ind,
                           const DT_* val,
                           const Index* row_ptr,
                           const Index size,
                           const Index offset)
        {
          dv_csr_generic(b, v, col_ind, val, row_ptr, size, offset);
        }

        template <typename DT_>
        static void dv_csr_generic(DT_* b,
                                   const DT_* v,
                                   const Index* col_ind,
                                   const DT_* val,
                                   const Index* row_ptr,
                                   const Index size,
                                   const Index offset);

      };

      template <>
      struct GatherPrim<Mem::CUDA>
      {
        template <typename DT_>
        static void dv_csr(DT_* b,
                           const DT_* v,
                           const Index* col_ind,
                           const DT_* val,
                           const Index* row_ptr,
                           const Index size,
                           const Index offset);
      };


    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/gather_prim_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_GATHER_PRIM_HPP
