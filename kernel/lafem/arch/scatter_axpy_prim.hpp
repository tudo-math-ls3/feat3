#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_HPP
#define KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_HPP 1

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
      struct ScatterAxpyPrim;

      template <>
      struct ScatterAxpyPrim<Mem::Main>
      {
        template <typename DT_, typename IT_>
        static void dv_csr(DT_* v,
                           const DT_* b,
                           const IT_* col_ind,
                           const DT_* val,
                           const IT_* row_ptr,
                           const DT_ alpha,
                           const Index size)
        {
          dv_csr_generic(v, b, col_ind, val, row_ptr, alpha, size);
        }

        template <typename DT_, typename IT_>
        static void dv_csr_generic(DT_* v,
                                   const DT_* b,
                                   const IT_* col_ind,
                                   const DT_* val,
                                   const IT_* row_ptr,
                                   const DT_ alpha,
                                   const Index size);

      };

      template <>
      struct ScatterAxpyPrim<Mem::CUDA>
      {
        template <typename DT_, typename IT_>
        static void dv_csr(DT_* v,
                           const DT_* b,
                           const IT_* col_ind,
                           const DT_* val,
                           const IT_* row_ptr,
                           const DT_ alpha,
                           const Index size);
      };


    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/scatter_axpy_prim_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_HPP
