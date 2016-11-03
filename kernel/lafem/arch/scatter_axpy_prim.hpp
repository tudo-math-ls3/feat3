#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_HPP
#define KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_HPP 1

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

        template <typename DT_, typename IT_, int BlockSize_>
        static void dvb_csr(DT_* v,
                           const DT_* b,
                           const IT_* col_ind,
                           const DT_* val,
                           const IT_* row_ptr,
                           const DT_ alpha,
                           const Index size)
        {
          dvb_csr_generic<DT_, IT_, BlockSize_>(v, b, col_ind, val, row_ptr, alpha, size);
        }

        template <typename DT_, typename IT_, int BlockSize_>
        static void dvb_csr_generic(DT_* v,
                                   const DT_* b,
                                   const IT_* col_ind,
                                   const DT_* val,
                                   const IT_* row_ptr,
                                   const DT_ alpha,
                                   const Index size);

      };

#ifdef FEAT_EICKT
      extern template void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(float*, const float*, const unsigned long*, const float*, const unsigned long*, const float alpha, const Index);
      extern template void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(double*, const double*, const unsigned long*, const double*, const unsigned long*, const double alpha, const Index);
      extern template void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(float*, const float*, const unsigned int*, const float*, const unsigned int*, const float alpha, const Index);
      extern template void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(double*, const double*, const unsigned int*, const double*, const unsigned int*, const double alpha, const Index);
#endif

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

        template <typename DT_, typename IT_, int BlockSize_>
        static void dvb_csr(DT_* v,
                           const DT_* b,
                           const IT_* col_ind,
                           const DT_* val,
                           const IT_* row_ptr,
                           const DT_ alpha,
                           const Index size);
      };


    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/scatter_axpy_prim_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_HPP
