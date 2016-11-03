#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {

      template <typename DT_, typename IT_>
      void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(DT_* v, const DT_* b, const IT_* col_ind, const DT_* val, const IT_* row_ptr, const DT_ alpha, const Index size)
      {
        // loop over all scatter-matrix rows
        for (Index row(0) ; row < size ; ++row)
        {
          // skip empty rows
          if(row_ptr[row] == row_ptr[row + 1])
            continue;

          DT_ sum(0);
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            sum += val[i] * (b[col_ind[i]]);
          }
          v[row] += alpha * sum;
        }
      }

      template <typename DT_, typename IT_, int BlockSize_>
      void ScatterAxpyPrim<Mem::Main>::dvb_csr_generic(DT_* v, const DT_* b, const IT_* col_ind, const DT_* val, const IT_* row_ptr, const DT_ alpha, const Index size)
      {
        Tiny::Vector<DT_, BlockSize_> * x(reinterpret_cast<Tiny::Vector<DT_, BlockSize_> *>(v));
        const Tiny::Vector<DT_, BlockSize_> * bbuffer(reinterpret_cast<const Tiny::Vector<DT_, BlockSize_> *>(b));

        for (Index row(0) ; row < size ; ++row)
        {
          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          Tiny::Vector<DT_, BlockSize_> bsum(0);
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i(row_ptr[row]) ; i < end ; ++i)
          {
            bsum += val[i] * bbuffer[col_ind[i]];
          }
          x[row] += bsum * alpha;
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_GENERIC_HPP
