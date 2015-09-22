#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {

      template <typename DT_>
      void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(DT_* v, const DT_* b, const Index* col_ind, const DT_* val, const Index* row_ptr, const DT_ alpha, const Index size, const Index offset)
      {
        // loop over all scatter-matrix rows
        for (Index row(0) ; row < size ; ++row)
        {
          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          DT_ sum(DT_(0));
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            sum += DT_(val[i]) * DT_(b[offset + col_ind[i]]);
          }
          v[row] += alpha * sum;
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_SCATTER_AXPY_PRIM_GENERIC_HPP
