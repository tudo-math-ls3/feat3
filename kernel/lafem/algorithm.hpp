#pragma once
#ifndef KERNEL_LAFEM_ALGORITHM_HPP
#define KERNEL_LAFEM_ALGORITHM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    template <typename Mem_, typename DT_>
    void copy(DenseVector<Mem_, DT_> & dest, const DenseVector<Mem_, DT_> & src)
    {
      if (dest.size() != src.size())
        throw InternalError("Vector size mismatch!");

      DT_ * pdest(dest.elements());
      const DT_ * psrc(src.elements());

      MemoryPool<Mem_>::template copy<DT_>(pdest, psrc, dest.size());
    }

    template <typename Mem_, typename Arch2_, typename DT_>
    void copy(DenseVector<Mem_, DT_> & dest, const DenseVector<Arch2_, DT_> & src)
    {
      if (dest.size() != src.size())
        throw InternalError("Vector size mismatch!");

      DenseVector<Mem_, DT_> temp(src);
      copy(dest, temp);
    }
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ALGORITHM_HPP
