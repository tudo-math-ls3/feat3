#pragma once
#ifndef KERNEL_LAFEM_ALGORITHM_HPP
#define KERNEL_LAFEM_ALGORITHM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Arch_, typename DT_>
    void copy(DenseVector<Arch_, DT_> & dest, const DenseVector<Arch_, DT_> & src)
    {
      if (dest.size() != src.size())
        throw InternalError("Vector size mismatch!");

      const Index size(dest.size());
      DT_ * pdest(dest.elements());
      const DT_ * psrc(src.elements());

      for (Index i(0) ; i < size ; ++i)
      {
        pdest[i] = psrc[i];
      }
    }

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ALGORITHM_HPP
