#pragma once
#ifndef KERNEL_LAFEM_SCALE_HPP
#define KERNEL_LAFEM_SCALE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Arch_, typename BType_>
    struct Scale
    {
    };

    template <>
    struct Scale<Archs::CPU, Archs::Generic>
    {
      template <typename DT_>
      static void value(DenseVector<Archs::CPU, DT_> & r, const DenseVector<Archs::CPU, DT_> & x, const DT_ s)
      {
        if (x.size() != r.size())
          throw InternalError("Vector size does not match!");

        const DT_ * xp(x.elements());
        DT_ * rp(r.elements());
        const Index size(r.size());

        for (Index i(0) ; i < size ; ++i)
        {
          rp[i] = xp[i] * s;
        }
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SCALE_HPP
