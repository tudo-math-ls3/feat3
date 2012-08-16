#pragma once
#ifndef KERNEL_HORNET_AXPY_HPP
#define KERNEL_HORNET_AXPY_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/hornet/dense_vector.hpp>



namespace FEAST
{
  template <typename Arch_, typename BType_>
  struct Axpy
  {
  };

  template <>
  struct Axpy <Archs::CPU, Archs::Generic>
  {
    template <typename DT_>
    static void value(DenseVector<Archs::CPU, DT_> & r, const DT_ a, const DenseVector<Archs::CPU, DT_> & x, const DenseVector<Archs::CPU, DT_> & y)
    {
      if (x.size() != y.size())
        throw InternalError("Vector size does not match!");
      if (x.size() != r.size())
        throw InternalError("Vector size does not match!");

      const DT_ * xp(x.elements());
      const DT_ * yp(y.elements());
      DT_ * rp(r.elements());
      const Index size(r.size());

      for (Index i(0) ; i < size ; ++i)
      {
        rp[i] = a * xp[i] + yp[i];
      }
    }

    template <typename DT_>
    static void value(DenseVector<Archs::CPU, DT_> & r, const DenseVector<Archs::CPU, DT_> & a, const DenseVector<Archs::CPU, DT_> & x, const DenseVector<Archs::CPU, DT_> & y)
    {
      if (x.size() != y.size())
        throw InternalError("Vector size does not match!");
      if (x.size() != r.size())
        throw InternalError("Vector size does not match!");

      const DT_ * ap(a.elements());
      const DT_ * xp(x.elements());
      const DT_ * yp(y.elements());
      DT_ * rp(r.elements());
      const Index size(r.size());

      for (Index i(0) ; i < size ; ++i)
      {
        rp[i] = ap[i] * xp[i] + yp[i];
      }
    }
  };

} // namespace FEAST

#endif // KERNEL_HORNET_AXPY_HPP
