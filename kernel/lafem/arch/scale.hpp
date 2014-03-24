#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCALE_HPP
#define KERNEL_LAFEM_ARCH_SCALE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_, typename Algo_>
      struct Scale;

      template <>
      struct Scale<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const DT_ s, const Index size)
        {
          if (x == r)
          {
            for (Index i(0) ; i < size ; ++i)
            {
              r[i] *= s;
            }
          }
          else
          {
            for (Index i(0) ; i < size ; ++i)
            {
              r[i] = x[i] * s;
            }
          }
        }
      };

      extern template void Scale<Mem::Main, Algo::Generic>::value(float *, const float * const, const float, const Index);
      extern template void Scale<Mem::Main, Algo::Generic>::value(double *, const double * const, const double, const Index);

      template <>
      struct Scale<Mem::Main, Algo::MKL>
      {
        static void value(float * r, const float * const x, const float, const Index size);
        static void value(double * r, const double * const x, const double, const Index size);
      };

      template <>
      struct Scale<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const DT_, const Index size);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_SCALE_HPP
