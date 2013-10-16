#pragma once
#ifndef KERNEL_LAFEM_NORM_HPP
#define KERNEL_LAFEM_NORM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>

#include <cmath>


namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct Norm2
    {
    };

    /**
     * \brief L2 norm calculations.
     *
     * This class calculates euclidean length (L2 norm) operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Norm2<Algo::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow \sqrt{ \sum\limits_i x_i \cdot x_i } \f$
       *
       * \param[in] x The input vector.
       *
       * \returns The L2 norm.
       */
      template <typename DT_>
      static DT_ value(const DenseVector<Mem::Main, DT_> & x)
      {
        const DT_ * xp(x.elements());
        DT_ r(0);
        const Index size(x.size());

        DT_ xpv;
        for (Index i(0) ; i < size ; ++i)
        {
          xpv = xp[i];
          r += xpv * xpv;
        }

        return (DT_)std::sqrt(r);
      }
    };

    template <>
    struct Norm2<Algo::MKL>
    {
      static float value(const DenseVector<Mem::Main, float> & x);
      static double value(const DenseVector<Mem::Main, double> & x);
    };

    template <>
    struct Norm2<Algo::CUDA>
    {
      template <typename DT_>
      static DT_ value(const DenseVector<Mem::CUDA, DT_> & x);
    };

    /**
     * \brief L2 norm calculations.
     *
     * This class calculates euclidean length (L2 norm) operations without the final square root.
     *
     * \author Dirk Ribbrock
     */
    template <typename Algo_>
    struct Norm2wosqrt
    {
      /**
       * \brief Calculate \f$r \leftarrow \sqrt{ \sum\limits_i x_i \cdot x_i } \f$
       *
       * \param[in] x The input vector.
       *
       * \returns The L2 norm.
       */
      template <typename Mem_, typename DT_>
      static DT_ value(const DenseVector<Mem_, DT_> & x)
      {
        DT_ result(Norm2<Algo_>::value(x));
        return result * result;
      }
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_NORM2_HPP
