#pragma once
#ifndef KERNEL_LAFEM_DOT_PRODUCT_HPP
#define KERNEL_LAFEM_DOT_PRODUCT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct DotProduct
    {
    };

    /**
     * \brief Inner product calculations.
     *
     * This class calculates dot product operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct DotProduct<Algo::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow x \cdot y\f$
       *
       * \param[out] r The dot product result.
       * \param[in] x.The first vector.
       * \param[in] y The second vector.
       */
      template <typename DT_>
      static DT_ value(const DenseVector<Mem::Main, DT_> & x, const DenseVector<Mem::Main, DT_> & y)
      {
        if (x.size() != y.size())
          throw InternalError("Vector size does not match!");

        const DT_ * xp(x.elements());
        const DT_ * yp(y.elements());
        DT_ r(0);
        const Index size(x.size());

        if(xp == yp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r += xp[i] * xp[i];
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r += xp[i] * yp[i];
          }
        }

        return r;
      }
    };

    template <>
    struct DotProduct<Algo::MKL>
    {
      static float value(const DenseVector<Mem::Main, float> & x, const DenseVector<Mem::Main, float> & y);
      static double value(const DenseVector<Mem::Main, double> & x, const DenseVector<Mem::Main, double> & y);
    };

    template <>
    struct DotProduct<Algo::CUDA>
    {
      template <typename DT_>
      static DT_ value(const DenseVector<Mem::CUDA, DT_> & x, const DenseVector<Mem::CUDA, DT_> & y);
    };


  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DOT_PRODUCT_HPP
