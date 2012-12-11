#pragma once
#ifndef KERNEL_LAFEM_ELEMENT_PRODUCT_HPP
#define KERNEL_LAFEM_ELEMENT_PRODUCT_HPP 1

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
    struct ComponentProduct
    {
    };

    /**
     * \brief Elementwise product calculations.
     *
     * This class calculates elementwise product operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct ComponentProduct<Algo::Generic>
    {
      /**
       * \brief Calculate \f$r_i \leftarrow x_i \cdot y_i\f$
       *
       * \param[out] r The elementwise product result.
       * \param[in] x.The first factor.
       * \param[in] y The second factor.
       */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const DenseVector<Mem::Main, DT_> & x, const DenseVector<Mem::Main, DT_> & y)
      {
        if (x.size() != y.size())
          throw InternalError("Vector size does not match!");
        if (x.size() != r.size())
          throw InternalError("Vector size does not match!");

        const DT_ * xp(x.elements());
        const DT_ * yp(y.elements());
        DT_ * rp(r.elements());
        const Index size(r.size());

        if (rp == xp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= yp[i];
          }
        }
        else if (rp == yp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= xp[i];
          }
        }
        else if (rp == xp && rp == yp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= rp[i];
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] = xp[i] * yp[i];
          }
        }
      }
    };

    template <>
    struct ComponentProduct<Algo::MKL>
    {
      static void value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const DenseVector<Mem::Main, float> & y);
      static void value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const DenseVector<Mem::Main, double> & y);
    };

    template <>
    struct ComponentProduct<Algo::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & x, const DenseVector<Mem::CUDA, DT_> & y);
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ELEMENT_PRODUCT_HPP
