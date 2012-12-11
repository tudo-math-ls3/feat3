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
    template <typename Algo_>
    struct Scale
    {
    };

    /**
     * \brief Scaling calculations.
     *
     * This class calculates scaling operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Scale<Algo::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow x \cdot s\f$
       *
       * \param[out] r The scaled vector.
       * \param[in] x The vector to be scaled.
       * \param[in] s A scalar to scale x with.
       */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const DenseVector<Mem::Main, DT_> & x, const DT_ s)
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

    template <>
    struct Scale<Algo::MKL>
    {
      static void value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const float s);
      static void value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const double s);
    };

    template <>
    struct Scale<Algo::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & x, const DT_ s);
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SCALE_HPP
