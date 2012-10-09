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

    /**
     * \brief Scaling calculations.
     *
     * This class calculates scaling operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Scale<Archs::CPU, Archs::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow x \cdot s\f$
       *
       * \param[out] r The scaled vector.
       * \param[in] x The vector to be scaled.
       * \param[in] s A scalar to scale x with.
       */
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

    template <>
    struct Scale<Archs::GPU, Archs::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Archs::GPU, DT_> & r, const DenseVector<Archs::GPU, DT_> & x, const DT_ s);
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SCALE_HPP
