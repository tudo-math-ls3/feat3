#pragma once
#ifndef KERNEL_LAFEM_AXPY_HPP
#define KERNEL_LAFEM_AXPY_HPP 1

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
    struct Axpy
    {
    };

    /**
     * \brief Axpy calculations.
     *
     * This class calculates axpy operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Axpy <Archs::CPU, Archs::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow ax + y\f$
       *
       * \param[out] r The axpy result.
       * \param[in] a A scalar to scale x with.
       * \param[in] x The vector to be scaled.
       * \param[in] y The other vector
       */
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

        if (rp == yp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] += a * xp[i];
          }
        }
        else if (rp == xp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= a;
            rp[i]+= yp[i];
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] = a * xp[i] + yp[i];
          }
        }
      }

      /**
       * \brief Calculate \f$r \leftarrow ax + y\f$
       *
       * \param[out] r The axpy result.
       * \param[in] a A vector to scale x with.
       * \param[in] x The vector to be scaled.
       * \param[in] y The other vector
       */
      template <typename DT_>
      static void value(DenseVector<Archs::CPU, DT_> & r, const DenseVector<Archs::CPU, DT_> & a, const DenseVector<Archs::CPU, DT_> & x, const DenseVector<Archs::CPU, DT_> & y)
      {
        if (x.size() != y.size())
          throw InternalError("Vector size does not match!");
        if (x.size() != r.size())
          throw InternalError("Vector size does not match!");
        if (a.size() != r.size())
          throw InternalError("Vector size does not match!");

        const DT_ * ap(a.elements());
        const DT_ * xp(x.elements());
        const DT_ * yp(y.elements());
        DT_ * rp(r.elements());
        const Index size(r.size());

        if (rp == yp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] += ap[i] * xp[i];
          }
        }
        else if(rp == xp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= ap[i];
            rp[i]+= yp[i];
          }
        }
        else if(rp == ap)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= xp[i];
            rp[i] += yp[i];
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] = ap[i] * xp[i] + yp[i];
          }
        }
      }
    };

    template <>
    struct Axpy <Archs::GPU, Archs::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Archs::GPU, DT_> & r, const DT_ a, const DenseVector<Archs::GPU, DT_> & x, const DenseVector<Archs::GPU, DT_> & y);

      template <typename DT_>
      static void value(DenseVector<Archs::GPU, DT_> & r, const DenseVector<Archs::GPU, DT_> & a, const DenseVector<Archs::GPU, DT_> & x, const DenseVector<Archs::GPU, DT_> & y);
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_AXPY_HPP
