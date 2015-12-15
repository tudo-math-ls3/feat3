#pragma once
#ifndef KERNEL_GLOBAL_GATE_HPP
#define KERNEL_GLOBAL_GATE_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>

namespace FEAST
{
  /**
   * \brief Global linear algebra namespace
   */
  namespace Global
  {
    /**
     * \brief Gate interface class template
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_>
    class Gate
    {
    public:
      typedef typename LocalVector_::MemType MemType;
      typedef typename LocalVector_::DataType DataType;
      typedef typename LocalVector_::IndexType IndexType;

    public:
      virtual ~Gate()
      {
      }

      /**
       * \brief Converts a type-1 vector into a type-0 vector.
       *
       * \param[inout] vector
       * On entry, the type-1 vector to be converted.
       * On exit, the converted type-0 vector.
       *
       * \note This function does not perform any synchronisation.
       */
      virtual void from_1_to_0(LocalVector_& vector) const = 0;

      /**
       * \brief Synchronises a type-0 vector, resulting in a type-1 vector.
       *
       * \param[inout] vector
       * On entry, the type-0 vector to be synchronised.\n
       * On exit, the synchronised type-1 vector.
       */
      virtual void sync_0(LocalVector_& vector) const = 0;

      /**
       * \brief Synchronises a type-1 vector, resulting in a type-1 vector.
       *
       * \param[inout] vector
       * On entry, the type-1 vector to be synchronised.\n
       * On exit, the synchronised type-1 vector.
       *
       * \note
       * This function effectively applies the from_1_to_0() and sync_0()
       * functions onto the input vector.
       */
      virtual void sync_1(LocalVector_& vector) const = 0;

      /**
       * \brief Computes a synchronised dot-product of two type-1 vectors.
       *
       * \param[in] x, y
       * The two type-1 vector whose dot-product is to be computed.
       *
       * \returns
       * The dot-product of \p x and \p y.
       */
      virtual DataType dot(const LocalVector_& x, const LocalVector_& y) const = 0;

      /**
       * \brief Computes a reduced sum over all processes.
       *
       * \param[in] x
       * The value that is to be summarised over all processes.
       *
       * \returns
       * The reduced sum of all \p x.
       */
      virtual DataType sum(DataType x) const = 0;

      /**
       * \brief Computes a reduced 2-norm over all processes.
       *
       * This function is equivalent to the call
       *    Math::sqrt(this->sum(x*x))
       *
       * \param[in] x
       * The value that is to be summarised over all processes.
       *
       * \returns
       * The reduced 2-norm of all \p x.
       */
      virtual DataType norm2(DataType x) const = 0;
    };
  } // namespace Global
} // namespace FEAST

#endif // KERNEL_GLOBAL_GATE_HPP
