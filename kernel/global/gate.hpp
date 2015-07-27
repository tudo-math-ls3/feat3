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

      virtual void sync_0(LocalVector_& vector) const = 0;
      virtual void sync_1(LocalVector_& vector) const = 0;
      virtual DataType dot(const LocalVector_& x, const LocalVector_& y) const = 0;
      virtual DataType sum(DataType x) const = 0;
      virtual DataType norm2(DataType x) const = 0;
    };
  } // namespace Global
} // namespace FEAST

#endif // KERNEL_GLOBAL_GATE_HPP
