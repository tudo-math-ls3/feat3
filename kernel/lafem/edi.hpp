#pragma once
#ifndef KERNEL_LAFEM_EDI_HPP
#define KERNEL_LAFEM_EDI_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/archs.hpp>


namespace FEAT
{
  /**
   * \brief LAFEM namespace
   */
  namespace LAFEM
  {

    /**
     * \brief EDI - Enhanced Datatype Item.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This is a wrapper for modifiable data items
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_>
    class EDI
    {
    private:
      DT_ _value;
      DT_ * _address;
      //bool _armed;

    public:
      EDI(DT_ value, DT_* address) :
        _value(value),
        _address(address)
        //_armed(false)
      {
      }

      ~EDI()
      {
        //if (_armed)
        {
          MemoryPool<Mem_>::set_memory(_address, _value);
        }
      }

      EDI(const EDI<Mem_, DT_> & other) :
        _value(other._value),
        _address(other._address)
        //_armed(true)
      {
        //if (other._armed)
        //  throw InternalError(__func__, __FILE__, __LINE__, "You may not use the EDI copy constructor on your own!");
      }

      EDI<Mem_, DT_> & operator=(DT_ value)
      {
        _value = value;
        return *this;
      }

      EDI<Mem_, DT_> & operator+=(DT_ value)
      {
        _value += value;
        return *this;
      }

      EDI<Mem_, DT_> & operator-=(DT_ value)
      {
        _value -= value;
        return *this;
      }

      EDI<Mem_, DT_> & operator*=(DT_ value)
      {
        _value *= value;
        return *this;
      }

      EDI<Mem_, DT_> & operator/=(DT_ value)
      {
        _value /= value;
        return *this;
      }

      EDI<Mem_, DT_> & operator%=(DT_ value)
      {
        _value %= value;
        return *this;
      }
    };
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_EDI_HPP
