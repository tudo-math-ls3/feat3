#pragma once
#ifndef KERNEL_HORNET_DENSE_VECTOR_HPP
#define KERNEL_HORNET_DENSE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/hornet/container.hpp>



namespace FEAST
{
  template <typename Arch_, typename DT_>
  class DenseVector : public Container<Arch_, DT_>
  {
    public:
      /// \todo replace std::vector by memory_arbiter calls
      DenseVector(Index size) :
        Container<Arch_, DT_>(size)
      {
        this->_size = size;
        this->_elements.push_back((DT_*)MemoryArbiter::instance()->allocate_memory(size * sizeof(DT_)));
      }

      DenseVector(Index size, DT_ value) :
        Container<Arch_, DT_>(size)
      {
        this->_size = size;
        this->_elements.push_back((DT_*)MemoryArbiter::instance()->allocate_memory(size * sizeof(DT_)));

        DT_ * elements(this->_elements.at(0));
        for (Index i(0) ; i < this->_size ; ++i)
        {
          elements[i] = value;
        }
      }

      DenseVector(const DenseVector<Arch_, DT_> & other) :
        Container<Arch_, DT_>(other)
      {
      }

      DT_ & operator[](Index index)
      {
        return (this->_elements.at(0))[index];
      }

      Index size()
      {
        return this->_size;
      }
  };
} // namespace FEAST

#endif // KERNEL_HORNET_DENSE_VECTOR_HPP
