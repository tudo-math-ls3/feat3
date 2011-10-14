#pragma once
#ifndef KERNEL_HORNET_CONTAINER_HPP
#define KERNEL_HORNET_CONTAINER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/hornet/memory_arbiter.hpp>

#include <vector>



namespace FEAST
{
  template <typename Arch_, typename DT_>
  class Container
  {
    protected:
      Index _size;
      std::vector<DT_*> _elements;
      std::vector<Index*> _indices;

    public:
      Container(Index size) :
        _size(size)
      {
      }

      ~Container()
      {
        for (Index i(0) ; i < _elements.size() ; ++i)
          MemoryArbiter::instance()->release_memory(_elements.at(i));
        for (Index i(0) ; i < _indices.size() ; ++i)
          MemoryArbiter::instance()->release_memory(_indices.at(i));
      }

      Container(const Container<Arch_, DT_> & other) :
        _size(other._size),
        _elements(other._elements),
        _indices(other._indices)
      {
        for (Index i(0) ; i < _elements.size() ; ++i)
          MemoryArbiter::instance()->increase_memory(_elements.at(i));
        for (Index i(0) ; i < _indices.size() ; ++i)
          MemoryArbiter::instance()->increase_memory(_indices.at(i));
      }

      DT_ * get_elements(Index i = 0)
      {
        return _elements.at(i);
      }

      Index * get_indices(Index i = 0)
      {
        return _indices.at(i);
      }
  };
} // namespace FEAST

#endif // KERNEL_HORNET_CONTAINER_HPP
