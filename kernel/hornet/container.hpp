#pragma once
#ifndef KERNEL_HORNET_CONTAINER_HPP
#define KERNEL_HORNET_CONTAINER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/hornet/memory_pool.hpp>
#include <kernel/archs.hpp>

#include <vector>
#include <limits>
#include <cmath>



namespace FEAST
{
  class ContainerBase
  {
    protected:
      Index _size;

      ContainerBase(Index size) :
        _size(size)
      {
      }

  };

  template <typename Arch_, typename DT_>
  class Container : public ContainerBase
  {
    protected:
      std::vector<DT_*> _elements;
      std::vector<Index*> _indices;

    public:
      Container(Index size) :
        ContainerBase(size)
      {
      }

      ~Container()
      {
        //TODO add arch
        for (Index i(0) ; i < _elements.size() ; ++i)
          MemoryPool<Arch_>::instance()->release_memory(_elements.at(i));
        for (Index i(0) ; i < _indices.size() ; ++i)
          MemoryPool<Arch_>::instance()->release_memory(_indices.at(i));
      }

      Container(const Container<Arch_, DT_> & other) :
        ContainerBase(other),
        _elements(other._elements),
        _indices(other._indices)
      {
        for (Index i(0) ; i < _elements.size() ; ++i)
          MemoryPool<Arch_>::instance()->increase_memory(_elements.at(i));
        for (Index i(0) ; i < _indices.size() ; ++i)
          MemoryPool<Arch_>::instance()->increase_memory(_indices.at(i));
      }

      template <typename Arch2_, typename DT2_>
      Container(const Container<Arch2_, DT2_> & other) :
        ContainerBase(other)
      {
        if (typeid(DT_) != typeid(DT2_))
            throw InternalError("type conversion not supported yet!");

        for (Index i(0) ; i < (other.get_elements()).size() ; ++i)
          MemoryPool<Arch_>::instance()->allocate_memory(this->_size * sizeof(DT_));
        for (Index i(0) ; i < other.get_indices().size() ; ++i)
          MemoryPool<Arch_>::instance()->allocate_memory(this->_size * sizeof(Index));

        //TODO copy memory from arch2 to arch

      }

      std::vector<DT_*> & get_elements()
      {
        return _elements;
      }

      const std::vector<DT_*> & get_elements() const
      {
        return _elements;
      }

      std::vector<Index*> & get_indices()
      {
        return _indices;
      }

      const std::vector<Index*> & get_indices() const
      {
        return _indices;
      }

      const Index & size() const
      {
        return this->_size;
      }
  };

  template <typename Arch_, typename DT_> bool operator== (const Container<Arch_, DT_> & a, const Container<Arch_, DT_> & b)
  {
    if (a.size() != b.size())
      return false;
    if (a.get_elements().size() != b.get_elements().size())
      return false;
    if (a.get_indices().size() != b.get_indices().size())
      return false;

    for (Index i(0) ; i < a.get_indices().size() ; ++i)
    {
      for (Index j(0) ; j < a.size() ; ++j)
        if (a.get_indices().at(i)[j] != b.get_indices().at(i)[j])
          return false;
    }

    if (typeid(DT_) == typeid(Index))
    {
      for (Index i(0) ; i < a.get_elements().size() ; ++i)
      {
        for (Index j(0) ; j < a.size() ; ++j)
          if (a.get_elements().at(i)[j] != b.get_elements().at(i)[j])
        return false;
      }
    }
    else
    {


      for (Index i(0) ; i < a.get_elements().size() ; ++i)
      {
        for (Index j(0) ; j < a.size() ; ++j)
          if (fabs(a.get_elements().at(i)[j] - b.get_elements().at(i)[j]) > std::numeric_limits<DT_>::epsilon())
            return false;
      }
    }

    return true;
  }

} // namespace FEAST

#endif // KERNEL_HORNET_CONTAINER_HPP
