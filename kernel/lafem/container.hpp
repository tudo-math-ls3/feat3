#pragma once
#ifndef KERNEL_LAFEM_CONTAINER_HPP
#define KERNEL_LAFEM_CONTAINER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/lafem/memory_pool.hpp>
#include <kernel/archs.hpp>

#include <vector>
#include <limits>
#include <cmath>
#include <typeinfo>


namespace FEAST
{
  namespace LAFEM
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
        std::vector<Index> _elements_size;
        std::vector<Index> _indices_size;

      public:
        Container(Index size) :
          ContainerBase(size)
      {
      }

      ~Container()
      {
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
          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(other.get_elements_size().at(i) * sizeof(DT_)));
        for (Index i(0) ; i < other.get_indices().size() ; ++i)
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(other.get_indices_size().at(i) * sizeof(Index)));

        for (Index i(0) ; i < (other.get_elements()).size() ; ++i)
        {
          Index src_size(other.get_elements_size().at(i) * sizeof(DT2_));
          Index dest_size(other.get_elements_size().at(i) * sizeof(DT_));
          void * temp(::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_elements().at(i), src_size);
          MemoryPool<Arch_>::upload(this->get_elements().at(i), temp, dest_size);
          ::free(temp);
        }
        for (Index i(0) ; i < other.get_indices().size() ; ++i)
        {
          Index src_size(other.get_indices_size().at(i) * sizeof(DT2_));
          Index dest_size(other.get_indices_size().at(i) * sizeof(DT_));
          void * temp(::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_indices().at(i), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(i), temp, dest_size);
          ::free(temp);
        }

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

      std::vector<Index> & get_elements_size()
      {
        return _elements_size;
      }

      const std::vector<Index> & get_elements_size() const
      {
        return _elements_size;
      }

      std::vector<Index> & get_indices_size()
      {
        return _indices_size;
      }

      const std::vector<Index> & get_indices_size() const
      {
        return _indices_size;
      }

      const Index & size() const
      {
        return this->_size;
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_CONTAINER_HPP
