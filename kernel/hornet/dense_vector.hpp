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
    private:
      DT_ * _pelements;

    public:
      DenseVector(Index size) :
        Container<Arch_, DT_>(size)
      {
        this->_size = size;
        this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(size * sizeof(DT_)));
        this->_pelements = this->_elements.at(0);
      }

      DenseVector(Index size, DT_ value) :
        Container<Arch_, DT_>(size)
      {
        this->_size = size;
        this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(size * sizeof(DT_)));
        this->_pelements = this->_elements.at(0);

        //TODO add arch, use memory arbiter set memory function
        for (Index i(0) ; i < this->_size ; ++i)
        {
          _pelements[i] = value;
        }
      }

      DenseVector(const DenseVector<Arch_, DT_> & other) :
        Container<Arch_, DT_>(other)
      {
        this->_pelements = this->_elements.at(0);
      }

      template <typename Arch2_, typename DT2_>
      DenseVector(const DenseVector<Arch2_, DT2_> & other) :
        Container<Arch_, DT_>(other)
      {
        this->_pelements = this->_elements.at(0);
      }

      DenseVector<Arch_, DT_> & operator= (DenseVector<Arch_, DT_> & other)
      {
        this->_size = other.size();

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

        this->_elements.clear();
        this->_indices.clear();

        std::vector<DT_ *> new_elements = other.get_elements();
        std::vector<unsigned long *> new_indices = other.get_indices();

        this->_elements.assign(new_elements.begin(), new_elements.end());
        this->_indices.assign(new_indices.begin(), new_indices.end());

        _pelements = this->_elements.at(0);

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

        return *this;
      }

      template <typename Arch2_, typename DT2_>
      DenseVector<Arch_, DT_> & operator= (DenseVector<Arch2_, DT2_> & other)
      {
        //TODO copy memory from arch2 to arch
        return *this;
      }

      const DT_ & operator()(Index index) const
      {
        return _pelements[index];
      }

      void operator()(Index index, DT_ value)
      {
        _pelements[index] = value;
      }

      Index size()
      {
        return this->_size;
      }
  };
} // namespace FEAST

#endif // KERNEL_HORNET_DENSE_VECTOR_HPP
