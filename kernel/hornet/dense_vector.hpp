#pragma once
#ifndef KERNEL_HORNET_DENSE_VECTOR_HPP
#define KERNEL_HORNET_DENSE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/hornet/container.hpp>
#include <kernel/hornet/absolute.hpp>


namespace FEAST
{
  template <typename Arch_, typename DT_>
  class DenseVector : public Container<Arch_, DT_>
  {
    private:
      DT_ * _pelements;

    public:
      typedef DT_ data_type;

      explicit DenseVector() :
        Container<Arch_, DT_> (0)
      {
      }

      explicit DenseVector(Index size) :
        Container<Arch_, DT_>(size)
      {
        this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(size * sizeof(DT_)));
        this->_pelements = this->_elements.at(0);
      }

      explicit DenseVector(Index size, DT_ value) :
        Container<Arch_, DT_>(size)
      {
        this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(size * sizeof(DT_)));
        this->_pelements = this->_elements.at(0);

        //TODO add arch, use memory arbiter set memory function
        for (Index i(0) ; i < this->_size ; ++i)
        {
          _pelements[i] = value;
        }
      }

      explicit DenseVector(Index size, DT_ * data) :
        Container<Arch_, DT_>(size)
      {
        this->_elements.push_back(data);
        this->_pelements = this->_elements.at(0);

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));
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
        if (this == &other)
          return *this;

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
        if (this == &other)
          return *this;

        //TODO copy memory from arch2 to arch
        return *this;
      }

      DT_ * elements()
      {
        return _pelements;
      }

      const DT_ * elements() const
      {
        return _pelements;
      }

      const DT_ & operator()(Index index) const
      {
        ASSERT(index < this->_size, "Error: " + stringify(index) + "exceeds dense vector size " + stringify(this->_size) + " !");
        return _pelements[index];
      }

      void operator()(Index index, DT_ value)
      {
        ASSERT(index < this->_size, "Error: " + stringify(index) + "exceeds dense vector size " + stringify(this->_size) + " !");
        _pelements[index] = value;
      }
  };

  template <typename Arch_, typename DT_> bool operator== (const DenseVector<Arch_, DT_> & a, const DenseVector<Arch_, DT_> & b)
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
          if (Absolute<DT_>::value(a.get_elements().at(i)[j] - b.get_elements().at(i)[j]) > std::numeric_limits<DT_>::epsilon())
            return false;
      }
    }

    return true;
  }


  template <typename Arch_, typename DT_>
  std::ostream &
  operator<< (std::ostream & lhs, const DenseVector<Arch_, DT_> & b)
  {
    lhs << "[";
    for (Index i(0) ; i < b.size() ; ++i)
    {
      lhs << "  " << b(i);
    }
    lhs << "]" << std::endl;

    return lhs;
  }

} // namespace FEAST

#endif // KERNEL_HORNET_DENSE_VECTOR_HPP
