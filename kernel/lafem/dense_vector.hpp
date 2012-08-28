#pragma once
#ifndef KERNEL_LAFEM_DENSE_VECTOR_HPP
#define KERNEL_LAFEM_DENSE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    template <typename Arch_, typename DT_>
    class DenseVector : public Container<Arch_, DT_>
    {
      private:
        DT_ * _pelements;

      public:
        typedef DT_ data_type;
        typedef Arch_ arch_type;

        explicit DenseVector() :
          Container<Arch_, DT_> (0)
        {
        }

        explicit DenseVector(Index size) :
          Container<Arch_, DT_>(size)
        {
          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(size * sizeof(DT_)));
          this->_elements_size.push_back(size);
          this->_pelements = this->_elements.at(0);
        }

        explicit DenseVector(Index size, DT_ value) :
          Container<Arch_, DT_>(size)
        {
          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(size * sizeof(DT_)));
          this->_elements_size.push_back(size);
          this->_pelements = this->_elements.at(0);

          MemoryPool<Arch_>::instance()->set_memory(_pelements, value, size);
        }

        explicit DenseVector(Index size, DT_ * data) :
          Container<Arch_, DT_>(size)
        {
          this->_elements.push_back(data);
          this->_elements_size.push_back(size);
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

        DenseVector<Arch_, DT_> & operator= (const DenseVector<Arch_, DT_> & other)
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
          this->_elements_size.clear();
          this->_indices_size.clear();

          std::vector<DT_ *> new_elements = other.get_elements();
          std::vector<unsigned long *> new_indices = other.get_indices();

          this->_elements.assign(new_elements.begin(), new_elements.end());
          this->_indices.assign(new_indices.begin(), new_indices.end());
          this->_elements_size.assign(other.get_elements_size().begin(), other.get_elements_size().end());
          this->_indices_size.assign(other.get_indices_size().begin(), other.get_indices_size().end());

          _pelements = this->_elements.at(0);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

          return *this;
        }

        template <typename Arch2_, typename DT2_>
        DenseVector<Arch_, DT_> & operator= (const DenseVector<Arch2_, DT2_> & other)
        {
          this->_size = other.size();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();


          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(other.size() * sizeof(DT_)));
          this->_elements_size.push_back(this->_size);
          this->_pelements = this->_elements.at(0);

          Index src_size(other.get_elements_size().at(0) * sizeof(DT2_));
          Index dest_size(other.get_elements_size().at(0) * sizeof(DT_));
          void * temp(::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_elements().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_elements().at(0), temp, dest_size);
          ::free(temp);

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

        const DT_ operator()(Index index) const
        {
          ASSERT(index < this->_size, "Error: " + stringify(index) + " exceeds dense vector size " + stringify(this->_size) + " !");
          return MemoryPool<Arch_>::get_element(_pelements, index);
        }

        void operator()(Index index, DT_ value)
        {
          ASSERT(index < this->_size, "Error: " + stringify(index) + " exceeds dense vector size " + stringify(this->_size) + " !");
          MemoryPool<Arch_>::modify_element(_pelements, index, value);
        }
    };

    template <typename Arch_, typename Arch2_, typename DT_> bool operator== (const DenseVector<Arch_, DT_> & a, const DenseVector<Arch2_, DT_> & b)
    {
      if (a.size() != b.size())
        return false;
      if (a.get_elements().size() != b.get_elements().size())
        return false;
      if (a.get_indices().size() != b.get_indices().size())
        return false;

      for (Index i(0) ; i < a.size() ; ++i)
        if (a(i) != b(i))
          return false;

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

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
