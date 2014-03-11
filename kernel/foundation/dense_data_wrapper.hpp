#pragma once
#ifndef KERNEL_FOUNDATION_DENSE_DATA_WRAPPER_HH
#define KERNEL_FOUNDATION_DENSE_DATA_WRAPPER_HH 1

#include <kernel/foundation/base.hpp>
#include <iterator>

namespace FEAST
{
  namespace Foundation
  {
    template<typename DT_>
    class DDVIterator : public std::iterator<std::random_access_iterator_tag, DT_>
    {
      public:
        DDVIterator(typename std::iterator<std::random_access_iterator_tag, DT_>::pointer x) :
          _iter(x)
        {
        }

        DDVIterator(const DDVIterator& i) :
          _iter(i._iter)
        {
        }

        DDVIterator() :
          _iter(nullptr)
        {
        }

        DDVIterator& operator++()
        {
          ++_iter;
          return *this;
        }

        DDVIterator operator++(int)
        {
          DDVIterator tmp(*this);
          operator++();
          return tmp;
        }

        DDVIterator& operator--()
        {
          --_iter;
          return *this;
        }

        DDVIterator operator--(int)
        {
          DDVIterator tmp(*this);
          operator--();
          return tmp;
        }

        bool operator==(const DDVIterator& rhs)
        {
          return _iter == rhs._iter;
        }

        bool operator!=(const DDVIterator& rhs)
        {
          return _iter != rhs._iter;
        }

        bool operator>(const DDVIterator& rhs)
        {
          return _iter > rhs._iter;
        }

        bool operator<(const DDVIterator& rhs)
        {
          return _iter < rhs._iter;
        }

        bool operator>=(const DDVIterator& rhs)
        {
          return _iter >= rhs._iter;
        }

        bool operator<=(const DDVIterator& rhs)
        {
          return _iter <= rhs._iter;
        }

        DDVIterator& operator=(const DDVIterator& r)
        {
          this->_iter = r._iter;
          return *this;
        }

        DDVIterator& operator+=(const DDVIterator& r)
        {
          for(Index i(0) ; i < r ; ++i)
            operator++();

          return *this;
        }

        DDVIterator& operator-=(const DDVIterator& r)
        {
          for(Index i(0) ; i < r ; ++i)
            operator--();

          return *this;
        }

        DT_ operator[](const Index i)
        {
          return _iter[i];
        }

        DT_ operator[](const Index i) const
        {
          return _iter[i];
        }

        DT_* operator->()
        {
          return _iter;
        }

        DT_& operator*()
        {
          return *_iter;
        }

        DT_* get()
        {
          return _iter;
        }

        DT_* get() const
        {
          return _iter;
        }

      private:

        DT_* _iter;
    };

    template<typename DT_>
    typename DDVIterator<DT_>::difference_type operator-(const DDVIterator<DT_>& l, const DDVIterator<DT_>& r)
    {
      return typename DDVIterator<DT_>::difference_type(l.get() - r.get());
    }

    template<typename DT_, typename IT_>
    DDVIterator<DT_> operator-(const DDVIterator<DT_>& l, const IT_& r)
    {
      return DDVIterator<DT_>((DT_*)((l.get() - r)));
    }

    template<typename DT_, typename IT_>
    DDVIterator<DT_> operator+(const DDVIterator<DT_>& l, const IT_& r)
    {
      return DDVIterator<DT_>(l.get() + r);
    }

    template<typename DT_, typename IT_>
    DDVIterator<DT_> operator+(const IT_& l, const DDVIterator<DT_>& r)
    {
      return DDVIterator<DT_>(l + r.get());
    }

    template<typename DT_, typename RT_>
    DDVIterator<DT_> operator*(const DDVIterator<DT_>& l, const RT_& r)
    {
      return DDVIterator<DT_>(&(l.get()[r]));
    }

    template<Index _i,
      typename Arch_,
      typename DT_,
      template<typename, typename, typename> class ContType_,
      typename IT_ = Index>
    class DenseDataWrapper
    {
      public:
        typedef DDVIterator<DT_> iterator;

        DenseDataWrapper() :
          _size(_i),
          _num_non_zeros(0),
          _data(_i)
        {
        }

        DenseDataWrapper(DenseDataWrapper && other) :
          _size(other._size),
          _num_non_zeros(other._num_non_zeros),
          _data(std::move(other._data))
        {
        }

        ~DenseDataWrapper()
        {
        }

        iterator begin()
        {
          return iterator(_data.elements());
        }

        iterator end()
        {
          return iterator(_data.elements() + _num_non_zeros);
        }

        IT_ size() const
        {
          return _num_non_zeros;
        }

        IT_ capacity() const
        {
          return _size - _num_non_zeros;
        }

        void push_back(DT_ value)
        {
          _data(_num_non_zeros, value);
          ++_num_non_zeros;
        }

        const DT_ at(IT_ i) const
        {
          //todo in non-zero range check
          return _data(i);
        }

        DT_ at(IT_ i)
        {
          //todo in non-zero range check
          return _data(i);
        }

        DT_ operator[](IT_ i) const
        {
          //todo in non-zero range check
          return _data(i);
        }

        DT_ operator[](IT_ i)
        {
          //todo in non-zero range check
          return _data(i);
        }

        ///TODO make use of () operators of LAFEM only
        void erase(iterator i)
        {
          for(; i < end() - 1; ++i)
            *i = *(i + 1);

          --_num_non_zeros;
        }

        ///TODO make use of () operators of LAFEM only
        void insert(iterator i, DT_ value)
        {
          for(iterator j(end()); j > i ; --j)
            *j = *(j - 1);

          *i = value;
          ++_num_non_zeros;
        }

        DenseDataWrapper& operator=(DenseDataWrapper&& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_size = rhs._size;
          this->_num_non_zeros = rhs._num_non_zeros;
          this->_data = std::move(rhs._data);

          return *this;
        }

      private:
        IT_ _size;
        IT_ _num_non_zeros;
        ContType_<Arch_, DT_, IT_> _data;
    };
  }
}
#endif
