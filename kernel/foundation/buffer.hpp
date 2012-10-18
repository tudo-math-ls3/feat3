#pragma once
#ifndef KERNEL_FOUNDATION_BUFFER_HH
#define KERNEL_FOUNDATION_BUFFER_HH 1

#include<kernel/base_header.hpp>
#include<kernel/util/cpp11_smart_pointer.hpp>
#include<vector>

namespace FEAST
{
  namespace Foundation
  {
    class SharedArrayBase
    {
      public:
        virtual Index size() const = 0;
        virtual ~SharedArrayBase()
        {
        }
    };

    template<typename T_>
    class BufferedSharedArray : public SharedArrayBase
    {
      public:
        typedef T_ data_type_;

        static std::shared_ptr<SharedArrayBase> create(Index i)
        {
          return std::shared_ptr<SharedArrayBase>(new BufferedSharedArray<T_>(i));
        }

        T_& operator[](const Index i) const
        {
          return _data[i];
        }

        T_& operator[](const Index i)
        {
          return _data[i];
        }

        Index size() const
        {
          return _size;
        }

        T_* get()
        {
          return _data;
        }

        const T_* get() const
        {
          return _data;
        }

        ~BufferedSharedArray()
        {
          delete[] _data;
        }

        BufferedSharedArray& operator=(const BufferedSharedArray& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_data = rhs._data;
          this->_size = rhs._size;
          return *this;
        }


        BufferedSharedArray(const BufferedSharedArray& other) :
          _data(other._data),
          _size(other._size)
        {
        }

      private:

        BufferedSharedArray(Index i) :
          _data(new T_[i]),
          _size(i)
        {
        }

        T_* _data;
        Index _size;
    };

    template<typename T_>
    struct Bufferable
    {
      ///must implement how buffers are created
      virtual T_ buffer(Index estimated_size_increase = 0) = 0;
      ///must implement how buffers are filled
      virtual void to_buffer(T_& buffer) = 0;
      ///must implement how data is transfered from buffers
      virtual void from_buffer(const T_& buffer) = 0;
    };

    template<template<typename, typename> class ST_ = std::vector>
    class BufferedData
    {
      public:
        typedef ST_<std::shared_ptr<SharedArrayBase>, std::allocator<std::shared_ptr<SharedArrayBase> > > compound_storage_type_;

        BufferedData() :
          _data()
        {
        }

        ST_<std::shared_ptr<SharedArrayBase>, std::allocator<std::shared_ptr<SharedArrayBase> > >& get()
        {
          return _data;
        }

        const ST_<std::shared_ptr<SharedArrayBase>, std::allocator<std::shared_ptr<SharedArrayBase> > >& get() const
        {
          return _data;
        }

      private:
        ST_<std::shared_ptr<SharedArrayBase>, std::allocator<std::shared_ptr<SharedArrayBase> > > _data;
    };
  }
}
#endif
