#pragma once
#ifndef KERNEL_FOUNDATION_BUFFER_HH
#define KERNEL_FOUNDATION_BUFFER_HH 1

#include<kernel/base_header.hpp>
#include<kernel/util/cpp11_smart_pointer.hpp>
#include<kernel/foundation/communication.hpp>
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

    ///does not need to be copy-assignable, since only shared_ptrs are going to be copied externally
    template<typename T_>
    class BufferedSharedArray : public SharedArrayBase
    {
      public:
        typedef T_ data_type_;

        static std::shared_ptr<SharedArrayBase> create(Index i)
        {
          return std::shared_ptr<SharedArrayBase>(new BufferedSharedArray<T_>(i));
        }

        T_& operator[](Index i) const
        {
          return _data[i];
        }

        T_& operator[](Index i)
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
      virtual T_ buffer() = 0;
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
          _data(),
          _sizes()
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

        ST_<Index, std::allocator<Index> >& get_sizes()
        {
          return _sizes;
        }

        const ST_<Index, std::allocator<Index> >& get_sizes() const
        {
          return _sizes;
        }

      private:
        ST_<std::shared_ptr<SharedArrayBase>, std::allocator<std::shared_ptr<SharedArrayBase> > > _data;
        ST_<Index, std::allocator<Index> > _sizes;
    };
  }
}
#endif
