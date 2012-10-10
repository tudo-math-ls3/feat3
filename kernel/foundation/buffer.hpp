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
    class BufferedSharedArray : public SharedArrayBase, public Communicateable<T_>
    {
      public:
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

        ~BufferedSharedArray()
        {
          delete[] _data;
        }

        ///Implementation of Communicateable interface
        void send_recv(int destrank,
                       T_& recvdata,
                       int sourcerank)
        {
          ///TODO change semantics of communicateable: THIS can be directly used as a buffer
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
      virtual T_ buffer() = 0;
      virtual void to_buffer(T_& buffer) = 0;
      virtual void from_buffer(const T_& buffer) = 0;
    };

    template<template<typename, typename> class ST_ = std::vector>
    class BufferedData
    {
      public:
        BufferedData() :
          _data()
        {
        }

        ST_<std::shared_ptr<SharedArrayBase>, std::allocator<std::shared_ptr<SharedArrayBase> > >& get()
        {
          return _data;
        }

        ST_<std::shared_ptr<SharedArrayBase>, std::allocator<std::shared_ptr<SharedArrayBase> > >& get() const
        {
          return _data;
        }

      private:
        ST_<std::shared_ptr<SharedArrayBase>, std::allocator<std::shared_ptr<SharedArrayBase> > > _data;
    };
  }
}
#endif
