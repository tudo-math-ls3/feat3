#pragma once
#ifndef KERNEL_FOUNDATION_ATTRIBUTE_HH
#define KERNEL_FOUNDATION_ATTRIBUTE_HH 1

#include<kernel/base_header.hpp>
#include <vector>
#include <kernel/foundation/attribute_error.hpp>
#include <kernel/foundation/communication.hpp>
#include<kernel/archs.hpp>

using namespace FEAST::Archs;

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief AttributeBase wraps STL containers in a polimorphic way
     *
     *
     * \tparam StorageType_
     * storage type must be STL conformal
     *
     *
     * \author Markus Geveler
     */
    template<template<typename, typename> class StorageType_ = std::vector>
    class AttributeBase
    {
      public:
        virtual Index size() = 0;
        virtual ~AttributeBase()
        {
        }
    };

    /**
     * \brief Attribute wraps STL containers in a polimorphic way
     *
     * \tparam DataType_
     * data type
     *
     * \tparam StorageType_
     * storage type must be STL conformal
     *
     *
     * \author Markus Geveler
     */
    template<typename DataType_, template<typename, typename> class StorageType_ = std::vector>
    class Attribute :
      public AttributeBase<StorageType_>,
      public Bufferable<BufferedData<StorageType_> >,
      public Communicateable<BufferedData<StorageType_>, com_send_receive >
    {
      public:
        Attribute() :
          _data()
        {
        }

        typedef DataType_ data_type_;
        typedef Attribute<DataType_, StorageType_> exact_type_;
        typedef BufferedData<StorageType_> buffer_type_;

        StorageType_<DataType_, std::allocator<DataType_> >& get_data()
        {
          return _data;
        }

        const StorageType_<DataType_, std::allocator<DataType_> >& get_data() const
        {
          return _data;
        }

        virtual Index size()
        {
          return Index(_data.size());
        }

        void push_back(const DataType_ d)
        {
          _data.push_back(d);
        }

        DataType_& at(const Index i)
        {
          return _data.at(i);
        }

        DataType_& at(const Index i) const
        {
          return _data.at(i);
        }

        DataType_& operator[](const Index i) const
        {
          return _data.at(i);
        }

        DataType_& operator[](const Index i)
        {
          return _data.at(i);
        }

        Attribute& operator=(const Attribute& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_data = rhs._data;
          return *this;
        }

        Attribute(const Attribute& other) :
          _data(other._data)
        {
        }

        ///implementation of Bufferable interface
        virtual BufferedData<StorageType_> buffer(Index estimated_size_increase = 0)
        {
          BufferedData<StorageType_> result;
          result.get().push_back(BufferedSharedArray<Index>::create(2));
          result.get().push_back(BufferedSharedArray<DataType_>::create(_data.size() + estimated_size_increase));

          (*(BufferedSharedArray<Index>*)((result.get().at(0).get())))[0] = 2;
          (*(BufferedSharedArray<Index>*)((result.get().at(0).get())))[1] = _data.size() + estimated_size_increase;

          return result;
        }

        virtual void to_buffer(BufferedData<StorageType_>& buffer)
        {
          for(Index i(0) ; i < _data.size() ; ++i)
          {
            (*(BufferedSharedArray<DataType_>*)((buffer.get().at(1).get())))[i] = _data.at(i);
          }
        }

        virtual void from_buffer(const BufferedData<StorageType_>& buffer)
        {
          _data.clear();
          for(Index i(0) ; i < (*(BufferedSharedArray<Index>*)((buffer.get().at(0).get())))[1] ; ++i)
          {
            _data.push_back( (DataType_)(*(BufferedSharedArray<DataType_>*)((buffer.get().at(1).get())))[i] );
          }
        }

        ///implementation of Communicateable interface
        void send_recv(BufferedData<StorageType_>& sendbuffers,
                       int destrank,
                       BufferedData<StorageType_>& recvbuffers,
                       int sourcerank)
        {
#ifndef SERIAL
          Comm<Parallel>::send_recv(((BufferedSharedArray<Index>*)sendbuffers.get().at(0).get())->get(),
              (*(BufferedSharedArray<Index>*)((sendbuffers.get().at(0).get())))[0],
              destrank,
              ((BufferedSharedArray<Index>*)recvbuffers.get().at(0).get())->get(),
              (*(BufferedSharedArray<Index>*)((recvbuffers.get().at(0).get())))[0],
              sourcerank);

          Comm<Parallel>::send_recv(((BufferedSharedArray<DataType_>*)sendbuffers.get().at(1).get())->get(),
              (*(BufferedSharedArray<Index>*)((sendbuffers.get().at(0).get())))[1],
              destrank,
              ((BufferedSharedArray<DataType_>*)recvbuffers.get().at(1).get())->get(),
              (*(BufferedSharedArray<Index>*)((recvbuffers.get().at(0).get())))[1],
              sourcerank);

#else
          Comm<Serial>::send_recv(((BufferedSharedArray<Index>*)sendbuffers.get().at(0).get())->get(),
              (*(BufferedSharedArray<Index>*)((sendbuffers.get().at(0).get())))[0],
              destrank,
              ((BufferedSharedArray<Index>*)recvbuffers.get().at(0).get())->get(),
              (*(BufferedSharedArray<Index>*)((recvbuffers.get().at(0).get())))[0],
              sourcerank);

          Comm<Serial>::send_recv(((BufferedSharedArray<DataType_>*)sendbuffers.get().at(1).get())->get(),
              (*(BufferedSharedArray<Index>*)((sendbuffers.get().at(0).get())))[1],
              destrank,
              ((BufferedSharedArray<DataType_>*)recvbuffers.get().at(1).get())->get(),
              (*(BufferedSharedArray<Index>*)((recvbuffers.get().at(0).get())))[1],
              sourcerank);

#endif
        }


      protected:
        StorageType_<DataType_, std::allocator<DataType_> > _data;
    };
  }
}
#endif
