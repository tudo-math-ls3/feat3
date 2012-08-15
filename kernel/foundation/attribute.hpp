#pragma once
#ifndef KERNEL_FOUNDATION_ATTRIBUTE_HH
#define KERNEL_FOUNDATION_ATTRIBUTE_HH 1

#include <vector>
#include <kernel/foundation/attribute_error.hpp>
#include<kernel/base_header.hpp>

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
    class Attribute : public AttributeBase<StorageType_>
    {
      public:
        Attribute() :
          _data()
        {
        }

        StorageType_<DataType_, std::allocator<DataType_> >& get_data()
        {
          return _data;
        }

        virtual Index size()
        {
          return _data.size();
        }

        void push_back(DataType_ d)
        {
            _data.push_back(d);
        }

        DataType_& at(Index i)
        {
          return _data.at(i);
        }

        Attribute& operator=(Attribute& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_data = rhs._data;
          return *this;
        }


      protected:
        StorageType_<DataType_, std::allocator<DataType_> > _data;
    };
  }
}
#endif
