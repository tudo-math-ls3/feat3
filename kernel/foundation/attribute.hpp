#pragma once
#ifndef KERNEL_FOUNDATION_ATTRIBUTE_HH
#define KERNEL_FOUNDATION_ATTRIBUTE_HH 1

#include <vector>
#include<kernel/base_header.hpp>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief bla
     *
     * bla
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

      protected:
        unsigned long _size;
    };

    template<typename DataType_, template<typename, typename> class StorageType_ = std::vector>
    class Attribute : public AttributeBase<StorageType_>
    {
      public:
        Attribute(unsigned long s) :
          _data(StorageType_<DataType_, std::allocator<DataType_> >())
        {
          this->_size = s;
        }

        StorageType_<DataType_, std::allocator<DataType_> >& get_data()
        {
          return _data;
        }

        virtual Index size()
        {
          return this->_size;
        }

      private:
        StorageType_<DataType_, std::allocator<DataType_> > _data;
    };

  }
}
#endif
