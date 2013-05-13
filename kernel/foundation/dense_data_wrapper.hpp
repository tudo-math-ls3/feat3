#pragma once
#ifndef KERNEL_FOUNDATION_DENSE_DATA_WRAPPER_HH
#define KERNEL_FOUNDATION_DENSE_DATA_WRAPPER_HH 1

#include <kernel/foundation/base.hpp>

namespace FEAST
{
  namespace Foundation
  {
    template<Index _i,
             typename Arch_,
             typename DT_,
             template<typename, typename> class ContType_>
      class DenseDataWrapper
      {
        public:
          DenseDataWrapper() :
            _size(_i),
            _num_non_zeros(0),
            _data(_i)
        {
        }

          ~DenseDataWrapper()
          {
          }

          Index size()
          {
            return _num_non_zeros;
          }

          Index capacity()
          {
            return _size - _num_non_zeros;
          }

          void push_back(DT_ value)
          {
            //todo capacity check
            _data(_num_non_zeros, value);
            ++_num_non_zeros;
          }

          const DT_ at(Index i)
          {
            //todo in non-zero range check
            return _data(i);
          }

          const DT_ operator[](Index i)
          {
            //todo in non-zero range check
            return _data(i);
          }

          DenseDataWrapper& operator=(const DenseDataWrapper& rhs)
          {
            if(this == &rhs)
              return *this;

            this->_size = rhs._size;
            this->_num_non_zeros = rhs._num_non_zeros;
            this->_data = rhs._data;

            return *this;
          }

        private:
          Index _size;
          Index _num_non_zeros;
          ContType_<Arch_, DT_> _data;
      };
  }
}
#endif
