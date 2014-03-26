#pragma once
#ifndef FOUNDATION_GUARD_HALO_FREQUENCIES_HPP
#define FOUNDATION_GUARD_HALO_FREQUENCIES_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/lafem/arch/scale.hpp>
#include<kernel/lafem/arch/component_product.hpp>
#include<kernel/lafem/arch/dot_product.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

///TODO add communicators
namespace FEAST
{
  namespace Foundation
  {
      template <typename Mem_, typename Algo_>
      struct HaloFrequencies
      {
      };

      template <>
      struct HaloFrequencies<Mem::Main, Algo::Generic>
      {
        public:

          template<typename VectorT_,
                   typename VectorMirrorT_,
                   template<typename, typename> class StorageT_>
          static VectorT_ value(const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                                StorageT_<VectorT_, std::allocator<VectorT_> >& mirror_buffers,
                                StorageT_<VectorT_, std::allocator<VectorT_> >& frequency_buffers)
          {
            if(frequency_buffers.size() == 0)
              return VectorT_();

            VectorT_ result(frequency_buffers.at(0).size(), typename VectorT_::DataType(1));

            //clear buffers
            for(auto& fb_i : frequency_buffers)
              for(Index i(0) ; i < fb_i.size() ; ++i)
                fb_i(i, typename VectorT_::DataType(0));

            for(auto& mb_i : mirror_buffers)
              for(Index i(0) ; i < mb_i.size() ; ++i)
                mb_i(i, typename VectorT_::DataType(1));

            typename VectorT_::DataType* result_data(result.elements());
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              const typename VectorT_::DataType* f_data(frequency_buffers.at(i).elements());
              mirrors.at(i).scatter_dual(frequency_buffers.at(i), mirror_buffers.at(i));
              for(Index j(0) ; j < result.size() ; ++j)
              {
                result_data[j] += f_data[j];
              }
            }

            return result;
          }
      };
  }
}


#endif
