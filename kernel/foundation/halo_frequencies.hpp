#pragma once
#ifndef FOUNDATION_GUARD_HALO_FREQUENCIES_HPP
#define FOUNDATION_GUARD_HALO_FREQUENCIES_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/lafem/arch/scale.hpp>
#include<kernel/lafem/arch/component_product.hpp>
#include<kernel/lafem/arch/dot_product.hpp>

namespace FEAST
{
  namespace Foundation
  {
      /// \todo add communicators
      template <typename Mem_, typename Algo_>
      struct HaloFrequencies
      {
        public:

          template<typename SystemVectorT_,
                   typename BufferVectorT_,
                   typename VectorMirrorT_,
                   template<typename, typename> class StorageT_>
          static SystemVectorT_ value(const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                                      StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >& mirror_buffers,
                                      SystemVectorT_& frequency_buffer)
          {
            if(mirrors.empty())
              return SystemVectorT_();

            frequency_buffer.format(typename SystemVectorT_::DataType(1));

            for(Index i(0); i < Index(mirrors.size()); ++i)
            {
              mirror_buffers.at(i).format(typename BufferVectorT_::DataType(1));
              mirrors.at(i).template scatter_axpy_dual<Algo_>(frequency_buffer, mirror_buffers.at(i));
            }

            // invert frequencies
            frequency_buffer.template component_invert<Algo_>(frequency_buffer);

            return frequency_buffer.clone();
          }
      };
  }
}


#endif
