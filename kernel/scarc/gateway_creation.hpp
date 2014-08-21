/**
 * \file
 * \brief FEAST milestone 2 gateway creation implementations
 * \author Markus Geveler
 * \date 2014
 *
 * See class documentation.
 *
 */

#pragma once
#ifndef SCARC_GUARD_SCARC_GATEWAY_CREATION_HPP
#define SCARC_GUARD_SCARC_GATEWAY_CREATION_HPP 1

#include<kernel/foundation/gateway.hpp>

using namespace FEAST;
using namespace FEAST::Foundation;

namespace FEAST
{
  namespace ScaRC
  {
    enum GlobalOpType
    {
      got_dot = 0,
      got_nrm2,
      got_nrm2sqr,
      got_synch_vec0,
      got_synch_vec1,
      got_product_mat0_vec1
    };

    template<GlobalOpType got_, typename Algo_, template<typename, typename> class StorageT_ = std::vector>
    struct GatewayCreation
    {
    };

    template<typename Algo_>
    struct GatewayCreation<got_dot, Algo_>
    {
      template<typename SynchedScaRCDataT_>
      static GlobalDotGateway<typename SynchedScaRCDataT_::mem_, Algo_, typename SynchedScaRCDataT_::vector_type_> value(SynchedScaRCDataT_& data)
      {
        return GlobalDotGateway<typename SynchedScaRCDataT_::mem_, Algo_, typename SynchedScaRCDataT_::vector_type_>(data.halo_frequencies());
      }
    };

    template<typename Algo_>
    struct GatewayCreation<got_nrm2, Algo_>
    {
      template<typename SynchedScaRCDataT_>
      static GlobalNorm2Gateway<typename SynchedScaRCDataT_::mem_, Algo_, typename SynchedScaRCDataT_::vector_type_> value(SynchedScaRCDataT_& data)
      {
        return GlobalNorm2Gateway<typename SynchedScaRCDataT_::mem_, Algo_, typename SynchedScaRCDataT_::vector_type_>(data.halo_frequencies());
      }
    };

    template<typename Algo_>
    struct GatewayCreation<got_nrm2sqr, Algo_>
    {
      template<typename SynchedScaRCDataT_>
      static GlobalNorm2SquaredGateway<typename SynchedScaRCDataT_::mem_, Algo_, typename SynchedScaRCDataT_::vector_type_> value(SynchedScaRCDataT_& data)
      {
        return GlobalNorm2SquaredGateway<typename SynchedScaRCDataT_::mem_, Algo_, typename SynchedScaRCDataT_::vector_type_>(data.halo_frequencies());
      }
    };

    template<typename Algo_, template<typename, typename> class StorageT_>
    struct GatewayCreation<got_synch_vec0, Algo_, StorageT_>
    {
      template<typename SynchedScaRCDataT_>
      static GlobalSynchVec0Gateway<typename SynchedScaRCDataT_::mem_,
                                    Algo_,
                                    typename SynchedScaRCDataT_::vector_type_,
                                    typename SynchedScaRCDataT_::vector_mirror_type_,
                                    StorageT_> value(SynchedScaRCDataT_& data)
      {
        return GlobalSynchVec0Gateway<typename SynchedScaRCDataT_::mem_,
                                      Algo_,
                                      typename SynchedScaRCDataT_::vector_type_,
                                      typename SynchedScaRCDataT_::vector_mirror_type_,
                                      StorageT_>(
                                                 data.vector_mirrors(),
                                                 data.dest_ranks(),
                                                 data.vector_mirror_sendbufs(),
                                                 data.vector_mirror_recvbufs(),
                                                 data.tags()
                                                );
      }
    };

    template<typename Algo_, template<typename, typename> class StorageT_>
    struct GatewayCreation<got_synch_vec1, Algo_, StorageT_>
    {
      template<typename SynchedScaRCDataT_>
      static GlobalSynchVec1Gateway<typename SynchedScaRCDataT_::mem_,
                                    Algo_,
                                    typename SynchedScaRCDataT_::vector_type_,
                                    typename SynchedScaRCDataT_::vector_mirror_type_,
                                    StorageT_> value(SynchedScaRCDataT_& data)
      {
        return GlobalSynchVec1Gateway<typename SynchedScaRCDataT_::mem_,
                                      Algo_,
                                      typename SynchedScaRCDataT_::vector_type_,
                                      typename SynchedScaRCDataT_::vector_mirror_type_,
                                      StorageT_>(
                                                 data.vector_mirrors(),
                                                 data.halo_frequencies(),
                                                 data.dest_ranks(),
                                                 data.vector_mirror_sendbufs(),
                                                 data.vector_mirror_recvbufs(),
                                                 data.tags()
                                                );
      }
    };

    template<typename Algo_, template<typename, typename> class StorageT_>
    struct GatewayCreation<got_product_mat0_vec1, Algo_, StorageT_>
    {
      template<typename SynchedScaRCDataT_>
      static GlobalProductMat0Vec1Gateway<typename SynchedScaRCDataT_::mem_,
                                    Algo_,
                                    typename SynchedScaRCDataT_::vector_type_,
                                    typename SynchedScaRCDataT_::matrix_type_,
                                    typename SynchedScaRCDataT_::vector_mirror_type_,
                                    StorageT_> value(SynchedScaRCDataT_& data)
      {
        return GlobalProductMat0Vec1Gateway<typename SynchedScaRCDataT_::mem_,
                                      Algo_,
                                      typename SynchedScaRCDataT_::vector_type_,
                                      typename SynchedScaRCDataT_::matrix_type_,
                                      typename SynchedScaRCDataT_::vector_mirror_type_,
                                      StorageT_>(
                                                 data.vector_mirrors(),
                                                 data.dest_ranks(),
                                                 data.vector_mirror_sendbufs(),
                                                 data.vector_mirror_recvbufs(),
                                                 data.tags()
                                                );
      }
    };
  }
}

#endif
