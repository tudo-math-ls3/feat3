#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_PRODUCT_MAT_VEC_HPP
#define FOUNDATION_GUARD_GLOBAL_PRODUCT_MAT_VEC_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/global_synch_vec.hpp>
#include<kernel/foundation/environment.hpp>

///TODO add communicators
namespace FEAST
{
  namespace Foundation
  {
      template <typename Mem_, typename Algo_>
      struct GlobalProductMat0Vec1
      {
        public:

#ifndef SERIAL
          template<typename MatrixT_, typename SystemVectorT_, typename VectorMirrorT_, typename BufferVectorT_, template<typename, typename> class StorageT_>
          static void exec(
                           SystemVectorT_& target,
                           const MatrixT_& A,
                           const SystemVectorT_& x,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >& sendbufs,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >& recvbufs,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Communicator communicator = Communicator(MPI_COMM_WORLD)
                           )
          {
            ///assumes type-1 vector (full entries at inner boundaries)
            ///assumes type-0 matrix (entry fractions at inner boundaries)

            A.template apply<Algo_>(target, x);
            GlobalSynchVec0<Mem_, Algo_>::exec(
                                               target,
                                               mirrors,
                                               other_ranks,
                                               sendbufs,
                                               recvbufs,
                                               tags,
                                               communicator);
          }
#else
          template<typename MatrixT_, typename SystemVectorT_, typename VectorMirrorT_, typename BufferVectorT_, template<typename, typename> class StorageT_>
          static void exec(
                           SystemVectorT_&,
                           const MatrixT_&,
                           const SystemVectorT_&,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >&,
                           StorageT_<Index, std::allocator<Index> >&,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >&,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >&,
                           const StorageT_<Index, std::allocator<Index> >&,
                           Communicator = Communicator(0)
                           )
          {
          }
#endif
      };
  }
}


#endif
