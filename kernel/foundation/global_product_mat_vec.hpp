#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_PRODUCT_MAT_VEC_HPP
#define FOUNDATION_GUARD_GLOBAL_PRODUCT_MAT_VEC_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/global_synch_vec.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

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
          template<typename MatrixT_, typename VectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static void exec(
                           VectorT_& target,
                           const MatrixT_& A,
                           const VectorT_& x,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                           StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                           Index tag = 0,
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
                                               tag,
                                               communicator);
          }
#else
          template<typename MatrixT_, typename VectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static void exec(
                           VectorT_&,
                           const MatrixT_&,
                           const VectorT_&,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >&,
                           StorageT_<Index, std::allocator<Index> >&,
                           StorageT_<VectorT_, std::allocator<VectorT_> >&,
                           StorageT_<VectorT_, std::allocator<VectorT_> >&,
                           Index = 0,
                           Communicator = Communicator(0)
                           )
          {
          }
#endif
      };
  }
}


#endif
