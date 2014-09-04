#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_DEFECT_HPP
#define FOUNDATION_GUARD_GLOBAL_DEFECT_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/global_product_mat_vec.hpp>
#include<kernel/lafem/arch/difference.hpp>
#include<kernel/foundation/environment.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

///TODO add communicators
namespace FEAST
{
  namespace Foundation
  {
      template <typename Mem_, typename Algo_>
      struct GlobalDefect
      {
        public:
#ifndef SERIAL
          template<typename MatrixT_, typename SystemVectorT_, typename VectorMirrorT_, typename BufferVectorT_, template<typename, typename> class StorageT_>
          static void exec(
                           SystemVectorT_& target,
                           const SystemVectorT_& b,
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

            GlobalProductMat0Vec1<Mem::Main, Algo::Generic>::exec(
                                                                  target,
                                                                  A,
                                                                  x,
                                                                  mirrors,
                                                                  other_ranks,
                                                                  sendbufs,
                                                                  recvbufs,
                                                                  tags,
                                                                  communicator);

            Difference<Mem_, Algo_>::value(target.elements(), b.elements(), target.elements(), target.size());
          }
#else
          template<typename MatrixT_, typename SystemVectorT_, typename VectorMirrorT_, typename BufferVectorT_, template<typename, typename> class StorageT_>
          static void exec(
                           SystemVectorT_& target,
                           const SystemVectorT_& b,
                           const MatrixT_& A,
                           const SystemVectorT_& x,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >&,
                           StorageT_<Index, std::allocator<Index> >&,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >&,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >&,
                           const StorageT_<Index, std::allocator<Index> >&,
                           Communicator = Communicator(0)
                           )
          {
            A.template apply<Algo_>(target, x);
            Difference<Mem_, Algo_>::value(target.elements(), b.elements(), target.elements(), target.size());
          }
#endif
      };
  }
}


#endif
