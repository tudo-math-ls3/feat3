#pragma once
#ifndef FOUNDATION_GUARD_GATEWAY_HPP
#define FOUNDATION_GUARD_GATEWAY_HPP 1

#include<kernel/foundation/global_dot.hpp>
#include<kernel/foundation/global_norm.hpp>
#include<kernel/foundation/global_synch_vec.hpp>
#include<kernel/foundation/global_product_mat_vec.hpp>

namespace FEAST
{
  namespace Foundation
  {

    template<typename Mem_, typename VectorT_>
    class GlobalDotGateway : public LAFEM::Arch::DotGatewayBase<Mem_, VectorT_>
    {
      public:
        GlobalDotGateway(VectorT_& frequencies) :
          _frequencies(frequencies)
        {
        }

        virtual typename VectorT_::DataType value(const VectorT_& x, const VectorT_& y) const override
        {
          typename VectorT_::DataType r = 0;
          GlobalDot<Mem_>::value(r,
                                        x,
                                        y,
                                        _frequencies);
          return r;
        }

      private:
        VectorT_& _frequencies;
    };

    template<typename Mem_, typename VectorT_>
    class GlobalNorm2Gateway : public LAFEM::Arch::Norm2GatewayBase<Mem_, VectorT_>
    {
      public:
        GlobalNorm2Gateway(VectorT_& frequencies) :
          _frequencies(frequencies)
        {
        }

        virtual typename VectorT_::DataType value(const VectorT_& x) const override
        {
          typename VectorT_::DataType r = 0;
          GlobalNorm2<Mem_>::value(r,
                                          x,
                                          _frequencies);
          return r;
        }

      private:
        VectorT_& _frequencies;
    };

    template<typename Mem_, typename VectorT_>
    class GlobalNorm2SquaredGateway : public LAFEM::Arch::Norm2SquaredGatewayBase<Mem_, VectorT_>
    {
      public:
        GlobalNorm2SquaredGateway(VectorT_& frequencies) :
          _frequencies(frequencies)
        {
        }

        virtual typename VectorT_::DataType value(const VectorT_& x) const override
        {
          typename VectorT_::DataType r = 0;
          GlobalNorm2Squared<Mem_>::value(r,
                                                 x,
                                                 _frequencies);
          return r;
        }

      private:
        VectorT_& _frequencies;
    };

    template<typename Mem_,
             typename VectorT_,
             typename VectorMirrorT_,
             template<typename, typename> class StorageT_ = std::vector>
    class GlobalSynchVec0Gateway : public LAFEM::Arch::SynchVec0GatewayBase<Mem_, VectorT_>
    {
      public:
#ifndef SERIAL
        GlobalSynchVec0Gateway(
                               const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                               StorageT_<Index, std::allocator<Index> >& other_ranks,
                               StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                               StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                               StorageT_<Index, std::allocator<Index> >& tags,
                               Communicator communicator = Communicator(MPI_COMM_WORLD)) :
          _mirrors(mirrors),
          _other_ranks(other_ranks),
          _sendbufs(sendbufs),
          _recvbufs(recvbufs),
          _tags(tags),
          _communicator(communicator)
        {
        }
#else
        GlobalSynchVec0Gateway(
                               const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                               StorageT_<Index, std::allocator<Index> >& other_ranks,
                               StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                               StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                               StorageT_<Index, std::allocator<Index> >& tags,
                               Communicator communicator = Communicator(0)) :
          _mirrors(mirrors),
          _other_ranks(other_ranks),
          _sendbufs(sendbufs),
          _recvbufs(recvbufs),
          _tags(tags),
          _communicator(communicator)
        {
        }
#endif

        virtual VectorT_& value(VectorT_& x) const override
        {
          GlobalSynchVec0<Mem_>::exec(
                                              x,
                                              _mirrors,
                                              _other_ranks,
                                              _sendbufs,
                                              _recvbufs,
                                              _tags,
                                              _communicator);
          return x;
        }

      private:
        const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& _mirrors;
        StorageT_<Index, std::allocator<Index> >& _other_ranks;
        StorageT_<VectorT_, std::allocator<VectorT_> >& _sendbufs;
        StorageT_<VectorT_, std::allocator<VectorT_> >& _recvbufs;
        StorageT_<Index, std::allocator<Index> >& _tags;
        Communicator _communicator;
    };

    template<typename Mem_,
             typename VectorT_,
             typename VectorMirrorT_,
             template<typename, typename> class StorageT_ = std::vector>
    class GlobalSynchVec1Gateway : public LAFEM::Arch::SynchVec1GatewayBase<Mem_, VectorT_>
    {
      public:
#ifndef SERIAL
        GlobalSynchVec1Gateway(
                               const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                               const VectorT_& frequencies,
                               StorageT_<Index, std::allocator<Index> >& other_ranks,
                               StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                               StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                               StorageT_<Index, std::allocator<Index> >& tags,
                               Communicator communicator = Communicator(MPI_COMM_WORLD)) :
          _mirrors(mirrors),
          _frequencies(frequencies),
          _other_ranks(other_ranks),
          _sendbufs(sendbufs),
          _recvbufs(recvbufs),
          _tags(tags),
          _communicator(communicator)
        {
        }
#else
        GlobalSynchVec1Gateway(
                               const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                               const VectorT_& frequencies,
                               StorageT_<Index, std::allocator<Index> >& other_ranks,
                               StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                               StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                               StorageT_<Index, std::allocator<Index> >& tags,
                               Communicator communicator = Communicator(0)) :
          _mirrors(mirrors),
          _frequencies(frequencies),
          _other_ranks(other_ranks),
          _sendbufs(sendbufs),
          _recvbufs(recvbufs),
          _tags(tags),
          _communicator(communicator)
        {
        }
#endif

        virtual VectorT_& value(VectorT_& x) const override
        {
          GlobalSynchVec1<Mem_>::exec(
                                              x,
                                              _mirrors,
                                              _frequencies,
                                              _other_ranks,
                                              _sendbufs,
                                              _recvbufs,
                                              _tags,
                                              _communicator);
          return x;
        }

      private:
        const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& _mirrors;
        const VectorT_& _frequencies;
        StorageT_<Index, std::allocator<Index> >& _other_ranks;
        StorageT_<VectorT_, std::allocator<VectorT_> >& _sendbufs;
        StorageT_<VectorT_, std::allocator<VectorT_> >& _recvbufs;
        StorageT_<Index, std::allocator<Index> >& _tags;
        Communicator _communicator;
    };

    template<typename Mem_,
             typename VectorT_,
             typename MatrixT_,
             typename VectorMirrorT_,
             template<typename, typename> class StorageT_ = std::vector>
    class GlobalProductMat0Vec1Gateway : public LAFEM::Arch::ProductMat0Vec1GatewayBase<Mem_, VectorT_, MatrixT_>
    {
      public:
#ifndef SERIAL
        GlobalProductMat0Vec1Gateway(
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                           StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                           StorageT_<Index, std::allocator<Index> >& tags,
                           Communicator communicator = Communicator(MPI_COMM_WORLD)
                               ) :
          _mirrors(mirrors),
          _other_ranks(other_ranks),
          _sendbufs(sendbufs),
          _recvbufs(recvbufs),
          _tags(tags),
          _communicator(communicator)
        {
        }
#else
        GlobalProductMat0Vec1Gateway(
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                           StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                           StorageT_<Index, std::allocator<Index> >& tags,
                           Communicator communicator = Communicator(0)
                               ) :
          _mirrors(mirrors),
          _other_ranks(other_ranks),
          _sendbufs(sendbufs),
          _recvbufs(recvbufs),
          _tags(tags),
          _communicator(communicator)
        {
        }
#endif

        virtual VectorT_& value(VectorT_& r, const MatrixT_& A, const VectorT_& x) override
        {
          GlobalProductMat0Vec1<Mem_>::exec(
                                                   r,
                                                   A,
                                                   x,
                                                   _mirrors,
                                                   _other_ranks,
                                                   _sendbufs,
                                                   _recvbufs,
                                                   _tags,
                                                   _communicator
                                                   );
          return r;
        }

      private:
        const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& _mirrors;
        StorageT_<Index, std::allocator<Index> >& _other_ranks;
        StorageT_<VectorT_, std::allocator<VectorT_> >& _sendbufs;
        StorageT_<VectorT_, std::allocator<VectorT_> >& _recvbufs;
        StorageT_<Index, std::allocator<Index> >& _tags;
        Communicator _communicator;
    };
  }
}

#endif
