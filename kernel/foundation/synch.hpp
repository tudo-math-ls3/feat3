#pragma once
#ifndef SCARC_GUARD_SYNCH_HH
#define SCARC_GUARD_SYNCH_HH 1

#include<kernel/foundation/communication.hpp>
#include<kernel/lafem/sum.hpp>
#include<kernel/lafem/scale.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;

///TODO add communicators
namespace FEAST
{
  namespace Foundation
  {
    template<typename Tag_, typename Arch_, Tier2CommModes cm_>
    struct SynchVec
    {
    };

    template<typename Tag_, typename Arch_>
    struct SynchVec<Tag_, Arch_, com_exchange>
    {
      //single target, single mirror
      template<typename VectorT_, typename VectorMirrorT_>
      static inline void execute(VectorT_& target,
                                 VectorMirrorT_& mirror,
                                 VectorT_& sendbuf,
                                 VectorT_& recvbuf,
                                 Index dest_rank,
                                 Index source_rank)
      {
        mirror.gather_dual(sendbuf, target);

        Comm<Arch_>::send_recv(sendbuf.elements(),
                               sendbuf.size(),
                               dest_rank,
                               recvbuf.elements(),
                               recvbuf.size(),
                               source_rank);

        mirror.scatter_dual(target, recvbuf);
      }

      //single target, multiple mirrors (stemming from multiple halos)
      template<template<typename, typename> class StorageT_, typename VectorT_, typename VectorMirrorT_>
      static inline void execute(VectorT_& target,
                                 StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                                 StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                                 StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                                 StorageT_<Index, std::allocator<Index> >& dest_ranks,
                                 StorageT_<Index, std::allocator<Index> >& source_ranks)
      {
        for(Index i(0) ; i < mirrors.size() ; ++i)
        {
          SynchVec<Tag_, Arch_, com_exchange>::execute(target,
                                                       mirrors.at(i),
                                                       sendbufs.at(i),
                                                       recvbufs.at(i),
                                                       dest_ranks.at(i),
                                                       source_ranks.at(i));
        }
      }
    };

    template<typename Tag_, typename Arch_>
    struct SynchVec<Tag_, Arch_, com_accumulate>
    {
      //single target, single mirror
      template<typename VectorT_, typename VectorMirrorT_>
      static inline void execute(VectorT_& target,
                                 VectorMirrorT_& mirror,
                                 VectorT_& sendbuf,
                                 VectorT_& recvbuf,
                                 Index dest_rank,
                                 Index source_rank)
      {
        mirror.gather_dual(sendbuf, target);

        Comm<Arch_>::send_recv(sendbuf.elements(),
                               sendbuf.size(),
                               dest_rank,
                               recvbuf.elements(),
                               recvbuf.size(),
                               source_rank);

        Sum<Tag_>::value(recvbuf, sendbuf, recvbuf);

        mirror.scatter_dual(target, recvbuf);
      }

      //single target, multiple mirrors (stemming from multiple halos)
      template<template<typename, typename> class StorageT_, typename VectorT_, typename VectorMirrorT_>
      static inline void execute(VectorT_& target,
                                 StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                                 StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                                 StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                                 StorageT_<Index, std::allocator<Index> >& dest_ranks,
                                 StorageT_<Index, std::allocator<Index> >& source_ranks)
      {
        for(Index i(0) ; i < mirrors.size() ; ++i)
        {
          SynchVec<Tag_, Arch_, com_accumulate>::execute(target,
                                                         mirrors.at(i),
                                                         sendbufs.at(i),
                                                         recvbufs.at(i),
                                                         dest_ranks.at(i),
                                                         source_ranks.at(i));
        }
      }
    };

    template<typename Tag_, typename Arch_>
    struct SynchVec<Tag_, Arch_, com_average>
    {
      //single target, single mirror
      template<typename VectorT_, typename VectorMirrorT_>
      static inline void execute(VectorT_& target,
                                 VectorMirrorT_& mirror,
                                 VectorT_& sendbuf,
                                 VectorT_& recvbuf,
                                 Index dest_rank,
                                 Index source_rank)
      {
        mirror.gather_dual(sendbuf, target);

        Comm<Arch_>::send_recv(sendbuf.elements(),
                               sendbuf.size(),
                               dest_rank,
                               recvbuf.elements(),
                               recvbuf.size(),
                               source_rank);

        Sum<Tag_>::value(recvbuf, sendbuf, recvbuf);
        Scale<Tag_>::value(recvbuf, recvbuf, typename VectorT_::DataType(0.5));

        mirror.scatter_dual(target, recvbuf);
      }

      //single target, multiple mirrors (stemming from multiple halos)
      template<template<typename, typename> class StorageT_, typename VectorT_, typename VectorMirrorT_>
      static inline void execute(VectorT_& target,
                                 StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                                 StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                                 StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                                 StorageT_<Index, std::allocator<Index> >& dest_ranks,
                                 StorageT_<Index, std::allocator<Index> >& source_ranks)
      {
        for(Index i(0) ; i < mirrors.size() ; ++i)
        {
          SynchVec<Tag_, Arch_, com_average>::execute(target,
                                                         mirrors.at(i),
                                                         sendbufs.at(i),
                                                         recvbufs.at(i),
                                                         dest_ranks.at(i),
                                                         source_ranks.at(i));
        }
      }
    };

    template<typename Arch_, Tier2CommModes cm_>
    struct SynchScal
    {
    };

    template<typename Arch_>
    struct SynchScal<Arch_, com_allreduce_sqrtsum>
    {
      //single target, single solver per process
      template<typename DataType_>
      static inline void execute(DataType_& target,
                                 DataType_& sendbuf,
                                 DataType_& recvbuf)
      {
        sendbuf = target;
        Comm<Arch_>::allreduce(&sendbuf, Index(1), &recvbuf);
        target = sqrt(recvbuf);
      }
    };

    template<>
    struct SynchScal<Serial, com_allreduce_sqrtsum>
    {
      //single target, single solver per process
      template<typename DataType_>
      static inline void execute(DataType_& target,
                                 DataType_& sendbuf,
                                 DataType_& recvbuf)
      {
        target = sqrt(target);
      }
    };
  }
}

#endif
