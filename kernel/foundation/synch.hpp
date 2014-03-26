#pragma once
#ifndef SCARC_GUARD_SYNCH_HH
#define SCARC_GUARD_SYNCH_HH 1

#include<kernel/foundation/comm_base.hpp>

#include<kernel/foundation/communication.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;

///TODO add communicators
namespace FEAST
{
  namespace Foundation
  {
    template<typename Tag_, Tier2CommModes cm_>
    struct SynchVec
    {
    };

    template<typename Tag_>
    struct SynchVec<Tag_, com_exchange>
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

        Status s;
        Comm::send_recv(sendbuf.elements(),
                        sendbuf.size(),
                        dest_rank,
                        recvbuf.elements(),
                        recvbuf.size(),
                        source_rank,
                        s);

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
          SynchVec<Tag_, com_exchange>::execute(target,
                                                       mirrors.at(i),
                                                       sendbufs.at(i),
                                                       recvbufs.at(i),
                                                       dest_ranks.at(i),
                                                       source_ranks.at(i));
        }
      }
    };

    template<typename Tag_>
    struct SynchVec<Tag_, com_accumulate>
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

        Status s;
        Comm::send_recv(sendbuf.elements(),
                               sendbuf.size(),
                               dest_rank,
                               recvbuf.elements(),
                               recvbuf.size(),
                               source_rank,
                               s);

        recvbuf.template axpy<Tag_>(sendbuf, recvbuf);

        mirror.scatter_dual(target, recvbuf);
      }

      //single target, multiple mirrors (stemming from multiple halos)
      template<template<typename, typename> class StorageT_, typename VectorT_, typename VectorMirrorT_>
      static inline void execute(VectorT_& target,
                                 StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                                 StorageT_<VectorT_, std::allocator<VectorT_> >& sendbufs,
                                 StorageT_<VectorT_, std::allocator<VectorT_> >& recvbufs,
                                 StorageT_<Index, std::allocator<Index> >& dest_ranks,
                                 StorageT_<Index, std::allocator<Index> >& source_ranks//,
                                 /*typename VectorT_::DataType weight_vertex = 0.25,
                                 typename VectorT_::DataType weight_edge = 0.5,
                                 typename VectorT_::DataType weight_face = 0.5,
                                 typename VectorT_::DataType weight_polyhedron = 0.5*/
                                 )
      {
        for(Index i(0) ; i < mirrors.size() ; ++i)
        {
          mirrors.at(i).gather_dual(sendbufs.at(i), target);

          Status s;
          Comm::send_recv(sendbufs.at(i).elements(),
              sendbufs.at(i).size(),
              dest_ranks.at(i),
              recvbufs.at(i).elements(),
              recvbufs.at(i).size(),
              source_ranks.at(i),
              s);
        }

        Comm::barrier();

        for(Index i(0) ; i < mirrors.size() ; ++i)
        {
          mirrors.at(i).gather_dual(sendbufs.at(i), target); //we dont need the sendbuf any more

          recvbufs.at(i).template axpy<Tag_>(sendbufs.at(i), recvbufs.at(i));

          mirrors.at(i).scatter_dual(target, recvbufs.at(i));
        }
      }
    };

    template<typename Tag_>
    struct SynchVec<Tag_, com_average>
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

        Comm::send_recv(sendbuf.elements(),
                               sendbuf.size(),
                               dest_rank,
                               recvbuf.elements(),
                               recvbuf.size(),
                               source_rank);

        recvbuf.template axpy<Tag_>(sendbuf, recvbuf);
        recvbuf.template scale<Tag_>(recvbuf, typename VectorT_::DataType(0.5));

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
        SynchVec<Tag_, com_accumulate>::execute(target, mirrors, sendbufs, recvbufs, dest_ranks, source_ranks);

        for(Index i(0) ; i < mirrors.size() ; ++i)
        {
          //average by two along non-vertex halos only (implies divide by four along vertex halos) TODO: generalize will in most cases not work
          if(mirrors.at(i).size() != 1)
          {
            mirrors.at(i).gather_dual(sendbufs.at(i), target);
            sendbufs.at(i).template scale<Tag_>(sendbufs.at(i), typename VectorT_::DataType(0.5));
            mirrors.at(i).scatter_dual(target, sendbufs.at(i));
          }
        }
      }
    };

    template<Tier2CommModes cm_>
    struct SynchScal
    {
    };

    template<>
    struct SynchScal<com_allreduce_sqrtsum>
    {
      //single target, single solver per process
      template<typename DataType_>
      static inline void execute(DataType_& target,
                                 DataType_& sendbuf,
                                 DataType_& recvbuf)
      {
        sendbuf = target;
        Comm::allreduce(&sendbuf, Index(1), &recvbuf);
        target = (DataType_)sqrt(recvbuf);
      }
    };
  }
}

#endif
