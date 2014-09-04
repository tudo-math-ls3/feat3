#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_SYNCH_VEC_HPP
#define FOUNDATION_GUARD_GLOBAL_SYNCH_VEC_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
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
      struct GlobalSynchVec0
      {
      };

      template <typename Algo_>
      struct GlobalSynchVec0<Mem::Main, Algo_>
      {
        public:

#ifndef SERIAL
          ///TODO implement version with active queue of requests
          template<typename TargetVectorT_, typename VectorMirrorT_, typename BufferVectorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_& target,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >& sendbufs,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >& recvbufs,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Communicator communicator = Communicator(MPI_COMM_WORLD)
                           )
          {
            if(mirrors.size() == 0)
              return;

            ///assumes type-0 vector (entry fractions at inner boundaries)

            ///start recv
            StorageT_<Request, std::allocator<Request> > recvrequests;
            StorageT_<Status, std::allocator<Status> > recvstatus;
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              Request rr;
              Status rs;

              recvrequests.push_back(std::move(rr));
              recvstatus.push_back(std::move(rs));

              Comm::irecv(recvbufs.at(i).elements(),
                          recvbufs.at(i).size(),
                          other_ranks.at(i),
                          recvrequests.at(i),
                          tags.at(i),
                          communicator
                         );

            }

            ///gather and start send
            StorageT_<Request, std::allocator<Request> > sendrequests;
            StorageT_<Status, std::allocator<Status> > sendstatus;
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              mirrors.at(i).template gather_dual<Algo_>(sendbufs.at(i), target);

              Request sr;
              Status ss;

              sendrequests.push_back(std::move(sr));
              sendstatus.push_back(std::move(ss));

              Comm::isend(sendbufs.at(i).elements(),
                          sendbufs.at(i).size(),
                          other_ranks.at(i),
                          sendrequests.at(i),
                          tags.at(i),
                          communicator
                          );
            }

            int* recvflags = new int[recvrequests.size()];
            int* taskflags = new int[recvrequests.size()];
            for(Index i(0) ; i < recvrequests.size() ; ++i)
            {
              recvflags[i] = 0;
              taskflags[i] = 0;
            }

            Index count(0);
            while(count != recvrequests.size())
            {
              //go through all requests round robin
              for(Index i(0) ; i < recvrequests.size() ; ++i)
              {
                if(taskflags[i] == 0)
                {
                  Comm::test(recvrequests.at(i), recvflags[i], recvstatus.at(i));
                  if(recvflags[i] != 0)
                  {
                    mirrors.at(i).template scatter_axpy_dual<Algo_>(target, recvbufs.at(i));
                    ++count;
                    taskflags[i] = 1;
                  }
                }
              }
            }

            for(Index i(0) ; i < sendrequests.size() ; ++i)
            {
              Status ws;
              Comm::wait(sendrequests.at(i), ws);
            }

            delete[] recvflags;
            delete[] taskflags;
          }
#else
          template<typename TargetVectorT_, typename VectorMirrorT_, typename BufferVectorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_&,
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

      template <typename Mem_, typename Algo_>
      struct GlobalSynchVec1
      {
      };

      template <typename Algo_>
      struct GlobalSynchVec1<Mem::Main, Algo_>
      {
        public:

#ifndef SERIAL
          ///TODO start all recvs in one sweep first
          template<typename TargetVectorT_, typename VectorMirrorT_, typename BufferVectorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_& target,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           const TargetVectorT_& frequencies,
                           StorageT_<Index, std::allocator<Index> > other_ranks,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >& sendbufs,
                           StorageT_<BufferVectorT_, std::allocator<BufferVectorT_> >& recvbufs,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Communicator communicator = Communicator(MPI_COMM_WORLD)
                           )
          {
            if(mirrors.size() == 0)
              return;

            ///assumes type-1 vector (full entries at inner boundaries)

            ///start recv
            StorageT_<Request, std::allocator<Request> > recvrequests;
            StorageT_<Status, std::allocator<Status> > recvstatus;
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              Request rr;
              Status rs;

              recvrequests.push_back(std::move(rr));
              recvstatus.push_back(std::move(rs));

              Comm::irecv(recvbufs.at(i).elements(),
                          recvbufs.at(i).size(),
                          other_ranks.at(i),
                          recvrequests.at(i),
                          tags.at(i),
                          communicator
                         );
            }

            ///gather and start send
            StorageT_<Request, std::allocator<Request> > sendrequests;
            StorageT_<Status, std::allocator<Status> > sendstatus;
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              mirrors.at(i).template gather_dual<Algo_>(sendbufs.at(i), target);

              Request sr;
              Status ss;

              sendrequests.push_back(std::move(sr));
              sendstatus.push_back(std::move(ss));

              Comm::isend(sendbufs.at(i).elements(),
                          sendbufs.at(i).size(),
                          other_ranks.at(i),
                          sendrequests.at(i),
                          tags.at(i),
                          communicator
                          );
            }

            int* recvflags = new int[recvrequests.size()];
            int* taskflags = new int[recvrequests.size()];
            for(Index i(0) ; i < recvrequests.size() ; ++i)
            {
              recvflags[i] = 0;
              taskflags[i] = 0;
            }


            Index count(0);

            ///handle receives
            while(count != recvrequests.size())
            {
              //go through all requests round robin
              for(Index i(0) ; i < recvrequests.size() ; ++i)
              {
                if(taskflags[i] == 0)
                {
                  Comm::test(recvrequests.at(i), recvflags[i], recvstatus.at(i));
                  if(recvflags[i] != 0)
                  {
                    mirrors.at(i).template scatter_axpy_dual<Algo_>(target, recvbufs.at(i));
                    ++count;
                    taskflags[i] = 1;
                  }
                }
              }
            }

            // scale target vector by frequencies
            target.template component_product<Algo_>(target, frequencies);

            for(Index i(0) ; i < sendrequests.size() ; ++i)
            {
              Status ws;
              Comm::wait(sendrequests.at(i), ws);
            }

            delete[] recvflags;
            delete[] taskflags;
          }
#else
          template<typename TargetVectorT_, typename VectorMirrorT_, typename BufferVectorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_&,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >&,
                           const TargetVectorT_&,
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
