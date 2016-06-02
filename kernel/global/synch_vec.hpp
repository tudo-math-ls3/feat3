#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_SYNCH_VEC_HPP
#define FOUNDATION_GUARD_GLOBAL_SYNCH_VEC_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/util/comm_base.hpp>
#include<kernel/util/time_stamp.hpp>
#include<kernel/util/statistics.hpp>
#include<kernel/global/ticket.hpp>

namespace FEAT
{
  namespace Global
  {
      /// \todo add communicators
      struct SynchVec0Async
      {
        public:

#ifdef FEAT_HAVE_MPI
          ///TODO implement version with active queue of requests
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static std::shared_ptr<VecTicket<TargetVectorT_, StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_>>>> exec(
                           TargetVectorT_& target,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           const TargetVectorT_& frequencies,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Util::Communicator communicator = Util::Communicator(MPI_COMM_WORLD)
                           )
          {
            TimeStamp ts_start;

            auto ticket = std::make_shared<VecTicket<TargetVectorT_, StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_>>>> (target, frequencies, mirrors);

            //start recv
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              Util::Comm::irecv(ticket->recvbufs.at(i).elements(),
                          ticket->recvbufs.at(i).size(),
                          other_ranks.at(i),
                          *ticket->recvrequests.at(i),
                          tags.at(i),
                          communicator
                         );

            }

            //gather and start send
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              mirrors.at(i).gather_dual(ticket->sendbufs.at(i), target);

              Util::Comm::isend(ticket->sendbufs.at(i).elements(),
                          ticket->sendbufs.at(i).size(),
                          other_ranks.at(i),
                          *ticket->sendrequests.at(i),
                          tags.at(i),
                          communicator
                          );
            }

            TimeStamp ts_stop;
            Statistics::add_time_mpi_execute(ts_stop.elapsed(ts_start));

            return ticket;
          }
#else
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static std::shared_ptr<VecTicket<TargetVectorT_, StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_>>>> exec(
                           TargetVectorT_&,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >&,
                           const TargetVectorT_&,
                           StorageT_<Index, std::allocator<Index> >&,
                           const StorageT_<Index, std::allocator<Index> >&,
                           Util::Communicator = Util::Communicator(0)
                           )
          {
            auto ticket = std::make_shared<VecTicket<TargetVectorT_, StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_>>>> ();
            return ticket;
          }
#endif
      };

      struct SynchVec1Async
      {
        public:

#ifdef FEAT_HAVE_MPI
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static std::shared_ptr<VecTicket<TargetVectorT_, StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_>>>>  exec(
                           TargetVectorT_& target,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           const TargetVectorT_& frequencies,
                           StorageT_<Index, std::allocator<Index> > other_ranks,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Util::Communicator communicator = Util::Communicator(MPI_COMM_WORLD)
                           )
          {
            auto ticket = SynchVec0Async::exec(target, mirrors, frequencies, other_ranks, tags, communicator);
            ticket->scale = true;
            return ticket;
          }
#else
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static std::shared_ptr<VecTicket<TargetVectorT_, StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_>>>>  exec(
                           TargetVectorT_&,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >&,
                           const TargetVectorT_&,
                           StorageT_<Index, std::allocator<Index> >&,
                           const StorageT_<Index, std::allocator<Index> >&,
                           Util::Communicator = Util::Communicator(0)
                           )
          {
            auto ticket = std::make_shared<VecTicket<TargetVectorT_, StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_>>>> ();
            return ticket;
          }
#endif
      };

      /// \todo add communicators
      struct SynchVec0
      {
        public:
#ifdef FEAT_HAVE_MPI
          ///TODO implement version with active queue of requests
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_& target,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           const TargetVectorT_& frequencies,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Util::Communicator communicator = Util::Communicator(MPI_COMM_WORLD)
                           )
          {
            auto ticket = SynchVec0Async::exec(target, mirrors, frequencies, other_ranks, tags, communicator);
            ticket->wait();
          }
#else
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_&,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >&,
                           const TargetVectorT_&,
                           StorageT_<Index, std::allocator<Index> >&,
                           const StorageT_<Index, std::allocator<Index> >&,
                           Util::Communicator = Util::Communicator(0)
                           )
          {
          }
#endif
      };

      /// \todo add communicators
      struct SynchVec1
      {
        public:
#ifdef FEAT_HAVE_MPI
          ///TODO implement version with active queue of requests
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_& target,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           const TargetVectorT_& frequencies,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Util::Communicator communicator = Util::Communicator(MPI_COMM_WORLD)
                           )
          {
            auto ticket = SynchVec0Async::exec(target, mirrors, frequencies, other_ranks, tags, communicator);
            ticket->scale = true;
            ticket->wait();
          }
#else
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_&,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >&,
                           const TargetVectorT_&,
                           StorageT_<Index, std::allocator<Index> >&,
                           const StorageT_<Index, std::allocator<Index> >&,
                           Util::Communicator = Util::Communicator(0)
                           )
          {
          }
#endif
      };


  }
}


#endif
