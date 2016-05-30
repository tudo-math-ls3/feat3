#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_SYNCH_VEC_HPP
#define FOUNDATION_GUARD_GLOBAL_SYNCH_VEC_HPP 1

#include<kernel/util/comm_base.hpp>
#include<kernel/util/time_stamp.hpp>
#include<kernel/util/statistics.hpp>

namespace FEAT
{
  namespace Global
  {
      /// \todo add communicators
      template <typename Mem_>
      struct SynchVec0
      {
        public:

#ifndef SERIAL
          ///TODO implement version with active queue of requests
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_& target,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           StorageT_<Index, std::allocator<Index> >& other_ranks,
                           StorageT_<LAFEM::DenseVector<Mem::Main, typename TargetVectorT_::DataType, typename TargetVectorT_::IndexType>, std::allocator<LAFEM::DenseVector<Mem::Main, typename TargetVectorT_::DataType, typename TargetVectorT_::IndexType> > >& sendbufs,
                           StorageT_<LAFEM::DenseVector<Mem::Main, typename TargetVectorT_::DataType, typename TargetVectorT_::IndexType>, std::allocator<LAFEM::DenseVector<Mem::Main, typename TargetVectorT_::DataType, typename TargetVectorT_::IndexType> > >& recvbufs,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Util::Communicator communicator = Util::Communicator(MPI_COMM_WORLD)
                           )
          {
            if(mirrors.size() == 0)
              return;

            TimeStamp ts_start;

            //start recv
            StorageT_<Util::CommRequest, std::allocator<Util::CommRequest> > recvrequests(mirrors.size());
            StorageT_<Util::CommStatus, std::allocator<Util::CommStatus> > recvstatus;
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              Util::CommStatus rs;

              recvstatus.push_back(std::move(rs));

              Util::Comm::irecv(recvbufs.at(i).elements(),
                          recvbufs.at(i).size(),
                          other_ranks.at(i),
                          recvrequests.at(i),
                          tags.at(i),
                          communicator
                         );

            }

            //gather and start send
            StorageT_<Util::CommRequest, std::allocator<Util::CommRequest> > sendrequests(mirrors.size());
            StorageT_<Util::CommStatus, std::allocator<Util::CommStatus> > sendstatus;
            for(Index i(0) ; i < mirrors.size() ; ++i)
            {
              mirrors.at(i).gather_dual(sendbufs.at(i), target);

              Util::CommStatus ss;

              sendstatus.push_back(std::move(ss));

              Util::Comm::isend(sendbufs.at(i).elements(),
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
            //handle receives
            while(count != recvrequests.size())
            {
              //go through all requests round robin
              for(Index i(0) ; i < recvrequests.size() ; ++i)
              {
                if(taskflags[i] == 0)
                {
                  Util::Comm::test(recvrequests.at(i), recvflags[i], recvstatus.at(i));
                  if(recvflags[i] != 0)
                  {
                    mirrors.at(i).scatter_axpy_dual(target, recvbufs.at(i));
                    ++count;
                    taskflags[i] = 1;
                  }
                }
              }
            }

            for(Index i(0) ; i < sendrequests.size() ; ++i)
            {
              Util::CommStatus ws;
              Util::Comm::wait(sendrequests.at(i), ws);
            }

            delete[] recvflags;
            delete[] taskflags;

            TimeStamp ts_stop;
            Statistics::add_time_mpi_execute(ts_stop.elapsed(ts_start));
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
                           Util::Communicator = Util::Communicator(0)
                           )
          {
          }
#endif
      };

      template <typename Mem_>
      struct SynchVec1
      {
        public:

#ifndef SERIAL
          template<typename TargetVectorT_, typename VectorMirrorT_, template<typename, typename> class StorageT_>
          static void exec(
                           TargetVectorT_& target,
                           const StorageT_<VectorMirrorT_, std::allocator<VectorMirrorT_> >& mirrors,
                           const TargetVectorT_& frequencies,
                           StorageT_<Index, std::allocator<Index> > other_ranks,
                           StorageT_<LAFEM::DenseVector<Mem::Main, typename TargetVectorT_::DataType, typename TargetVectorT_::IndexType>, std::allocator<LAFEM::DenseVector<Mem::Main, typename TargetVectorT_::DataType, typename TargetVectorT_::IndexType> > >& sendbufs,
                           StorageT_<LAFEM::DenseVector<Mem::Main, typename TargetVectorT_::DataType, typename TargetVectorT_::IndexType>, std::allocator<LAFEM::DenseVector<Mem::Main, typename TargetVectorT_::DataType, typename TargetVectorT_::IndexType> > >& recvbufs,
                           const StorageT_<Index, std::allocator<Index> >& tags,
                           Util::Communicator communicator = Util::Communicator(MPI_COMM_WORLD)
                           )
          {
            if(mirrors.size() == 0)
              return;

            SynchVec0<Mem_>::exec(target, mirrors, other_ranks, sendbufs, recvbufs, tags, communicator);

            // scale target vector by frequencies
            target.component_product(target, frequencies);
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
                           Util::Communicator = Util::Communicator(0)
                           )
          {
          }
#endif
      };
  }
}


#endif
