#pragma once
#ifndef KERNEL_FOUNDATION_COMM_BASE_HH
#define KERNEL_FOUNDATION_COMM_BASE_HH 1

#ifndef SERIAL
#include<mpi.h>
#include<memory>
#include<kernel/foundation/base.hpp>
#include<kernel/util/time_stamp.hpp>
#include<kernel/util/statistics.hpp>

namespace FEAST
{
  namespace Foundation
  {
    template <typename DT_>
      class MPIType
      {
      };

    template <>
      class MPIType<float>
      {
        public:
          static inline MPI_Datatype value()
          {
            return MPI_FLOAT;
          }
      };

    template <>
      class MPIType<double>
      {
        public:
          static inline MPI_Datatype value()
          {
            return MPI_DOUBLE;
          }
      };

    template <>
      class MPIType<unsigned long>
      {
        public:
          static inline MPI_Datatype value()
          {
            return MPI_UNSIGNED_LONG;
          }
      };

    template <>
      class MPIType<unsigned>
      {
        public:
          static inline MPI_Datatype value()
          {
            return MPI_UNSIGNED;
          }
      };

    template <>
      class MPIType<int>
      {
        public:
          static inline MPI_Datatype value()
          {
            return MPI_INT;
          }
      };

    template <>
      class MPIType<char>
      {
        public:
          static inline MPI_Datatype value()
          {
            return MPI_CHAR;
          }
      };

    class Communicator
    {
      public:
        Communicator(MPI_Comm comm) :
          _comm(comm)
      {
      }

      MPI_Comm mpi_comm()
      {
        return _comm;
      }

      MPI_Comm* get_mpi_comm()
      {
        return &_comm;
      }

      private:
        MPI_Comm _comm;
    };

    class Operation
    {
      public:
        Operation(MPI_Op op) :
          _op(op)
      {
      }

      MPI_Op mpi_op()
      {
        return _op;
      }

      private:
        MPI_Op _op;
    };

    class Request
    {
      public:
        Request() :
          _r(MPI_Request())
        {
        }

        MPI_Request& mpi_request()
        {
          return _r;
        }

        virtual ~Request()
        {
          if(_r != MPI_REQUEST_NULL)
          {
            MPI_Request_free(&_r);
          }
        }

      private:
        MPI_Request _r;
    };

    template<template<typename, typename> class ST_>
    class RequestSeq
    {
      public:
        RequestSeq(Index size) :
          _data(ST_<MPI_Request, std::allocator<MPI_Request> >(size))
        {
        }

        MPI_Request* mpi_requests()
        {
          return _data.data();
        }

        virtual ~RequestSeq()
        {
          for(auto& r : _data)
            if(r != MPI_REQUEST_NULL)
            {
              MPI_Request_free(&r);
            }
        }

      private:
        ST_<MPI_Request, std::allocator<MPI_Request> > _data;
    };

    class Status
    {
      public:
        Status() :
          _s(MPI_Status())
        {
        }

        MPI_Status& mpi_status()
        {
          return _s;
        }

      private:
        MPI_Status _s;
    };

    template<template<typename, typename> class ST_>
    class StatusSeq
    {
      public:
        StatusSeq(Index size) :
          _data(ST_<MPI_Status, std::allocator<MPI_Status> >(size))
        {
        }

        MPI_Status* mpi_statuses()
        {
          return _data.data();
        }

      private:
        ST_<MPI_Status, std::allocator<MPI_Status> > _data;
    };

      class Comm
      {
        public:
          template<typename DataType1_, typename DataType2_>
            static inline void send_recv(DataType1_ * sendbuf,
                                         Index num_elements_to_send,
                                         Index dest_rank,
                                         DataType2_* recvbuf,
                                         Index num_elements_to_recv,
                                         Index source_rank,
                                         Status& s,
                                         Index send_tag = 0,
                                         Index recv_tag = 0,
                                         Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Sendrecv(sendbuf,
                           (int)num_elements_to_send,
                           MPIType<DataType1_>::value(),
                           (int)dest_rank,
                           (int)send_tag,
                           recvbuf,
                           (int)num_elements_to_recv,
                           MPIType<DataType2_>::value(),
                           (int)source_rank,
                           (int)recv_tag,
                           communicator.mpi_comm(),
                           &(s.mpi_status()));
            }

          template<typename DataType_>
            static inline void send(DataType_ * sendbuf,
                                    Index num_elements_to_send,
                                    Index dest_rank,
                                    Index send_tag = 0,
                                    Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Send(sendbuf,
                        (int)num_elements_to_send,
                        MPIType<DataType_>::value(),
                        (int)dest_rank,
                        (int)send_tag,
                        communicator.mpi_comm());
            }

          template<typename DataType_>
          static inline void isend(DataType_ * sendbuf,
                                   Index num_elements_to_send,
                                   Index dest_rank,
                                   Request& r,
                                   Index send_tag = 0,
                                   Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Isend(sendbuf,
                        (int)num_elements_to_send,
                        MPIType<DataType_>::value(),
                        (int)dest_rank,
                        (int)send_tag,
                        communicator.mpi_comm(),
                        &(r.mpi_request()));
            }

          template<typename DataType_>
          static inline void irecv(DataType_ * recvbuf,
                                  Index num_elements_to_recv,
                                  Index src_rank,
                                  Request& r,
                                  Index recv_tag = 0,
                                  Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Irecv(recvbuf,
                        (int)num_elements_to_recv,
                        MPIType<DataType_>::value(),
                        (int)src_rank,
                        (int)recv_tag,
                        communicator.mpi_comm(),
                        &(r.mpi_request()));
            }

          template<typename DataType_>
            static inline void recv(DataType_ * recvbuf,
                Index num_elements_to_recv,
                Index src_rank,
                Status& s,
                Index recv_tag = 0,
                Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Recv(recvbuf,
                       (int)num_elements_to_recv,
                       MPIType<DataType_>::value(),
                       (int)src_rank,
                       (int)recv_tag,
                       communicator.mpi_comm(),
                       &(s.mpi_status()));
            }

          static inline void wait(Request& r, Status& s)
          {
            TimeStamp ts_start;
            MPI_Wait(&(r.mpi_request()), &(s.mpi_status()));
            TimeStamp ts_stop;
            Statistics::add_time_mpi_wait(ts_stop.elapsed(ts_start));
          }

          template<template<typename, typename> class ST_>
          static inline void waitall(RequestSeq<ST_>& r, StatusSeq<ST_>& s)
          {
            TimeStamp ts_start;
            MPI_Waitall(r.size(), r.mpi_requests(), s.mpi_statuses());
            TimeStamp ts_stop;
            Statistics::add_time_mpi_wait(ts_stop.elapsed(ts_start));
          }

          static inline void test(Request& r, int& flag, Status& s)
          {
            MPI_Test(&(r.mpi_request()), &flag, &(s.mpi_status()));
          }

          static inline void barrier(Communicator communicator = Communicator(MPI_COMM_WORLD))
          {
            TimeStamp ts_start;
            MPI_Barrier(communicator.mpi_comm());
            TimeStamp ts_stop;
            Statistics::add_time_mpi_wait(ts_stop.elapsed(ts_start));
          }

          template<typename DataType_>
            static inline void bcast(DataType_ * buf,
                Index num_elements,
                Index root,
                Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Bcast(buf, (int)num_elements, MPIType<DataType_>::value(), (int)root, communicator.mpi_comm());
            }

          template<typename DataType1_, typename DataType2_>
            static inline void scatter(DataType1_ * sendbuf,
                Index num_elements_to_send,
                DataType2_ * recvbuf,
                Index num_elements_to_recv,
                Index root,
                Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Scatter(sendbuf, (int)num_elements_to_send, MPIType<DataType1_>::value(), recvbuf, (int)num_elements_to_recv,
                  MPIType<DataType2_>::value(), (int)root, communicator.mpi_comm());
            }

          template<typename DataType1_, typename DataType2_>
            static inline void gather(DataType1_ * sendbuf,
                Index num_elements_to_send,
                DataType2_ * recvbuf,
                Index num_elements_to_recv,
                Index root,
                Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Gather(sendbuf, (int)num_elements_to_send, MPIType<DataType1_>::value(), recvbuf, (int)num_elements_to_recv,
                  MPIType<DataType2_>::value(), (int)root, communicator.mpi_comm());
            }

          template<typename DataType_>
            static inline void reduce(DataType_ * sendbuf,
                DataType_ * recvbuf,
                Index num_elements_to_send,
                Operation op,
                Index root,
                Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Reduce(sendbuf, recvbuf, (int)num_elements_to_send, MPIType<DataType_>::value(), op.mpi_op(), (int)root, communicator.mpi_comm());
            }

          template<typename DataType1_, typename DataType2_>
            static inline void allgather(DataType1_ * sendbuf,
                Index num_elements_to_send,
                DataType2_ * recvbuf,
                Index num_elements_to_recv,
                Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Allgather(sendbuf, (int)num_elements_to_send, MPIType<DataType1_>::value(),
                  recvbuf, (int)num_elements_to_recv, MPIType<DataType2_>::value(), communicator.mpi_comm());
            }

          ///TODO delegate Op to an MPI_Op resolver as in MPI_Type resolver
          template<typename DataType1_>
            static inline void allreduce(DataType1_ * sendbuf,
                                         Index num_elements_to_send_and_receive,
                                         DataType1_ * recvbuf,
                                         Operation op = Operation(MPI_SUM),
                                         Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Allreduce(sendbuf,
                            recvbuf,
                            (int)num_elements_to_send_and_receive,
                            MPIType<DataType1_>::value(),
                            op.mpi_op(),
                            communicator.mpi_comm());
            }

            template<typename DataType_>
            static inline int alltoallv(const DataType_ * sendbuf,
                                        const int sendcounts[],
                                        const int sdispls[],
                                        DataType_* recvbuf,
                                        const int recvcounts[],
                                        const int rdispls[],
                                        Communicator communicator = Communicator(MPI_COMM_WORLD),
                                        bool in_place = false)
            {
              if(!in_place)
              {
                return MPI_Alltoallv(sendbuf,
                                     sendcounts,
                                     sdispls,
                                     MPIType<DataType_>::value(),
                                     recvbuf,
                                     recvcounts,
                                     rdispls,
                                     MPIType<DataType_>::value(),
                                     communicator.mpi_comm());
              }
              else
              {
                return MPI_Alltoallv(MPI_IN_PLACE,
                                     sendcounts,
                                     sdispls,
                                     MPIType<DataType_>::value(),
                                     recvbuf,
                                     recvcounts,
                                     rdispls,
                                     MPIType<DataType_>::value(),
                                     communicator.mpi_comm());
              }
            }
          template<typename DataType1_>
            static inline void iallreduce(DataType1_ * sendbuf,
                                          Index num_elements_to_send_and_receive,
                                          DataType1_ * recvbuf,
                                          Request& r,
                                          Operation op = Operation(MPI_SUM),
                                          Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Iallreduce(sendbuf,
                            recvbuf,
                            (int)num_elements_to_send_and_receive,
                            MPIType<DataType1_>::value(),
                            op.mpi_op(),
                            communicator.mpi_comm(),
                            &(r.mpi_request()));
            }

#ifndef MSMPI_VER
          template<typename DataType_>
            static inline void ialltoallv(const DataType_ * sendbuf,
                                          const int sendcounts[],
                                          const int sdispls[],
                                          DataType_* recvbuf,
                                          const int recvcounts[],
                                          const int rdispls[],
                                          Request& r,
                                          Communicator communicator = Communicator(MPI_COMM_WORLD))
            {
              MPI_Ialltoallv(sendbuf,
                             sendcounts,
                             sdispls,
                             MPIType<DataType_>::value(),
                             recvbuf,
                             recvcounts,
                             rdispls,
                             MPIType<DataType_>::value(),
                             communicator.mpi_comm(),
                             &(r.mpi_request()));
            }
#endif

          static inline Index rank(Communicator c = Communicator(MPI_COMM_WORLD))
          {
            int r;

            MPI_Comm_rank(c.mpi_comm(), &r);

            return Index(r);
          }

          static inline Index size(Communicator c = Communicator(MPI_COMM_WORLD))
          {
            int r;

            MPI_Comm_size(c.mpi_comm(), &r);

            return Index(r);
          }

          //TODO
      };
  }
}
#else //SERIAL MODE
namespace FEAST
{
  namespace Foundation
  {
    class Communicator
    {
      public:
        explicit Communicator(Index comm) :
          _comm(comm)
      {
      }

      private:
        Index _comm;
    };

    class Operation
    {
      public:
        explicit Operation(Index op) :
          _op(op)
      {
      }

      private:
        Index _op;
    };

    class Request
    {
      public:
        Request() :
          _r(0)
        {
        }

        Index mpi_request()
        {
          return _r;
        }

      private:
        Index _r;
    };

    class Status
    {
      public:
        Status() :
          _s(0)
        {
        }

        Index mpi_status()
        {
          return _s;
        }

      private:
        Index _s;
    };

    class Comm
    {
      public:
        template<typename DataType1_, typename DataType2_>
          static inline void send_recv(DataType1_ *,
              Index,
              Index,
              DataType2_*,
              Index,
              Index,
              Status,
              Index = 0,
              Communicator = Communicator(0))
          {
          }

        template<typename DataType_>
          static inline void send(DataType_,
              Index,
              Index,
              Index,
              Communicator)
          {
          }

        template<typename DataType_>
          static inline void recv(DataType_ *,
              Index,
              Index,
              Index,
              Communicator)
          {
          }

        static inline void barrier(Communicator = Communicator(0))
        {
        }

        template<typename DataType_>
          static inline void bcast(DataType_*,
              Index,
              Index,
              Communicator)
          {
          }

        template<typename DataType1_, typename DataType2_>
          static inline void scatter(DataType1_*,
              Index,
              DataType2_*,
              Index,
              Index,
              Communicator)
          {
          }

        template<typename DataType1_, typename DataType2_>
          static inline void gather(DataType1_*,
              Index,
              DataType2_*,
              Index,
              Index,
              Communicator)
          {
          }

        template<typename DataType_>
          static inline void reduce(DataType_*,
              DataType_*,
              Index,
              Operation,
              Index,
              Communicator)
          {
          }

        template<typename DataType1_, typename DataType2_>
          static inline void allgather(DataType1_*,
              Index,
              DataType2_*,
              Index,
              Communicator)
          {
          }

        ///TODO delegate Op to an MPI_Op resolver as in MPI_Type resolver
        template<typename DataType1_>
          static inline void allreduce(DataType1_*,
              Index,
              DataType1_*,
              Operation = Operation(0),
              Communicator = Communicator(0))
          {
          }

          template<typename DataType1_>
            static inline void iallreduce(DataType1_ *,
                                          Index,
                                          DataType1_ *,
                                          Request&,
                                          Operation,
                                          Communicator)
            {
            }

          static inline void wait(Request&, Status&)
          {
          }

          static inline void test(Request&, int&, Status&)
          {
          }

          static inline Index rank(Communicator = Communicator(0))
          {
            return Index(0);
          }

          static inline Index size(Communicator = Communicator(0))
          {
            return Index(1);
          }
        //TODO
    };
  }
}
#endif//ifndef SERIAL
#endif //GUARD
