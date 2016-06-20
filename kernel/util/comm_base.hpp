#pragma once
#ifndef KERNEL_UTIL_COMM_BASE_HPP
#define KERNEL_UTIL_COMM_BASE_HPP 1

#include<kernel/base_header.hpp>

#include<iostream> // for std::istream
#include<cstring>  // for std::strcpy

#ifdef FEAT_HAVE_MPI
#include<mpi.h>
#include<memory>
#include<kernel/util/time_stamp.hpp>
#include<kernel/util/statistics.hpp>

namespace FEAT
{
  namespace Util
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

        Index size()
        {
          int r;

          MPI_Comm_size(_comm, &r);

          return Index(r);
        }


      private:
        MPI_Comm _comm;
    };

    class CommOperation
    {
      public:
        CommOperation(MPI_Op op) :
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

    class CommOperationSum : public CommOperation
    {
      public:
        CommOperationSum() :
          CommOperation(MPI_SUM)
      {
      }

        MPI_Op mpi_op()
        {
          return _op;
        }

      private:
        MPI_Op _op;
    };

    class CommOperationMax : public CommOperation
    {
      public:
        CommOperationMax() :
          CommOperation(MPI_MAX)
      {
      }

        MPI_Op mpi_op()
        {
          return _op;
        }

      private:
        MPI_Op _op;
    };

    class CommRequest
    {
      public:
        CommRequest() :
          _r(MPI_Request())
      {
      }

        MPI_Request& mpi_request()
        {
          return _r;
        }

        virtual ~CommRequest()
        {
          if(_r != MPI_REQUEST_NULL)
          {
            MPI_Request_free(&_r);
          }
        }

        /// Unwanted copy constructor: Do not implement!
        CommRequest(const CommRequest &) = delete;
        /// Unwanted copy assignment operator: Do not implement!
        CommRequest & operator=(const CommRequest &) = delete;

      private:
        MPI_Request _r;
    };

    template<template<typename, typename> class ST_>
      class CommRequestSeq
      {
        public:
          CommRequestSeq(Index size) :
            _data(ST_<MPI_Request, std::allocator<MPI_Request> >(size))
        {
        }

          MPI_Request* mpi_requests()
          {
            return _data.data();
          }

          virtual ~CommRequestSeq()
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

    class CommStatus
    {
      public:
        CommStatus() :
          _s(MPI_Status())
      {
      }

        MPI_Status * mpi_status()
        {
          return &_s;
        }

      private:
        MPI_Status _s;
    };

    class CommStatusIgnore : public CommStatus
    {
      public:
        CommStatusIgnore()
      {
      }

        MPI_Status * mpi_status()
        {
          return MPI_STATUSES_IGNORE;
        }
    };

    template<template<typename, typename> class ST_>
      class CommStatusSeq
      {
        public:
          CommStatusSeq(Index size) :
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
              CommStatus& s,
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
                s.mpi_status());
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
              CommRequest& r,
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
              CommRequest& r,
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
              CommStatus& s,
              Index recv_tag = 0,
              Communicator communicator = Communicator(MPI_COMM_WORLD))
          {
            MPI_Recv(recvbuf,
                (int)num_elements_to_recv,
                MPIType<DataType_>::value(),
                (int)src_rank,
                (int)recv_tag,
                communicator.mpi_comm(),
                s.mpi_status());
          }

        static inline void wait(CommRequest& r, CommStatus& s)
        {
          TimeStamp ts_start;
          MPI_Wait(&(r.mpi_request()), s.mpi_status());
          TimeStamp ts_stop;
        }

        template<template<typename, typename> class ST_>
          static inline void waitall(CommRequestSeq<ST_>& r, CommStatusSeq<ST_>& s)
          {
            TimeStamp ts_start;
            MPI_Waitall(r.size(), r.mpi_requests(), s.mpi_statuses());
            TimeStamp ts_stop;
          }

        static inline void test(CommRequest& r, int& flag, CommStatus& s)
        {
          MPI_Test(&(r.mpi_request()), &flag, s.mpi_status());
        }

        static inline void barrier(Communicator communicator = Communicator(MPI_COMM_WORLD))
        {
          TimeStamp ts_start;
          MPI_Barrier(communicator.mpi_comm());
          TimeStamp ts_stop;
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
              CommOperation op,
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
              CommOperation op = CommOperationSum(),
              Communicator communicator = Communicator(MPI_COMM_WORLD))
          {
            if (sendbuf == recvbuf)
            {
              MPI_Allreduce(MPI_IN_PLACE,
              recvbuf,
              (int)num_elements_to_send_and_receive,
              MPIType<DataType1_>::value(),
              op.mpi_op(),
              communicator.mpi_comm());
            }
            else
            {
              MPI_Allreduce(sendbuf,
              recvbuf,
              (int)num_elements_to_send_and_receive,
              MPIType<DataType1_>::value(),
              op.mpi_op(),
              communicator.mpi_comm());
            }
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
        template<typename DataType_>
          static inline int allgatherv(const DataType_ * sendbuf,
              const int sendcount,
              DataType_* recvbuf,
              const int recvcounts[],
              const int rdispls[],
              Communicator communicator = Communicator(MPI_COMM_WORLD),
              bool in_place = false)
          {
            if(!in_place)
            {
              return MPI_Allgatherv(sendbuf,
                  sendcount,
                  MPIType<DataType_>::value(),
                  recvbuf,
                  recvcounts,
                  rdispls,
                  MPIType<DataType_>::value(),
                  communicator.mpi_comm());
            }
            else
            {
              return MPI_Allgatherv(MPI_IN_PLACE,
                  sendcount,
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
              CommRequest& r,
              CommOperation op = CommOperationSum(),
              Communicator communicator = Communicator(MPI_COMM_WORLD))
          {
            if (sendbuf == recvbuf)
            {
              MPI_Iallreduce(MPI_IN_PLACE,
              recvbuf,
              (int)num_elements_to_send_and_receive,
              MPIType<DataType1_>::value(),
              op.mpi_op(),
              communicator.mpi_comm(),
              &(r.mpi_request()));
            }
            else
            {
              MPI_Iallreduce(sendbuf,
              recvbuf,
              (int)num_elements_to_send_and_receive,
              MPIType<DataType1_>::value(),
              op.mpi_op(),
              communicator.mpi_comm(),
              &(r.mpi_request()));
            }
          }

#ifndef MSMPI_VER
        template<typename DataType_>
          static inline void ialltoallv(const DataType_ * sendbuf,
              const int sendcounts[],
              const int sdispls[],
              DataType_* recvbuf,
              const int recvcounts[],
              const int rdispls[],
              CommRequest& r,
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

        /**
         * \brief Synchronises a stringstream from rank 0 to all ranks in the communicator
         *
         * \param[in,out] iss
         * The std::stringstream to synchronise, gets overwritten by the output.
         *
         * \param[in] comm
         * Communicator for the synch operation, defaults to MPI_COMM_WORLD.
         *
         */
        static void synch_stringstream(std::stringstream& iss, const Communicator comm = Communicator(MPI_COMM_WORLD))
        {
          Index my_rank(Util::Comm::rank(comm));
          Index size;
          String str;

          // Get size to broadcast
          if(my_rank == 0)
          {
            str = (iss.str());
            size = Index(str.length());
          }
          // synchronize length
          Util::Comm::bcast(&size, 1, 0, comm);

          if(size > Index(0))
          {
            // allocate
            char* buf = new char[size + 1];

            //fill
            if(my_rank == 0) //master
            {
              std::strcpy(buf, str.c_str());
            }

            //bcast data
            Util::Comm::bcast(buf, size, 0, comm);

            //convert
            if(my_rank != 0)
            {
              String res_str(buf, size);
              iss << res_str;
            }

            delete[] buf;
          }
        }

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
namespace FEAT
{
  namespace Util
  {
    class Communicator
    {
      public:
        explicit Communicator(Index comm) :
          _comm(comm)
        {
        }

        Index size()
        {
          return Index(1);
        }

      private:
        Index _comm;
    };

    class CommOperation
    {
      public:
        explicit CommOperation(Index op) :
          _op(op)
      {
      }

      private:
        Index _op;
    };

    class CommOperationSum : public CommOperation
    {
      public:
        CommOperationSum() :
          CommOperation(0)
      {
      }
    };

    class CommOperationMax : public CommOperation
    {
      public:
        CommOperationMax() :
          CommOperation(0)
      {
      }
    };

    class CommRequest
    {
      public:
        CommRequest() :
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

    class CommStatus
    {
      public:
        CommStatus() :
          _s(0)
      {
      }

        Index * mpi_status()
        {
          return &_s;
        }

      private:
        Index _s;
    };

    class CommStatusIgnore : public CommStatus
    {
      public:
        CommStatusIgnore() :
          _s(0)
      {
      }

        Index * mpi_status()
        {
          return &_s;
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
              CommStatus,
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
              CommOperation,
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
              CommOperation = CommOperation(0),
              Communicator = Communicator(0))
          {
          }

        template<typename DataType1_>
          static inline void iallreduce(DataType1_ *,
              Index,
              DataType1_ *,
              CommRequest&,
              CommOperation,
              Communicator)
          {
          }

        static inline void wait(CommRequest&, CommStatus&)
        {
        }

        static inline void test(CommRequest&, int&, CommStatus&)
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

        static void synch_stringstream(std::stringstream& DOXY(iss), Communicator DOXY(comm) = Communicator(0))
        {
        }
        //TODO
    };
  }
}
#endif//ifndef SERIAL
#endif //GUARD
