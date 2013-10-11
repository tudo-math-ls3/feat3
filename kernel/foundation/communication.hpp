#pragma once
#ifndef KERNEL_FOUNDATION_COMMUNICATION_HH
#define KERNEL_FOUNDATION_COMMUNICATION_HH 1

#ifndef SERIAL
#include<mpi.h>
#endif
#include <kernel/foundation/base.hpp>
#include <kernel/archs.hpp>
#include <kernel/foundation/communication_error.hpp>
#include <kernel/foundation/buffer.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

#include <vector>

using namespace FEAST::Archs;
using namespace FEAST::LAFEM;


namespace FEAST
{
  namespace Foundation
  {
    ///Communication modes used in Comm implementation or pass-over to backends
    enum Tier0CommModes
    {
      com_send = 0,
      com_receive,
      com_send_receive,
      com_gather,
      com_scatter,
      com_reduce_min
        //TODO...
    };

    enum Tier2CommModes
    {
      com_send_replace = 0,
      com_recv_replace,
      com_exchange,
      com_average,
      com_accumulate,
      com_allreduce_sqrtsum,
      com_min,
      com_max
    };

#ifndef SERIAL
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
#endif

    /**
     * \brief Tier-0 Communication implementation or backend pass-over
     *
     * Foundation Tier-0 Comm protocols involve simple pointers as buffers.
     *
     * See specialisations.
     *
     * \tparam Tag_
     * backend specifier
     *
     * \author Markus Geveler
     */
    template<typename Tag_>
      class Comm
      {
      };

    ///example shared-mem exchange
    template<>
      class Comm<Archs::Serial>
      {
        public:
          template<typename DataType1_, typename DataType2_>
            static inline void send_recv(DataType1_ * /*sendbuf*/,
                                         Index /*num_elements_to_send*/,
                                         Index /*dest_rank*/,
                                         DataType2_* /*recvbuf*/,
                                         Index /*num_elements_to_recv*/,
                                         Index /*source_rank*/,
                                         Index /*send_tag*/ = 0,
                                         Index /*recv_tag*/ = 0,
                                         Index /*communicator*/ = 0)
            {
              /*const Index send_end(num_elements_to_send);
              const Index recv_end(num_elements_to_recv);
              DataType1_ bufsend(0);
              DataType2_ bufrecv(0);

              for(Index i(0) ; i < send_end ; ++i)
              {
                bufsend = (DataType1_)recvbuf[i];
                recvbuf[i] = (DataType2_)sendbuf[i];
                recvbuf[i] = bufsend;
              }
              for(Index i(0) ; i < recv_end ; ++i)
              {
                bufrecv = (DataType2_)sendbuf[i];
                sendbuf[i] = (DataType1_)recvbuf[i];
                recvbuf[i] = bufrecv;
              }*/
            }

          template<typename DataType1_>
            static inline void allreduce(DataType1_ * /*sendbuf*/,
                                         Index /*num_elements_to_send_and_receive*/,
                                         DataType1_ * /*recvbuf*/)
            {
            }

          //TODO
      };

#ifndef SERIAL
    template<>
      class Comm<Archs::Parallel>
      {
        public:
          template<typename DataType1_, typename DataType2_>
            static inline void send_recv(DataType1_ * sendbuf,
                                         Index num_elements_to_send,
                                         Index dest_rank,
                                         DataType2_* recvbuf,
                                         Index num_elements_to_recv,
                                         Index source_rank,
                                         Index send_tag = 0,
                                         Index recv_tag = 0,
                                         MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Status status;

              MPI_Sendrecv(sendbuf,
                           (unsigned)num_elements_to_send,
                           MPIType<DataType1_>::value(),
                           (unsigned)dest_rank,
                           (unsigned)send_tag,
                           recvbuf,
                           (unsigned)num_elements_to_recv,
                           MPIType<DataType2_>::value(),
                           (unsigned)source_rank,
                           (unsigned)recv_tag,
                           communicator,
                           &status);
            }

          template<typename DataType_>
            static inline void send(DataType_ * sendbuf,
                Index num_elements_to_send,
                Index dest_rank,
                Index send_tag = 0,
                MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Send(sendbuf,
                           (unsigned)num_elements_to_send,
                           MPIType<DataType_>::value(),
                           (unsigned)dest_rank,
                           (unsigned)send_tag,
                           communicator);
            }

          template<typename DataType_>
            static inline void recv(DataType_ * recvbuf,
                Index num_elements_to_recv,
                Index src_rank,
                Index recv_tag = 0,
                MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Status status;

              MPI_Recv(recvbuf,
                           (unsigned)num_elements_to_recv,
                           MPIType<DataType_>::value(),
                           (unsigned)src_rank,
                           (unsigned)recv_tag,
                           communicator,
                           &status);
            }

          static inline void barrier(MPI_Comm communicator = MPI_COMM_WORLD)
          {
            MPI_Barrier(communicator);
          }

          template<typename DataType_>
            static inline void bcast(DataType_ * buf,
                Index num_elements,
                Index root,
                MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Bcast(buf, (unsigned)num_elements, MPIType<DataType_>::value(), (unsigned)root, communicator);
            }

          template<typename DataType1_, typename DataType2_>
            static inline void scatter(DataType1_ * sendbuf,
                Index num_elements_to_send,
                DataType2_ * recvbuf,
                Index num_elements_to_recv,
                Index root,
                MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Scatter(sendbuf, (unsigned)num_elements_to_send, MPIType<DataType1_>::value(), recvbuf, (unsigned)num_elements_to_recv,
                  MPIType<DataType2_>::value(), (unsigned)root, communicator);
            }

          template<typename DataType1_, typename DataType2_>
            static inline void gather(DataType1_ * sendbuf,
                Index num_elements_to_send,
                DataType2_ * recvbuf,
                Index num_elements_to_recv,
                Index root,
                MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Gather(sendbuf, (unsigned)num_elements_to_send, MPIType<DataType1_>::value(), recvbuf, (unsigned)num_elements_to_recv,
                  MPIType<DataType2_>::value(), (unsigned)root, communicator);
            }

          template<typename DataType_>
            static inline void reduce(DataType_ * sendbuf,
                DataType_ * recvbuf,
                Index num_elements_to_send,
                MPI_Op op,
                Index root,
                MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Reduce(sendbuf, recvbuf, (unsigned)num_elements_to_send, MPIType<DataType_>::value(), op, (unsigned)root, communicator);
            }

          template<typename DataType1_, typename DataType2_>
            static inline void allgather(DataType1_ * sendbuf,
                Index num_elements_to_send,
                DataType2_ * recvbuf,
                Index num_elements_to_recv,
                MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Allgather(sendbuf, (unsigned)num_elements_to_send, MPIType<DataType1_>::value(),
                  recvbuf, (unsigned)num_elements_to_recv, MPIType<DataType2_>::value(), communicator);
            }

          ///TODO delegate Op to an MPI_Op resolver as in MPI_Type resolver
          template<typename DataType1_>
            static inline void allreduce(DataType1_ * sendbuf,
                                         Index num_elements_to_send_and_receive,
                                         DataType1_ * recvbuf,
                                         MPI_Op op = MPI_SUM,
                                         MPI_Comm communicator = MPI_COMM_WORLD)
            {
              MPI_Allreduce(sendbuf,
                            recvbuf,
                            (unsigned)num_elements_to_send_and_receive,
                            MPIType<DataType1_>::value(),
                            op,
                            communicator);
            }

          //TODO
      };
#endif

    ///Tier-1 implementation: Foundation Tier-1 Comm protocols use Tier-0 protocols and define how Foundation datastructures can be communicated.
    template<typename B_, Tier0CommModes c_>
    class Communicateable
    {
    };

    ///implemented by Bufferable Foundation datastructures that can be communicated
    template<typename BufferType_>
    class Communicateable<BufferType_, com_send_receive>
    {
      public:
        virtual void send_recv(BufferType_& senddata,
                               int destrank,
                               BufferType_& recvdata,
                               int sourcerank) = 0;
    };

    template<typename T_, Tier0CommModes c_>
    class CommunicateableByAggregates
    {
    };

    ///inherited by Foundation datastructures that can be communicated but don't need to be buffered because all their aggregates already are
    template<typename T_>
    class CommunicateableByAggregates<T_, com_send_receive>
    {
      public:
        template<typename AggregateStorageType_>
        void send_recv(AggregateStorageType_& aggregates_to_communicate,
                       int destrank,
                       int sourcerank,
                       Index estimated_size_increase = 0)
        {
          for(Index i(0) ; i < aggregates_to_communicate.size() ; ++i)
          {
            typename T_::buffer_type_ sendbuf(aggregates_to_communicate.at(i).buffer());
            typename T_::buffer_type_ recvbuf(aggregates_to_communicate.at(i).buffer(estimated_size_increase));

            aggregates_to_communicate.at(i).to_buffer(sendbuf);

            aggregates_to_communicate.at(i).send_recv(sendbuf, destrank, recvbuf, sourcerank);

            aggregates_to_communicate.at(i).from_buffer(recvbuf);
          }
        }
    };

    /**
     * \brief Tier-2 Communication implementation or backend pass-over
     *
     * Foundation Tier-2 Comm protocols use Foundation datastructures to orchestrate communication.
     *
     * See specialisations.
     *
     * \author Markus Geveler
     */
    template<typename ContainerBackend_, Tier2CommModes op_>
    class InterfacedComm
    {
    };

    template<typename ContainerBackend_>
    class InterfacedComm<ContainerBackend_, com_exchange>
    {
      public:
        ///TODO think Halo<...pl_vertex...> implementations

        ///implementation for Foundation datastructures (Attribute<.>)
        template<
          unsigned delta_,
          typename a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
                   typename,
                   typename,
                   template<typename, typename> class>
           class HaloType_,
           typename AT_>
         static void execute(const HaloType_<delta_, a_, b_, c_>& interface, AT_& fct)
         {
           //acquire buffers
           std::shared_ptr<SharedArrayBase > sendbuf(BufferedSharedArray<typename AT_::data_type_>::create(interface.size()));
           std::shared_ptr<SharedArrayBase > recvbuf(BufferedSharedArray<typename AT_::data_type_>::create(interface.size()));

           //collect data
           for(Index i(0) ; i < interface.size() ; ++i)
           {
             (*((BufferedSharedArray<typename AT_::data_type_>*)(sendbuf.get())))[i] = fct.at(interface.get_element(i));
           }

           //post send_recv
           ///TODO validate via mesh reference, that polytope level is correct
#ifndef SERIAL
          Comm<Parallel>::send_recv(((BufferedSharedArray<typename AT_::data_type_>*)(sendbuf.get()))->get(),
              interface.size(),
              interface.get_other(),
              ((BufferedSharedArray<typename AT_::data_type_>*)(recvbuf.get()))->get(),
              interface.size(),
              interface.get_other());
#else
          Comm<Serial>::send_recv(((BufferedSharedArray<typename AT_::data_type_>*)(sendbuf.get()))->get(),
              interface.size(),
              interface.get_other(),
              ((BufferedSharedArray<typename AT_::data_type_>*)(recvbuf.get()))->get(),
              interface.size(),
              interface.get_other());
#endif
          //reset values
           for(Index i(0) ; i < interface.size() ; ++i)
           {
              fct.at(interface.get_element(i)) = (*((BufferedSharedArray<typename AT_::data_type_>*)(recvbuf.get())))[i];
           }

           //buffers are destroyed automatically
         }

        ///implementation for lafem dv
        template<
          unsigned delta_,
          typename a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
                   typename,
                   typename,
                   template<typename, typename> class>
           class HaloType_,
           typename DT_>
         static void execute(const HaloType_<delta_, a_, b_, c_>& interface, DenseVector<ContainerBackend_, DT_>& fct)
         {
           ///TODO: does assume, attribute for level a_ is stored in fct
           //acquire buffers
           DenseVector<ContainerBackend_, DT_> sendbuf(interface.size());
           DenseVector<ContainerBackend_, DT_> recvbuf(interface.size());

           //collect data
           for(Index i(0) ; i < interface.size() ; ++i)
           {
             sendbuf(i, fct(interface.get_element(i)));
           }

           //post send_recv
           ///TODO validate via mesh reference, that polytope level is correct
#ifndef SERIAL
          Comm<Parallel>::send_recv(sendbuf.elements(),
              interface.size(),
              interface.get_other(),
              recvbuf.elements(),
              interface.size(),
              interface.get_other());
#else
          Comm<Serial>::send_recv(sendbuf.elements(),
              interface.size(),
              interface.get_other(),
              recvbuf.elements(),
              interface.size(),
              interface.get_other());
#endif
          //reset values
           for(Index i(0) ; i < interface.size() ; ++i)
           {
              fct(interface.get_element(i), recvbuf(i));
           }
         }

        ///implementation for lafem smcsr
        template<
          unsigned delta_,
          typename a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
                   typename,
                   typename,
                   template<typename, typename> class>
           class HaloType_,
           typename DT_>
         static void execute(const HaloType_<delta_, a_, b_, c_>& interface, SparseMatrixCSR<ContainerBackend_, DT_>& mat)
         {
           DT_* val(mat.val());
           Index* row_ptr(mat.row_ptr());
           Index* row_ptr_end(mat.row_ptr_end());
           //Index* col_ind(mat.col_ind());

           ///TODO: does assume, matrix rows are associated with polytope level a_
           ///TODO: does assume, interface uses same index set (numbering!!) as matrix
           //directly replace row i
           for(Index i(0); i < interface.size() ; ++i)
           {
#ifndef SERIAL
             Comm<Parallel>::send_recv(&val[row_ptr[interface.get_element(i)]],
                 row_ptr_end[interface.get_element(i)] - row_ptr[interface.get_element(i)],
                 interface.get_other(),
                 &val[row_ptr[interface.get_element(i)]],
                 row_ptr_end[interface.get_element(i)] - row_ptr[interface.get_element(i)],
                 interface.get_other());
#else
             Comm<Serial>::send_recv(&val[row_ptr[interface.get_element(i)]],
                 row_ptr_end[interface.get_element(i)] - row_ptr[interface.get_element(i)],
                 interface.get_other(),
                 &val[row_ptr[interface.get_element(i)]],
                 row_ptr_end[interface.get_element(i)] - row_ptr[interface.get_element(i)],
                 interface.get_other());
#endif
           }
         }
    };
  }
}
#endif
