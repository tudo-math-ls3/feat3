#pragma once
#ifndef KERNEL_FOUNDATION_COMMUNICATION_HH
#define KERNEL_FOUNDATION_COMMUNICATION_HH 1

#ifndef FEAST_SERIAL_MODE
#include<mpi.h>
#endif
#include <kernel/archs.hpp>
#include <kernel/foundation/communication_error.hpp>

#include <vector>


namespace FEAST
{
  namespace Foundation
  {
    ///Communication modes used in Comm implementation or pass-over to backends
    enum CommModes
    {
      com_send = 0,
      com_receive,
      com_send_receive,
      com_average
        //TODO...
    };

    template<typename B_, CommModes c_>
    class Communicateable
    {
    };

    ///TODO send_recv->tparam CommOp_ ?
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

    template<typename T_, CommModes c_>
    class CommunicateableByAggregates
    {
    };

    ///implemented by Foundation datastructures that can be communicated but don't need to be buffered because all their aggregates already are
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

#ifndef FEAST_SERIAL_MODE
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

    template<typename TopologyType_>
      struct CommStructures
      {
        CommStructures(const TopologyType_& n, const TopologyType_& p) :
          network(n),
          patch_mesh(p),
          patch_process_map(TopologyType_())
        {
        }

        const TopologyType_& network;
        const TopologyType_& patch_mesh;
        TopologyType_ patch_process_map;
        TopologyType_ process_patch_map;
      };


    /**
     * \brief Communication implementation or backend pass-over
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
            static inline void send_recv(DataType1_ * sendbuf,
                                         Index num_elements_to_send,
                                         Index dest_rank,
                                         DataType2_* recvbuf,
                                         Index num_elements_to_recv,
                                         Index source_rank)
            {
              const Index send_end(num_elements_to_send);
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
              }
            }

          //TODO
      };

#ifndef FEAST_SERIAL_MODE
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
                                         Index source_rank)
            {
              MPI_Status status;

              MPI_Sendrecv(sendbuf,
                           num_elements_to_send,
                           MPIType<DataType1_>::value(),
                           dest_rank,
                           999,
                           recvbuf,
                           num_elements_to_recv,
                           MPIType<DataType2_>::value(),
                           source_rank,
                           999,
                           MPI_COMM_WORLD,
                           &status);
            }

          //TODO
      };
#endif

    /**
     * \brief Communication implementation or backend pass-over
     *
     * See specialisations.
     *
     * \tparam i_
     * overlap of the patches
     *
     * \tparam j_
     * communication mode
     *
     * \tparam AttributeType_
     * data type to be transferred
     *
     * \tparam Tag_
     * backend specifier
     *
     * \author Markus Geveler
     */
    /*template<unsigned i_ = 1, CommModes j_ = com_send_receive, typename AttributeType_ = double, typename Tag_ = Nil> //overlap is a compile-time decision now, if not feasible, move to inner function template
      class Communication
      {
        public:
          //example: Halo-based
          template<typename HaloType_, typename MeshType_>
            static void execute(
                HaloType_ & halo,
                unsigned attr_index,
                MeshType_ & other_mesh,
                Index other_rank) //TODO other_mesh resolved by process mesh list (mesh by id), needs to be a local thing
            {
#ifdef FOUNDATION_DEBUG
              if(i_ != halo.get_overlap())
                throw CommunicationHaloOverlapMismatch(i_, halo.get_overlap());
#endif
              //temp example
              switch(j_)
              {
                case com_send_receive:
                  {
                    //temp until some sort of patch-internal buffer or master bufferpool available. Emulates Pack(). if rank=.. states here
                    AttributeType_* sendbuf(new AttributeType_[halo.size()]);
                    AttributeType_* recvbuf(new AttributeType_[halo.size()]);
                    for(Index i(0) ; i < halo.size() ; ++i)
                    {
                      sendbuf[i] = ((Attribute<AttributeType_>*)(halo.get_mesh().get_attributes()->at(attr_index).get()))->get_data().at(halo.get_element(i));
                      recvbuf[i] = ((Attribute<AttributeType_>*)(other_mesh.get_attributes()->at(attr_index).get()))->get_data().at(halo.get_element_counterpart(i));
                    }
                    //'post'
                    Foundation::Comm<Tag_>::send_recv(sendbuf,
                                                      halo.size(),
                                                      other_rank,
                                                      recvbuf,
                                                      halo.size(),
                                                      halo.get_mesh().get_pp_rank());
                    for(Index i(0) ; i < halo.size() ; ++i)
                    {
                      ((Attribute<AttributeType_>*)(halo.get_mesh().get_attributes()->at(attr_index).get()))->get_data().at(halo.get_element(i)) = sendbuf[i];
                      ((Attribute<AttributeType_>*)(other_mesh.get_attributes()->at(attr_index).get()))->get_data().at(halo.get_element_counterpart(i)) = recvbuf[i];
                    }

                    delete[] sendbuf;
                    delete[] recvbuf;
                  }
              }
            }
          //TODO
      };*/
  }
}
#endif
