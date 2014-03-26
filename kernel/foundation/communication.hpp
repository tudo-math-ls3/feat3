#pragma once
#ifndef KERNEL_FOUNDATION_COMMUNICATION_HH
#define KERNEL_FOUNDATION_COMMUNICATION_HH 1

#include<kernel/foundation/comm_base.hpp>

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
                               Index destrank,
                               BufferType_& recvdata,
                               Index sourcerank) = 0;
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
                       Index destrank,
                       Index sourcerank,
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
                   typename,
                   template<typename, typename> class>
           class HaloType_,
           typename AT_>
         static void execute(const HaloType_<delta_, a_, b_, typename AT_::data_type_, c_>& interface, AT_& fct)
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
           Status s;
           Comm::send_recv(((BufferedSharedArray<typename AT_::data_type_>*)(sendbuf.get()))->get(),
              interface.size(),
              interface.get_other(),
              ((BufferedSharedArray<typename AT_::data_type_>*)(recvbuf.get()))->get(),
              interface.size(),
              interface.get_other(),
              s);

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
                   typename,
                   template<typename, typename> class>
           class HaloType_,
           typename DT_ = double>
         static void execute(const HaloType_<delta_, a_, b_, DT_, c_>& interface, DenseVector<ContainerBackend_, DT_>& fct)
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
           Status s;
           Comm::send_recv(sendbuf.elements(),
              interface.size(),
              interface.get_other(),
              recvbuf.elements(),
              interface.size(),
              interface.get_other(),
              s);

          //reset values
           for(Index i(0) ; i < interface.size() ; ++i)
           {
              fct(interface.get_element(i), recvbuf(i));
           }
         }

        ///implementation for lafem smcsr
        //TODO fails since openmpi v.1.8.0
        /*template<
          unsigned delta_,
          typename a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
                   typename,
                   typename,
                   typename,
                   template<typename, typename> class>
           class HaloType_,
           typename DT_ = double>
         static void execute(const HaloType_<delta_, a_, b_, DT_, c_>& interface, SparseMatrixCSR<ContainerBackend_, DT_>& mat)
         {
           DT_* val(mat.val());
           Index* row_ptr(mat.row_ptr());
           //Index* col_ind(mat.col_ind());

           ///TODO: does assume, matrix rows are associated with polytope level a_
           ///TODO: does assume, interface uses same index set (numbering!!) as matrix
           //directly replace row i
           for(Index i(0); i < interface.size() ; ++i)
           {
             Status s;
             Comm::send_recv(&val[row_ptr[interface.get_element(i)]],
                 row_ptr[interface.get_element(i) + 1] - row_ptr[interface.get_element(i)],
                 interface.get_other(),
                 &val[row_ptr[interface.get_element(i)]],
                 row_ptr[interface.get_element(i) + 1] - row_ptr[interface.get_element(i)],
                 interface.get_other(),
                 s,
                 i,
                 i);
           }
         }*/
    };
  }
}
#endif
