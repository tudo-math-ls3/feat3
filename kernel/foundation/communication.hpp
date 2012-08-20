#pragma once
#ifndef KERNEL_FOUNDATION_COMMUNICATION_HH
#define KERNEL_FOUNDATION_COMMUNICATION_HH 1

#include <kernel/archs.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/foundation/attribute.hpp>
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
    class Comm<Archs::None>
    {
      public:
        template<typename DataType_>
        static inline void send_recv(DataType_ * sendbuf, Index /*dest_rank*/, Index num_elements, DataType_* recvbuf, Index /*source_rank*/)
          {
            DataType_ buf;
            const Index i_end(num_elements);
            for(Index i(0) ; i < i_end ; ++i)
            {
              buf = sendbuf[i];
              sendbuf[i] = recvbuf[i];
              recvbuf[i] = buf;
            }
          }

        //TODO
    };

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
    template<unsigned i_ = 1, CommModes j_ = com_send_receive, typename AttributeType_ = double, typename Tag_ = Nil> //overlap is a compile-time decision now, if not feasible, move to inner function template
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
                  Foundation::Comm<Tag_>::send_recv(sendbuf, other_rank, halo.size(), recvbuf, halo.get_mesh().get_pp_rank());
                  for(Index i(0) ; i < halo.size() ; ++i)
                  {
                    ((Attribute<AttributeType_>*)(halo.get_mesh().get_attributes()->at(attr_index).get()))->get_data().at(halo.get_element(i)) = sendbuf[i];
                    ((Attribute<AttributeType_>*)(other_mesh.get_attributes()->at(attr_index).get()))->get_data().at(halo.get_element_counterpart(i)) = recvbuf[i];
                  }
                }
            }
          }
        //TODO
    };
  }
}
#endif
