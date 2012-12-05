#pragma once
#ifndef SCARC_GUARD_SYNCH_HH
#define SCARC_GUARD_SYNCH_HH 1

#include<kernel/foundation/communication.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;

namespace FEAST
{
  namespace Foundation
  {
    template<typename Arch_, Tier2CommModes cm_>
    struct Synch
    {
    };

    template<typename Arch_>
    struct Synch<Arch_, com_exchange>
    {
      template<typename VectorT_, typename VectorMirrorT_>
      static void execute(VectorT_& target,
                          const VectorMirrorT_& mirror,
                          VectorT_& sendbuf,
                          VectorT_& recvbuf,
                          Index dest_rank,
                          Index source_rank)
      {
        mirror.gather_dual(sendbuf, target);

        Comm<Arch_>::send_recv(sendbuf.elements(),
                               sendbuf.size(),
                               dest_rank,
                               recvbuf.elements(),
                               recvbuf.size(),
                               source_rank);

        mirror.scatter_dual(target, recvbuf);
      }
    };
  }
}

#endif
