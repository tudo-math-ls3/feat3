#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_SYNCH_SCAL_HPP
#define FOUNDATION_GUARD_GLOBAL_SYNCH_SCAL_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/lafem/arch/scale.hpp>
#include<kernel/lafem/arch/component_product.hpp>
#include<kernel/lafem/arch/dot_product.hpp>
#include<kernel/util/time_stamp.hpp>
#include<kernel/util/statistics.hpp>

namespace FEAST
{
  namespace Foundation
  {
      /// \todo add communicators
      template <typename Mem_>
      struct GlobalSynchScal0
      {
        public:
#ifndef SERIAL
          template<typename DT_>
          static DT_ value(DT_& r,
                           DT_& x)
          {
            TimeStamp ts_start;

            DT_ sendbuf(x), recvbuf;

            Status stat;

// The current MS-MPI implementation does not offer the MPI_Iallreduce function...
/// \todo Remove this workaround once MS does its homework.
#ifdef MSMPI_VER
            Comm::allreduce(&sendbuf, Index(1), &recvbuf);
#else
            Request req;

            Comm::iallreduce(&sendbuf, Index(1), &recvbuf, req);

            Comm::wait(req, stat);
#endif // MSMPI_VER
            r = recvbuf;

            TimeStamp ts_stop;
            Statistics::add_time_mpi_execute(ts_stop.elapsed(ts_start));
            return r;
          }
#else
          template<typename DT_>
          static DT_ value(DT_& r,
                           DT_& x)
          {
            r = x;
            return r;
          }
#endif
      };
  }
}


#endif
