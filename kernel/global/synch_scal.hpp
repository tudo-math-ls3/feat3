#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_SYNCH_SCAL_HPP
#define FOUNDATION_GUARD_GLOBAL_SYNCH_SCAL_HPP 1

#include<kernel/util/comm_base.hpp>
#include<kernel/lafem/arch/scale.hpp>
#include<kernel/lafem/arch/component_product.hpp>
#include<kernel/lafem/arch/dot_product.hpp>
#include<kernel/util/time_stamp.hpp>
#include<kernel/util/statistics.hpp>
#include<kernel/global/ticket.hpp>

namespace FEAT
{
  namespace Global
  {
      /// \todo add communicators
      struct SynchScal0Async
      {
        public:
#ifdef FEAT_HAVE_MPI
          template<typename DT_>
          static std::shared_ptr<ScalTicket<DT_>> value(DT_ & x, Util::CommOperation op = Util::CommOperationSum())
          {
            TimeStamp ts_start;

            auto ticket = std::make_shared<ScalTicket<DT_>>(x);

            Util::Comm::iallreduce(&(ticket->x), &(ticket->r), 1, *ticket->req, op);

            TimeStamp ts_stop;
            Statistics::add_time_mpi_execute(ts_stop.elapsed(ts_start));

            return ticket;
          }
#else
          template<typename DT_>
          static std::shared_ptr<ScalTicket<DT_>> value(DT_& x, Util::CommOperation /*op*/ = Util::CommOperationSum())
          {
            auto ticket = std::make_shared<ScalTicket<DT_>>(x);
            ticket->r = x;
            return ticket;
          }
#endif
      };

      /// \todo add communicators
      struct SynchScal0
      {
        public:
#ifdef FEAT_HAVE_MPI
          template<typename DT_>
          static DT_ value(DT_& x, Util::CommOperation op = Util::CommOperationSum())
          {
            return SynchScal0Async::value(x, op)->wait();
          }
#else
          template<typename DT_>
          static DT_ value(DT_& x, Util::CommOperation /*op*/ = Util::CommOperationSum())
          {
            return x;
          }
#endif
      };
  }
}


#endif
