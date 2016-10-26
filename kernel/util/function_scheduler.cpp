#include <kernel/util/function_scheduler.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/dist.hpp>

using namespace FEAT;
using namespace FEAT::Util;

void FEAT::Util::schedule_function(std::function<void (void)> func, ScheduleMode mode)
{
  Dist::Comm comm(Dist::Comm::world());

  int rank(comm.rank());
  int ranks(comm.size());

  switch (mode)
  {
    case ScheduleMode::sequential:
      {
        for (int i(0) ; i < ranks ; ++i)
        {
          if (i == rank)
          {
            func();
            comm.barrier();
          }
          else
          {
            comm.barrier();
          }
        }
      }
      break;

    case ScheduleMode::clustered:
      {
        int concurrent = 100;  /// \todo replace hardcoded number by entry from feat.ini
        int sweeps = (ranks / concurrent) + 1;

        for (int i(0) ; i < sweeps && i < ranks ; ++i)
        {
          if (rank >= i * concurrent && rank < (i+1) * concurrent)
          {
            func();
            comm.barrier();
          }
          else
          {
            comm.barrier();
          }
        }
      }
      break;

    default:
      throw InternalError(__func__, __FILE__, __LINE__, "ScheduleMode not supported!");
  }
}
