#include<kernel/util/function_scheduler.hpp>

using namespace FEAST;
using namespace FEAST::Util;

void FEAST::Util::schedule_function(std::function<void (void)> func, ScheduleMode mode)
{
  Index rank(Comm::rank());
  Index ranks(Comm::size());

  switch (mode)
  {
    case ScheduleMode::sequential:
      {
        for (Index i(0) ; i < ranks ; ++i)
        {
          if (i == rank)
          {
            func();
            Comm::barrier();
          }
          else
          {
            Comm::barrier();
          }
        }
      }
      break;

    case ScheduleMode::clustered:
      {
        Index concurrent = 100;  /// \todo replace hardcoded number by entry from feast.ini
        Index sweeps = (ranks / concurrent) + 1;

        for (Index i(0) ; i < sweeps ; ++i)
        {
          if (rank >= i * concurrent && rank < (i+1) * concurrent)
          {
            func();
            Comm::barrier();
          }
          else
          {
            Comm::barrier();
          }
        }
      }
      break;

    default:
      throw InternalError(__func__, __FILE__, __LINE__, "ScheduleMode not supported!");
  }
}
