#include <kernel/foundation/process.hpp>
#ifdef PARALLEL

namespace FEAST
{

  // define static member variables
  int Process::rank = MPI_PROC_NULL;
  int Process::rank_master = MPI_PROC_NULL;
  unsigned int Process::num_processes = 0;
  bool Process::is_master = false;

} // namespace FEAST

#endif // PARALLEL
