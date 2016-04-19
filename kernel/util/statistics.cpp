#include <kernel/util/statistics.hpp>
#include <kernel/util/function_scheduler.hpp>

using namespace FEAST;

// static member initialisation
Index Statistics::_flops = Index(0);
KahanAccumulation Statistics::_time_reduction;
KahanAccumulation Statistics::_time_spmv;
KahanAccumulation Statistics::_time_axpy;
KahanAccumulation Statistics::_time_precon;
KahanAccumulation Statistics::_time_mpi_execute;
KahanAccumulation Statistics::_time_mpi_wait;
std::map<FEAST::String, SolverStatistics> Statistics::_solver_statistics;

void Statistics::write_out_solver_statistics_scheduled(Index rank, size_t la_bytes, size_t domain_bytes, size_t mpi_bytes, String filename)
{
  auto func = [&] () { write_out_solver_statistics(rank, la_bytes, domain_bytes, mpi_bytes, filename); };
  Util::schedule_function(func, Util::ScheduleMode::clustered);
}
