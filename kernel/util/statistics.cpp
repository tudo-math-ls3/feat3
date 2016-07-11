#include <kernel/util/statistics.hpp>
#include <kernel/util/function_scheduler.hpp>

using namespace FEAT;

// static member initialisation
Index Statistics::_flops = Index(0);
KahanAccumulation Statistics::_time_reduction;
KahanAccumulation Statistics::_time_spmv;
KahanAccumulation Statistics::_time_axpy;
KahanAccumulation Statistics::_time_precon;
KahanAccumulation Statistics::_time_mpi_execute;
KahanAccumulation Statistics::_time_mpi_wait;
std::map<FEAT::String, SolverStatistics> Statistics::_solver_statistics;
double Statistics::toe_partition;
double Statistics::toe_assembly;
double Statistics::toe_solve;

void Statistics::write_out_solver_statistics_scheduled(Index rank, size_t la_bytes, size_t domain_bytes, size_t mpi_bytes, Index cells, Index dofs, Index nzes, String filename)
{
  auto func = [&] () { write_out_solver_statistics(rank, la_bytes, domain_bytes, mpi_bytes, cells, dofs, nzes, filename); };
  Util::schedule_function(func, Util::ScheduleMode::clustered);
}

String Statistics::get_formatted_times(double total_time)
{
  String result = "Total time: " + stringify(total_time) + "s";
  if (total_time == 0.)
    return result;

  KahanAccumulation measured_time;
  measured_time = KahanSum(measured_time, get_time_reduction());
  measured_time = KahanSum(measured_time, get_time_spmv());
  measured_time = KahanSum(measured_time, get_time_axpy());
  measured_time = KahanSum(measured_time, get_time_precon());
  measured_time = KahanSum(measured_time, get_time_mpi_execute());
  measured_time = KahanSum(measured_time, get_time_mpi_wait());

  if (measured_time.sum > total_time)
    throw InternalError("Accumulated op time (" + stringify(measured_time.sum) + ") is greater as the provided total execution time (" + stringify(total_time) + ") !");

  result += "\n";
  result += "Accumulated op time: " + stringify(measured_time.sum) + "\n";

  result += "\n";

  double t_local = get_time_reduction();
  double t_max;
  double t_min;
  Util::Comm::allreduce(&t_local, 1, &t_max, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, 1, &t_min, Util::CommOperationMin());
  result += String("Reductions:").pad_back(17) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_axpy();
  Util::Comm::allreduce(&t_local, 1, &t_max, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, 1, &t_min, Util::CommOperationMin());
  result += String("Blas-1:").pad_back(17) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_spmv();
  Util::Comm::allreduce(&t_local, 1, &t_max, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, 1, &t_min, Util::CommOperationMin());
  result += String("Blas-2:").pad_back(17) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_precon();
  Util::Comm::allreduce(&t_local, 1, &t_max, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, 1, &t_min, Util::CommOperationMin());
  result += String("Precon Kernels:").pad_back(17) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_mpi_execute();
  Util::Comm::allreduce(&t_local, 1, &t_max, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, 1, &t_min, Util::CommOperationMin());
  result += String("MPI Execution:").pad_back(17) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_mpi_wait();
  Util::Comm::allreduce(&t_local, 1, &t_max, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, 1, &t_min, Util::CommOperationMin());
  result += String("MPI Wait:").pad_back(17) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = total_time - measured_time.sum;
  Util::Comm::allreduce(&t_local, 1, &t_max, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, 1, &t_min, Util::CommOperationMin());
  result += String("Not covered:").pad_back(17) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";
  return result;
}
