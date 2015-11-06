#include <kernel/util/statistics.hpp>

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

KahanAccumulation FEAST::KahanSum(KahanAccumulation accumulation, double value)
{
  KahanAccumulation result;
  double y = value - accumulation.correction;
  double t = accumulation.sum + y;
  result.correction = (t - accumulation.sum) - y;
  result.sum = t;
  return result;
};
