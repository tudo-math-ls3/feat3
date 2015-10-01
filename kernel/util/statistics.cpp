#include <kernel/util/statistics.hpp>

using namespace FEAST;

// static member initialisation
Index Statistics::_flops = Index(0);
double Statistics::_time_reduction = 0.;
double Statistics::_time_spmv = 0.;
double Statistics::_time_axpy = 0.;
double Statistics::_time_precon = 0.;
double Statistics::_time_mpi_execute = 0.;
double Statistics::_time_mpi_wait = 0.;
std::map<FEAST::String, SolverStatistics> Statistics::_solver_statistics;
