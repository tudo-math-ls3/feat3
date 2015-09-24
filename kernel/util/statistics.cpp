#include <kernel/util/statistics.hpp>
#include <kernel/util/exception.hpp>

using namespace FEAST;

// static member initialisation
Index Statistics::_flops = Index(0);
double Statistics::_time_reduction = 0.;
double Statistics::_time_spmv = 0.;
double Statistics::_time_axpy = 0.;
double Statistics::_time_precon = 0.;
double Statistics::_time_mpi_execute = 0.;
double Statistics::_time_mpi_wait = 0.;

void Statistics::add_flops(Index flops)
{
  _flops += flops;
}

Index Statistics::get_flops()
{
  return _flops;
}

String Statistics::get_formated_flops(double seconds)
{
  double flops((double)_flops);
  flops /= seconds;
  flops /= 1000.; // kilo
  flops /= 1000.; // mega
  flops /= 1000.; // giga
  return stringify(flops) + " GFlop/s";
}

void Statistics::reset_flops()
{
  _flops = Index(0);
}

void Statistics::add_time_reduction(double seconds)
{
  _time_reduction += seconds;
}
void Statistics::add_time_spmv(double seconds)
{
  _time_spmv += seconds;
}
void Statistics::add_time_axpy(double seconds)
{
  _time_axpy += seconds;
}
void Statistics::add_time_precon(double seconds)
{
  _time_precon += seconds;
}
void Statistics::add_time_mpi_execute(double seconds)
{
  _time_mpi_execute += seconds;
}
void Statistics::add_time_mpi_wait(double seconds)
{
  _time_mpi_wait += seconds;
}

double Statistics::get_time_reduction()
{
  return _time_reduction;
}
double Statistics::get_time_spmv()
{
  return _time_spmv;
}
double Statistics::get_time_axpy()
{
  return _time_axpy;
}
double Statistics::get_time_precon()
{
  return _time_precon;
}
double Statistics::get_time_mpi_execute()
{
  return _time_mpi_execute;
}
double Statistics::get_time_mpi_wait()
{
  return _time_mpi_wait;
}

String Statistics::get_formated_times()
{
  double total_time = _time_reduction + _time_spmv + _time_axpy + _time_precon + _time_mpi_execute + _time_mpi_wait;
  return get_formated_times(total_time);
}

String Statistics::get_formated_times(double total_time)
{
  String result = "Total time: " + stringify(total_time) + "s";
  if (total_time == 0.)
    return result;

  double measured_time = _time_reduction + _time_spmv + _time_axpy + _time_precon + _time_mpi_execute + _time_mpi_wait;
  if (measured_time > total_time)
    throw InternalError("Accumulated op time (" + stringify(measured_time) + ") is greater as the provided total execution time (" + stringify(total_time) + ") !");

  result += "\n";
  result += "Reductions: " + stringify(_time_reduction / total_time * 100.) + "%\n";
  result += "Blas-1: " + stringify(_time_axpy / total_time * 100.) + "%\n";
  result += "Blas-2: " + stringify(_time_spmv / total_time * 100.) + "%\n";
  result += "Precon Kernels: " + stringify(_time_precon / total_time * 100.) + "%\n";
  result += "MPI Execution: " + stringify(_time_mpi_execute / total_time * 100.) + "%\n";
  result += "MPI Wait: " + stringify(_time_mpi_wait / total_time * 100.) + "%\n";
  result += "Non-Flop: " + stringify( (total_time - measured_time) / total_time * 100.) + "%";
  return result;
}

void Statistics::reset_times()
{
  _time_reduction = 0.;
  _time_spmv = 0.;
  _time_axpy = 0.;
  _time_precon = 0.;
  _time_mpi_execute = 0.;
  _time_mpi_wait = 0.;
}
