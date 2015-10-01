#pragma once
#ifndef KERNEL_STATISTICS_HPP
#define KERNEL_STATISTICS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>

#include <list>
#include <map>

namespace FEAST
{
  struct SolverStatistics
  {
    ///contains def_init and def_cur per iteration, implicitly contains the iteration count, too.
    /// \todo use proper datatype via template parameter, eg quad for quad based solvers
    std::vector<double> defect;

    /// Time of each iteration in seconds
    std::vector<double> toe;
  };

  /**
   * \brief Statistics collection class
   *
   * The class Statistics encapsulates various hardware counters,
   * that collects e.g. the accumulated flop count of an linear solver.
   */
  class Statistics
  {
    private:

      /// global flop counter
      static Index _flops;

      /// global time of execution for reduction type operations
      static double _time_reduction;

      /// global time of execution for blas-2 type operations
      static double _time_spmv;

      /// global time of execution for blas-1 type operations
      static double _time_axpy;

      /// global time of execution for special preconditioner kernel type operations
      static double _time_precon;

      /// global time of execution for mpi related operations, e.g. send/recv or gather/scatter
      static double _time_mpi_execute;

      /// global time of execution for mpi related idle/wait tasks
      static double _time_mpi_wait;

      static std::map<FEAST::String, SolverStatistics> _solver_statistics;

      static String _format_solver_statistics(String branch, SolverStatistics & stat)
      {
        String result;

        result += "\n" + branch + "\n";
        result +="Iteration Defects:\n";

        for (Index i(0) ; i < stat.defect.size() ; ++i)
        {
          if (i == 0 && stat.defect.at(i) == double(-1))
            continue;
          if (stat.defect.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify(stat.defect.at(i)) + "\n";
        }

        result +="Iteration Timings:\n";
        for (Index i(0) ; i < stat.toe.size() ; ++i)
        {
          if (i == 0 && stat.toe.at(i) == double(-1))
            continue;
          if (stat.toe.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify(stat.toe.at(i)) + "\n";
        }

        return result;
      }

    public:

      /// Add an amount of flops to the global flop counter
      static void add_flops(Index flops)
      {
        _flops += flops;
      }

      /// Retrieve current global flop counter
      static Index get_flops()
      {
        return _flops;
      }

      /// Retrieve formated flops per second string
      static String get_formated_flops(double seconds)
      {
        double flops((double)_flops);
        flops /= seconds;
        flops /= 1000.; // kilo
        flops /= 1000.; // mega
        flops /= 1000.; // giga
        return stringify(flops) + " GFlop/s";
      }

      /// Reset global flop counter
      static void reset_flops()
      {
        _flops = Index(0);
      }

      inline static void add_time_reduction(double seconds)
      {
        _time_reduction += seconds;
      }
      inline static void add_time_spmv(double seconds)
      {
        _time_spmv += seconds;
      }
      inline static void add_time_axpy(double seconds)
      {
        _time_axpy += seconds;
      }
      inline static void add_time_precon(double seconds)
      {
        _time_precon += seconds;
      }
      inline static void add_time_mpi_execute(double seconds)
      {
        _time_mpi_execute += seconds;
      }
      inline static void add_time_mpi_wait(double seconds)
      {
        _time_mpi_wait += seconds;
      }

      inline static double get_time_reduction()
      {
        return _time_reduction;
      }
      inline static double get_time_spmv()
      {
        return _time_spmv;
      }
      inline static double get_time_axpy()
      {
        return _time_axpy;
      }
      inline static double get_time_precon()
      {
        return _time_precon;
      }
      inline static double get_time_mpi_execute()
      {
        return _time_mpi_execute;
      }
      inline static double get_time_mpi_wait()
      {
        return _time_mpi_wait;
      }

      inline static void add_solver_toe(String solver, double seconds)
      {
        auto it = _solver_statistics.find(solver);

        if (it != _solver_statistics.end())
        {
          it->second.toe.push_back(seconds);
        }
        else
        {
          SolverStatistics temp;
          temp.toe.push_back(seconds);
          _solver_statistics[solver] = temp;
        }
      }

      inline static void add_solver_defect(String solver, double defect)
      {
        auto it = _solver_statistics.find(solver);
        if (it != _solver_statistics.end())
        {
          it->second.defect.push_back(defect);
        }
        else
        {
          SolverStatistics temp;
          temp.defect.push_back(defect);
          _solver_statistics[solver] = temp;
        }
      }

      inline static void reset_solver_statistics()
      {
        _solver_statistics.clear();
      }

      /**
       * \brief Generate output with detailed solver statistics
       *
       * The generated String contains detailed detailed defect norms and execution times
       * for each used solver.
       *
       * The output per solver contains the following informations:
       * - solver name (branch name)
       * - defect norm of each iteration
       * - time of exection of each iteration
       *
       * \note All FEAST::Solver::IterativeSolver instances contain the initial defect norm as the first iteration defect norm.
       *
       * \note The FEAST::Solver::BasicVCycle prints its TOE from coarse to fine level.
       */
      static String get_formated_solvers()
      {
        String result;

        for (auto stat : _solver_statistics)
        {
          result += _format_solver_statistics(stat.first, stat.second);
          result += "\n";
        }

        return result;
      }

      /// \copydoc get_formated_solvers()
      static String get_formated_solver(String branch)
      {
        auto it = _solver_statistics.find(branch);
        if (it == _solver_statistics.end())
          throw InternalError("solver with branch " + branch + " not found!");
        else
          return _format_solver_statistics(branch, it->second);
      }

      /// Retrieve formated time consumption overview in percent
      static String get_formated_times()
      {
        double total_time = _time_reduction + _time_spmv + _time_axpy + _time_precon + _time_mpi_execute + _time_mpi_wait;
        return get_formated_times(total_time);
      }

      /// Retrieve formated time consumption overview in percent relative to some provided total time
      static String get_formated_times(double total_time)
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

      /// Reset all global timer counters
      static void reset_times()
      {
        _time_reduction = 0.;
        _time_spmv = 0.;
        _time_axpy = 0.;
        _time_precon = 0.;
        _time_mpi_execute = 0.;
        _time_mpi_wait = 0.;
      }

  };
} // namespace FEAST

#endif // KERNEL_STATISTICS_HPP
