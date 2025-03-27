// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/kahan_accumulator.hpp>
#include <kernel/solver/expression.hpp>

#include <list>
#include <map>
#include <time.h>
#include <iostream>
#include <fstream>
#include <memory>

namespace FEAT
{
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
      static KahanAccumulator<double> _time_reduction;

      /// global time of execution for blas-2 type operations
      static KahanAccumulator<double> _time_blas2;

      /// global time of execution for blas-3 type operations
      static KahanAccumulator<double> _time_blas3;

      /// global time of execution for blas-1 type operations
      static KahanAccumulator<double> _time_axpy;

      /// global time of execution for special preconditioner kernel type operations
      static KahanAccumulator<double> _time_precon;

      /// global time of execution for mpi related idle/wait tasks of (scalar) reduction operations
      static KahanAccumulator<double> _time_mpi_execute_reduction;

      /// global time of execution for mpi related idle/wait tasks of blas-2 operations
      static KahanAccumulator<double> _time_mpi_execute_blas2;

      /// global time of execution for mpi related idle/wait tasks of blas-3 operations
      static KahanAccumulator<double> _time_mpi_execute_blas3;

      /// global time of execution for mpi related idle/wait tasks of collective operations (without scalar reduction)
      static KahanAccumulator<double> _time_mpi_execute_collective;

      /// global time of wait execution for mpi related idle/wait tasks of (scalar) reduction operations
      static KahanAccumulator<double> _time_mpi_wait_reduction;

      /// global time of wait execution for mpi related idle/wait tasks of blas2 operations
      static KahanAccumulator<double> _time_mpi_wait_blas2;

      /// global time of wait execution for mpi related idle/wait tasks of blas3 operations
      static KahanAccumulator<double> _time_mpi_wait_blas3;

      /// global time of wait execution for mpi related idle/wait tasks of collective operations (without scalar reduction)
      static KahanAccumulator<double> _time_mpi_wait_collective;

      /// a consecutive list of all solver actions
      static std::map<String, std::list<std::shared_ptr<Solver::ExpressionBase>>> _solver_expressions;

      /// mapping of solver target name to formatted solver tree string
      static std::map<String, String> _formatted_solver_trees;

      /// overall time per reset call per solver name string.
      static std::map<String, std::list<double>> _overall_toe;
      static std::map<String, std::list<Index>> _overall_iters;
      static std::map<String, std::list<double>> _overall_mpi_execute_reduction;
      static std::map<String, std::list<double>> _overall_mpi_execute_blas2;
      static std::map<String, std::list<double>> _overall_mpi_execute_blas3;
      static std::map<String, std::list<double>> _overall_mpi_execute_collective;
      static std::map<String, std::list<double>> _overall_mpi_wait_reduction;
      static std::map<String, std::list<double>> _overall_mpi_wait_blas2;
      static std::map<String, std::list<double>> _overall_mpi_wait_blas3;
      static std::map<String, std::list<double>> _overall_mpi_wait_collective;
      /// mapping of solver name to list of outer multigrid level timings. each std::vector holds a complete level hierarchy of timings.
      static std::map<String, std::list<std::vector<double>>> _outer_mg_toe;
      static std::map<String, std::list<std::vector<double>>> _outer_mg_mpi_execute_reduction;
      static std::map<String, std::list<std::vector<double>>> _outer_mg_mpi_execute_blas2;
      static std::map<String, std::list<std::vector<double>>> _outer_mg_mpi_execute_blas3;
      static std::map<String, std::list<std::vector<double>>> _outer_mg_mpi_execute_collective;
      static std::map<String, std::list<std::vector<double>>> _outer_mg_mpi_wait_reduction;
      static std::map<String, std::list<std::vector<double>>> _outer_mg_mpi_wait_blas2;
      static std::map<String, std::list<std::vector<double>>> _outer_mg_mpi_wait_blas3;
      static std::map<String, std::list<std::vector<double>>> _outer_mg_mpi_wait_collective;
      /// overall time of outer schwarz preconditioners internal solver
      static std::map<String, std::list<double>> _outer_schwarz_toe;
      /// overall iterations of outer schwarz preconditioners internal solver
      static std::map<String, std::list<Index>> _outer_schwarz_iters;

      /*static String _format_solver_statistics(String branch, SolverStatistics & stat)
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
            result += stringify_fp_sci(stat.defect.at(i)) + "\n";
        }

        result +="Iteration Timings [s]:\n";
        for (Index i(0) ; i < stat.toe.size() ; ++i)
        {
          if (i == 0 && stat.toe.at(i) == double(-1))
            continue;
          if (stat.toe.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify_fp_sci(stat.toe.at(i)) + "\n";
        }

        result +="MPI Execution Iteration Timings [s]:\n";
        for (Index i(0) ; i < stat.mpi_execute.size() ; ++i)
        {
          if (i == 0 && stat.mpi_execute.at(i) == double(-1))
            continue;
          if (stat.mpi_execute.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify_fp_sci(stat.mpi_execute.at(i)) + "\n";
        }

        result +="MPI Wait Iteration Timings [s]:\n";
        for (Index i(0) ; i < stat.mpi_wait.size() ; ++i)
        {
          if (i == 0 && stat.mpi_wait.at(i) == double(-1))
            continue;
          if (stat.mpi_wait.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify_fp_sci(stat.mpi_wait.at(i)) + "\n";
        }

        return result;
      }*/

      static String _generate_formatted_solver_tree(String target);

    public:

      /// specifies whether collection of solver expressions is to be enabled
      static bool enable_solver_expressions;
      /// time of partitioning in seconds, needs initialization
      static double toe_partition;
      /// time of assembly in seconds, needs initialization
      static double toe_assembly;
      /// time of solution in seconds, needs initialization
      static double toe_solve;
      /// the current solver's descriptive string
      static String expression_target;


      /// Reset all internal counters (flops/times/solver_statistics)
      static void reset()
      {
        reset_flops();
        reset_times();
        _solver_expressions.clear();
        _overall_toe.clear();
        _overall_iters.clear();
        _overall_mpi_execute_reduction.clear();
        _overall_mpi_execute_blas2.clear();
        _overall_mpi_execute_blas3.clear();
        _overall_mpi_execute_collective.clear();
        _overall_mpi_wait_reduction.clear();
        _overall_mpi_wait_blas2.clear();
        _overall_mpi_wait_blas3.clear();
        _overall_mpi_wait_collective.clear();
        _outer_mg_toe.clear();
        _outer_mg_mpi_execute_reduction.clear();
        _outer_mg_mpi_execute_blas2.clear();
        _outer_mg_mpi_execute_blas3.clear();
        _outer_mg_mpi_execute_collective.clear();
        _outer_mg_mpi_wait_reduction.clear();
        _outer_mg_mpi_wait_blas2.clear();
        _outer_mg_mpi_wait_blas3.clear();
        _outer_mg_mpi_wait_collective.clear();
        _outer_schwarz_toe.clear();
        _outer_schwarz_iters.clear();
      }

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

      /// Retrieve formatted flops per second string
      static String get_formatted_flops(double seconds, int ranks = 1)
      {
        double flops((double)_flops);
        flops /= seconds;
        flops /= 1000.; // kilo
        flops /= 1000.; // mega
        flops /= 1000.; // giga
        return stringify(flops) + " GFlop/s [" + stringify(flops * double(ranks)) + " GFlops/s]";
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
      inline static void add_time_blas2(double seconds)
      {
        _time_blas2 += seconds;
      }
      inline static void add_time_blas3(double seconds)
      {
        _time_blas3 += seconds;
      }
      inline static void add_time_axpy(double seconds)
      {
        _time_axpy += seconds;
      }
      inline static void add_time_precon(double seconds)
      {
        _time_precon += seconds;
      }
      inline static void add_time_mpi_execute_reduction(double seconds)
      {
        _time_mpi_execute_reduction += seconds;
      }
      inline static void add_time_mpi_execute_blas2(double seconds)
      {
        _time_mpi_execute_blas2 += seconds;
      }
      inline static void add_time_mpi_execute_blas3(double seconds)
      {
        _time_mpi_execute_blas3 += seconds;
      }
      inline static void add_time_mpi_execute_collective(double seconds)
      {
        _time_mpi_execute_collective += seconds;
      }
      inline static void add_time_mpi_wait_reduction(double seconds)
      {
        _time_mpi_wait_reduction += seconds;
      }
      inline static void add_time_mpi_wait_blas2(double seconds)
      {
        _time_mpi_wait_blas2 += seconds;
      }
      inline static void add_time_mpi_wait_blas3(double seconds)
      {
        _time_mpi_wait_blas3 += seconds;
      }
      inline static void add_time_mpi_wait_collective(double seconds)
      {
        _time_mpi_wait_collective += seconds;
      }

      inline static double get_time_reduction()
      {
        return _time_reduction;
      }
      inline static double get_time_blas2()
      {
        return _time_blas2.value;
      }
      inline static double get_time_blas3()
      {
        return _time_blas3.value;
      }
      inline static double get_time_axpy()
      {
        return _time_axpy.value;
      }
      inline static double get_time_precon()
      {
        return _time_precon.value;
      }
      inline static double get_time_mpi_execute_reduction()
      {
        return _time_mpi_execute_reduction.value;
      }
      inline static double get_time_mpi_execute_blas2()
      {
        return _time_mpi_execute_blas2.value;
      }
      inline static double get_time_mpi_execute_blas3()
      {
        return _time_mpi_execute_blas3.value;
      }
      inline static double get_time_mpi_execute_collective()
      {
        return _time_mpi_execute_collective.value;
      }
      inline static double get_time_mpi_wait_reduction()
      {
        return _time_mpi_wait_reduction.value;
      }
      inline static double get_time_mpi_wait_blas2()
      {
        return _time_mpi_wait_blas2.value;
      }
      inline static double get_time_mpi_wait_blas3()
      {
        return _time_mpi_wait_blas3.value;
      }
      inline static double get_time_mpi_wait_collective()
      {
        return _time_mpi_wait_collective.value;
      }

      inline static void add_solver_expression(std::shared_ptr<Solver::ExpressionBase> expression)
      {
        if(enable_solver_expressions)
          _solver_expressions[expression_target].push_back(expression);
      }

      static const std::list<std::shared_ptr<Solver::ExpressionBase>> & get_solver_expressions()
      {
        return _solver_expressions.at(expression_target);
      }

      /**
       * \brief Returns a descriptive string of the complete solver tree.
       *
       * Create and format a string describing the complete solver tree.
       *
       * \note This method makes some simplifications, e.g. stating only one smoother
       * for the complete FEAT::Solver::BasicVCycle.
       *
       * \note The solver must have been executed (successfully) at least one time, to make the solver tree available.
       */
      static String get_formatted_solver_tree(String target = "default")
      {
        if (_formatted_solver_trees.count(target) == 0)
          compress_solver_expressions();
        if (_formatted_solver_trees.count(target) == 0)
          return String("no solver tree for target '") + target + "' found";
        else
          return _formatted_solver_trees.at(target);
      }

      /// print out the complete solver expression list
      static void print_solver_expressions();

      ///compress solver statistics (toe / defect norm / mpi timings) from previous calls
      static void compress_solver_expressions();

      /// retrieve list of all overall solver toe entries
      static inline std::list<double> & get_time_toe(String target)
      {
        return _overall_toe.at(target);
      }

      /// retrieve list of all overall solver iteration entries
      static inline std::list<Index> & get_iters(String target)
      {
        return _overall_iters.at(target);
      }

      /// retrieve list of all overall solver mpi execute reduction toe entries
      static inline std::list<double> & get_time_mpi_execute_reduction(String target)
      {
        return _overall_mpi_execute_reduction.at(target);
      }

      /// retrieve list of all overall solver mpi execute blas2 toe entries
      static inline std::list<double> & get_time_mpi_execute_blas2(String target)
      {
        return _overall_mpi_execute_blas2.at(target);
      }

      /// retrieve list of all overall solver mpi execute blas3 toe entries
      static inline std::list<double> & get_time_mpi_execute_blas3(String target)
      {
        return _overall_mpi_execute_blas3.at(target);
      }

      /// retrieve list of all overall solver mpi execute collective toe entries
      static inline std::list<double> & get_time_mpi_execute_collective(String target)
      {
        return _overall_mpi_execute_collective.at(target);
      }

      /// retrieve list of all overall solver mpi reduction wait toe entries
      static inline std::list<double> & get_time_mpi_wait_reduction(String target)
      {
        return _overall_mpi_wait_reduction.at(target);
      }

      /// retrieve list of all overall solver mpi blas2 wait toe entries
      static inline std::list<double> & get_time_mpi_wait_blas2(String target)
      {
        return _overall_mpi_wait_blas2.at(target);
      }

      /// retrieve list of all overall solver mpi blas3 wait toe entries
      static inline std::list<double> & get_time_mpi_wait_blas3(String target)
      {
        return _overall_mpi_wait_blas3.at(target);
      }

      /// retrieve list of all overall solver mpi collective wait toe entries
      static inline std::list<double> & get_time_mpi_wait_collective(String target)
      {
        return _overall_mpi_wait_collective.at(target);
      }

      /// retrieve list of all overall solver toe entries per mg level
      static inline std::list<std::vector<double>> & get_time_mg(String target)
      {
        return _outer_mg_toe.at(target);
      }

      /// retrieve list of all overall solver mpi execute reduction toe entries per level
      static inline std::list<std::vector<double>> & get_time_mg_mpi_execute_reduction(String target)
      {
        return _outer_mg_mpi_execute_reduction.at(target);
      }

      /// retrieve list of all overall solver mpi execute blas2 toe entries per level
      static inline std::list<std::vector<double>> & get_time_mg_mpi_execute_blas2(String target)
      {
        return _outer_mg_mpi_execute_blas2.at(target);
      }

      /// retrieve list of all overall solver mpi execute blas3 toe entries per level
      static inline std::list<std::vector<double>> & get_time_mg_mpi_execute_blas3(String target)
      {
        return _outer_mg_mpi_execute_blas3.at(target);
      }

      /// retrieve list of all overall solver mpi execute collective toe entries per level
      static inline std::list<std::vector<double>> & get_time_mg_mpi_execute_collective(String target)
      {
        return _outer_mg_mpi_execute_collective.at(target);
      }

      /// retrieve list of all overall solver mpi reduction wait toe entries per level
      static inline std::list<std::vector<double>> & get_time_mg_mpi_wait_reduction(String target)
      {
        return _outer_mg_mpi_wait_reduction.at(target);
      }

      /// retrieve list of all overall solver mpi blas2 wait toe entries per level
      static inline std::list<std::vector<double>> & get_time_mg_mpi_wait_blas2(String target)
      {
        return _outer_mg_mpi_wait_blas2.at(target);
      }

      /// retrieve list of all overall solver mpi blas3 wait toe entries per level
      static inline std::list<std::vector<double>> & get_time_mg_mpi_wait_blas3(String target)
      {
        return _outer_mg_mpi_wait_blas3.at(target);
      }

      /// retrieve list of all overall solver mpi collective wait toe entries per level
      static inline std::list<std::vector<double>> & get_time_mg_mpi_wait_collective(String target)
      {
        return _outer_mg_mpi_wait_collective.at(target);
      }

      /// retrieve list of all outer schwarz solver call toe entries
      static inline std::list<double> & get_time_schwarz(String target)
      {
        return _outer_schwarz_toe.at(target);
      }

      /// retrieve list of all outer schwarz solver call iteration count entries
      static inline std::list<Index> & get_iters_schwarz(String target)
      {
        return _outer_schwarz_iters.at(target);
      }

      /**
       * \brief Retrieve formatted time consumption overview in percent relative to some provided total time
       *
       * \note This method uses mpi collectives and thus needs to be called by all ranks, even if you don't use the result on every rank on your own.
       */
      static String get_formatted_times(double total_time);

      /// Retrieve formatted timings and iteration counts of internal solver structures for the provided solver target
      static String get_formatted_solver_internals(String target = "default");

      /// Reset all global timer counters
      static void reset_times()
      {
        _time_reduction.clear();
        _time_blas2.clear();
        _time_blas3.clear();
        _time_axpy.clear();
        _time_precon.clear();
        _time_mpi_execute_reduction.clear();
        _time_mpi_execute_blas2.clear();
        _time_mpi_execute_blas3.clear();
        _time_mpi_execute_collective.clear();
        _time_mpi_wait_reduction.clear();
        _time_mpi_wait_blas2.clear();
        _time_mpi_wait_blas3.clear();
        _time_mpi_wait_collective.clear();
      }

      /*
      static void write_out_solver_statistics_scheduled(Index rank, size_t la_bytes, size_t domain_bytes, size_t mpi_bytes, Index cells, Index dofs, Index nzes, String filename = "solver_stats");

      static void write_out_solver_statistics(Index rank, size_t la_bytes, size_t domain_bytes, size_t mpi_bytes, Index cells, Index dofs, Index nzes, String filename = "solver_stats")
      {
        filename += ".";
        filename += stringify(rank);

        std::ofstream file(filename.c_str(), std::ofstream::out);
        if (! file.is_open())
          XABORTM("Unable to open statistics file " + filename);

        {
          time_t t = time(nullptr);
          struct tm * ts = localtime(&t);

          file << "timestamp " << stringify(ts->tm_mon + 1) << "/" << stringify(ts->tm_mday) << "/" << stringify(ts->tm_year + 1900) << " " << stringify(ts->tm_hour)
            << ":" << stringify(ts->tm_min) << ":" << stringify(ts->tm_sec) << "\n";
        }

        file << "\n";

        file << "#global time\n";
        file << "partition " << stringify(toe_partition) << "\n";
        file << "assembly " << stringify(toe_assembly) << "\n";
        file << "solve " << stringify(toe_solve) << "\n";

        file << "\n";

        file << "#solver time\n";
        file << "reduction " << stringify(get_time_reduction()) << "\n";
        file << "axpy " << stringify(get_time_axpy()) << "\n";
        file << "blas2 " << stringify(get_time_blas2()) << "\n";
        file << "precon " << stringify(get_time_precon()) << "\n";
        file << "mpi " << stringify(get_time_mpi_execute()) << "\n";
        file << "wait reduction" << stringify(get_time_mpi_wait_reduction()) << "\n";
        file << "wait blas2" << stringify(get_time_mpi_wait_blas2()) << "\n";

        file << "\n";

        file << "#global size\n";
        file << "domain " << stringify(domain_bytes) << "\n";
        file << "mpi " << stringify(mpi_bytes) << "\n";
        file << "la " << stringify(la_bytes) << "\n";

        file << "\n";

        file <<"#problem size\n";
        file << "cells " << stringify(cells) << "\n";
        file << "dofs " << stringify(dofs) << "\n";
        file << "nzes " << stringify(nzes) << "\n";

        file << "\n";

        file << "#solver statistics\n";
        for (auto stat : _solver_statistics)
        {
          file << stat.first << "\n";

          file << "#defects\n";
          for (Index i(0) ; i < stat.second.defect.size() ; ++i)
          {
            if (i == 0 && stat.second.defect.at(i) == double(-1))
              continue;
            if (stat.second.defect.at(i) == double(-1))
              file << "-\n";
            else
              file << stringify_fp_sci(stat.second.defect.at(i)) << "\n";
          }

          file << "#toe\n";
          for (Index i(0) ; i < stat.second.toe.size() ; ++i)
          {
            if (i == 0 && stat.second.toe.at(i) == double(-1))
              continue;
            if (stat.second.toe.at(i) == double(-1))
              file << "-\n";
            else
              file << stringify_fp_sci(stat.second.toe.at(i)) << "\n";
          }

          file << "#mpi_execute\n";
          for (Index i(0) ; i < stat.second.mpi_execute.size() ; ++i)
          {
            if (i == 0 && stat.second.mpi_execute.at(i) == double(-1))
              continue;
            if (stat.second.mpi_execute.at(i) == double(-1))
              file << "-\n";
            else
              file << stringify_fp_sci(stat.second.mpi_execute.at(i)) << "\n";
          }

          file << "#mpi_wait\n";
          for (Index i(0) ; i < stat.second.mpi_wait.size() ; ++i)
          {
            if (i == 0 && stat.second.mpi_wait.at(i) == double(-1))
              continue;
            if (stat.second.mpi_wait.at(i) == double(-1))
              file << "-\n";
            else
              file << stringify_fp_sci(stat.second.mpi_wait.at(i)) << "\n";
          }
        }

        file.close();
      }*/
  };
} // namespace FEAT
