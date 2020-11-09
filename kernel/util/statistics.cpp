// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/statistics.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/solver/base.hpp>

#include <queue>
#include <numeric>

using namespace FEAT;

// static member initialization
Index Statistics::_flops = Index(0);
KahanAccumulation Statistics::_time_reduction;
KahanAccumulation Statistics::_time_blas2;
KahanAccumulation Statistics::_time_blas3;
KahanAccumulation Statistics::_time_axpy;
KahanAccumulation Statistics::_time_precon;
KahanAccumulation Statistics::_time_mpi_execute_reduction;
KahanAccumulation Statistics::_time_mpi_execute_blas2;
KahanAccumulation Statistics::_time_mpi_execute_blas3;
KahanAccumulation Statistics::_time_mpi_execute_collective;
KahanAccumulation Statistics::_time_mpi_wait_reduction;
KahanAccumulation Statistics::_time_mpi_wait_blas2;
KahanAccumulation Statistics::_time_mpi_wait_blas3;
KahanAccumulation Statistics::_time_mpi_wait_collective;
std::map<String, std::list<std::shared_ptr<Solver::ExpressionBase>>> Statistics::_solver_expressions;
std::map<String, String> Statistics::_formatted_solver_trees;
std::map<String, std::list<double>> Statistics::_overall_toe;
std::map<String, std::list<Index>> Statistics::_overall_iters;
std::map<String, std::list<double>> Statistics::_overall_mpi_execute_reduction;
std::map<String, std::list<double>> Statistics::_overall_mpi_execute_blas2;
std::map<String, std::list<double>> Statistics::_overall_mpi_execute_blas3;
std::map<String, std::list<double>> Statistics::_overall_mpi_execute_collective;
std::map<String, std::list<double>> Statistics::_overall_mpi_wait_reduction;
std::map<String, std::list<double>> Statistics::_overall_mpi_wait_blas2;
std::map<String, std::list<double>> Statistics::_overall_mpi_wait_blas3;
std::map<String, std::list<double>> Statistics::_overall_mpi_wait_collective;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_toe;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_execute_reduction;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_execute_blas2;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_execute_blas3;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_execute_collective;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_wait_reduction;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_wait_blas2;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_wait_blas3;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_wait_collective;
std::map<String, std::list<double>> Statistics::_outer_schwarz_toe;
std::map<String, std::list<Index>> Statistics::_outer_schwarz_iters;
String Statistics::expression_target = "default";
double Statistics::toe_partition;
double Statistics::toe_assembly;
double Statistics::toe_solve;

String Statistics::_generate_formatted_solver_tree(String target)
{
  std::list<String> names;
  std::list<int> found; //number of found preconds / smoothers / coarse solvers. 0 = nothing found, 1 = smoother / s found, 2 = coarse solver / a found, 3 = all found

  if(_solver_expressions.count(target) == 0)
  {
    XABORTM("target "+target+" not present in _solver_expressions");
  }

  if (_solver_expressions[target].front()->get_type() != Solver::ExpressionType::start_solve)
  {
    XABORTM("Should never happen - _solver_expressions list did not start with start solve expression!");
  }

  // process the very first entry, e.g. the outer most solver
  auto it = _solver_expressions[target].begin();
  String tree((*it)->solver_name);
  names.push_back((*it)->solver_name);
  found.push_back(0);

  // process current last element in the list, until no element needs to be processed
  while (names.size() > 0 && it != _solver_expressions[target].end())
  {
    ++it;

    // the current solver is of multigrid type, search for smoother and coarse solver (smoother always comes first).
    // if smoother and coarse solver have been found, skip everything until solver end statement has been found
    if (names.back().starts_with("MultiGrid") || names.back().starts_with("VCycle") || names.back().starts_with("ScaRCMultiGrid"))
    {
      while (it != _solver_expressions[target].end())
      {
        auto expression = *it;

        if (expression->get_type() == Solver::ExpressionType::call_smoother && expression->solver_name == names.back() && (found.back() == 0 || found.back() == 2))
        {
          found.back() += 1;
          auto j = it;
          ++j;
          // smoother call found, that is a solver on its own. break processing of current solver and dive on step deeper into smoother solver processing
          if ((*j)->get_type() == Solver::ExpressionType::start_solve)
          {
            ++it; //shift over to solver start expression
            tree += " ( S: " + (*it)->solver_name;
            names.push_back((*it)->solver_name);
            found.push_back(0);
            break;
          }
          // smoother is no solver on its own, we can continue with the current solver's end statement search
          else
          {
            auto t = dynamic_cast<Solver::ExpressionCallSmoother*>(expression.get());
            tree += " ( S: " + t->smoother_name;
          }
        }

        if (expression->get_type() == Solver::ExpressionType::call_coarse_solver && expression->solver_name == names.back() && (found.back() == 0 || found.back() == 1))
        {
          found.back() += 2;
          auto j = it;
          ++j;
          // coarse solver call found, that is a solver on its own. break processing of current solver and dive on step deeper into coarse solver processing
          if ((*j)->get_type() == Solver::ExpressionType::start_solve)
          {
            ++it; //shift over to solver start expression
            tree += " / C: " + (*it)->solver_name;
            names.push_back((*it)->solver_name);
            found.push_back(0);
            break;
          }
          // coarse solver is no solver on its own, we can continue with the current solver's end statement search
          else
          {
            auto t = dynamic_cast<Solver::ExpressionCallCoarseSolver*>(expression.get());
            tree += " / C: " + t->coarse_solver_name;
          }
        }

        if (expression->get_type() == Solver::ExpressionType::end_solve && expression->solver_name == names.back())
        {
          tree += " )";
          names.pop_back();
          found.pop_back();
          break;
        }
        ++it;
      }
    }

    // the current solver is of uzawa complement type, search for a and s solver
    // if both have been found, skip everything until solver end statement has been found
    else if (names.back().starts_with("Uzawa"))
    {
      while (it != _solver_expressions[target].end())
      {
        auto expression = *it;

        if (expression->get_type() == Solver::ExpressionType::call_uzawa_s && expression->solver_name == names.back() && (found.back() == 0 || found.back() == 2))
        {
          found.back() += 1;
          auto j = it;
          ++j;
          // s solver call found, that is a solver on its own. break processing of current solver and dive on step deeper into s solver processing
          if ((*j)->get_type() == Solver::ExpressionType::start_solve)
          {
            ++it; //shift over to solver start expression
            if (found.back() == 1)
              tree += " ( S: " + (*it)->solver_name;
            else
            tree += " / S: " + (*it)->solver_name;
            names.push_back((*it)->solver_name);
            found.push_back(0);
            break;
          }
          // s solver is no solver on its own, we can continue with the current solver's end statement search
          else
          {
            auto t = dynamic_cast<Solver::ExpressionCallUzawaS*>(expression.get());
            if (found.back() == 1)
              tree += " ( S: " + t->solver_s_name;
            else
              tree += " / S: " + t->solver_s_name;
          }
        }

        if (expression->get_type() == Solver::ExpressionType::call_uzawa_a && expression->solver_name == names.back() && (found.back() == 0 || found.back() == 1))
        {
          found.back() += 2;
          auto j = it;
          ++j;
          // a solver call found, that is a solver on its own. break processing of current solver and dive on step deeper into a solver processing
          if ((*j)->get_type() == Solver::ExpressionType::start_solve)
          {
            ++it; //shift over to solver start expression
            if (found.back() == 2)
              tree += " ( A: " + (*it)->solver_name;
            else
              tree += " / A: " + (*it)->solver_name;
            names.push_back((*it)->solver_name);
            found.push_back(0);
            break;
          }
          // a solver is no solver on its own, we can continue with the current solver's end statement search
          else
          {
            auto t = dynamic_cast<Solver::ExpressionCallUzawaA*>(expression.get());
            if (found.back() == 2)
              tree += " ( A: " + t->solver_a_name;
            else
              tree += " / A: " + t->solver_a_name;
          }
        }

        if (expression->get_type() == Solver::ExpressionType::end_solve && expression->solver_name == names.back())
        {
          tree += " )";
          names.pop_back();
          found.pop_back();
          break;
        }
        ++it;
      }
    }

    // the current solver uses l and r preconditioners, search for both of them
    // if both have been found, skip everything until solver end statement has been found
    else if (names.back().starts_with("PCGNR"))
    {
      while (it != _solver_expressions[target].end())
      {
        auto expression = *it;

        if (expression->get_type() == Solver::ExpressionType::call_precond_l && expression->solver_name == names.back() && (found.back() == 0 || found.back() == 2))
        {
          found.back() += 1;
          auto j = it;
          ++j;
          // l solver call found, that is a solver on its own. break processing of current solver and dive on step deeper into l solver processing
          if ((*j)->get_type() == Solver::ExpressionType::start_solve)
          {
            ++it; //shift over to solver start expression
            if (found.back() == 1)
              tree += " ( L: " + (*it)->solver_name;
            else
            tree += " / L: " + (*it)->solver_name;
            names.push_back((*it)->solver_name);
            found.push_back(0);
            break;
          }
          // l solver is no solver on its own, we can continue with the current solver's end statement search
          else
          {
            auto t = dynamic_cast<Solver::ExpressionCallPrecondL*>(expression.get());
            if (found.back() == 1)
              tree += " ( L: " + t->precond_name;
            else
              tree += " / L: " + t->precond_name;
          }
        }

        if (expression->get_type() == Solver::ExpressionType::call_precond_r && expression->solver_name == names.back() && (found.back() == 0 || found.back() == 1))
        {
          found.back() += 2;
          auto j = it;
          ++j;
          // r solver call found, that is a solver on its own. break processing of current solver and dive on step deeper into r solver processing
          if ((*j)->get_type() == Solver::ExpressionType::start_solve)
          {
            ++it; //shift over to solver start expression
            if (found.back() == 2)
              tree += " ( R: " + (*it)->solver_name;
            else
              tree += " / R: " + (*it)->solver_name;
            names.push_back((*it)->solver_name);
            found.push_back(0);
            break;
          }
          // r solver is no solver on its own, we can continue with the current solver's end statement search
          else
          {
            auto t = dynamic_cast<Solver::ExpressionCallPrecondR*>(expression.get());
            if (found.back() == 2)
              tree += " ( R: " + t->precond_name;
            else
              tree += " / R: " + t->precond_name;
          }
        }

        if (expression->get_type() == Solver::ExpressionType::end_solve && expression->solver_name == names.back())
        {
          tree += " )";
          names.pop_back();
          found.pop_back();
          break;
        }
        ++it;
      }
    }

    else
    {
      // the current solver is not of multigrid or uzawa type, i.e. is uses at most one preconditioner, search for its call or the solvers end.
      while (it != _solver_expressions[target].end())
      {
        auto expression = *it;

        if (expression->get_type() == Solver::ExpressionType::call_precond && found.back() == 0)
        {
          found.back() += 1;
          auto j = it;
          ++j;
          // preconditioner call found, that is a solver on its own. break processing of current solver and dive on step deeper into preconditioner solver processing
          if ((*j)->get_type() == Solver::ExpressionType::start_solve)
          {
            ++it; //shift over to solver start expression
            tree += " ( " + (*it)->solver_name;
            names.push_back((*it)->solver_name);
            found.push_back(0);
            break;
          }
          // preconditioner is no solver on its own, we can continue with the current solver's end statement search
          else
          {
            auto t = dynamic_cast<Solver::ExpressionCallPrecond*>(expression.get());
            tree += " ( " + t->precond_name;
          }
        }

        if (expression->get_type() == Solver::ExpressionType::end_solve && expression->solver_name == names.back())
        {
          if (found.back() == 0)
            tree += " ( none";
          tree += " )";
          names.pop_back();
          found.pop_back();
          break;
        }
        ++it;
      }
    }
  }

  if (names.size() > 0)
  {
    XABORTM("Should never happen - not all solver calls were parsed to the end!");
  }

  return tree;
}

void Statistics::print_solver_expressions()
{
  size_t padding(0);

  for (auto expression : _solver_expressions[expression_target])
  {
    switch (expression->get_type())
    {
      case Solver::ExpressionType::start_solve:
        {
          String s = stringify(expression->get_type()) + "[" + expression->solver_name + "]";
          std::cout<<String(padding, ' ') << s << std::endl;
          padding += 2;
          break;
        }
      case Solver::ExpressionType::end_solve:
        {
          auto t = dynamic_cast<Solver::ExpressionEndSolve*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + stringify(t->status) + " / " + stringify(t->iters) + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          padding -= 2;
          break;
        }
      case Solver::ExpressionType::defect:
        {
          auto t = dynamic_cast<Solver::ExpressionDefect*>(expression.get());
          String s = stringify(t->get_type()) + "[" + stringify(t->solver_name) + "] (" + stringify(t->def) + " / " + stringify(t->iter) + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::call_precond:
        {
          auto t = dynamic_cast<Solver::ExpressionCallPrecond*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + t->precond_name + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::call_precond_l:
        {
          auto t = dynamic_cast<Solver::ExpressionCallPrecondL*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + t->precond_name + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::call_precond_r:
        {
          auto t = dynamic_cast<Solver::ExpressionCallPrecondR*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + t->precond_name + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::call_smoother:
        {
          auto t = dynamic_cast<Solver::ExpressionCallSmoother*>(expression.get());
          String s = stringify(t->get_type()) + "[" + stringify(t->solver_name) + "] (" + t->smoother_name + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::call_coarse_solver:
        {
          auto t = dynamic_cast<Solver::ExpressionCallCoarseSolver*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + t->coarse_solver_name + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::prol:
        {
          auto t = dynamic_cast<Solver::ExpressionProlongation*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + stringify(t->level) + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::rest:
        {
          auto t = dynamic_cast<Solver::ExpressionRestriction*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + stringify(t->level) + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::timings:
        {
          auto t = dynamic_cast<Solver::ExpressionTimings*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + stringify(t->solver_toe) + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::level_timings:
        {
          auto t = dynamic_cast<Solver::ExpressionLevelTimings*>(expression.get());
          String s = stringify(t->get_type()) + "[" + t->solver_name + "] (" + stringify(t->level) + " / " + stringify(t->level_toe) + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::call_uzawa_s:
        {
          auto t = dynamic_cast<Solver::ExpressionCallUzawaS*>(expression.get());
          String s = stringify(t->get_type()) + "[" + stringify(t->solver_name) + "] (" + t->solver_s_name + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::call_uzawa_a:
        {
          auto t = dynamic_cast<Solver::ExpressionCallUzawaA*>(expression.get());
          String s = stringify(t->get_type()) + "[" + stringify(t->solver_name) + "] (" + t->solver_a_name + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      default:
        {
          String s = stringify(expression->get_type()) + "[" + expression->solver_name + "]";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
    }
  }
  std::cout<<std::endl;
}

String Statistics::get_formatted_times(double total_time)
{
  String result = "Total time: " + stringify(total_time) + "s";
  if (total_time == 0.)
    return result;

  KahanAccumulation measured_time;
  measured_time = KahanSum(measured_time, get_time_reduction());
  measured_time = KahanSum(measured_time, get_time_blas2());
  measured_time = KahanSum(measured_time, get_time_blas3());
  measured_time = KahanSum(measured_time, get_time_axpy());
  measured_time = KahanSum(measured_time, get_time_precon());
  measured_time = KahanSum(measured_time, get_time_mpi_execute_reduction());
  measured_time = KahanSum(measured_time, get_time_mpi_execute_blas2());
  measured_time = KahanSum(measured_time, get_time_mpi_execute_blas3());
  measured_time = KahanSum(measured_time, get_time_mpi_execute_collective());
  measured_time = KahanSum(measured_time, get_time_mpi_wait_reduction());
  measured_time = KahanSum(measured_time, get_time_mpi_wait_blas2());
  measured_time = KahanSum(measured_time, get_time_mpi_wait_blas3());
  measured_time = KahanSum(measured_time, get_time_mpi_wait_collective());

  result += "\n";
  result += "Accumulated op time: " + stringify(measured_time.sum) + "\n";

  result += "\n";

  Dist::Comm comm(Dist::Comm::world());

  double t_max[14];
  double t_min[14];
  double t_local[14];
  t_local[0] = total_time - measured_time.sum;
  t_local[1] = get_time_reduction();
  t_local[2] = get_time_axpy();
  t_local[3] = get_time_blas2();
  t_local[4] = get_time_blas3();
  t_local[5] = get_time_precon();
  t_local[6] = get_time_mpi_execute_reduction();
  t_local[7] = get_time_mpi_execute_blas2();
  t_local[8] = get_time_mpi_execute_blas3();
  t_local[9] = get_time_mpi_execute_collective();
  t_local[10] = get_time_mpi_wait_reduction();
  t_local[11] = get_time_mpi_wait_blas2();
  t_local[12] = get_time_mpi_wait_blas3();
  t_local[13] = get_time_mpi_wait_collective();

  comm.allreduce(t_local, t_max, std::size_t(14), Dist::op_max);
  comm.allreduce(t_local, t_min, std::size_t(14), Dist::op_min);

  result += String("Reductions:").pad_back(22) + "max: " + stringify(t_max[1]) + ", min: " + stringify(t_min[1]) + ", local: " + stringify(t_local[1]) + "\n";

  result += String("Blas-1:").pad_back(22) + "max: " + stringify(t_max[2]) + ", min: " + stringify(t_min[2]) + ", local: " + stringify(t_local[2]) + "\n";

  result += String("Blas-2:").pad_back(22) + "max: " + stringify(t_max[3]) + ", min: " + stringify(t_min[3]) + ", local: " + stringify(t_local[3]) + "\n";

  result += String("Blas-3:").pad_back(22) + "max: " + stringify(t_max[4]) + ", min: " + stringify(t_min[4]) + ", local: " + stringify(t_local[4]) + "\n";

  result += String("Precon Kernels:").pad_back(22) + "max: " + stringify(t_max[5]) + ", min: " + stringify(t_min[5]) + ", local: " + stringify(t_local[5]) + "\n";

  result += String("MPI Exec Reduction:").pad_back(22) + "max: " + stringify(t_max[6]) + ", min: " + stringify(t_min[6]) + ", local: " + stringify(t_local[6]) + "\n";

  result += String("MPI Exec Blas-2:").pad_back(22) + "max: " + stringify(t_max[7]) + ", min: " + stringify(t_min[7]) + ", local: " + stringify(t_local[7]) + "\n";

  result += String("MPI Exec Blas-3:").pad_back(22) + "max: " + stringify(t_max[8]) + ", min: " + stringify(t_min[8]) + ", local: " + stringify(t_local[8]) + "\n";

  result += String("MPI Exec Collective:").pad_back(22) + "max: " + stringify(t_max[9]) + ", min: " + stringify(t_min[9]) + ", local: " + stringify(t_local[9]) + "\n";

  result += String("MPI Wait Reduction:").pad_back(22) + "max: " + stringify(t_max[10]) + ", min: " + stringify(t_min[10]) + ", local: " + stringify(t_local[10]) + "\n";

  result += String("MPI Wait Blas-2:").pad_back(22) + "max: " + stringify(t_max[11]) + ", min: " + stringify(t_min[11]) + ", local: " + stringify(t_local[11]) + "\n";

  result += String("MPI Wait Blas-3:").pad_back(22) + "max: " + stringify(t_max[12]) + ", min: " + stringify(t_min[12]) + ", local: " + stringify(t_local[12]) + "\n";

  result += String("MPI Wait Collective:").pad_back(22) + "max: " + stringify(t_max[13]) + ", min: " + stringify(t_min[13]) + ", local: " + stringify(t_local[13]) + "\n";

  result += String("Not covered:").pad_back(22) + "max: " + stringify(t_max[0]) + ", min: " + stringify(t_min[0]) + ", local: " + stringify(t_local[0]) + "\n";

  if (t_min[0] < 0.0)
  {
    // total_time < measured_time.sum for at least one process
    result += "WARNING: Accumulated op time is greater than the provided total execution time !\n";
  }
  return result;
}

String Statistics::get_formatted_solver_internals(String target)
{
  Dist::Comm comm(Dist::Comm::world());

  if (_formatted_solver_trees.count(target) == 0)
    compress_solver_expressions();

  auto solver_time_mg = FEAT::Statistics::get_time_mg(target);
  auto solver_time_mg_mpi_execute_reduction = FEAT::Statistics::get_time_mg_mpi_execute_reduction(target);
  auto solver_time_mg_mpi_execute_blas2 = FEAT::Statistics::get_time_mg_mpi_execute_blas2(target);
  auto solver_time_mg_mpi_execute_blas3 = FEAT::Statistics::get_time_mg_mpi_execute_blas3(target);
  auto solver_time_mg_mpi_execute_collective = FEAT::Statistics::get_time_mg_mpi_execute_collective(target);
  auto solver_time_mg_mpi_wait_reduction = FEAT::Statistics::get_time_mg_mpi_wait_reduction(target);
  auto solver_time_mg_mpi_wait_blas2 = FEAT::Statistics::get_time_mg_mpi_wait_blas2(target);
  auto solver_time_mg_mpi_wait_blas3 = FEAT::Statistics::get_time_mg_mpi_wait_blas3(target);
  auto solver_time_mg_mpi_wait_collective = FEAT::Statistics::get_time_mg_mpi_wait_collective(target);

  Index item_count(Index(12) + Index(solver_time_mg.front().size() + solver_time_mg_mpi_execute_reduction.front().size() + solver_time_mg_mpi_execute_blas2.front().size() +
        solver_time_mg_mpi_execute_blas3.front().size() + solver_time_mg_mpi_execute_collective.front().size() +
        solver_time_mg_mpi_wait_reduction.front().size() + solver_time_mg_mpi_wait_blas2.front().size() + solver_time_mg_mpi_wait_blas3.front().size() +
        solver_time_mg_mpi_wait_collective.front().size()));

  double * t_local = new double[item_count];
  double * t_max = new double[item_count];
  double * t_min = new double[item_count];

  /*
   * array to value mapping:
   * 0 solver_toe
   * 1 solver_iters
   * 2 solver_mpi_execute_reduction
   * 3 solver_mpi_execute_blas2
   * 4 solver_mpi_execute_blas3
   * 5 solver_mpi_execute_collective
   * 6 solver_mpi_wait_reduction
   * 7 solver_mpi_wait_blas2
   * 8 solver_mpi_wait_blas3
   * 9 solver_mpi_wait_collective
   * 10 solver_schwarz_toe
   * 11 solver_schwarz_iters
   * n solver_mg_toe
   * n solver_mg_mpi_execute_reduction
   * n solver_mg_mpi_execute_blas2
   * n solver_mg_mpi_execute_blas3
   * n solver_mg_mpi_execute_collective
   * n solver_mg_mpi_wait_reduction
   * n solver_mg_mpi_wait_blas2
   * n solver_mg_mpi_wait_blas3
   * n solver_mg_mpi_wait_collective
   */

  t_local[0] = std::accumulate(FEAT::Statistics::get_time_toe(target).begin(),
        FEAT::Statistics::get_time_toe(target).end(), 0.);
  t_local[1] = double(std::accumulate(FEAT::Statistics::get_iters(target).begin(),
        FEAT::Statistics::get_iters(target).end(), Index(0)));
  t_local[2] = std::accumulate(FEAT::Statistics::get_time_mpi_execute_reduction(target).begin(),
        FEAT::Statistics::get_time_mpi_execute_reduction(target).end(), 0.);
  t_local[3] = std::accumulate(FEAT::Statistics::get_time_mpi_execute_blas2(target).begin(),
        FEAT::Statistics::get_time_mpi_execute_blas2(target).end(), 0.);
  t_local[4] = std::accumulate(FEAT::Statistics::get_time_mpi_execute_blas3(target).begin(),
        FEAT::Statistics::get_time_mpi_execute_blas3(target).end(), 0.);
  t_local[5] = std::accumulate(FEAT::Statistics::get_time_mpi_execute_collective(target).begin(),
        FEAT::Statistics::get_time_mpi_execute_collective(target).end(), 0.);
  t_local[6] = std::accumulate(FEAT::Statistics::get_time_mpi_wait_reduction(target).begin(),
        FEAT::Statistics::get_time_mpi_wait_reduction(target).end(), 0.);
  t_local[7] = std::accumulate(FEAT::Statistics::get_time_mpi_wait_blas2(target).begin(),
        FEAT::Statistics::get_time_mpi_wait_blas2(target).end(), 0.);
  t_local[8] = std::accumulate(FEAT::Statistics::get_time_mpi_wait_blas3(target).begin(),
        FEAT::Statistics::get_time_mpi_wait_blas3(target).end(), 0.);
  t_local[9] = std::accumulate(FEAT::Statistics::get_time_mpi_wait_collective(target).begin(),
        FEAT::Statistics::get_time_mpi_wait_collective(target).end(), 0.);
  t_local[10] = std::accumulate(FEAT::Statistics::get_time_schwarz(target).begin(),
        FEAT::Statistics::get_time_schwarz(target).end(), 0.);
  t_local[11] = double(std::accumulate(FEAT::Statistics::get_iters_schwarz(target).begin(),
        FEAT::Statistics::get_iters_schwarz(target).end(), Index(0)));
  t_local[12] /= double(FEAT::Statistics::get_iters_schwarz(target).size());

  Index offset(12);
  Index levels(Index(solver_time_mg.front().size()));

  for (auto& step : solver_time_mg)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  offset += levels;

  for (auto& step : solver_time_mg_mpi_execute_reduction)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  offset += levels;

  for (auto& step : solver_time_mg_mpi_execute_blas2)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  offset += levels;

  for (auto& step : solver_time_mg_mpi_execute_blas3)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  offset += levels;

  for (auto& step : solver_time_mg_mpi_execute_collective)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  offset += levels;

  for (auto& step : solver_time_mg_mpi_wait_reduction)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  offset += levels;

  for (auto& step : solver_time_mg_mpi_wait_blas2)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  offset += levels;

  for (auto& step : solver_time_mg_mpi_wait_blas3)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  offset += levels;

  for (auto& step : solver_time_mg_mpi_wait_collective)
  {
    XASSERT(step.size() == levels);
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] = double(0);
    }
    for (Index i(0) ; i < step.size() ; ++i)
    {
      t_local[offset + i] += step.at(i);
    }
  }

  comm.allreduce(t_local, t_max, std::size_t(item_count), Dist::op_max);
  comm.allreduce(t_local, t_min, std::size_t(item_count), Dist::op_min);

  String result = target + "\n";
  result += String("toe:").pad_back(27) + "max: " + stringify(t_max[0]) + ", min: " + stringify(t_min[0]) + ", local: " +
      stringify(t_local[0]) + "\n";
  result += String("iters:").pad_back(27) + "max: " + stringify(Index(t_max[1])) + ", min: " + stringify(Index(t_min[1])) + ", local: " +
      stringify(Index(t_local[1])) + "\n";
  result += String("mpi exe reduction:").pad_back(27) + "max: " + stringify(t_max[2]) + ", min: " + stringify(t_min[2]) + ", local: " +
      stringify(t_local[2]) + "\n";
  result += String("mpi exec blas2:").pad_back(27) + "max: " + stringify(t_max[3]) + ", min: " + stringify(t_min[3]) + ", local: " +
      stringify(t_local[3]) + "\n";
  result += String("mpi exec blas3:").pad_back(27) + "max: " + stringify(t_max[3]) + ", min: " + stringify(t_min[3]) + ", local: " +
      stringify(t_local[3]) + "\n";
  result += String("mpi exec collective:").pad_back(27) + "max: " + stringify(t_max[4]) + ", min: " + stringify(t_min[4]) + ", local: " +
      stringify(t_local[4]) + "\n";
  result += String("mpi wait reduction:").pad_back(27) + "max: " + stringify(t_max[5]) + ", min: " + stringify(t_min[5]) + ", local: " +
      stringify(t_local[5]) + "\n";
  result += String("mpi wait blas2:").pad_back(27) + "max: " + stringify(t_max[6]) + ", min: " + stringify(t_min[6]) + ", local: " +
      stringify(t_local[6]) + "\n";
  result += String("mpi wait blas3:").pad_back(27) + "max: " + stringify(t_max[6]) + ", min: " + stringify(t_min[6]) + ", local: " +
      stringify(t_local[6]) + "\n";
  result += String("mpi wait collective:").pad_back(27) + "max: " + stringify(t_max[7]) + ", min: " + stringify(t_min[7]) + ", local: " +
      stringify(t_local[7]) + "\n";
  result += String("schwarz toe:").pad_back(27) + "max: " + stringify(t_max[8]) + ", min: " + stringify(t_min[8]) + ", local: " +
      stringify(t_local[8]) + "\n";
  result += String("schwarz iters:").pad_back(27) + "max: " + stringify(Index(t_max[9])) + ", min: " + stringify(Index(t_min[9])) + ", local: " +
      stringify(Index(t_local[9])) + "\n";

  offset = 12;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("toe lvl ") + stringify(i).pad_back(19) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) + ", local: " +
        stringify(t_local[offset + i]) + "\n";
  }
  offset += levels;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("mpi exec reduction lvl ") + stringify(i).pad_back(4) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) +
        ", local: " + stringify(t_local[offset + i]) + "\n";
  }
  offset += levels;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("mpi exec blas2 lvl ") + stringify(i).pad_back(9) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) +
        ", local: " + stringify(t_local[offset + i]) + "\n";
  }
  offset += levels;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("mpi exec blas3 lvl ") + stringify(i).pad_back(9) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) +
        ", local: " + stringify(t_local[offset + i]) + "\n";
  }
  offset += levels;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("mpi exec collective lvl ") + stringify(i).pad_back(3) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) +
        ", local: " + stringify(t_local[offset + i]) + "\n";
  }
  offset += levels;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("mpi wait red lvl ") + stringify(i).pad_back(10) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) +
        ", local: " + stringify(t_local[offset + i]) + "\n";
  }
  offset += levels;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("mpi wait blas2 lvl ") + stringify(i).pad_back(9) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) +
        ", local: " + stringify(t_local[offset + i]) + "\n";
  }
  offset += levels;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("mpi wait blas3 lvl ") + stringify(i).pad_back(9) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) +
        ", local: " + stringify(t_local[offset + i]) + "\n";
  }
  offset += levels;
  for (Index i(0) ; i < levels ; ++i)
  {
    result += String("mpi wait collective lvl ") + stringify(i).pad_back(3) + "max: " + stringify(t_max[offset + i]) + ", min: " + stringify(t_min[offset + i]) +
        ", local: " + stringify(t_local[offset + i]);
    if (i != levels - 1)
      result += "\n";
  }

  delete[] t_local;
  delete[] t_max;
  delete[] t_min;

  return result;
}

void Statistics::compress_solver_expressions()
{
  for (auto itarget = _solver_expressions.begin() ; itarget != _solver_expressions.end() ; ++itarget)
  {
    auto target = itarget->first;

    // generate solver tree string for target, if it was not created earlier
    if (_formatted_solver_trees.count(target) == 0)
      _formatted_solver_trees[target] = _generate_formatted_solver_tree(target);

    // explore solver tree, gather top level timings
    std::vector<String> names;

    if(_solver_expressions.count(target) == 0)
    {
      XABORTM("target "+target+" not present in _solver_expressions");
    }

    if (_solver_expressions[target].front()->get_type() != Solver::ExpressionType::start_solve)
    {
      XABORTM("Should never happen - _solver_expressions list did not start with start solve expression!");
    }

    _overall_toe[target].push_back(0.);
    _overall_iters[target].push_back(Index(0));
    _overall_mpi_execute_reduction[target].push_back(0.);
    _overall_mpi_execute_blas2[target].push_back(0.);
    _overall_mpi_execute_blas3[target].push_back(0.);
    _overall_mpi_execute_collective[target].push_back(0.);
    _overall_mpi_wait_reduction[target].push_back(0.);
    _overall_mpi_wait_blas2[target].push_back(0.);
    _overall_mpi_wait_blas3[target].push_back(0.);
    _overall_mpi_wait_collective[target].push_back(0.);
    _outer_mg_toe[target].emplace_back();
    _outer_mg_mpi_execute_reduction[target].emplace_back();
    _outer_mg_mpi_execute_blas2[target].emplace_back();
    _outer_mg_mpi_execute_blas3[target].emplace_back();
    _outer_mg_mpi_execute_collective[target].emplace_back();
    _outer_mg_mpi_wait_reduction[target].emplace_back();
    _outer_mg_mpi_wait_blas2[target].emplace_back();
    _outer_mg_mpi_wait_blas3[target].emplace_back();
    _outer_mg_mpi_wait_collective[target].emplace_back();
    _outer_schwarz_toe[target].push_back(0.);
    _outer_schwarz_iters[target].push_back(Index(0));

    Index outer_mg_depth(0);
    Index outer_schwarz_depth(0);

    for (auto it = _solver_expressions.at(target).begin() ; it != _solver_expressions.at(target).end() ; ++it)
    {
      auto expression = *it;

      if (expression->get_type() == Solver::ExpressionType::start_solve)
      {
        names.push_back((*it)->solver_name);

        // set outest mg depth to first mg found in solver tree while descending
        if (names.size() > 0 && (names.back().starts_with("MultiGrid") || names.back().starts_with("VCycle") || names.back().starts_with("ScaRCMultiGrid")) && outer_mg_depth == 0)
          outer_mg_depth = (Index)names.size();

        // set depth of schwarz preconditioner in solver tree while descending
        if (names.size() > 0 && names.back().starts_with("Schwarz") && outer_schwarz_depth == 0)
          outer_schwarz_depth = (Index)names.size();
      }

      //fetch iters from top-lvl (lying insided of schwarz or global)
      if (expression->get_type() == Solver::ExpressionType::end_solve && expression->solver_name == names.back())
      {
        if (names.size() > 1 && names.size() == outer_schwarz_depth + 1 && names.at(names.size() - 2) == "Schwarz")
        {
          auto t = dynamic_cast<Solver::ExpressionEndSolve*>(expression.get());
          _outer_schwarz_iters[target].back() += t->iters;
        }
        if (names.size() < 2)
        {
          auto t = dynamic_cast<Solver::ExpressionEndSolve*>(expression.get());

          _overall_iters[target].back() += t->iters;
        }

        names.pop_back();
        continue;
      }

      if ((names.size() < 2) && expression->get_type() == FEAT::Solver::ExpressionType::timings)
      {
        auto t = dynamic_cast<Solver::ExpressionTimings*>(expression.get());

        _overall_toe[target].back() += t->solver_toe;
        _overall_mpi_execute_reduction[target].back() += t->mpi_execute_reduction;
        _overall_mpi_execute_blas2[target].back() += t->mpi_execute_blas2;
        _overall_mpi_execute_blas3[target].back() += t->mpi_execute_blas3;
        _overall_mpi_execute_collective[target].back() += t->mpi_execute_collective;
        _overall_mpi_wait_reduction[target].back() += t->mpi_wait_reduction;
        _overall_mpi_wait_blas2[target].back() += t->mpi_wait_blas2;
        _overall_mpi_wait_blas3[target].back() += t->mpi_wait_blas3;
        _overall_mpi_wait_collective[target].back() += t->mpi_wait_collective;
      }

      if (names.size() == outer_mg_depth && expression->get_type() == FEAT::Solver::ExpressionType::level_timings)
      {
        auto t = dynamic_cast<Solver::ExpressionLevelTimings*>(expression.get());
        // add new vector entry for current level if vector does not already contain an entry for this level
        if (_outer_mg_toe[target].back().size() <= t->level)
        {
          _outer_mg_toe[target].back().push_back(t->level_toe);
          _outer_mg_mpi_execute_reduction[target].back().push_back(t->mpi_execute_reduction);
          _outer_mg_mpi_execute_blas2[target].back().push_back(t->mpi_execute_blas2);
          _outer_mg_mpi_execute_blas3[target].back().push_back(t->mpi_execute_blas3);
          _outer_mg_mpi_execute_collective[target].back().push_back(t->mpi_execute_collective);
          _outer_mg_mpi_wait_reduction[target].back().push_back(t->mpi_wait_reduction);
          _outer_mg_mpi_wait_blas2[target].back().push_back(t->mpi_wait_blas2);
          _outer_mg_mpi_wait_blas3[target].back().push_back(t->mpi_wait_blas3);
          _outer_mg_mpi_wait_collective[target].back().push_back(t->mpi_wait_collective);
        }
        else
        {
          _outer_mg_toe[target].back().at(t->level) += t->level_toe;
          _outer_mg_mpi_execute_reduction[target].back().at(t->level) += t->mpi_execute_reduction;
          _outer_mg_mpi_execute_blas2[target].back().at(t->level) += t->mpi_execute_blas2;
          _outer_mg_mpi_execute_blas3[target].back().at(t->level) += t->mpi_execute_blas3;
          _outer_mg_mpi_execute_collective[target].back().at(t->level) += t->mpi_execute_collective;
          _outer_mg_mpi_wait_reduction[target].back().at(t->level) += t->mpi_wait_reduction;
          _outer_mg_mpi_wait_blas2[target].back().at(t->level) += t->mpi_wait_blas2;
          _outer_mg_mpi_wait_blas3[target].back().at(t->level) += t->mpi_wait_blas3;
          _outer_mg_mpi_wait_collective[target].back().at(t->level) += t->mpi_wait_collective;
        }
      }

      //grep the solver which lies in the schwarz solver to get its toe
      if (names.size() > 1 && names.size() == outer_schwarz_depth + 1 && names.at(names.size() - 2) == "Schwarz" && expression->get_type() == FEAT::Solver::ExpressionType::timings)
      {
        auto t = dynamic_cast<Solver::ExpressionTimings*>(expression.get());
        _outer_schwarz_toe[target].back() += t->solver_toe;
      }
    }
  }

  // clear raw solver expressions
  _solver_expressions.clear();
}
