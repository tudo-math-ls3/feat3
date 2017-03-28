#include <kernel/util/statistics.hpp>
#include <kernel/util/function_scheduler.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/solver/base.hpp>

#include <queue>
#include <numeric>

using namespace FEAT;

// static member initialisation
Index Statistics::_flops = Index(0);
KahanAccumulation Statistics::_time_reduction;
KahanAccumulation Statistics::_time_spmv;
KahanAccumulation Statistics::_time_axpy;
KahanAccumulation Statistics::_time_precon;
KahanAccumulation Statistics::_time_mpi_execute;
KahanAccumulation Statistics::_time_mpi_wait_reduction;
KahanAccumulation Statistics::_time_mpi_wait_spmv;
std::map<String, std::list<std::shared_ptr<Solver::ExpressionBase>>> Statistics::_solver_expressions;
std::map<String, String> Statistics::_formatted_solver_trees;
std::map<String, std::list<double>> Statistics::_overall_toe;
std::map<String, std::list<Index>> Statistics::_overall_iters;
std::map<String, std::list<double>> Statistics::_overall_mpi_execute;
std::map<String, std::list<double>> Statistics::_overall_mpi_wait_reduction;
std::map<String, std::list<double>> Statistics::_overall_mpi_wait_spmv;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_toe;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_execute;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_wait_reduction;
std::map<String, std::list<std::vector<double>>> Statistics::_outer_mg_mpi_wait_spmv;
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
    throw InternalError(__func__, __FILE__, __LINE__,
    "target "+target+" not present in _solver_expressions");
  }

  if (_solver_expressions[target].front()->get_type() != Solver::ExpressionType::start_solve)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "Should never happen - _solver_expressions list did not start with start solve expression!");
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

    // the current solver is of schur complement type, search for a and s solver
    // if both have been found, skip everything until solver end statement has been found
    else if (names.back().starts_with("Schur"))
    {
      while (it != _solver_expressions[target].end())
      {
        auto expression = *it;

        if (expression->get_type() == Solver::ExpressionType::call_schur_s && expression->solver_name == names.back() && (found.back() == 0 || found.back() == 2))
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
            auto t = dynamic_cast<Solver::ExpressionCallSchurS*>(expression.get());
            if (found.back() == 1)
              tree += " ( S: " + t->solver_s_name;
            else
              tree += " / S: " + t->solver_s_name;
          }
        }

        if (expression->get_type() == Solver::ExpressionType::call_schur_a && expression->solver_name == names.back() && (found.back() == 0 || found.back() == 1))
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
            auto t = dynamic_cast<Solver::ExpressionCallSchurA*>(expression.get());
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
      // the current solver is not of multigrid or schur type, i.e. is uses at most one preconditioner, search for its call or the solvers end.
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
    throw InternalError(__func__, __FILE__, __LINE__, "Should never happen - not all solver calls were parsed to the end!");
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
      case Solver::ExpressionType::call_schur_s:
        {
          auto t = dynamic_cast<Solver::ExpressionCallSchurS*>(expression.get());
          String s = stringify(t->get_type()) + "[" + stringify(t->solver_name) + "] (" + t->solver_s_name + ")";
          std::cout<<String(padding, ' ') << s << std::endl;
          break;
        }
      case Solver::ExpressionType::call_schur_a:
        {
          auto t = dynamic_cast<Solver::ExpressionCallSchurA*>(expression.get());
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

/*
void Statistics::write_out_solver_statistics_scheduled(Index rank, size_t la_bytes, size_t domain_bytes, size_t mpi_bytes, Index cells, Index dofs, Index nzes, String filename)
{
  auto func = [&] () { write_out_solver_statistics(rank, la_bytes, domain_bytes, mpi_bytes, cells, dofs, nzes, filename); };
  Util::schedule_function(func, Util::ScheduleMode::clustered);
}*/

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
  measured_time = KahanSum(measured_time, get_time_mpi_wait_reduction());
  measured_time = KahanSum(measured_time, get_time_mpi_wait_spmv());

  if (measured_time.sum > total_time)
    throw InternalError("Accumulated op time (" + stringify(measured_time.sum) + ") is greater as the provided total execution time (" + stringify(total_time) + ") !");

  result += "\n";
  result += "Accumulated op time: " + stringify(measured_time.sum) + "\n";

  result += "\n";

  Dist::Comm comm(Dist::Comm::world());

  double t_local = get_time_reduction();
  double t_max;
  double t_min;
  comm.allreduce(&t_local, &t_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&t_local, &t_min, std::size_t(1), Dist::op_min);
  result += String("Reductions:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_axpy();
  comm.allreduce(&t_local, &t_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&t_local, &t_min, std::size_t(1), Dist::op_min);
  result += String("Blas-1:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_spmv();
  comm.allreduce(&t_local, &t_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&t_local, &t_min, std::size_t(1), Dist::op_min);
  result += String("Blas-2:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_precon();
  comm.allreduce(&t_local, &t_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&t_local, &t_min, std::size_t(1), Dist::op_min);
  result += String("Precon Kernels:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_mpi_execute();
  comm.allreduce(&t_local, &t_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&t_local, &t_min, std::size_t(1), Dist::op_min);
  result += String("MPI Execution:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_mpi_wait_reduction();
  comm.allreduce(&t_local, &t_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&t_local, &t_min, std::size_t(1), Dist::op_min);
  result += String("MPI Wait Reduction:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_mpi_wait_spmv();
  comm.allreduce(&t_local, &t_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&t_local, &t_min, std::size_t(1), Dist::op_min);
  result += String("MPI Wait Blas-2:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = total_time - measured_time.sum;
  comm.allreduce(&t_local, &t_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&t_local, &t_min, std::size_t(1), Dist::op_min);
  result += String("Not covered:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";
  return result;
}

String Statistics::get_formatted_solver_internals(String target)
{
  Dist::Comm comm(Dist::Comm::world());

  double solver_toe(std::accumulate(FEAT::Statistics::get_time_toe(target).begin(),
        FEAT::Statistics::get_time_toe(target).end(), 0.));
  double solver_toe_max;
  double solver_toe_min;
  comm.allreduce(&solver_toe, &solver_toe_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&solver_toe, &solver_toe_min, std::size_t(1), Dist::op_min);
  Index solver_iters(std::accumulate(FEAT::Statistics::get_iters(target).begin(),
        FEAT::Statistics::get_iters(target).end(), Index(0)));
  Index solver_iters_max;
  Index solver_iters_min;
  comm.allreduce(&solver_iters, &solver_iters_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&solver_iters, &solver_iters_min, std::size_t(1), Dist::op_min);
  double solver_mpi_execute(std::accumulate(FEAT::Statistics::get_time_mpi_execute(target).begin(),
        FEAT::Statistics::get_time_mpi_execute(target).end(), 0.));
  double solver_mpi_execute_max;
  double solver_mpi_execute_min;
  comm.allreduce(&solver_mpi_execute, &solver_mpi_execute_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&solver_mpi_execute, &solver_mpi_execute_min, std::size_t(1), Dist::op_min);
  double solver_mpi_wait_reduction(std::accumulate(FEAT::Statistics::get_time_mpi_wait_reduction(target).begin(),
        FEAT::Statistics::get_time_mpi_wait_reduction(target).end(), 0.));
  double solver_mpi_wait_reduction_max;
  double solver_mpi_wait_reduction_min;
  comm.allreduce(&solver_mpi_wait_reduction, &solver_mpi_wait_reduction_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&solver_mpi_wait_reduction, &solver_mpi_wait_reduction_min, std::size_t(1), Dist::op_min);
  double solver_mpi_wait_spmv(std::accumulate(FEAT::Statistics::get_time_mpi_wait_spmv(target).begin(),
        FEAT::Statistics::get_time_mpi_wait_spmv(target).end(), 0.));
  double solver_mpi_wait_spmv_max;
  double solver_mpi_wait_spmv_min;
  comm.allreduce(&solver_mpi_wait_spmv, &solver_mpi_wait_spmv_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&solver_mpi_wait_spmv, &solver_mpi_wait_spmv_min, std::size_t(1), Dist::op_min);

  auto solver_time_mg = FEAT::Statistics::get_time_mg(target);
  std::vector<double> solver_mg_toe(solver_time_mg.front().size(), double(0.));
  for (auto& step : solver_time_mg)
  {
    for (Index i(0) ; i < step.size() ; ++i)
    {
      solver_mg_toe.at(i) += step.at(i);
    }
  }
  std::vector<double> solver_mg_toe_max(solver_mg_toe.size(), double(0.));
  std::vector<double> solver_mg_toe_min(solver_mg_toe.size(), double(0.));
  for (Index i(0) ; i < solver_mg_toe.size() ; ++i)
  {
    comm.allreduce(&solver_mg_toe.at(i), &solver_mg_toe_max.at(i), std::size_t(1), Dist::op_max);
    comm.allreduce(&solver_mg_toe.at(i), &solver_mg_toe_min.at(i), std::size_t(1), Dist::op_min);
  }

  auto solver_time_mg_mpi_execute = FEAT::Statistics::get_time_mg_mpi_execute(target);
  std::vector<double> solver_mg_mpi_execute(solver_time_mg_mpi_execute.front().size(), double(0.));
  for (auto& step : solver_time_mg_mpi_execute)
  {
    for (Index i(0) ; i < step.size() ; ++i)
    {
      solver_mg_mpi_execute.at(i) += step.at(i);
    }
  }
  std::vector<double> solver_mg_mpi_execute_max(solver_mg_mpi_execute.size(), double(0.));
  std::vector<double> solver_mg_mpi_execute_min(solver_mg_mpi_execute.size(), double(0.));
  for (Index i(0) ; i < solver_mg_mpi_execute.size() ; ++i)
  {
    comm.allreduce(&solver_mg_mpi_execute.at(i), &solver_mg_mpi_execute_max.at(i), std::size_t(1), Dist::op_max);
    comm.allreduce(&solver_mg_mpi_execute.at(i), &solver_mg_mpi_execute_min.at(i), std::size_t(1), Dist::op_min);
  }

  auto solver_time_mg_mpi_wait_reduction = FEAT::Statistics::get_time_mg_mpi_wait_reduction(target);
  std::vector<double> solver_mg_mpi_wait_reduction(solver_time_mg_mpi_wait_reduction.front().size(), double(0.));
  for (auto& step : solver_time_mg_mpi_wait_reduction)
  {
    for (Index i(0) ; i < step.size() ; ++i)
    {
      solver_mg_mpi_wait_reduction.at(i) += step.at(i);
    }
  }
  std::vector<double> solver_mg_mpi_wait_reduction_max(solver_mg_mpi_wait_reduction.size(), double(0.));
  std::vector<double> solver_mg_mpi_wait_reduction_min(solver_mg_mpi_wait_reduction.size(), double(0.));
  for (Index i(0) ; i < solver_mg_mpi_wait_reduction.size() ; ++i)
  {
    comm.allreduce(&solver_mg_mpi_wait_reduction.at(i), &solver_mg_mpi_wait_reduction_max.at(i), std::size_t(1), Dist::op_max);
    comm.allreduce(&solver_mg_mpi_wait_reduction.at(i), &solver_mg_mpi_wait_reduction_min.at(i), std::size_t(1), Dist::op_min);
  }

  auto solver_time_mg_mpi_wait_spmv = FEAT::Statistics::get_time_mg_mpi_wait_spmv(target);
  std::vector<double> solver_mg_mpi_wait_spmv(solver_time_mg_mpi_wait_spmv.front().size(), double(0.));
  for (auto& step : solver_time_mg_mpi_wait_spmv)
  {
    for (Index i(0) ; i < step.size() ; ++i)
    {
      solver_mg_mpi_wait_spmv.at(i) += step.at(i);
    }
  }
  std::vector<double> solver_mg_mpi_wait_spmv_max(solver_mg_mpi_wait_spmv.size(), double(0.));
  std::vector<double> solver_mg_mpi_wait_spmv_min(solver_mg_mpi_wait_spmv.size(), double(0.));
  for (Index i(0) ; i < solver_mg_mpi_wait_spmv.size() ; ++i)
  {
    comm.allreduce(&solver_mg_mpi_wait_spmv.at(i), &solver_mg_mpi_wait_spmv_max.at(i), std::size_t(1), Dist::op_max);
    comm.allreduce(&solver_mg_mpi_wait_spmv.at(i), &solver_mg_mpi_wait_spmv_min.at(i), std::size_t(1), Dist::op_min);
  }

  double solver_schwarz_toe(std::accumulate(FEAT::Statistics::get_time_schwarz(target).begin(),
        FEAT::Statistics::get_time_schwarz(target).end(), 0.));
  double solver_schwarz_toe_max;
  double solver_schwarz_toe_min;
  comm.allreduce(&solver_schwarz_toe, &solver_schwarz_toe_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&solver_schwarz_toe, &solver_schwarz_toe_min, std::size_t(1), Dist::op_min);

  Index solver_schwarz_iters(std::accumulate(FEAT::Statistics::get_iters_schwarz(target).begin(),
        FEAT::Statistics::get_iters_schwarz(target).end(), Index(0)));
  solver_schwarz_iters = solver_schwarz_iters / (Index)FEAT::Statistics::get_iters_schwarz(target).size();
  Index solver_schwarz_iters_max;
  Index solver_schwarz_iters_min;
  comm.allreduce(&solver_schwarz_iters, &solver_schwarz_iters_max, std::size_t(1), Dist::op_max);
  comm.allreduce(&solver_schwarz_iters, &solver_schwarz_iters_min, std::size_t(1), Dist::op_min);


  String result = target + "\n";
  result += String("toe:").pad_back(20) + "max: " + stringify(solver_toe_max) + ", min: " + stringify(solver_toe_min) + ", local: " +
      stringify(solver_toe) + "\n";
  result += String("ites:").pad_back(20) + "max: " + stringify(solver_iters_max) + ", min: " + stringify(solver_iters_min) + ", local: " +
      stringify(solver_iters) + "\n";
  result += String("mpi execute:").pad_back(20) + "max: " + stringify(solver_mpi_execute_max) + ", min: " + stringify(solver_mpi_execute_min) + ", local: " +
      stringify(solver_mpi_execute) + "\n";
  result += String("mpi wait reduction:").pad_back(20) + "max: " + stringify(solver_mpi_wait_reduction_max) + ", min: " + stringify(solver_mpi_wait_reduction_min) + ", local: " +
      stringify(solver_mpi_wait_reduction) + "\n";
  result += String("mpi wait spmv:").pad_back(20) + "max: " + stringify(solver_mpi_wait_spmv_max) + ", min: " + stringify(solver_mpi_wait_spmv_min) + ", local: " +
      stringify(solver_mpi_wait_spmv) + "\n";
  for (Index i(0) ; i < solver_mg_toe.size() ; ++i)
  {
    result += String("toe lvl ") + stringify(i).pad_back(12) + "max: " + stringify(solver_mg_toe_max.at(i)) + ", min: " + stringify(solver_mg_toe_min.at(i)) + ", local: " +
        stringify(solver_mg_toe.at(i)) + "\n";
  }
  for (Index i(0) ; i < solver_mg_toe.size() ; ++i)
  {
    result += String("mpi execute lvl ") + stringify(i).pad_back(4) + "max: " + stringify(solver_mg_mpi_execute_max.at(i)) + ", min: " + stringify(solver_mg_mpi_execute_min.at(i)) +
        ", local: " + stringify(solver_mg_mpi_execute.at(i)) + "\n";
  }
  for (Index i(0) ; i < solver_mg_toe.size() ; ++i)
  {
    result += String("mpi wait red lvl ") + stringify(i).pad_back(3) + "max: " + stringify(solver_mg_mpi_wait_reduction_max.at(i)) + ", min: " + stringify(solver_mg_mpi_wait_reduction_min.at(i)) +
        ", local: " + stringify(solver_mg_mpi_wait_reduction.at(i)) + "\n";
  }
  for (Index i(0) ; i < solver_mg_toe.size() ; ++i)
  {
    result += String("mpi wait spmv lvl ") + stringify(i).pad_back(2) + "max: " + stringify(solver_mg_mpi_wait_spmv_max.at(i)) + ", min: " + stringify(solver_mg_mpi_wait_spmv_min.at(i)) +
        ", local: " + stringify(solver_mg_mpi_wait_spmv.at(i)) + "\n";
  }
  result += String("schwarz toe:").pad_back(20) + "max: " + stringify(solver_schwarz_toe_max) + ", min: " + stringify(solver_schwarz_toe_min) + ", local: " +
      stringify(solver_schwarz_toe) + "\n";
  result += String("schwarz iters:").pad_back(20) + "max: " + stringify(solver_schwarz_iters_max) + ", min: " + stringify(solver_schwarz_iters_min) + ", local: " +
      stringify(solver_schwarz_iters);

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
      throw InternalError(__func__, __FILE__, __LINE__,
          "target "+target+" not present in _solver_expressions");
    }

    if (_solver_expressions[target].front()->get_type() != Solver::ExpressionType::start_solve)
    {
      throw InternalError(__func__, __FILE__, __LINE__, "Should never happen - _solver_expressions list did not start with start solve expression!");
    }

    _overall_toe[target].push_back(0.);
    _overall_iters[target].push_back(Index(0));
    _overall_mpi_execute[target].push_back(0.);
    _overall_mpi_wait_reduction[target].push_back(0.);
    _overall_mpi_wait_spmv[target].push_back(0.);
    _outer_mg_toe[target].emplace_back();
    _outer_mg_mpi_execute[target].emplace_back();
    _outer_mg_mpi_wait_reduction[target].emplace_back();
    _outer_mg_mpi_wait_spmv[target].emplace_back();
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
        _overall_mpi_execute[target].back() += t->mpi_execute;
        _overall_mpi_wait_reduction[target].back() += t->mpi_wait_reduction;
        _overall_mpi_wait_spmv[target].back() += t->mpi_wait_spmv;
      }

      if (names.size() == outer_mg_depth && expression->get_type() == FEAT::Solver::ExpressionType::level_timings)
      {
        auto t = dynamic_cast<Solver::ExpressionLevelTimings*>(expression.get());
        // add new vector entry for current level if vector does not already contain an entry for this level
        if (_outer_mg_toe[target].back().size() <= t->level)
        {
          _outer_mg_toe[target].back().push_back(t->level_toe);
          _outer_mg_mpi_execute[target].back().push_back(t->mpi_execute);
          _outer_mg_mpi_wait_reduction[target].back().push_back(t->mpi_wait_reduction);
          _outer_mg_mpi_wait_spmv[target].back().push_back(t->mpi_wait_spmv);
        }
        else
        {
          _outer_mg_toe[target].back().at(t->level) += t->level_toe;
          _outer_mg_mpi_execute[target].back().at(t->level) += t->mpi_execute;
          _outer_mg_mpi_wait_reduction[target].back().at(t->level) += t->mpi_wait_reduction;
          _outer_mg_mpi_wait_spmv[target].back().at(t->level) += t->mpi_wait_spmv;
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
