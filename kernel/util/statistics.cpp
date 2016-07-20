#include <kernel/util/statistics.hpp>
#include <kernel/util/function_scheduler.hpp>
#include <kernel/solver/base.hpp>

#include <queue>

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
std::list<std::shared_ptr<Solver::ExpressionBase>> Statistics::_solver_expressions;
double Statistics::toe_partition;
double Statistics::toe_assembly;
double Statistics::toe_solve;

String Statistics::get_formatted_solver_tree()
{
  std::list<String> names;
  std::list<int> found; //number of found preconds / smoothers / coarse solvers. 0 = nothing found, 1 = smoother / s found, 2 = coarse solver / a found, 3 = all found

  if (_solver_expressions.front()->get_type() != Solver::ExpressionType::start_solve)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "Should never happen - _solver_expressions list did not start with start solve expression!");
  }

  // process the very first entry, e.g. the outer most solver
  auto it = _solver_expressions.begin();
  String tree((*it)->solver_name);
  names.push_back((*it)->solver_name);
  found.push_back(0);

  // process current last element in the list, until no element needs to be processed
  while (names.size() > 0 && it != _solver_expressions.end())
  {
    ++it;

    // the current solver is of multigrid type, search for smoother and coarse solver (smoother allways comes first).
    // if smoother and coarse solver have been found, skip everything until solver end statement has been found
    if (names.back().starts_with("MultiGrid") || names.back().starts_with("VCycle") || names.back().starts_with("ScaRCMultiGrid"))
    {
      while (it != _solver_expressions.end())
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
      while (it != _solver_expressions.end())
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
      while (it != _solver_expressions.end())
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
      while (it != _solver_expressions.end())
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

  for (auto expression : _solver_expressions)
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

  double t_local = get_time_reduction();
  double t_max;
  double t_min;
  Util::Comm::allreduce(&t_local, &t_max, 1, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, &t_min, 1, Util::CommOperationMin());
  result += String("Reductions:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_axpy();
  Util::Comm::allreduce(&t_local, &t_max, 1, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, &t_min, 1, Util::CommOperationMin());
  result += String("Blas-1:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_spmv();
  Util::Comm::allreduce(&t_local, &t_max, 1, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, &t_min, 1, Util::CommOperationMin());
  result += String("Blas-2:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_precon();
  Util::Comm::allreduce(&t_local, &t_max, 1, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, &t_min, 1, Util::CommOperationMin());
  result += String("Precon Kernels:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_mpi_execute();
  Util::Comm::allreduce(&t_local, &t_max, 1, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, &t_min, 1, Util::CommOperationMin());
  result += String("MPI Execution:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_mpi_wait_reduction();
  Util::Comm::allreduce(&t_local, &t_max, 1, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, &t_min, 1, Util::CommOperationMin());
  result += String("MPI Wait Reduction:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = get_time_mpi_wait_spmv();
  Util::Comm::allreduce(&t_local, &t_max, 1, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, &t_min, 1, Util::CommOperationMin());
  result += String("MPI Wait Blas-2:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";

  t_local = total_time - measured_time.sum;
  Util::Comm::allreduce(&t_local, &t_max, 1, Util::CommOperationMax());
  Util::Comm::allreduce(&t_local, &t_min, 1, Util::CommOperationMin());
  result += String("Not covered:").pad_back(20) + "max: " + stringify(t_max) + ", min: " + stringify(t_min) + ", local: " + stringify(t_local) + "\n";
  return result;
}
