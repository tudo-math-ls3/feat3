#pragma once
#ifndef KERNEL_SOLVER_EXPRESSION_HPP
#define KERNEL_SOLVER_EXPRESSION_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/string.hpp>

namespace FEAT
{
  namespace Solver
  {
    // forward declaration
    enum class Status;

    /// Expression Type enumeration
    enum class ExpressionType
    {
      /// start new solve process
      start_solve = 0,
      /// end last solve process
      end_solve,
      /// call preconditioner
      call_precond,
      /// call L preconditioner
      call_precond_l,
      /// call R preconditioner
      call_precond_r,
      /// call smoother
      call_smoother,
      /// call coarse grid solver
      call_coarse_solver,
      /// annotate iterations defect
      defect,
      /// annotate iterations timings
      timings,
      /// annotate level timings
      level_timings,
      /// prolonation (multigrid)
      prol,
      /// restriction (multigrid)
      rest,
      /// call S matrix solver (schur complement)
      call_schur_s,
      /// call A matrix solver (schur complement)
      call_schur_a
    };

    /// \cond internal
    inline std::ostream& operator<<(std::ostream& os, ExpressionType type)
    {
      switch(type)
      {
        case ExpressionType::start_solve:
        return os << "start";
        case ExpressionType::end_solve:
        return os << "end";
        case ExpressionType::call_precond:
        return os << "precond";
        case ExpressionType::call_precond_l:
        return os << "precond_l";
        case ExpressionType::call_precond_r:
        return os << "precond_r";
        case ExpressionType::call_smoother:
        return os << "smoother";
        case ExpressionType::call_coarse_solver:
        return os << "coarse";
        case ExpressionType::defect:
        return os << "defect";
        case ExpressionType::timings:
        return os << "timings";
        case ExpressionType::level_timings:
        return os << "level_timings";
        case ExpressionType::prol:
        return os << "prol";
        case ExpressionType::rest:
        return os << "rest";
        case ExpressionType::call_schur_s:
        return os << "schur_s";
        case ExpressionType::call_schur_a:
        return os << "schur_a";
      default:
        return os << "unknown";
      }
    }
    /// \endcond

    class ExpressionBase
    {
      public:
        String solver_name;

        explicit ExpressionBase(String name) :
          solver_name(name)
        {
        }

        virtual ~ExpressionBase()
        {
        }

        virtual ExpressionType get_type() = 0;
    };

    class ExpressionStartSolve : public ExpressionBase
    {
      public:
        explicit ExpressionStartSolve(String name) :
          ExpressionBase(name)
        {
        }

        virtual ~ExpressionStartSolve()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::start_solve;
        }
    };

    class ExpressionEndSolve : public ExpressionBase
    {
      public:
        /// the status result of the solver call
        Status status;
        /// the iteration count needed in this solve process
        Index iters;

        explicit ExpressionEndSolve(String name, Status end_status, Index iterations) :
          ExpressionBase(name),
          status(end_status),
          iters(iterations)
        {
        }
        virtual ~ExpressionEndSolve()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::end_solve;
        }
    };

    class ExpressionCallPrecond : public ExpressionBase
    {
      public:
        String precond_name;

        explicit ExpressionCallPrecond(String name, String precond_name_in) :
          ExpressionBase(name),
          precond_name(precond_name_in)
        {
        }

        virtual ~ExpressionCallPrecond()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::call_precond;
        }
    };

    class ExpressionCallPrecondL : public ExpressionBase
    {
      public:
        String precond_name;

        explicit ExpressionCallPrecondL(String name, String precond_name_in) :
          ExpressionBase(name),
          precond_name(precond_name_in)
        {
        }

        virtual ~ExpressionCallPrecondL()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::call_precond_l;
        }
    };

    class ExpressionCallPrecondR : public ExpressionBase
    {
      public:
        String precond_name;

        explicit ExpressionCallPrecondR(String name, String precond_name_in) :
          ExpressionBase(name),
          precond_name(precond_name_in)
        {
        }

        virtual ~ExpressionCallPrecondR()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::call_precond_r;
        }
    };

    class ExpressionCallSmoother : public ExpressionBase
    {
      public:
        String smoother_name;

        explicit ExpressionCallSmoother(String name, String smoother_name_in) :
          ExpressionBase(name),
          smoother_name(smoother_name_in)
        {
        }

        virtual ~ExpressionCallSmoother()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::call_smoother;
        }
    };

    class ExpressionCallCoarseSolver : public ExpressionBase
    {
      public:
        String coarse_solver_name;

        explicit ExpressionCallCoarseSolver(String name, String coarse_solver_name_in) :
          ExpressionBase(name),
          coarse_solver_name(coarse_solver_name_in)
        {
        }

        virtual ~ExpressionCallCoarseSolver()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::call_coarse_solver;
        }
    };

    class ExpressionProlongation : public ExpressionBase
    {
      public:
        Index level;

        explicit ExpressionProlongation(String name, Index level_in) :
          ExpressionBase(name),
          level(level_in)
        {
        }

        virtual ~ExpressionProlongation()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::prol;
        }
    };

    class ExpressionRestriction : public ExpressionBase
    {
      public:
        Index level;

        explicit ExpressionRestriction(String name, Index level_in) :
          ExpressionBase(name),
          level(level_in)
        {
        }

        virtual ~ExpressionRestriction()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::rest;
        }
    };

    class ExpressionDefect : public ExpressionBase
    {
      public:
        double def;
        Index iter;

        explicit ExpressionDefect(String name, double defect, Index iter_in) :
          ExpressionBase(name),
          def(defect),
          iter(iter_in)
        {
        }

        virtual ~ExpressionDefect()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::defect;
        }
    };

    class ExpressionTimings : public ExpressionBase
    {
      public:
        double solver_toe, mpi_execute, mpi_wait_reduction, mpi_wait_spmv;

        explicit ExpressionTimings(String name, double solver_toe_in, double mpi_execute_in, double mpi_wait_reduction_in, double mpi_wait_spmv_in) :
          ExpressionBase(name),
          solver_toe(solver_toe_in),
          mpi_execute(mpi_execute_in),
          mpi_wait_reduction(mpi_wait_reduction_in),
          mpi_wait_spmv(mpi_wait_spmv_in)
        {
        }

        virtual ~ExpressionTimings()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::timings;
        }
    };

    class ExpressionLevelTimings : public ExpressionBase
    {
      public:
        Index level;
        double level_toe, mpi_execute, mpi_wait_reduction, mpi_wait_spmv;

        explicit ExpressionLevelTimings(String name, Index level_in, double level_toe_in, double mpi_execute_in, double mpi_wait_reduction_in, double mpi_wait_spmv_in) :
          ExpressionBase(name),
          level(level_in),
          level_toe(level_toe_in),
          mpi_execute(mpi_execute_in),
          mpi_wait_reduction(mpi_wait_reduction_in),
          mpi_wait_spmv(mpi_wait_spmv_in)
        {
        }

        virtual ~ExpressionLevelTimings()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::level_timings;
        }
    };

    class ExpressionCallSchurS : public ExpressionBase
    {
      public:
        String solver_s_name;

        explicit ExpressionCallSchurS(String name, String solver_s_name_in) :
          ExpressionBase(name),
          solver_s_name(solver_s_name_in)
        {
        }

        virtual ~ExpressionCallSchurS()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::call_schur_s;
        }
    };

    class ExpressionCallSchurA : public ExpressionBase
    {
      public:
        String solver_a_name;

        explicit ExpressionCallSchurA(String name, String solver_a_name_in) :
          ExpressionBase(name),
          solver_a_name(solver_a_name_in)
        {
        }

        virtual ~ExpressionCallSchurA()
        {
        }

        virtual ExpressionType get_type() override
        {
          return ExpressionType::call_schur_a;
        }
    };
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_EXPRESSION_HPP