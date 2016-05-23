#pragma once
#ifndef CONTROL_MESHOPT_MESHOPT_SOLVER_FACTORY_HPP
#define CONTROL_MESHOPT_MESHOPT_SOLVER_FACTORY_HPP 1

#include <kernel/solver/alglib_wrapper.hpp>
#include <kernel/solver/basic_vcycle.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/solver/nlsd.hpp>
#include <kernel/solver/nlopt_precond.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>

#include <deque>

namespace FEAST
{
  namespace Control
  {
    struct MeshoptSolverFactory
    {
      private:
        /// \todo make internal
        template<typename VectorType_>
        static void configure_iterative_solver(PropertyMap * section, std::shared_ptr<Solver::IterativeSolver<VectorType_> > solver)
        {
          auto plot_p = section->get_entry("plot");
          if (plot_p.second)
          {
            Index plot(std::stoul(plot_p.first));
            if (plot == 0)
            {
              solver->set_plot(false);
            }
            else if (plot == 1)
            {
              solver->set_plot(Util::Comm::rank()== 0);
            }
            else
            {
              throw InternalError(__func__, __FILE__, __LINE__, "plot value " + stringify(plot) + " unknown!");
            }
          }

          auto tol_rel_p = section->get_entry("tol_rel");
          if (tol_rel_p.second)
            solver->set_tol_rel(std::stod(tol_rel_p.first));

          auto tol_abs_p = section->get_entry("tol_abs");
          if (tol_abs_p.second)
            solver->set_tol_abs(std::stod(tol_abs_p.first));

          auto max_iter_p = section->get_entry("max_iter");
          if (max_iter_p.second)
            solver->set_max_iter(std::stoul(max_iter_p.first));

          auto min_iter_p = section->get_entry("min_iter");
          if (min_iter_p.second)
            solver->set_min_iter(std::stoul(min_iter_p.first));
        }

        template <typename VectorType_>
        static void configure_iterative_solver(PropertyMap * section, std::shared_ptr<Solver::PreconditionedIterativeSolver<VectorType_> > solver)
        {
          auto plot_p = section->get_entry("plot");
          if (plot_p.second)
          {
            Index plot(std::stoul(plot_p.first));
            if (plot == 0)
            {
              solver->set_plot(false);
            }
            else if (plot == 1)
            {
              solver->set_plot(Util::Comm::rank() == 0);
            }
            else
            {
              throw InternalError(__func__, __FILE__, __LINE__, "plot value " + stringify(plot) + " unknown!");
            }
          }

          auto tol_rel_p = section->get_entry("tol_rel");
          if (tol_rel_p.second)
            solver->set_tol_rel(std::stod(tol_rel_p.first));

          auto tol_abs_p = section->get_entry("tol_abs");
          if (tol_abs_p.second)
            solver->set_tol_abs(std::stod(tol_abs_p.first));

          auto max_iter_p = section->get_entry("max_iter");
          if (max_iter_p.second)
            solver->set_max_iter(std::stoul(max_iter_p.first));

          auto min_iter_p = section->get_entry("min_iter");
          if (min_iter_p.second)
            solver->set_min_iter(std::stoul(min_iter_p.first));
        }

        /// \todo make internal
        static String get_section_path(PropertyMap * base, PropertyMap * section, String path, String name)
        {
          if (section->get_sub_section(name) != nullptr)
          {
            return path + "/" + name;
          }
          else if (base->get_sub_section(name) != nullptr)
          {
            return name;
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "section " + name + " not found!");
        }

        /// returns the object, if T_ has a GateType, i.e. is a GlobalVector - SFINAE at its best
        template <typename Evaluator_, typename T_>
        static T_ & derefer(T_ & object, typename Evaluator_::GateType *)
        {
          return object;
        }

        /// returns the dereferenced object, if T_ has no GateType, i.e. is a LocalVector - SFINAE at its best
        template <typename Evaluator_, typename T_>
        static auto derefer(T_ & object, ...) -> decltype(*object)
        {
          return *object;
        }

        template <typename Evaluator_, typename SystemLevelType_, typename TransferLevelType_>
        static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVectorR> >
        create_schwarz_precon(std::deque<SystemLevelType_*> & system_levels, std::deque<TransferLevelType_*> & transfer_levels, PropertyMap * base, String solver_name, PropertyMap * section, typename Evaluator_::GateType *)
        {
          typedef typename SystemLevelType_::GlobalSystemVectorR SolverVectorType;
          std::shared_ptr<Solver::SolverBase<typename SolverVectorType::LocalVectorType> > precon_schwarz;
          auto schwarz_p = section->query("solver");
          if (schwarz_p.second)
          {
            precon_schwarz = create_linear_solver<SystemLevelType_, TransferLevelType_, typename SolverVectorType::LocalVectorType>(system_levels, transfer_levels, base, get_section_path(base, section, solver_name, schwarz_p.first));
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "Schwarz precon section without solver key is not allowed!");

          return Solver::new_schwarz_precond(precon_schwarz, system_levels.back()->filter_sys);
        }

        template <typename Evaluator_, typename SystemLevelType_, typename TransferLevelType_>
        static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVectorR::LocalVectorType> >
        create_schwarz_precon(std::deque<SystemLevelType_*> & /*system_levels*/, std::deque<TransferLevelType_*> & /*transfer_levels*/, PropertyMap * /*base*/, String /*solver_name*/,
        PropertyMap * /*section*/, ...)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Schwarz precon section is only allowed in global context! Maybe you have two in one solver branch?");
          return nullptr;
        }

        //template <typename Evaluator_, typename SystemLevelType_>
        //static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVectorR> >
        //create_ilu_precon(std::deque<SystemLevelType_*> & /*system_levels*/, typename Evaluator_::GateType *)
        //{
        //  throw InternalError(__func__, __FILE__, __LINE__, "ilu precon section is only allowed in local context!");
        //  return nullptr;
        //}

        //template <typename Evaluator_, typename SystemLevelType_>
        //static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVectorR::LocalVectorType> >
        //create_ilu_precon(std::deque<SystemLevelType_*> & system_levels, ...)
        //{
        //  auto result = Solver::new_ilu_precond(*system_levels.back()->op_sys, *system_levels.back()->filter_sys, 0ul);
        //  return result;
        //}

        //template <typename Evaluator_, typename SystemLevelType_>
        //static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVectorR> >
        //create_ssor_precon(std::deque<SystemLevelType_*> & /*system_levels*/, PropertyMap * /*section*/, typename Evaluator_::GateType *)
        //{
        //  throw InternalError(__func__, __FILE__, __LINE__, "ssor precon section is only allowed in local context!");
        //  return nullptr;
        //}

        //template <typename Evaluator_, typename SystemLevelType_>
        //static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVectorR::LocalVectorType> >
        //create_ssor_precon(std::deque<SystemLevelType_*> & system_levels, PropertyMap * section, ...)
        //{
        //  auto omega_p = section->get_entry("omega");
        //  if (omega_p.second)
        //  {
        //    auto result = Solver::new_ssor_precond(*system_levels.back()->op_sys, *system_levels.back()->filter_sys, std::stod(omega_p.first));
        //    return result;
        //  }
        //  else
        //  {
        //    auto result = Solver::new_ssor_precond(*system_levels.back()->op_sys, *system_levels.back()->filter_sys);
        //    return result;
        //  }
        //}

      public:

        /**
         * \brief Create solver tree based on PropertyMap
         *
         * \param[in] system_levels SystemLevel hierarchy
         * \param[in] transfer_levels TransferLevel hierarchy
         * \param[in] base A pointer to the PropertyMap that contains all solver related informations
         * \param[in] solver_name The name of the solver tree's root section
         */
        template <typename SystemLevelType_, typename TransferLevelType_, typename SolverVectorType_ = typename SystemLevelType_::GlobalSystemVectorR>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_linear_solver(std::deque<SystemLevelType_*> & system_levels,
        std::deque<TransferLevelType_*> & transfer_levels, PropertyMap * base, String solver_name)
        {
          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > result;

          auto section = base->query_section(solver_name);
          if(section == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,"No solver section "+solver_name+" found");

          auto solver_p = section->query("type");
          if (!solver_p.second)
            throw InternalError(__func__, __FILE__, __LINE__, "no type key found in property map: " + solver_name + "!");
          String solver_type = solver_p.first;

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon = nullptr;
          auto precon_p = section->query("precon");
          if (precon_p.second)
          {
            if (precon_p.first == "none")
              precon = nullptr;
            else
              precon = create_linear_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(system_levels, transfer_levels, base, get_section_path(base, section, solver_name, precon_p.first));
          }

          if (solver_type == "pcg")
          {
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_pcg(derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "bicgstab")
          {
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_bicgstab(derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "fgmres")
          {
            auto krylov_dim_p = section->get_entry("krylov_dim");
            Index krylov_dim;
            if (krylov_dim_p.second)
              krylov_dim = std::stoul(krylov_dim_p.first);
            else
              throw InternalError(__func__, __FILE__, __LINE__, "no krylov_dim key found in fgmres section!");

            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_fgmres(derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), krylov_dim, 0.0, precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "richardson")
          {
            auto omega_p = section->get_entry("omega");
            double omega;
            if (omega_p.second)
              omega = stod(omega_p.first);
            else
              omega = 1.0;

            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_richardson(derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), omega, precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "jac")
          {
            auto omega_p = section->get_entry("omega");
            if (omega_p.second)
            {
              result = Solver::new_jacobi_precond(derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), std::stod(omega_p.first));
            }
            else
            {
              result = Solver::new_jacobi_precond(derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr));
            }
          }
          else if (solver_type == "scale")
          {
            auto omega_p = section->get_entry("omega");
            if (omega_p.second)
            {
              result = Solver::new_scale_precond(derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), std::stod(omega_p.first));
            }
            else
              throw InternalError(__func__, __FILE__, __LINE__, "no omega key found in scale section!");
          }
          //else if (solver_type == "ilu")
          //{
          //  result = create_ilu_precon<SolverVectorType_>(system_levels, nullptr);
          //}
          //else if (solver_type == "ssor")
          //{
          //  result = create_ssor_precon<SolverVectorType_>(system_levels, section, nullptr);
          //}
          else if (solver_type == "mgv")
          {
            typename std::remove_reference<decltype(derefer<SolverVectorType_>(system_levels.front()->op_sys, nullptr))>::type dummy_sys;
            typename std::remove_reference<decltype(derefer<SolverVectorType_>(system_levels.front()->filter_sys, nullptr))>::type dummy_filter;
            typename std::remove_reference<decltype(derefer<SolverVectorType_>(transfer_levels.front()->prol_sys, nullptr))>::type dummy_prol;
            auto mgv = std::make_shared<
              Solver::BasicVCycle<
              decltype(dummy_sys),
              decltype(dummy_filter),
              decltype(dummy_prol)
            > >();

            std::shared_ptr<Solver::SolverBase<SolverVectorType_> > coarse_solver;
            auto coarse_solver_p = section->query("coarse");
            if (coarse_solver_p.second)
            {
              //artificial deque containing the coarsest level only
              typename std::remove_reference<decltype(system_levels)>::type coarse_system_level(system_levels.begin(), ++system_levels.begin());
              typename std::remove_reference<decltype(transfer_levels)>::type coarse_transfer_level(transfer_levels.begin(), ++transfer_levels.begin());
              coarse_solver = create_linear_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(coarse_system_level, coarse_transfer_level, base, get_section_path(base, section, solver_name, coarse_solver_p.first));
            }
            else
              throw InternalError(__func__, __FILE__, __LINE__, "mgv precon section without coarse key is not allowed!");

            mgv->set_coarse_level(derefer<SolverVectorType_>(system_levels.front()->op_sys, nullptr), derefer<SolverVectorType_>(system_levels.front()->filter_sys, nullptr), coarse_solver);

            auto jt = transfer_levels.begin();
            for (auto it = ++system_levels.begin(); it != system_levels.end(); ++it, ++jt)
            {
              std::shared_ptr<Solver::SolverBase<SolverVectorType_> > smoother;
              auto smoother_p = section->query("smoother");
              if (smoother_p.second)
              {
                //artificial deque containing all levels up to the current level, that shall be smoothed
                typename std::remove_reference<decltype(system_levels)>::type smoother_system_levels(system_levels.begin(), it+1);
                typename std::remove_reference<decltype(transfer_levels)>::type smoother_transfer_levels(transfer_levels.begin(), jt+1);
                smoother = create_linear_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(smoother_system_levels, smoother_transfer_levels, base, get_section_path(base, section, solver_name, smoother_p.first));
              }
              else
                throw InternalError(__func__, __FILE__, __LINE__, "mgv precon section without smoother key is not allowed!");

              mgv->push_level(derefer<SolverVectorType_>((*it)->op_sys, nullptr), derefer<SolverVectorType_>((*it)->filter_sys, nullptr), derefer<SolverVectorType_>((*jt)->prol_sys, nullptr),
              derefer<SolverVectorType_>((*jt)->rest_sys, nullptr), smoother, smoother);
            }

            result = mgv;
          }
          else if (solver_type == "schwarz")
          {
            result = create_schwarz_precon<SolverVectorType_>(system_levels, transfer_levels, base, solver_name, section, nullptr);
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "solver with type " + solver_type + " unkown!");

          return result;
        }

        template
        <
          typename SystemLevelType_,
          /*typename TransferLevelType_,*/
          typename SolverVectorType_ = typename SystemLevelType_::GlobalSystemVectorR
        >
        static std::shared_ptr<Solver::Linesearch<typename SystemLevelType_::GlobalQualityFunctional, typename SystemLevelType_::GlobalSystemFilter> >
        create_linesearch( std::deque<SystemLevelType_*> & system_levels,
        /* std::deque<TransferLevelType_*>& transfer_levels, */
        PropertyMap* base, String solver_name)
        {
          typedef typename SystemLevelType_::GlobalQualityFunctional OperatorType;
          typedef typename SystemLevelType_::GlobalSystemFilter FilterType;
          typedef typename SolverVectorType_::DataType DataType;

          std::shared_ptr<Solver::Linesearch<OperatorType, FilterType> > result;

          auto section = base->query_section(solver_name);

          auto solver_p = section->query("type");
          if (!solver_p.second)
            throw InternalError(__func__, __FILE__, __LINE__, "no type key found in property map: " + solver_name + "!");

          // \todo: NewtonRaphsonLinesearch requires the operator to compute hessians, which has to be caught at
          // runtime
          String solver_type = solver_p.first;
          /*if(solver_type == "NewtonRaphsonLinesearch")
            {
            result = Solver::new_newton_raphson_linesearch(
              derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr),
              derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr));
              }
              else */
          if(solver_type == "SecantLinesearch")
          {
            DataType initial_step(Solver::SecantLinesearch<OperatorType, FilterType>::initial_step_default);
            // Check if the initial step was specified
            auto initial_step_p = section->query("initial_step");
            if(initial_step_p.second)
              initial_step = DataType(std::stod(initial_step_p.first));

            bool keep_iterates(false);
            // Check if we have to keep the iterates
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second && std::stoul(keep_iterates_p.first) == 1)
              keep_iterates = true;

            result = Solver::new_secant_linesearch(
              derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr),
              derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr),
              initial_step,
              keep_iterates);
          }
          else if(solver_type == "StrongWolfeLinesearch")
          {
            bool keep_iterates(false);
            // Check if we have to keep the iterates
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second && std::stoul(keep_iterates_p.first) == 1)
              keep_iterates = true;

            // Set default function value decrease tolerance for the strong Wolfe condition
            DataType tol_decrease(Solver::StrongWolfeLinesearch<OperatorType, FilterType>::tol_decrease_default);
            // Check if something was specified
            auto tol_decrease_p = section->query("tol_decrease");
            if(tol_decrease_p.second)
              tol_decrease = DataType(std::stod(tol_decrease_p.first));

            // Set default curvature tolerance for the strong Wolfe condition
            DataType tol_curvature(Solver::StrongWolfeLinesearch<OperatorType, FilterType>::tol_curvature_default);
            // Check if something was specified
            auto tol_curvature_p = section->query("tol_curvature");
            if(tol_curvature_p.second)
              tol_curvature = DataType(std::stod(tol_curvature_p.first));

            result = Solver::new_strong_wolfe_linesearch(
              derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr),
              derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr),
              keep_iterates, tol_decrease, tol_curvature);
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "Unknown linesearch type " + solver_name + "!");

          // Some compilers (i.e. clang 3.7.1) cannot deduce the appropriate template parameter for the base class
          // so we need to pass result through a base class pointer
          std::shared_ptr<Solver::IterativeSolver<SolverVectorType_>> tmp = result;
          configure_iterative_solver(section, tmp);

          return result;

        }

        template
        <
          typename SystemLevelType_,
          typename TransferLevelType_,
          typename SolverVectorType_ = typename SystemLevelType_::GlobalSystemVectorR
        >
        static std::shared_ptr<Solver::IterativeSolver<SolverVectorType_>>
        create_nonlinear_optimiser( std::deque<SystemLevelType_*> & system_levels,
        std::deque<TransferLevelType_*>& DOXY(transfer_levels),
        PropertyMap* base, String solver_name,
        std::shared_ptr
        <
          Solver::NLOptPrecond<SolverVectorType_, typename SystemLevelType_::GlobalSystemFilter>
        >
        precon = nullptr)
        {
          typedef typename SystemLevelType_::GlobalQualityFunctional OperatorType;
          typedef typename SystemLevelType_::GlobalSystemFilter FilterType;

          std::shared_ptr<Solver::IterativeSolver<SolverVectorType_> > result;

          auto section = base->query_section(solver_name);

          auto solver_p = section->query("type");
          if (!solver_p.second)
            throw InternalError(__func__, __FILE__, __LINE__, "no type key found in property map: " + solver_name + "!");
          String solver_type = solver_p.first;

          if (solver_type == "ALGLIBMinLBFGS")
          {
#ifndef FEAST_HAVE_ALGLIB
            throw InternalError(__func__,__FILE__,__LINE__,
            "ALGLIBMinLBFGS is only available if FEAST was built with the alglib token in the buildid.");
#else

            if( Util::Comm::size() != 1)
              throw InternalError(__func__, __FILE__, __LINE__, "ALGLIBMinLBFGS is only available with 1 process!");

            // Default LBFGS update dimension is 0 so the constructor can choose is automatically.
            alglib::ae_int_t lbfgs_dim(0);
            auto lbfgs_dim_p = section->query("lbfgs_dim");
            // Get LBFGS update dimension.
            if(lbfgs_dim_p.second)
              lbfgs_dim = alglib::ae_int_t(std::stoul(lbfgs_dim_p.first));

            // Check if we have to keep the iterates
            bool keep_iterates(false);
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second && std::stoul(keep_iterates_p.first) == 1)
              keep_iterates = true;

            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_alglib_minlbfgs(
              derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr),
              derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr),
              lbfgs_dim,
              keep_iterates);
            configure_iterative_solver(section, solver);
            result = solver;
#endif // FEAST_HAVE_ALGLIB
          }
          else if (solver_type == "ALGLIBMinCG")
          {
#ifndef FEAST_HAVE_ALGLIB
            throw InternalError(__func__,__FILE__,__LINE__,
            "ALGLIBMinCG is only available if FEAST was built with the alglib token in the buildid.");
#else
            if( Util::Comm::size() != 1)
              throw InternalError(__func__, __FILE__, __LINE__, "ALGLIBMinCG is only available with 1 process!");

            // Get default direction update
            Solver::NLCGDirectionUpdate my_update(
              Solver::ALGLIBMinCG<OperatorType, FilterType>::direction_update_default);

            // Get direction update
            auto direction_update_p = section->query("direction_update");
            if(direction_update_p.second)
              my_update << direction_update_p.first;

            // Check if we have to keep the iterates
            bool keep_iterates(false);
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second && std::stoul(keep_iterates_p.first) == 1)
              keep_iterates = true;

            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_alglib_mincg(
              derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr),
              derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr),
              my_update,
              keep_iterates);
            configure_iterative_solver(section, solver);
            result = solver;
#endif // FEAST_HAVE_ALGLIB
          }
          else if (solver_type == "NLCG")
          {
            // Get default direction update
            Solver::NLCGDirectionUpdate my_update(Solver::NLCG<OperatorType, FilterType>::direction_update_default);
            // Get direction update
            auto direction_update_p = section->query("direction_update");
            if(direction_update_p.second)
              my_update << direction_update_p.first;

            // Check if we have to keep the iterates
            bool keep_iterates(false);
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second && std::stoul(keep_iterates_p.first) == 1)
              keep_iterates = true;

            std::shared_ptr<Solver::Linesearch
            <typename SystemLevelType_::GlobalQualityFunctional,
            typename SystemLevelType_::GlobalSystemFilter>> my_linesearch;

            // Default linesearch is StrongWolfeLinesearch
            String linesearch_name("StrongWolfeLinesearch");
            auto linesearch_p = section->query("linesearch");
            if(linesearch_p.second)
              linesearch_name = linesearch_p.first;

            my_linesearch = create_linesearch(system_levels, /* transfer_levels, */ base, linesearch_name);

            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_nlcg(
              derefer<SolverVectorType_>(system_levels.back()->op_sys, nullptr),
              derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr),
              my_linesearch,
              my_update,
              keep_iterates,
              precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else
            throw InternalError(__func__,__FILE__,__LINE__,"Solver type key "+stringify(solver_type)+" unknown.");
          return result;
        } // create_nonlinear_optimiser

    }; // MeshoptSolverFactory
  } // namespace Control
} // namespace FEAST

#endif // CONTROL_MESHOPT_MESHOPT_SOLVER_FACTORY_HPP
