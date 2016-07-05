#pragma once
#ifndef CONTROL_SCALAR_SOLVER_HPP
#define CONTROL_SCALAR_SOLVER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/basic_vcycle.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>

namespace FEAT
{
  namespace Control
  {
    struct SolverFactory
    {
      private:

        template <typename VectorType_>
        static void configure_iterative_solver(PropertyMap * section, std::shared_ptr<Solver::PreconditionedIterativeSolver<VectorType_> > solver)
        {
          Index rank = Util::Comm::rank();

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
              solver->set_plot(rank == 0);
            }
            else
            {
              throw InternalError(__func__, __FILE__, __LINE__, "plot value " + stringify(plot) + " unknown!");
            }
          }

          auto tol_rel_p = section->get_entry("tol_rel");
          if (tol_rel_p.second)
            solver->set_tol_rel(std::stod(tol_rel_p.first));

          auto max_iter_p = section->get_entry("max_iter");
          if (max_iter_p.second)
            solver->set_max_iter(std::stoul(max_iter_p.first));

          auto min_iter_p = section->get_entry("min_iter");
          if (min_iter_p.second)
            solver->set_min_iter(std::stoul(min_iter_p.first));
        }

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
        static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector> >
        create_schwarz_precon(std::deque<SystemLevelType_*> & system_levels, std::deque<TransferLevelType_*> & transfer_levels, PropertyMap * base, String solver_name, PropertyMap * section, typename Evaluator_::GateType *)
        {
          typedef typename SystemLevelType_::GlobalSystemVector SolverVectorType;
          std::shared_ptr<Solver::SolverBase<typename SolverVectorType::LocalVectorType> > precon_schwarz;
          auto schwarz_p = section->query("solver");
          if (schwarz_p.second)
          {
            precon_schwarz = create_scalar_solver<SystemLevelType_, TransferLevelType_, typename SolverVectorType::LocalVectorType>(system_levels, transfer_levels, base, get_section_path(base, section, solver_name, schwarz_p.first));
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "Schwarz precon section without solver key is not allowed!");

          return Solver::new_schwarz_precond(precon_schwarz, system_levels.back()->filter_sys);
        }

        template <typename Evaluator_, typename SystemLevelType_, typename TransferLevelType_>
        static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector::LocalVectorType> >
        create_schwarz_precon(std::deque<SystemLevelType_*> & /*system_levels*/, std::deque<TransferLevelType_*> & /*transfer_levels*/, PropertyMap * /*base*/, String /*solver_name*/,
        PropertyMap * /*section*/, ...)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Schwarz precon section is only allowed in global context! Maybe you have two in one solver branch?");
          return nullptr;
        }

        template <typename Evaluator_, typename SystemLevelType_>
        static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector> >
        create_ilu_precon(std::deque<SystemLevelType_*> & /*system_levels*/, typename Evaluator_::GateType *)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "ilu precon section is only allowed in local context!");
          return nullptr;
        }

        template <typename Evaluator_, typename SystemLevelType_>
        static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector::LocalVectorType> >
        create_ilu_precon(std::deque<SystemLevelType_*> & system_levels, ...)
        {
          auto result = Solver::new_ilu_precond(*system_levels.back()->matrix_sys, *system_levels.back()->filter_sys, 0ul);
          return result;
        }

        template <typename Evaluator_, typename SystemLevelType_>
        static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector> >
        create_ssor_precon(std::deque<SystemLevelType_*> & /*system_levels*/, PropertyMap * /*section*/, typename Evaluator_::GateType *)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "ssor precon section is only allowed in local context!");
          return nullptr;
        }

        template <typename Evaluator_, typename SystemLevelType_>
        static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector::LocalVectorType> >
        create_ssor_precon(std::deque<SystemLevelType_*> & system_levels, PropertyMap * section, ...)
        {
          auto omega_p = section->get_entry("omega");
          if (omega_p.second)
          {
            auto result = Solver::new_ssor_precond(*system_levels.back()->matrix_sys, *system_levels.back()->filter_sys, std::stod(omega_p.first));
            return result;
          }
          else
          {
            auto result = Solver::new_ssor_precond(*system_levels.back()->matrix_sys, *system_levels.back()->filter_sys);
            return result;
          }
        }

      public:

        /**
         * \brief Create solver tree based on PropertyMap
         *
         * \param[in] system_levels SystemLevel hierarchy
         * \param[in] transfer_levels TransferLevel hierarchy
         * \param[in] base A pointer to the PropertyMap that contains all solver related informations
         * \param[in] solver_name The name of the solver tree's root section
         */
        template <typename SystemLevelType_, typename TransferLevelType_, typename SolverVectorType_ = typename SystemLevelType_::GlobalSystemVector>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_scalar_solver(std::deque<SystemLevelType_*> & system_levels, std::deque<TransferLevelType_*> & transfer_levels, PropertyMap * base, String solver_name)
        {
          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > result;

          auto section = base->query_section(solver_name);

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
              precon = create_scalar_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(system_levels, transfer_levels, base, get_section_path(base, section, solver_name, precon_p.first));
          }

          if (solver_type == "pcg")
          {
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_pcg(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "bicgstab")
          {
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_bicgstab(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), precon);
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
            solver = Solver::new_fgmres(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), krylov_dim, 0.0, precon);
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
            solver = Solver::new_richardson(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), omega, precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "jac")
          {
            auto omega_p = section->get_entry("omega");
            if (omega_p.second)
            {
              result = Solver::new_jacobi_precond(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), std::stod(omega_p.first));
            }
            else
            {
              result = Solver::new_jacobi_precond(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr));
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
          else if (solver_type == "ilu")
          {
            result = create_ilu_precon<SolverVectorType_>(system_levels, nullptr);
          }
          else if (solver_type == "ssor")
          {
            result = create_ssor_precon<SolverVectorType_>(system_levels, section, nullptr);
          }
          else if (solver_type == "mgv")
          {
            typename std::remove_reference<decltype(derefer<SolverVectorType_>(system_levels.front()->matrix_sys, nullptr))>::type dummy_sys;
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
              coarse_solver = create_scalar_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(coarse_system_level, coarse_transfer_level, base, get_section_path(base, section, solver_name, coarse_solver_p.first));
            }
            else
              throw InternalError(__func__, __FILE__, __LINE__, "mgv precon section without coarse key is not allowed!");

            mgv->set_coarse_level(derefer<SolverVectorType_>(system_levels.front()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.front()->filter_sys, nullptr), coarse_solver);

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
                smoother = create_scalar_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(smoother_system_levels, smoother_transfer_levels, base, get_section_path(base, section, solver_name, smoother_p.first));
              }
              else
                throw InternalError(__func__, __FILE__, __LINE__, "mgv precon section without smoother key is not allowed!");

              mgv->push_level(derefer<SolverVectorType_>((*it)->matrix_sys, nullptr), derefer<SolverVectorType_>((*it)->filter_sys, nullptr), derefer<SolverVectorType_>((*jt)->prol_sys, nullptr),
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
    }; // SolverFactory

    /**
     * \brief High level linear system solver wrapper
     *
     * This control class holds the complete linear system hierarchy and optionally a corresponding
     * linear solver.
     *
     * The linear system mimics the provided MatrixTypeSolve_ type with regard to sparse matrix format,
     * memory-, data- and index type.
     *
     * The optionally created linear solver will be based on the provided ParamSection.
     * If no ParamSection is provided, the solver will read as follows:
     *
     * PCG ( VCycle ( S: Richardson ( Jacobi )  / C: Richardson ( Jacobi )  )  )
     *
     * with 4 pre- and postmoothing steps each and 4 jacobi calls as the coarse level 'solver'.
     */
    template <typename MatrixTypeSolve_, typename SystemLevelType_, typename TransferLevelType_>
    class ScalarSolver
    {
      public:
        /// The solvers lgs scalar matrix type
        using MatrixTypeSolve = MatrixTypeSolve_;
        /// The solvers lgs memory type
        using MemTypeSolve = typename MatrixTypeSolve_::MemType;
        /// The solvers lgs data type
        using DataTypeSolve = typename MatrixTypeSolve_::DataType;
        /// The solvers lgs index type
        using IndexTypeSolve = typename MatrixTypeSolve_::IndexType;
        /// The solvers lgs system level type
        using SystemLevelTypeSolve =  typename SystemLevelType_::template BaseType<MemTypeSolve, DataTypeSolve, IndexTypeSolve, MatrixTypeSolve_>;
        /// The solvers lgs transfer level type
        using TransferLevelTypeSolve = typename TransferLevelType_::template BaseType<SystemLevelTypeSolve>;
        /// The solvers lgs global system vector type
        using GlobalSystemVectorSolve = typename SystemLevelTypeSolve::GlobalSystemVector;
        /// The solvers lgs global system matrix type
        using GlobalSystemMatrixSolve = typename SystemLevelTypeSolve::GlobalSystemMatrix;
        /// The solvers lgs global system filter type
        using GlobalSystemFilterSolve = typename SystemLevelTypeSolve::GlobalSystemFilter;


      private:
        std::shared_ptr<Solver::SolverBase<GlobalSystemVectorSolve> > _solver;
        std::deque<SystemLevelTypeSolve*> _system_levels_solve;
        std::deque<TransferLevelTypeSolve*> _transfer_levels_solve;

        GlobalSystemVectorSolve _vec_sol_solve;
        GlobalSystemVectorSolve _vec_rhs_solve;

      public:

        /**
         * \brief Create ScalarSolver controll object
         *
         * This creates a new set of system levels, transfer levels and vectors from the containers provided,
         * that fit the MatrixTypeSolve_ type.
         *
         * Optionally a solver can be created to solve this system.
         *
         * \param[in] system_levels The assembled system_levels
         * \param[in] transfer_levels The assembled transfer_levels
         * \param[in] vec_sol The initialised solution vector
         * \param[in] vec_rhs The assembled right hand side
         */
        explicit ScalarSolver(std::deque<SystemLevelType_*> & system_levels, std::deque<TransferLevelType_*> & transfer_levels,
          typename SystemLevelType_::GlobalSystemVector & vec_sol, typename SystemLevelType_::GlobalSystemVector & vec_rhs) :
          _solver(nullptr)
        {
          for (Index i(0); i < system_levels.size(); ++i)
          {
            //system levels must be converted first, because transfer levels use their converted gates
            _system_levels_solve.push_back(new SystemLevelTypeSolve());
            _system_levels_solve.back()->convert(*system_levels.at(i));
            if (i > 0)
            {
              _transfer_levels_solve.push_back(new TransferLevelTypeSolve());
              _transfer_levels_solve.back()->convert(*_system_levels_solve.at(i-1), *_system_levels_solve.at(i), *transfer_levels.at(i-1));
            }
          }

          _vec_sol_solve.convert(&_system_levels_solve.back()->gate_sys, vec_sol);
          _vec_rhs_solve.convert(&_system_levels_solve.back()->gate_sys, vec_rhs);
        }

        /// DTOR
        ~ScalarSolver()
        {
          while (!_transfer_levels_solve.empty())
          {
            delete _transfer_levels_solve.back();
            _transfer_levels_solve.pop_back();
          }
          while (!_system_levels_solve.empty())
          {
            delete _system_levels_solve.back();
            _system_levels_solve.pop_back();
          }
        }

        /*/// Create vector in solve memory from given source vector
        GlobalSystemVectorSolve create_solve_vector(const typename SystemLevelType_::GlobalSystemVector & other)
        {
          GlobalSystemVectorSolve t;
          t.convert(&_system_levels_solve.back()->gate_sys, other);
          return t;
        }*/

        /// Returns an already created solver, or nullptr
        std::shared_ptr<Solver::SolverBase<GlobalSystemVectorSolve> > operator*()
        {
          return _solver;
        }

        /**
         * \brief Retrieve our solver pointer
         *
         * \note The solver will be created the first time this method is called.
         * \note The user is responsible to call init and done methods of the solver object where necessary.
         */
        std::shared_ptr<Solver::IterativeSolver<GlobalSystemVectorSolve> > create_default_solver()
        {
          if (_solver == nullptr)
          {
              //PCG ( VCycle ( S: Richardson ( Jacobi )  / C: Richardson ( Jacobi )  )  )
              auto mgv = std::make_shared<
                Solver::BasicVCycle<
                GlobalSystemMatrixSolve,
                GlobalSystemFilterSolve,
                typename TransferLevelTypeSolve::GlobalSystemTransferMatrix
              > >();

              auto coarse_precond = Solver::new_jacobi_precond(_system_levels_solve.front()->matrix_sys, _system_levels_solve.front()->filter_sys, 0.7);
              auto coarse_solver = Solver::new_richardson(_system_levels_solve.front()->matrix_sys, _system_levels_solve.front()->filter_sys, 1.0, coarse_precond);
              coarse_solver->set_min_iter(4);
              coarse_solver->set_max_iter(4);
              mgv->set_coarse_level(_system_levels_solve.front()->matrix_sys, _system_levels_solve.front()->filter_sys, coarse_solver);

              auto jt = _transfer_levels_solve.begin();
              for (auto it = ++_system_levels_solve.begin(); it != _system_levels_solve.end(); ++it, ++jt)
              {
                auto jac_smoother = Solver::new_jacobi_precond((*it)->matrix_sys, (*it)->filter_sys, 0.7);
                auto smoother = Solver::new_richardson((*it)->matrix_sys, (*it)->filter_sys, 1.0, jac_smoother);
                smoother->set_min_iter(4);
                smoother->set_max_iter(4);
                mgv->push_level((*it)->matrix_sys, (*it)->filter_sys, (*jt)->prol_sys, (*jt)->rest_sys, smoother, smoother);
              }
              auto solver = Solver::new_pcg(_system_levels_solve.back()->matrix_sys, _system_levels_solve.back()->filter_sys, mgv);
              _solver = solver;
              return solver;
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "cannot create solver object twice!");
        }

        /**
         * \brief Create solver tree based on PropertyMap
         *
         * \param[in] pmbase A pointer to the PropertyMap that contains all solver related informations
         * \param[in] pmbase_name The name of the solver tree's root section
         */
        std::shared_ptr<Solver::SolverBase<GlobalSystemVectorSolve> > create_solver(PropertyMap * pmbase, String pmbase_name)
        {
          _solver = Control::SolverFactory::create_scalar_solver(_system_levels_solve, _transfer_levels_solve, pmbase, pmbase_name);
          return _solver;
        }

        /// Retrieve system_levels_solve
        std::deque<SystemLevelTypeSolve*> & get_system_levels_solve()
        {
          return _system_levels_solve;
        }

        /// Retrieve transfer_levels_solve
        std::deque<TransferLevelTypeSolve*> & get_trasnfer_levels_solve()
        {
          return _transfer_levels_solve;
        }

        /// Retrieve vec_sol_solve
        GlobalSystemVectorSolve & get_vec_sol_solve()
        {
          return _vec_sol_solve;
        }

        /// Retrieve vec_rhs_solve
        GlobalSystemVectorSolve & get_vec_rhs_solve()
        {
          return _vec_rhs_solve;
        }

        /// execute the instances own solver
        Solver::Status solve()
        {
          if (_solver)
            return Solver::solve(*_solver, _vec_sol_solve, _vec_rhs_solve, _system_levels_solve.back()->matrix_sys, _system_levels_solve.back()->filter_sys);
          else
            throw InternalError(__func__, __FILE__, __LINE__, "calling solve with no solver allocated!");
        }

    }; // ScalarSolver
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_SCALAR_SOLVER_HPP
