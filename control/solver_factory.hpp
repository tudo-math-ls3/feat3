#pragma once
#ifndef CONTROL_SOLVER_FACTORY_HPP
#define CONTROL_SOLVER_FACTORY_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/basic_vcycle.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/pmr.hpp>
#include <kernel/solver/pcr.hpp>
#include <kernel/solver/psd.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/sor_precond.hpp>
#include <kernel/solver/spai_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/convert_precond.hpp>
#include <kernel/solver/matrix_stock.hpp>

#include <kernel/solver/alglib_wrapper.hpp>
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/solver/nlsd.hpp>
#include <kernel/solver/nloptls.hpp>
#include <kernel/solver/nlopt_precond.hpp>
#include <kernel/solver/qpenalty.hpp>

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
          using DataType = typename VectorType_::DataType;

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
            solver->set_tol_rel((DataType)std::stod(tol_rel_p.first));

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

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_schwarz_precon(MST_ & matrix_stock, PropertyMap * base, String solver_name, PropertyMap * section, size_t back_level, typename SolverVectorType_::GateType *)
        {
          using SolverVectorType = SolverVectorType_;
          std::shared_ptr<Solver::SolverBase<typename SolverVectorType::LocalVectorType> > precon_schwarz;
          auto schwarz_p = section->query("solver");
          if (schwarz_p.second)
          {
            precon_schwarz = create_scalar_solver_by_section<MST_, typename SolverVectorType::LocalVectorType>(matrix_stock, base, get_section_path(base, section, solver_name, schwarz_p.first),
                back_level);
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "Schwarz precon section without solver key is not allowed!");

          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          return Solver::new_schwarz_precond(precon_schwarz, filters.at(back_level));
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_schwarz_precon(MST_ &, PropertyMap *, String, PropertyMap *, size_t, ...)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Schwarz precon section is only allowed in global context! Maybe you have two in one solver branch?");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_ilu_precon(MST_ & , size_t, typename SolverVectorType_::GateType *)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "ilu precon section is only allowed in local context!");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_ilu_precon(MST_ & matrix_stock, size_t back_level, ...)
        {
          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto result = Solver::new_ilu_precond(systems.at(back_level), filters.at(back_level), 0ul);
          return result;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_spai_precon(MST_ & , size_t, typename SolverVectorType_::GateType *)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "spai precon section is only allowed in local context!");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_spai_precon(MST_ & matrix_stock, size_t back_level, ...)
        {
          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto result = Solver::new_spai_precond(systems.at(back_level), filters.at(back_level));
          return result;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_sor_precon(MST_ & , PropertyMap * /*section*/, size_t, typename SolverVectorType_::GateType *)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "sor precon section is only allowed in local context!");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_sor_precon(MST_ & matrix_stock, PropertyMap * section, size_t back_level, ...)
        {
          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);

          auto omega_p = section->get_entry("omega");
          if (omega_p.second)
          {
            auto result = Solver::new_sor_precond(systems.at(back_level), filters.at(back_level), (typename SolverVectorType_::DataType)std::stod(omega_p.first));
            return result;
          }
          else
          {
            auto result = Solver::new_sor_precond(systems.at(back_level), filters.at(back_level));
            return result;
          }
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_ssor_precon(MST_ & , PropertyMap * /*section*/, size_t, typename SolverVectorType_::GateType *)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "ssor precon section is only allowed in local context!");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_ssor_precon(MST_ & matrix_stock, PropertyMap * section, size_t back_level, ...)
        {
          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);

          auto omega_p = section->get_entry("omega");
          if (omega_p.second)
          {
            auto result = Solver::new_ssor_precond(systems.at(back_level), filters.at(back_level), (typename SolverVectorType_::DataType)std::stod(omega_p.first));
            return result;
          }
          else
          {
            auto result = Solver::new_ssor_precond(systems.at(back_level), filters.at(back_level));
            return result;
          }
        }

        template <typename MST_, typename SolverVectorType_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> > create_scalar_solver_by_section(MST_ & matrix_stock, PropertyMap * base,
            String precon_section_path, size_t back_level)
        {
          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon = nullptr;
          auto precon_section = base->query_section(precon_section_path);
          auto precon_memory = precon_section->query("memory", "main");
          auto precon_datatype = precon_section->query("datatype", "double");
          auto precon_indextype = precon_section->query("indextype", "unsigned long");

          if (precon_memory == SolverVectorType_::MemType::name() &&
              precon_datatype == Type::Traits<typename SolverVectorType_::DataType>::name() &&
              precon_indextype == Type::Traits<typename SolverVectorType_::IndexType>::name())
          {
            precon = create_scalar_solver<MST_, SolverVectorType_>(matrix_stock, base, precon_section_path, back_level);
          }
          else if (precon_memory == "main" && precon_datatype == "float" && precon_indextype == "unsigned long")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::Main, float, unsigned long>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, back_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
          else if (precon_memory == "main" && precon_datatype == "double" && precon_indextype == "unsigned long")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::Main, double, unsigned long>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, back_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#ifdef FEAT_HAVE_CUDA
          else if (precon_memory == "cuda" && precon_datatype == "float" && precon_indextype == "unsigned long")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::CUDA, float, unsigned long>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, back_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
          else if (precon_memory == "cuda" && precon_datatype == "double" && precon_indextype == "unsigned long")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::CUDA, double, unsigned long>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, back_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#endif
          else if (precon_memory == "main" && precon_datatype == "float" && precon_indextype == "unsigned int")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::Main, float, unsigned int>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, back_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
          else if (precon_memory == "main" && precon_datatype == "double" && precon_indextype == "unsigned int")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::Main, double, unsigned int>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, back_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#ifdef FEAT_HAVE_CUDA
          else if (precon_memory == "cuda" && precon_datatype == "float" && precon_indextype == "unsigned int")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::CUDA, float, unsigned int>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, back_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
          else if (precon_memory == "cuda" && precon_datatype == "double" && precon_indextype == "unsigned int")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::CUDA, double, unsigned int>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, back_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#endif
          else
          {
            throw InternalError(__func__, __FILE__, __LINE__, "memory/datatype/indextype combination unknown!");
          }

          return precon;
        }

      public:

        /**
         * \brief Create solver tree based on PropertyMap
         *
         * \param[in] matrix_stock A MatrixStock object, initialised with Systemlevels etc
         * \param[in] base A pointer to the PropertyMap that contains all solver related informations
         * \param[in] solver_name The name of the solver tree's root section
         */
        template <typename MST_, typename SolverVectorType_ = typename MST_::VectorType>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_scalar_solver(MST_ & matrix_stock, PropertyMap * base, String solver_name, size_t back_level = std::numeric_limits<size_t>::max())
        {
          if (back_level == std::numeric_limits<size_t>::max())
            back_level = matrix_stock.systems.size() - 1;

          using MemType = typename SolverVectorType_::MemType;
          using DataType = typename SolverVectorType_::DataType;
          using IndexType = typename SolverVectorType_::IndexType;

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > result;

          auto section = base->query_section(solver_name);
          if (section == nullptr)
            throw InternalError(__func__, __FILE__, __LINE__, "section not found in property map: " + solver_name + "!");


          auto solver_p = section->query("type");
          if (!solver_p.second)
            throw InternalError(__func__, __FILE__, __LINE__, "no type key found in property map section: " + solver_name + "!");
          String solver_type = solver_p.first;

          auto section_memory = section->query("memory", "main");
          auto section_datatype = section->query("datatype", "double");
          auto section_indextype = section->query("indextype", "unsigned long");
          XASSERT(section_memory == MemType::name());
          XASSERT(section_datatype == Type::Traits<DataType>::name());
          XASSERT(section_indextype == Type::Traits<IndexType>::name());

          matrix_stock.template compile_systems<SolverVectorType_>(nullptr);

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon = nullptr;

          auto precon_p = section->query("precon");
          if (precon_p.second)
          {
            if (precon_p.first == "none")
              precon = nullptr;
            else
            {
              auto precon_section_path = get_section_path(base, section, solver_name, precon_p.first);
              precon = create_scalar_solver_by_section<MST_, SolverVectorType_>(matrix_stock, base, precon_section_path, back_level);
            }
          }

          if (solver_type == "pcg")
          {
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            solver = Solver::new_pcg(systems.at(back_level), filters.at(back_level), precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "bicgstab")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_bicgstab(systems.at(back_level), filters.at(back_level), precon);
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

            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_fgmres(systems.at(back_level), filters.at(back_level), krylov_dim, 0.0, precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "richardson")
          {
            auto omega_p = section->get_entry("omega");
            DataType omega;
            if (omega_p.second)
              omega = (DataType)stod(omega_p.first);
            else
              omega = 1.0;

            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            solver = Solver::new_richardson(systems.at(back_level), filters.at(back_level), omega, precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "pmr")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_pmr(systems.at(back_level), filters.at(back_level), precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "pcr")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_pcr(systems.at(back_level), filters.at(back_level), precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "psd")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
            solver = Solver::new_psd(systems.at(back_level), filters.at(back_level), precon);
            configure_iterative_solver(section, solver);
            result = solver;
          }
          else if (solver_type == "jac")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto omega_p = section->get_entry("omega");
            if (omega_p.second)
            {
              result = Solver::new_jacobi_precond(systems.at(back_level), filters.at(back_level),
                  (DataType)std::stod(omega_p.first));
            }
            else
            {
              result = Solver::new_jacobi_precond(systems.at(back_level), filters.at(back_level));
            }
          }
          else if (solver_type == "scale")
          {
            auto omega_p = section->get_entry("omega");
            if (omega_p.second)
            {
              auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
              result = Solver::new_scale_precond(filters.at(back_level), (DataType)std::stod(omega_p.first));
            }
            else
              throw InternalError(__func__, __FILE__, __LINE__, "no omega key found in scale section!");
          }
          else if (solver_type == "ilu")
          {
            result = create_ilu_precon<SolverVectorType_>(matrix_stock, back_level, nullptr);
          }
          else if (solver_type == "spai")
          {
            result = create_spai_precon<SolverVectorType_>(matrix_stock, back_level, nullptr);
          }
          else if (solver_type == "sor")
          {
            result = create_sor_precon<SolverVectorType_>(matrix_stock, section, back_level, nullptr);
          }
          else if (solver_type == "ssor")
          {
            result = create_ssor_precon<SolverVectorType_>(matrix_stock, section, back_level, nullptr);
          }
          else if (solver_type == "mg" || solver_type == "scarcmg")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& transfers = matrix_stock.template get_transfers<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& hmap = matrix_stock.template get_hierarchy_map<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);

            typename std::remove_reference<decltype(hmap)>::type::mapped_type hierarchy; //shared pointer to our hierarchy

            auto hierarchy_p = section->get_entry("hierarchy");
            if (!hierarchy_p.second)
              throw InternalError(__func__, __FILE__, __LINE__, "no hierarchy key found in mg section!");

            if (hmap.count(hierarchy_p.first) > 0)
            {
              hierarchy = hmap.at(hierarchy_p.first);
            }
            else
            {
              hierarchy = std::make_shared<typename decltype(hierarchy)::element_type>();
              hmap[hierarchy_p.first] = hierarchy;

              auto hierarchy_section_path = get_section_path(base, section, solver_name, hierarchy_p.first);
              auto hierarchy_section = base->query_section(hierarchy_section_path);

              auto coarse_solver_p = hierarchy_section->query("coarse");
              if (!coarse_solver_p.second)
                throw InternalError(__func__, __FILE__, __LINE__, "mg section without coarse key is not allowed!");
              auto coarse_solver_section_path = get_section_path(base, section, solver_name, coarse_solver_p.first);

              auto smoother_p = hierarchy_section->query("smoother");
              if (!smoother_p.second)
                throw InternalError(__func__, __FILE__, __LINE__, "mg section without smoother key is not allowed!");
              auto smoother_section_path = get_section_path(base, section, solver_name, smoother_p.first);

              //for(Index level(0) ; level < matrix_stock.systems.size() ; ++level)
              for(Index level(matrix_stock.systems.size()) ; level > Index(0) ; )
              {
                --level;
                auto coarse_solver = create_scalar_solver_by_section<MST_, SolverVectorType_>(matrix_stock, base, coarse_solver_section_path, level);
                if (level == 0)
                {
                  hierarchy->push_level(systems.at(level), filters.at(level), coarse_solver);
                }
                else
                {
                  auto smoother = create_scalar_solver_by_section<MST_, SolverVectorType_>(matrix_stock, base, smoother_section_path, level);
                  hierarchy->push_level(systems.at(level), filters.at(level), transfers.at(level), smoother, smoother, smoother, coarse_solver);
                }
              }
            }

            auto cycle_p = section->query("cycle");
            if (!cycle_p.second)
              throw InternalError(__func__, __FILE__, __LINE__, "mg section without cycle key is not allowed!");
            auto lvl_min_s = section->query("lvl_min", "0");
            auto lvl_min = std::stoi(lvl_min_s);
            auto lvl_max_s = section->query("lvl_max", "-1");
            auto lvl_max = std::stoi(lvl_max_s);
            if (solver_type == "scarcmg" && lvl_max != -1)
                throw InternalError(__func__, __FILE__, __LINE__, "You really should not manually set the max lvl for a scarc mg solver!");
            lvl_max = Math::max(lvl_max, (int)back_level);

            if (solver_type == "mg")
            {
              if (cycle_p.first == "v")
              {
                auto mgv = Solver::new_multigrid(hierarchy, Solver::MultiGridCycle::V, lvl_max, lvl_min);
                result = mgv;
              }
              else if (cycle_p.first == "w")
              {
                auto mgv = Solver::new_multigrid(hierarchy, Solver::MultiGridCycle::W, lvl_max, lvl_min);
                result = mgv;
              }
              else if (cycle_p.first == "f")
              {
                auto mgv = Solver::new_multigrid(hierarchy, Solver::MultiGridCycle::F, lvl_max, lvl_min);
                result = mgv;
              }
              else
                throw InternalError(__func__, __FILE__, __LINE__, "mg cycle " + cycle_p.first + " unknown!");
            }
            else if (solver_type == "scarcmg")
            {
              if (cycle_p.first == "v")
              {
                auto mgv = Solver::new_scarcmultigrid(hierarchy, Solver::MultiGridCycle::V, lvl_max, lvl_min);
                result = mgv;
              }
              else if (cycle_p.first == "w")
              {
                auto mgv = Solver::new_scarcmultigrid(hierarchy, Solver::MultiGridCycle::W, lvl_max, lvl_min);
                result = mgv;
              }
              else if (cycle_p.first == "f")
              {
                auto mgv = Solver::new_scarcmultigrid(hierarchy, Solver::MultiGridCycle::F, lvl_max, lvl_min);
                result = mgv;
              }
              else
                throw InternalError(__func__, __FILE__, __LINE__, "mg cycle " + cycle_p.first + " unknown!");
            }
            else
              throw InternalError(__func__, __FILE__, __LINE__, "mg type " + solver_type + " unknown!");
          }
          else if (solver_type == "schwarz")
          {
            result = create_schwarz_precon<SolverVectorType_>(matrix_stock, base, solver_name, section, back_level, nullptr);
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "solver with type " + solver_type + " unkown!");

          return result;
        }

        /**
         * \brief Creates a new Solver::Linesearch object according to a configuration
         *
         * \param[in] base
         * The PropertyMap containing the configuration of this solver
         *
         * \param[in] solver_name
         * The name to identify our section in the PropertyMap
         *
         * \returns An std::shared_ptr to the new solver object.
         */
        template
        <
          typename Functional_,
          typename Filter_
        >
        static std::shared_ptr<Solver::Linesearch<Functional_, Filter_>> create_linesearch(
          Functional_& functional, Filter_& filter, PropertyMap* base, const String& solver_name)
        {

          typedef typename Functional_::VectorTypeR VectorTypeR;

          std::shared_ptr<Solver::Linesearch<Functional_, Filter_> > result;

          // Get the section where our solver is configured
          auto section = base->query_section(solver_name);
          if(section == nullptr)
          {
            throw InternalError(__func__, __FILE__, __LINE__,
            "could not find section "+solver_name+" in PropertyMap!");
          }

          // Get the required type String
          auto solver_p = section->query("type");
          if (!solver_p.second)
          {
            throw InternalError(__func__, __FILE__, __LINE__,
            "No type key found in PropertyMap section with name " + solver_name + "!");
          }

          String solver_type = solver_p.first;

          // \todo: NewtonRaphsonLinesearch requires the operator to compute hessians, which has to be caught at
          // runtime
          /*if(solver_type == "NewtonRaphsonLinesearch")
            {
            result = Solver::new_newton_raphson_linesearch(
              derefer<VectorTypeR>(system_levels.at(back_level)->op_sys, nullptr),
              derefer<VectorTypeR>(system_levels.at(back_level)->filter_sys, nullptr));
              }
              else */
          if(solver_type == "SecantLinesearch")
          {
            result = Solver::new_secant_linesearch(
              derefer<VectorTypeR>(functional, nullptr),
              derefer<VectorTypeR>(filter, nullptr));
          }
          else if(solver_type == "StrongWolfeLinesearch")
          {
            result = Solver::new_strong_wolfe_linesearch(
              derefer<VectorTypeR>(functional, nullptr),
              derefer<VectorTypeR>(filter, nullptr));
          }
          else
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Unknown linesearch type " + solver_name + "!");
          }

          // Now read the other configuration options from this section
          result->read_config(section);

          return result;
        }

        /**
         *
         */
        template
        <
          typename Functional_,
          typename Filter_
        >
        static std::shared_ptr<Solver::IterativeSolver<typename Functional_::VectorTypeR>>
        create_nonlinear_optimiser(Functional_& functional, Filter_& filter, PropertyMap* base, String solver_name,
        std::shared_ptr
        <
          Solver::NLOptPrecond<typename Functional_::VectorTypeR, Filter_>
        >
        precon = nullptr)
        {
          typedef typename Functional_::DataType DataType;
          typedef typename Functional_::VectorTypeR VectorTypeR;

          // At the end, we return this guy
          std::shared_ptr<Solver::IterativeSolver<VectorTypeR>> result;

          auto section = base->query_section(solver_name);

          auto solver_p = section->query("type");
          if (!solver_p.second)
          {
            throw InternalError(__func__, __FILE__, __LINE__,
            "no type key found in property map: " + solver_name + "!");
          }

          String solver_type = solver_p.first;

          if(solver_type == "QPenalty")
          {
            // Create inner solver first
            auto inner_solver_p = section->query("inner_solver");
            // Safety catches for recursions
            XASSERTM(inner_solver_p.second, "QPenalty solver section is missing mandatory inner_solver key.");
            XASSERTM(inner_solver_p.first != "QPenalty", "QPenalty cannot be the inner solver for QPenalty.");

            std::shared_ptr<Solver::IterativeSolver<VectorTypeR>> inner_solver;
            inner_solver = create_nonlinear_optimiser(functional, filter, base, inner_solver_p.first, precon);

            DataType initial_penalty_param(1);
            result = Solver::new_qpenalty(derefer<VectorTypeR>(functional, nullptr), inner_solver,
            initial_penalty_param);

          }
          else if(solver_type == "ALGLIBMinLBFGS")
          {
#ifndef FEAT_HAVE_ALGLIB
            throw InternalError(__func__,__FILE__,__LINE__,
            "ALGLIBMinLBFGS is only available if FEAT was built with the alglib token in the buildid.");
#else

            if( Util::Comm::size() != 1)
            {
              throw InternalError(__func__, __FILE__, __LINE__, "ALGLIBMinLBFGS is only available with 1 process!");
            }

            // Get some default parameters for the constructor
            alglib::ae_int_t lbfgs_dim(0);
            bool keep_iterates(false);

            auto solver = Solver::new_alglib_minlbfgs(
              derefer<VectorTypeR>(functional, nullptr),
              derefer<VectorTypeR>(filter, nullptr),
              lbfgs_dim,
              keep_iterates);

            result = solver;
#endif // FEAT_HAVE_ALGLIB
          }
          else if (solver_type == "ALGLIBMinCG")
          {
#ifndef FEAT_HAVE_ALGLIB
            throw InternalError(__func__,__FILE__,__LINE__,
            "ALGLIBMinCG is only available if FEAT was built with the alglib token in the buildid.");
#else
            if( Util::Comm::size() != 1)
            {
              throw InternalError(__func__, __FILE__, __LINE__, "ALGLIBMinCG is only available with 1 process!");
            }

            // Get default direction update
            Solver::NLCGDirectionUpdate my_update(
              Solver::ALGLIBMinCG<Functional_, Filter_>::direction_update_default);

            // By default, do not keep the iterates
            bool keep_iterates(false);

            auto solver = Solver::new_alglib_mincg(
              derefer<VectorTypeR>(functional, nullptr),
              derefer<VectorTypeR>(filter, nullptr),
              my_update,
              keep_iterates);

            result = solver;
#endif // FEAT_HAVE_ALGLIB
          }
          else if (solver_type == "NLCG")
          {
            // Get default direction update
            Solver::NLCGDirectionUpdate my_update(Solver::NLCG<Functional_, Filter_>::direction_update_default);
            // Set this to false for the constructor
            bool keep_iterates(false);

            std::shared_ptr<Solver::Linesearch<Functional_, Filter_>> my_linesearch;

            // Default linesearch is StrongWolfeLinesearch
            String linesearch_name("StrongWolfeLinesearch");
            auto linesearch_p = section->query("linesearch");
            if(linesearch_p.second)
            {
              linesearch_name = linesearch_p.first;
            }

            my_linesearch = create_linesearch(functional, filter, base, linesearch_name);

            auto solver = Solver::new_nlcg(
              derefer<VectorTypeR>(functional, nullptr),
              derefer<VectorTypeR>(filter, nullptr),
              my_linesearch,
              my_update,
              keep_iterates,
              precon);
            result = solver;
          }
          else
          {
            throw InternalError(__func__,__FILE__,__LINE__,"Solver type key "+stringify(solver_type)+" unknown.");
          }

          result->read_config(section);

          return result;
        } // create_nonlinear_optimiser

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

    }; // SolverFactory

  } // namespace Control
} // namespace FEAT

#endif // CONTROL_SOLVER_FACTORY_HPP
