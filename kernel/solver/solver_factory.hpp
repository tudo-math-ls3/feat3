#pragma once
#ifndef KERNEL_SOLVER_SOLVER_FACTORY_HPP
#define KERNEL_SOLVER_SOLVER_FACTORY_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/pmr.hpp>
#include <kernel/solver/pcr.hpp>
#include <kernel/solver/psd.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/bicgstabl.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/chebyshev.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/rgcr.hpp>
#include <kernel/solver/pipepcg.hpp>
#include <kernel/solver/gropppcg.hpp>
#include <kernel/solver/rbicgstab.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/sor_precond.hpp>
#include <kernel/solver/spai_precond.hpp>
#include <kernel/solver/polynomial_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/convert_precond.hpp>
#include <kernel/solver/matrix_stock.hpp>
#include <kernel/solver/idrs.hpp>

// Includes for nonlinear optimisers
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/fixed_step_linesearch.hpp>
#include <kernel/solver/newton_raphson_linesearch.hpp>
#include <kernel/solver/secant_linesearch.hpp>
#include <kernel/solver/mqc_linesearch.hpp>
#include <kernel/solver/alglib_wrapper.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/solver/nlsd.hpp>
#include <kernel/solver/nloptls.hpp>
#include <kernel/solver/nlopt_precond.hpp>
#include <kernel/solver/qpenalty.hpp>

#include <kernel/util/dist.hpp>

namespace FEAT
{
  namespace Solver
  {
    struct SolverFactory
    {
      private:
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
        create_schwarz_precon(MST_& matrix_stock, PropertyMap* base, String section_name, PropertyMap* section,
            size_t solver_level, typename SolverVectorType_::GateType*)
        {
          using SolverVectorType = SolverVectorType_;
          std::shared_ptr<Solver::SolverBase<typename SolverVectorType::LocalVectorType> > precon_schwarz;
          auto schwarz_p = section->query("solver");
          if (schwarz_p.second)
          {
            precon_schwarz = create_scalar_solver_by_section<MST_, typename SolverVectorType::LocalVectorType>(matrix_stock, base, get_section_path(base, section, section_name, schwarz_p.first),
                solver_level);
          }
          else
          {
            throw InternalError(__func__, __FILE__, __LINE__,
            "Schwarz precon section without solver key is not allowed!");
          }

          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          return Solver::new_schwarz_precond(section_name, section, precon_schwarz, filters.at(solver_level));
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_schwarz_precon(MST_ &, PropertyMap *, String, PropertyMap *, size_t, ...)
        {
          throw InternalError(__func__, __FILE__, __LINE__,
          "Schwarz precon section is only allowed in global context! Maybe you have two in one solver branch?");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_pipepcg(MST_& matrix_stock, PropertyMap* base, String section_name, PropertyMap* section,
            size_t solver_level, typename SolverVectorType_::GateType*)
        {

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon(nullptr);
          auto precon_p = section->query("precon");
          if (precon_p.second)
          {
            if (precon_p.first == "none")
            {
              precon = nullptr;
            }
            else
            {
              auto precon_section_path = get_section_path(base, section, section_name, precon_p.first);
              precon = create_scalar_solver_by_section<MST_, SolverVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            }
          }

          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          return Solver::new_pipepcg(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_pipepcg(MST_ &, PropertyMap *, String, PropertyMap *, size_t, ...)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "pipepcg solver section is only allowed in global context!");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_gropppcg(MST_& matrix_stock, PropertyMap* base, String section_name, PropertyMap* section,
            size_t solver_level, typename SolverVectorType_::GateType*)
        {

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon(nullptr);
          auto precon_p = section->query("precon");
          if (precon_p.second)
          {
            if (precon_p.first == "none")
            {
              precon = nullptr;
            }
            else
            {
              auto precon_section_path = get_section_path(base, section, section_name, precon_p.first);
              precon = create_scalar_solver_by_section<MST_, SolverVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            }
          }

          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          return Solver::new_gropppcg(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_gropppcg(MST_ &, PropertyMap *, String, PropertyMap *, size_t, ...)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "gropppcg solver section is only allowed in global context!");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_rbicgstab(MST_& matrix_stock, PropertyMap* base, String section_name, PropertyMap* section,
            size_t solver_level, typename SolverVectorType_::GateType*)
        {

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon(nullptr);
          auto precon_p = section->query("precon");
          if (precon_p.second)
          {
            if (precon_p.first == "none")
            {
              precon = nullptr;
            }
            else
            {
              auto precon_section_path = get_section_path(base, section, section_name, precon_p.first);
              precon = create_scalar_solver_by_section<MST_, SolverVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            }
          }

          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          return Solver::new_rbicgstab(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_rbicgstab(MST_ &, PropertyMap *, String, PropertyMap *, size_t, ...)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "rbicgstab solver section is only allowed in global context!");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_ilu_precon(MST_&, const String&, PropertyMap*, size_t, typename SolverVectorType_::GateType *)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "ilu precon section is only allowed in local context!");
          return nullptr;
        }

        template <typename SolverVectorType_, typename MST_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_ilu_precon(MST_ & matrix_stock, const String& section_name, PropertyMap* section, size_t solver_level, ...)
        {
          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto result = Solver::new_ilu_precond(section_name, section, systems.at(solver_level), filters.at(solver_level));
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
        create_spai_precon(MST_ & matrix_stock, size_t solver_level, ...)
        {
          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto result = Solver::new_spai_precond(systems.at(solver_level), filters.at(solver_level));
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
        create_sor_precon(MST_ & matrix_stock, const String& section_name, PropertyMap* section, size_t solver_level,
        ...)
        {
          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);

          auto result = Solver::new_sor_precond(
            section_name, section, systems.at(solver_level), filters.at(solver_level));
          return result;
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
        create_ssor_precon(MST_ & matrix_stock, const String& section_name, PropertyMap* section, size_t solver_level,
        ...)
        {
          auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
          auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);

          auto result = Solver::new_ssor_precond(
            section_name, section, systems.at(solver_level), filters.at(solver_level));

          return result;
        }

        template <typename MST_, typename SolverVectorType_>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> > create_scalar_solver_by_section(
          MST_& matrix_stock, PropertyMap* base, const String& precon_section_path, size_t solver_level)
        {
          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon = nullptr;
          auto precon_section = base->query_section(precon_section_path);
          auto precon_memorytype = precon_section->query("memorytype", "main");
          auto precon_datatype = precon_section->query("datatype", "double");
          String precon_indextype;
          if (precon_memorytype == "cuda")
              precon_indextype = precon_section->query("indextype", "unsigned int");
          else
              precon_indextype = precon_section->query("indextype", "unsigned long");

          if (precon_memorytype == SolverVectorType_::MemType::name() &&
              precon_datatype == Type::Traits<typename SolverVectorType_::DataType>::name() &&
              precon_indextype == Type::Traits<typename SolverVectorType_::IndexType>::name())
          {
            precon = create_scalar_solver<MST_, SolverVectorType_>(matrix_stock, base, precon_section_path, solver_level);
          }
#ifdef FEAT_SF_ESOTERIC
          else if (precon_memorytype == "main" && precon_datatype == "float" && precon_indextype == "unsigned long")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::Main, float, unsigned long>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#endif
          else if (precon_memorytype == "main" && precon_datatype == "double" && precon_indextype == "unsigned long")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::Main, double, unsigned long>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#ifdef FEAT_SF_ESOTERIC
#ifdef FEAT_HAVE_CUDA
          else if (precon_memorytype == "cuda" && precon_datatype == "float" && precon_indextype == "unsigned long")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::CUDA, float, unsigned long>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
          else if (precon_memorytype == "cuda" && precon_datatype == "double" && precon_indextype == "unsigned long")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::CUDA, double, unsigned long>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#endif
#endif
#ifdef FEAT_SF_ESOTERIC
          else if (precon_memorytype == "main" && precon_datatype == "float" && precon_indextype == "unsigned int")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::Main, float, unsigned int>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
          else if (precon_memorytype == "main" && precon_datatype == "double" && precon_indextype == "unsigned int")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::Main, double, unsigned int>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#endif
#ifdef FEAT_HAVE_CUDA
#ifdef FEAT_SF_ESOTERIC
          else if (precon_memorytype == "cuda" && precon_datatype == "float" && precon_indextype == "unsigned int")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::CUDA, float, unsigned int>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#endif
          else if (precon_memorytype == "cuda" && precon_datatype == "double" && precon_indextype == "unsigned int")
          {
            using NextVectorType_ = typename SolverVectorType_::template ContainerTypeByMDI<Mem::CUDA, double, unsigned int>;
            auto precon_next = create_scalar_solver<MST_, NextVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            precon = std::make_shared<Solver::ConvertPrecond<SolverVectorType_, NextVectorType_>>(precon_next);
          }
#endif
          else
          {
            throw InternalError(__func__, __FILE__, __LINE__, "memorytype/datatype/indextype combination unknown!\n Did you try configure with the --sf_esoteric flag?");
          }

          return precon;
        }

      public:

        /**
         * \brief Create solver tree based on PropertyMap
         *
         * \param[in] matrix_stock A MatrixStock object, initialised with Systemlevels etc
         * \param[in] base A pointer to the PropertyMap that contains all solver related informations
         * \param[in] section_name The name of the solver tree's root section
         */
        template <typename MST_, typename SolverVectorType_ = typename MST_::VectorType>
        static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
        create_scalar_solver(MST_ & matrix_stock, PropertyMap* base, const String& section_name, std::size_t solver_level = std::size_t(0))
        {
          using MemType = typename SolverVectorType_::MemType;
          using DataType = typename SolverVectorType_::DataType;
          using IndexType = typename SolverVectorType_::IndexType;

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > result;

          auto section = base->query_section(section_name);
          if (section == nullptr)
          {
            throw InternalError(__func__, __FILE__, __LINE__,
            "section not found in property map: " + section_name + "!");
          }

          auto solver_p = section->query("type");
          if (!solver_p.second)
            throw InternalError(__func__, __FILE__, __LINE__, "no type key found in property map section: " + section_name + "!");
          String solver_type = solver_p.first;

          auto section_memorytype = section->query("memorytype", "main");
          auto section_datatype = section->query("datatype", "double");
          String section_indextype;
          if (section_memorytype == "cuda")
            section_indextype = section->query("indextype", "unsigned int");
          else
            section_indextype = section->query("indextype", "unsigned long");

          XASSERT(section_memorytype == MemType::name());
          XASSERT(section_datatype == Type::Traits<DataType>::name());
          XASSERT(section_indextype == Type::Traits<IndexType>::name());

          matrix_stock.template compile_systems<SolverVectorType_>(nullptr);

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon = nullptr;

          auto precon_p = section->query("precon");
          if (precon_p.second)
          {
            if (precon_p.first == "none")
            {
              precon = nullptr;
            }
            else
            {
              auto precon_section_path = get_section_path(base, section, section_name, precon_p.first);
              precon = create_scalar_solver_by_section<MST_, SolverVectorType_>(matrix_stock, base, precon_section_path, solver_level);
            }
          }

          if (solver_type == "pcg")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_pcg(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "bicgstab")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_bicgstab(
              section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "bicgstabl")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_bicgstabl(
              section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "fgmres")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_fgmres(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "richardson")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_richardson(
              section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "chebyshev")
          {
            XASSERTM(precon == nullptr, "Chebyshev solver accepts no preconditioner!");
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_chebyshev(
              section_name, section, systems.at(solver_level), filters.at(solver_level));
          }
          else if (solver_type == "pmr")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_pmr(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "pcr")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_pcr(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "psd")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_psd(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "rgcr")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_rgcr(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "idrs")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_idrs(section_name, section, systems.at(solver_level), filters.at(solver_level), precon);
          }
          else if (solver_type == "pipepcg")
          {
            result = create_pipepcg<SolverVectorType_>(matrix_stock, base, section_name, section, solver_level, nullptr);
          }
          else if (solver_type == "gropppcg")
          {
            result = create_gropppcg<SolverVectorType_>(matrix_stock, base, section_name, section, solver_level, nullptr);
          }
          else if (solver_type == "rbicgstab")
          {
            result = create_rbicgstab<SolverVectorType_>(matrix_stock, base, section_name, section, solver_level, nullptr);
          }
          else if (solver_type == "jacobi")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_jacobi_precond(section_name, section, systems.at(solver_level), filters.at(solver_level));
          }
          else if (solver_type == "polynomial")
          {
            auto& systems = matrix_stock.template get_systems<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_polynomial_precond(section_name, section, systems.at(solver_level), filters.at(solver_level));
          }
          else if (solver_type == "scale")
          {
            auto& filters = matrix_stock.template get_filters<SolverVectorType_>(nullptr, nullptr, nullptr, nullptr);
            result = Solver::new_scale_precond(section_name, section, filters.at(solver_level));
          }
          else if (solver_type == "ilu")
          {
            result = create_ilu_precon<SolverVectorType_>(matrix_stock, section_name, section, solver_level, nullptr);
          }
          else if (solver_type == "spai")
          {
            result = create_spai_precon<SolverVectorType_>(matrix_stock, solver_level, nullptr);
          }
          else if (solver_type == "sor")
          {
            result = create_sor_precon<SolverVectorType_>(matrix_stock, section_name, section, solver_level, nullptr);
          }
          else if (solver_type == "ssor")
          {
            result = create_ssor_precon<SolverVectorType_>(matrix_stock, section_name, section, solver_level, nullptr);
          }
          else if (solver_type == "mg")
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
              hierarchy = std::make_shared<typename decltype(hierarchy)::element_type>(matrix_stock.size_virtual);
              hmap[hierarchy_p.first] = hierarchy;

              auto hierarchy_section_path = get_section_path(base, section, section_name, hierarchy_p.first);
              auto hierarchy_section = base->query_section(hierarchy_section_path);

              auto hierarchy_type_p = hierarchy_section->query("type");
              if (!hierarchy_type_p.second)
                throw InternalError(__func__, __FILE__, __LINE__, "hierarchy section without type key is not allowed!");
              XASSERTM(hierarchy_type_p.first == "hierarchy", "hierarchy key type must match string 'hierarchy'!");

              auto coarse_solver_p = hierarchy_section->query("coarse");
              if (!coarse_solver_p.second)
                throw InternalError(__func__, __FILE__, __LINE__, "hierarchy section without coarse key is not allowed!");
              auto coarse_solver_section_path = get_section_path(base, section, section_name, coarse_solver_p.first);

              auto smoother_p = hierarchy_section->query("smoother");
              if (!smoother_p.second)
                throw InternalError(__func__, __FILE__, __LINE__, "hierarchy section without smoother key is not allowed!");
              auto smoother_section_path = get_section_path(base, section, section_name, smoother_p.first);

              XASSERT(matrix_stock.systems.size() <= matrix_stock.size_virtual);
              for(Index level(0) ; level < matrix_stock.systems.size() ; ++level)
              {
                auto coarse_solver = create_scalar_solver_by_section<MST_, SolverVectorType_>(matrix_stock, base, coarse_solver_section_path, level);
                if ((level+1) == matrix_stock.systems.size())
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

            // These parameters have to be known before calling the constructor, so we have to parse them here
            auto cycle_p = section->query("cycle");
            if (!cycle_p.second)
              throw InternalError(__func__, __FILE__, __LINE__, "mg section without cycle key is not allowed!");
            auto lvl_min_s = section->query("lvl_min", "-1");
            auto lvl_min = std::stoi(lvl_min_s);
            auto lvl_max_s = section->query("lvl_max", "0");
            auto lvl_max = std::stoi(lvl_max_s);
            lvl_max = Math::max(lvl_max, int(solver_level));
            lvl_min = Math::min(lvl_min, int(matrix_stock.size_virtual)-1);

            Solver::MultiGridCycle cycle(Solver::MultiGridCycle::V);

            if(!cycle_p.first.parse(cycle))
            {
              throw InternalError(__func__, __FILE__, __LINE__, "mg cycle " + cycle_p.first + " unknown!");
            }

            auto mgv = Solver::new_multigrid(hierarchy, cycle, lvl_max, lvl_min);

            auto adapt_cgc_p = section->query("adapt_cgc");
            if(adapt_cgc_p.second)
            {
              if(adapt_cgc_p.first == "fixed")
                mgv->set_adapt_cgc(Solver::MultiGridAdaptCGC::Fixed);
              else if(adapt_cgc_p.first == "min_energy")
                mgv->set_adapt_cgc(Solver::MultiGridAdaptCGC::MinEnergy);
              else if(adapt_cgc_p.first == "min_defect")
                mgv->set_adapt_cgc(Solver::MultiGridAdaptCGC::MinDefect);
              else
                throw InternalError(__func__, __FILE__, __LINE__, "unknown coarse grid correction adaptivity mode: " + adapt_cgc_p.first);
            }
            result = mgv;
          }
          else if (solver_type == "schwarz")
          {
            result = create_schwarz_precon<SolverVectorType_>(matrix_stock, base, section_name, section, solver_level, nullptr);
          }
          else
          {
            throw InternalError(__func__, __FILE__, __LINE__, "solver with type " + solver_type + " unkown!");
          }

          return result;
        }

        /**
         * \brief Creates a new Solver::Linesearch object according to a configuration
         *
         * \param[in] base
         * The PropertyMap containing the configuration of this solver
         *
         * \param[in] functional
         * The nonlinear functional the line search is used for
         *
         * \param[in] filter
         * The filter for the nonlinear functional's gradient
         *
         * \param[in] section_name
         * The name to identify our section in the PropertyMap
         *
         * \returns An std::shared_ptr to the new solver object.
         */
        template<typename Functional_, typename Filter_>
        static std::shared_ptr<Solver::Linesearch<Functional_, Filter_>> create_linesearch(
          Functional_& functional, Filter_& filter, PropertyMap* base, const String& section_name)
        {

          // This is the type for the state vector used in the functional
          typedef typename Functional_::VectorTypeR VectorTypeR;

          // In the end, we return this guy
          std::shared_ptr<Solver::Linesearch<Functional_, Filter_> > result;

          // Get the section where our solver is configured
          auto section = base->query_section(section_name);
          if(section == nullptr)
          {
            throw InternalError(__func__, __FILE__, __LINE__,
            "could not find section "+section_name+" in PropertyMap!");
          }

          // Get the required type String
          auto solver_p = section->query("type");
          if (!solver_p.second)
          {
            throw InternalError(__func__, __FILE__, __LINE__,
            "No type key found in PropertyMap section with name " + section_name + "!");
          }
          String solver_type(solver_p.first);

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
              section_name, section,derefer<VectorTypeR>(functional, nullptr), derefer<VectorTypeR>(filter, nullptr));
          }
          else if(solver_type == "MQCLinesearch")
          {
            result = Solver::new_mqc_linesearch(
              section_name, section,derefer<VectorTypeR>(functional, nullptr), derefer<VectorTypeR>(filter, nullptr));
          }
          else
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Unknown linesearch type " + section_name + "!");
          }

          return result;
        }

        /**
         * \brief Creates a new nonlinear optimiser according to a configuration
         *
         * \param[in] base
         * The PropertyMap containing the configuration of this solver
         *
         * \param[in] functional
         * The nonlinear functional the line search is used for
         *
         * \param[in] filter
         * The filter for the nonlinear functional's gradient
         *
         * \param[in] section_name
         * The name to identify our section in the PropertyMap
         *
         * \param[in] precon
         * The preconditioner
         *
         * \returns An IterativeSolver std::shared_ptr to the new solver object.
         */
        template<typename Functional_, typename Filter_>
        static std::shared_ptr<Solver::IterativeSolver<typename Functional_::VectorTypeR>>
        create_nonlinear_optimiser(Functional_& functional, Filter_& filter, PropertyMap* base, String section_name,
        std::shared_ptr<Solver::NLOptPrecond<typename Functional_::VectorTypeR, Filter_>> precon = nullptr)
        {
          typedef typename Functional_::VectorTypeR VectorTypeR;

          // At the end, we return this guy
          std::shared_ptr<Solver::IterativeSolver<VectorTypeR>> result;

          auto section = base->query_section(section_name);

          // Get the String identifying the solver type
          auto solver_p = section->query("type");
          if (!solver_p.second)
          {
            throw InternalError(__func__, __FILE__, __LINE__,
            "no type key found in property map: " + section_name + "!");
          }
          String solver_type(solver_p.first);

          if(solver_type == "QPenalty")
          {
            // Create inner solver first
            auto inner_solver_p = section->query("inner_solver");
            // Safety catches for recursions
            XASSERTM(inner_solver_p.second, "QPenalty solver section is missing mandatory inner_solver key.");
            XASSERTM(inner_solver_p.first != "QPenalty", "QPenalty cannot be the inner solver for QPenalty.");

            std::shared_ptr<Solver::IterativeSolver<VectorTypeR>> inner_solver;
            inner_solver = create_nonlinear_optimiser(functional, filter, base, inner_solver_p.first, precon);

            result = Solver::new_qpenalty<Functional_>(
              section_name, section, derefer<VectorTypeR>(functional, nullptr), inner_solver);

          }
          else if(solver_type == "ALGLIBMinLBFGS")
          {
#ifndef FEAT_HAVE_ALGLIB
            throw InternalError(__func__,__FILE__,__LINE__,
            "ALGLIBMinLBFGS is only available if FEAT was built with the alglib token in the buildid.");
#else

            Dist::Comm comm_world(Dist::Comm::world());
            if(comm_world.size() != 1)
            {
              throw InternalError(__func__, __FILE__, __LINE__, "ALGLIBMinLBFGS is only available with 1 process!");
            }

            result = Solver::new_alglib_minlbfgs<Functional_, Filter_>(
              section_name, section,
              derefer<VectorTypeR>(functional, nullptr), derefer<VectorTypeR>(filter, nullptr));

#endif // FEAT_HAVE_ALGLIB
          }
          else if (solver_type == "ALGLIBMinCG")
          {
#ifndef FEAT_HAVE_ALGLIB
            throw InternalError(__func__,__FILE__,__LINE__,
            "ALGLIBMinCG is only available if FEAT was built with the alglib token in the buildid.");
#else
            Dist::Comm comm_world(Dist::Comm::world());
            if(comm_world.size() != 1)
            {
              throw InternalError(__func__, __FILE__, __LINE__, "ALGLIBMinLBFGS is only available with 1 process!");
            }

            result = Solver::new_alglib_mincg<Functional_, Filter_>(
              section_name, section, derefer<VectorTypeR>(functional, nullptr), derefer<VectorTypeR>(filter, nullptr));

#endif // FEAT_HAVE_ALGLIB
          }
          else if (solver_type == "NLCG")
          {
            std::shared_ptr<Solver::Linesearch<Functional_, Filter_>> my_linesearch;

            String linesearch_name("");
            auto linesearch_p = section->query("linesearch");
            if(linesearch_p.second)
            {
              linesearch_name = linesearch_p.first;
            }
            else
            {
              throw InternalError(__func__,__FILE__,__LINE__,
              "NLCG config section "+section_name+" is missing the mandatory linesearch key!");
            }

            my_linesearch = create_linesearch(functional, filter, base, linesearch_name);

            result = Solver::new_nlcg(
              section_name, section, derefer<VectorTypeR>(functional, nullptr), derefer<VectorTypeR>(filter, nullptr),
              my_linesearch, precon);

          }
          else
          {
            throw InternalError(__func__,__FILE__,__LINE__,"Solver type key "+stringify(solver_type)+" unknown.");
          }

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

#ifdef FEAT_EICKT
    using MST1_ = Solver::MatrixStock<
      Global::Matrix<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Filter<LAFEM::UnitFilter<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Transfer<LAFEM::Transfer<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>>, LAFEM::VectorMirror<Mem::Main, double, Index>>>;

    extern template std::shared_ptr<Solver::SolverBase<MST1_::VectorType>>  SolverFactory::create_scalar_solver<
      MST1_,
      MST1_::VectorType>
      (MST1_&, PropertyMap*, const String&, std::size_t);

    extern template std::shared_ptr<Solver::SolverBase<MST1_::VectorType::LocalVectorType>>  SolverFactory::create_scalar_solver<
      MST1_,
      MST1_::VectorType::LocalVectorType>
      (MST1_&, PropertyMap*, const String&, std::size_t);

    /// \compilerhack ICC < 18 fail to link this with undefined reference to `__must_be_linked_with_icc_or_xild' error
#if not (defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1800))
    extern template std::shared_ptr<Solver::SolverBase<MST1_::VectorType>> SolverFactory::create_scalar_solver_by_section<
      MST1_,
      MST1_::VectorType>
      (MST1_&, PropertyMap*, const String&, size_t);

    extern template std::shared_ptr<Solver::SolverBase<MST1_::VectorType::LocalVectorType>> SolverFactory::create_scalar_solver_by_section<
      MST1_,
      MST1_::VectorType::LocalVectorType>
      (MST1_&, PropertyMap*, const String&, size_t);
#endif //FEAT_COMPILER_INTEL
#endif

  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_SOLVER_FACTORY_HPP
