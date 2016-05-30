#pragma once
#ifndef CONTROL_MESHOPT_MESHOPT_PRECOND_FACTORY_HPP
#define CONTROL_MESHOPT_MESHOPT_PRECOND_FACTORY_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/solver/nlopt_precond.hpp>

#include <control/meshopt/dudv_functional_control.hpp>

#include <deque>

namespace FEAT
{
  namespace Control
  {
    namespace Meshopt
    {
      struct MeshoptPrecondFactory
      {
        template
        <
          typename MeshoptCtrl_,
          typename DomCtrl_ = typename MeshoptCtrl_::DomainControlType,
          typename SolverVectorType_ = typename MeshoptCtrl_::SystemLevelType::GlobalSystemVectorR,
          typename FilterType_ = typename MeshoptCtrl_::SystemLevelType::GlobalSystemFilter
        >
        static std::shared_ptr<Solver::NLOptPrecond<SolverVectorType_, FilterType_>>
        create_nlopt_precond(MeshoptCtrl_& my_ctrl, DomCtrl_& dom_ctrl)
        {
          std::shared_ptr<Solver::NLOptPrecond<SolverVectorType_, FilterType_> > result;

          auto& solver_config = my_ctrl.solver_config;

          // Check the solver section to see if we have to create a precon_control
          auto solver_section = solver_config.query_section(my_ctrl.solver_name);
          if(solver_section == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,"Solver section "+my_ctrl.solver_name+" not found!");

          auto precon_p = solver_section->query("precon");

          if(precon_p.second && precon_p.first != "none")
          {
            auto precon_section = solver_config.query_section(precon_p.first);
            if(precon_section == nullptr)
              throw InternalError(__func__,__FILE__,__LINE__,"precon_section "+precon_p.first+" not found!");

            auto precon_type_p = precon_section->query("type");
            if(!precon_type_p.second)
              throw InternalError(__func__,__FILE__,__LINE__,
              "Section "+precon_p.first+" does not have a type key!");

            if(precon_type_p.first == "DuDvPrecon")
            {
              auto precon_solver_p = precon_section->query("linear_solver");
              if(!precon_solver_p.second)
                throw InternalError(__func__,__FILE__,__LINE__,
                "precon_section "+precon_type_p.first+" has no linear_solver key!");

              std::deque<String> precon_dirichlet_list;
              std::deque<String> precon_slip_list;

              auto dirichlet_list_p = precon_section->query("dirichlet_boundaries");
              dirichlet_list_p.first.split_by_charset(precon_dirichlet_list, " ");

              auto slip_list_p = precon_section->query("slip_boundaries");
              slip_list_p.first.split_by_charset(precon_slip_list, " ");

              typedef DuDvFunctionalControl
              <
                typename MeshoptCtrl_::MemType,
                typename MeshoptCtrl_::DataType,
                typename MeshoptCtrl_::IndexType,
                DomCtrl_,
                typename MeshoptCtrl_::TrafoType
              > PreconControlType;

              result = Solver::new_nonlinear_operator_precond_wrapper<PreconControlType>
                (dom_ctrl, precon_dirichlet_list, precon_slip_list, precon_solver_p.first, solver_config);
            }
            else if(precon_p.first != "none")
              throw InternalError(__func__,__FILE__,__LINE__,
              "Unsupport nonlinear optimiser precon: "+precon_p.first);
          }

          return result;
        } // create_nonlinear_optimiser_precon
      };
    } // namespace Meshopt
  } // namespace Control
} //namespace FEAT
#endif // CONTROL_MESHOPT_MESHOPT_PRECOND_FACTORY_HPP
