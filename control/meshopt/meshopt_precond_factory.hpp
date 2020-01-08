// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
      /**
       * \brief Factory class for specialised preconditioners
       */
      struct MeshoptPrecondFactory
      {
        /**
         * \brief Creates a preconditioner for nonlinear mesh optimisation problems
         *
         * \tparam MeshoptCtrl_
         * The type of the control object for mesh optimisation
         *
         * \tparam DomCtrl_
         * The type of the control object for the domain.
         *
         * \tparam SolverVectorType_
         * The input vector type for the preconditioner.
         *
         * \tparam FilterType_
         * The filter type for the preconditioner.
         *
         * \param[in] my_ctrl
         * The mesh optimisation control object.
         *
         * \param[in] dom_ctrl
         * The domain control object.
         *
         * \param[in] current_section
         * The configuration section of the solver this preconditioner is built for.
         *
         * \returns An std::shared_ptr to the new NLOptPrecond object.
         */
        template
        <
          typename MeshoptCtrl_,
          typename DomCtrl_ = typename MeshoptCtrl_::DomainControlType,
          typename SolverVectorType_ = typename MeshoptCtrl_::SystemLevelType::GlobalSystemVectorL,
          typename FilterType_ = typename MeshoptCtrl_::SystemLevelType::GlobalSystemFilter
        >
        static std::shared_ptr<Solver::NLOptPrecond<SolverVectorType_, FilterType_>>
        create_nlopt_precond(MeshoptCtrl_& my_ctrl, DomCtrl_& dom_ctrl, PropertyMap* current_section)
        {
          std::shared_ptr<Solver::NLOptPrecond<SolverVectorType_, FilterType_> > result;

          auto& solver_config = my_ctrl.solver_config;

          auto precon_p = current_section->query("precon");

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

              auto dirichlet_list_p = precon_section->query("dirichlet_boundaries");
              std::deque<String> precon_dirichlet_list = dirichlet_list_p.first.split_by_whitespaces();

              auto slip_list_p = precon_section->query("slip_boundaries");
              std::deque<String> precon_slip_list = slip_list_p.first.split_by_whitespaces();

              bool fixed_reference_domain(false);
              auto fixed_reference_domain_p = precon_section->query("fixed_reference_domain");
              if(fixed_reference_domain_p.second)
              {
                fixed_reference_domain = (std::stoi(fixed_reference_domain_p.first) == 1);
              }

              int level(my_ctrl.meshopt_lvl);

              typedef DuDvFunctionalControl
              <
                typename MeshoptCtrl_::MemType,
                typename MeshoptCtrl_::DataType,
                typename MeshoptCtrl_::IndexType,
                DomCtrl_
              > PreconControlType;

              result = Solver::new_nonlinear_operator_precond_wrapper<PreconControlType>
                (dom_ctrl, level, precon_dirichlet_list, precon_slip_list,
                precon_solver_p.first, solver_config, fixed_reference_domain);
            }
            else if(precon_p.first != "none")
            {
              throw InternalError(__func__,__FILE__,__LINE__,
              "Unsupport nonlinear optimiser precon: "+precon_p.first);
            }
          }
          else
          {
            auto inner_solver_p = current_section->query("inner_solver");
            if(inner_solver_p.second)
            {
              auto inner_solver_section = solver_config.query_section(inner_solver_p.first);
              result = create_nlopt_precond(my_ctrl, dom_ctrl, inner_solver_section);
            }
          }

          return result;
        } // create_nlopt_precond
      };
    } // namespace Meshopt
  } // namespace Control
} //namespace FEAT
#endif // CONTROL_MESHOPT_MESHOPT_PRECOND_FACTORY_HPP
