#pragma once
#ifndef FEAT_CONTROL_MESHOPT_MESHOPT_CONTROL_FACTORY_HPP
#define FEAT_CONTROL_MESHOPT_MESHOPT_CONTROL_FACTORY_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/meshopt/hyperelasticity_functional.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1split.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>

#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_solver_factory.hpp>
#include <control/meshopt/dudv_functional_control.hpp>
#include <control/meshopt/hyperelasticity_functional_control.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Meshopt
    {
      template
      <
        template<typename, typename, typename, typename, typename, typename> class Functional_,
        typename LocalFunctional_
      >
      class SetFunctional
      {
        public:
          template<typename A, typename B, typename C, typename D>
          using Functional = Functional_<A, B, C, D, LocalFunctional_, FEAT::Meshopt::RumpfTrafo<D, typename D::CoordType>>;
      };

      template<typename Mem_, typename DT_, typename IT_, typename Trafo_>
      struct ControlFactory
      {
        typedef typename Trafo_::MeshType MeshType;

        template<typename DomCtrl_>
        static std::shared_ptr <Control::Meshopt::MeshoptControlBase<DomCtrl_, Trafo_>>
        create_dudv_control(DomCtrl_& dom_ctrl, const String& section_key, PropertyMap* meshopt_config, PropertyMap* solver_config)
        {
          std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_, Trafo_>> result(nullptr);
          std::deque<String> dirichlet_list;
          std::deque<String> slip_list;

          // Get Meshopt configuration section
          auto meshopt_section = meshopt_config->query_section(section_key);
          if(meshopt_section == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,
            "Application config is missing the mandatory MeshOptimiser section!");

          // Verify the type
          auto type_p = meshopt_section->query("type");
          if(!type_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "MeshOptimiser section is missing the mandatory type!");
          else if (type_p.first != "DuDv")
            throw InternalError(__func__,__FILE__,__LINE__,
            "Invalid type "+type_p.first);

          // Get list of boundary conditions
          auto dirichlet_list_p = meshopt_section->query("dirichlet_boundaries");
          dirichlet_list_p.first.split_by_charset(dirichlet_list, " ");

          // Get list of boundary conditions
          auto slip_list_p = meshopt_section->query("slip_boundaries");
          slip_list_p.first.split_by_charset(slip_list, " ");

          auto config_section_p = meshopt_section->query("config_section");
          if(!config_section_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "MeshOptimiser config section is missing config_section entry!");

          auto dudv_config_section = meshopt_config->query_section(config_section_p.first);
          if(dudv_config_section == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,
            "config_section "+config_section_p.first+" not found!");

          auto solver_p = dudv_config_section->query("solver_config");
          if(!solver_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "DuDv config section is missing solver entry!");

          typedef Control::Meshopt::DuDvFunctionalControl<Mem_, DT_, IT_, DomCtrl_, Trafo_> DuDvCtrl;
          result = std::make_shared<DuDvCtrl>(dom_ctrl, dirichlet_list, slip_list, solver_p.first, *solver_config);

          return result;
        } // create_dudv_control

        template<typename DomCtrl_>
        static std::shared_ptr <Control::Meshopt::MeshoptControlBase<DomCtrl_, Trafo_>>
        create_hyperelasticity_control(DomCtrl_& dom_ctrl, const String& section_key, PropertyMap* meshopt_config,
        PropertyMap* solver_config)
        {
          std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_, Trafo_>> result(nullptr);

          std::deque<String> dirichlet_list;
          std::deque<String> slip_list;

          DT_ fac_norm(0);
          DT_ fac_det(0);
          DT_ fac_cof(0);
          DT_ fac_reg(0);

          bool split_hypercubes(false);

          // Get Meshopt configuration section
          auto meshopt_section = meshopt_config->query_section(section_key);
          if(meshopt_section == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,
            "Application config is missing the mandatory MeshOptimiser section!");

          // Verify the type
          auto type_p = meshopt_section->query("type");
          if(!type_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "MeshOptimiser section is missing the mandatory type!");
          else if (type_p.first != "Hyperelasticity")
            throw InternalError(__func__,__FILE__,__LINE__,
            "Invalid type "+type_p.first);

          // Get list of boundary conditions
          auto dirichlet_list_p = meshopt_section->query("dirichlet_boundaries");
          dirichlet_list_p.first.split_by_charset(dirichlet_list, " ");

          // Get list of boundary conditions
          auto slip_list_p = meshopt_section->query("slip_boundaries");
          slip_list_p.first.split_by_charset(slip_list, " ");

          // Get the name of the MeshOptimiser configuration section
          auto config_section_p = meshopt_section->query("config_section");
          if(!config_section_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "config_section "+config_section_p.first+" not found");

          // Get the MeshOptimiser configuration section
          auto hyperelasticity_config_section = meshopt_config->query_section(config_section_p.first);

          // Get fac_norm
          auto fac_norm_p = hyperelasticity_config_section->query("fac_norm");
          if(!fac_norm_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "config_section "+config_section_p.first+" is missing fac_norm entry!");
          fac_norm = DT_(std::stod(fac_norm_p.first));

          // Get fac_det
          auto fac_det_p = hyperelasticity_config_section->query("fac_det");
          if(!fac_det_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "config_section "+config_section_p.first+" is missing fac_det entry!");
          fac_det = DT_(std::stod(fac_det_p.first));

          // Get fac_cof
          auto fac_cof_p = hyperelasticity_config_section->query("fac_cof");
          if(!fac_cof_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "config_section "+config_section_p.first+" is missing fac_cof entry!");
          fac_cof = DT_(std::stod(fac_cof_p.first));

          // Get fac_reg
          auto fac_reg_p = hyperelasticity_config_section->query("fac_reg");
          if(!fac_reg_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "config_section "+config_section_p.first+" is missing fac_reg entry!");
          fac_reg = DT_(std::stod(fac_reg_p.first));

          // Get the local functional
          auto local_functional_p = hyperelasticity_config_section->query("local_functional");
          if(!local_functional_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "config_section "+config_section_p.first+" is missing local_functional entry!");

          // Check if we want to split hypercubes into simplices for the evaluation of the local functional
          auto split_hc_p = hyperelasticity_config_section->query("split_hypercubes");
          if(split_hc_p.second)
          split_hypercubes = (std::stoi(split_hc_p.first) == 1);

          // Get the local functional
          if(local_functional_p.first == "RumpfFunctional")
          {
            if(split_hypercubes)
            {
              // The underlying functional is RumpfFunctional, and this is passed to the Split functional as parameter
              typedef FEAT::Meshopt::RumpfFunctionalQ1Split
              <
                DT_,
                typename MeshType::ShapeType,
                FEAT::Meshopt::RumpfFunctional
              > FunctionalType;

              std::shared_ptr<FunctionalType> my_functional = std::make_shared<FunctionalType>
                (fac_norm, fac_det, fac_cof, fac_reg);
              result = create_hyperelasticity_control_with_functional
                (dom_ctrl, hyperelasticity_config_section, solver_config, my_functional, dirichlet_list, slip_list);
            }
            else
            {
              typedef FEAT::Meshopt::RumpfFunctional<DT_, typename MeshType::ShapeType> FunctionalType;
              std::shared_ptr<FunctionalType> my_functional = std::make_shared<FunctionalType>
                (fac_norm, fac_det, fac_cof, fac_reg);

              result = create_hyperelasticity_control_with_functional
                (dom_ctrl, hyperelasticity_config_section, solver_config, my_functional, dirichlet_list, slip_list);
            }
          }
          else if(local_functional_p.first == "RumpfFunctional_D2")
          {
            if(split_hypercubes)
            {
              // The underlying functional is RumpfFunctional_D2, and this is passed to the Split functional as parameter
              typedef FEAT::Meshopt::RumpfFunctionalQ1Split
              <
                DT_,
                typename MeshType::ShapeType,
                FEAT::Meshopt::RumpfFunctional_D2
              > FunctionalType;

              std::shared_ptr<FunctionalType> my_functional = std::make_shared<FunctionalType>
                (fac_norm, fac_det, fac_cof, fac_reg);
              result = create_hyperelasticity_control_with_functional(dom_ctrl, hyperelasticity_config_section,
              solver_config, my_functional, dirichlet_list, slip_list);
            }
            else
            {
              typedef FEAT::Meshopt::RumpfFunctional_D2<DT_, typename MeshType::ShapeType> FunctionalType;
              std::shared_ptr<FunctionalType> my_functional = std::make_shared<FunctionalType>
                (fac_norm, fac_det, fac_cof, fac_reg);
              result = create_hyperelasticity_control_with_functional(dom_ctrl, hyperelasticity_config_section,
              solver_config, my_functional, dirichlet_list, slip_list);
            }
          }
          else
            throw InternalError(__func__,__FILE__,__LINE__,"Unhandled local_functional "+local_functional_p.first);
          return result;
        } // create_hyperelasticity_control

        template<typename DomCtrl_, typename FunctionalType_>
        static std::shared_ptr <Control::Meshopt::MeshoptControlBase<DomCtrl_, Trafo_>>
        create_hyperelasticity_control_with_functional( DomCtrl_& dom_ctrl,
        PropertyMap* hyperelasticity_config_section, PropertyMap* solver_config,
        std::shared_ptr<FunctionalType_> my_functional,
        const std::deque<String>& dirichlet_list, const std::deque<String>& slip_list)
        {
          std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_, Trafo_>> result(nullptr);

          // Get the global mesh quality functional
          auto global_functional_p = hyperelasticity_config_section->query("global_functional");
          if(!global_functional_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "Hyperelasticity config section is missing global_functional entry!");

          // Get scale computation type, default is once_uniform
          FEAT::Meshopt::ScaleComputation scale_computation(FEAT::Meshopt::ScaleComputation::once_uniform);
          auto scale_computation_p = hyperelasticity_config_section->query("scale_computation");
          if(scale_computation_p.second)
            scale_computation << scale_computation_p.first;

          std::deque<String> distance_charts;
          auto distance_charts_p = hyperelasticity_config_section->query("distance_charts");
          if(distance_charts_p.second)
            distance_charts_p.first.split_by_charset(distance_charts, " ");

          // Get the solver config section
          auto solver_p = hyperelasticity_config_section->query("solver_config");
          if(!solver_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "Hyperelasticity config section is missing solver entry!");

          if(global_functional_p.first == "HyperelasticityFunctional")
          {
            result = std::make_shared<Control::Meshopt::HyperelasticityFunctionalControl
            <Mem_, DT_, IT_, DomCtrl_, Trafo_,
            SetFunctional<FEAT::Meshopt::HyperelasticityFunctional, FunctionalType_>::template Functional>>
              (dom_ctrl, dirichlet_list, slip_list, solver_p.first, *solver_config, my_functional,
              scale_computation, distance_charts);
          }
          else
            throw InternalError(__func__,__FILE__,__LINE__,
            "Unknown global_functional "+global_functional_p.first);

          return result;
        } // create_hyperelasticity_control_with_functional

        template<typename DomCtrl_>
        static std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_, Trafo_>>
        create_meshopt_control(DomCtrl_& dom_ctrl, const String& section_key,
        PropertyMap* meshopt_config, PropertyMap* solver_config)
        {
          std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_, Trafo_>> result(nullptr);

          String type("");
          std::deque<String> dirichlet_list;
          std::deque<String> slip_list;

          auto meshopt_section = meshopt_config->query_section(section_key);
          if(meshopt_section == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,
            "Meshopt config is missing the ["+section_key+"] section!");
          else
          {
            // Get mandatory quality functional entry
            auto type_p = meshopt_section->query("type");
            if(!type_p.second)
              throw InternalError(__func__,__FILE__,__LINE__,
              "MeshOptimiser section is missing the mandatory type!");
            else
              type = type_p.first;

            auto dirichlet_list_p = meshopt_section->query("dirichlet_boundaries");
            dirichlet_list_p.first.split_by_charset(dirichlet_list, " ");

            auto slip_list_p = meshopt_section->query("slip_boundaries");
            slip_list_p.first.split_by_charset(slip_list, " ");

          }

          if(type == "DuDv")
          {
            result = create_dudv_control(dom_ctrl, section_key, meshopt_config, solver_config);
          }
          if(type == "Hyperelasticity")
          {
            result = create_hyperelasticity_control(dom_ctrl, section_key, meshopt_config, solver_config);
          }

          //if(result == nullptr)
          //{
          //  throw InternalError(__func__,__FILE__,__LINE__,
          //  "MeshOptimiser section has unhandled type "+type);
          //}

          return result;
        }
      }; // struct ControlFactory

    } // namespace Meshopt
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_MESHOPT_MESHOPT_SOLVER_FACTORY_HPP