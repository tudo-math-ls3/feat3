// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include <kernel/base_header.hpp>

#include <kernel/meshopt/rumpf_functionals/p1.hpp>
#include <kernel/meshopt/rumpf_functionals/q1.hpp>

//#include <kernel/meshopt/rumpf_functionals/2d_q1_unrolled.hpp>
//#include <kernel/meshopt/rumpf_functionals/2d_p1_unrolled.hpp>
//#include <kernel/meshopt/rumpf_functionals/3d_p1_unrolled.hpp>

#include <kernel/meshopt/hyperelasticity_functional.hpp>
#include <kernel/meshopt/mesh_concentration_function.hpp>

#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/dudv_functional_control.hpp>
#include <control/meshopt/hyperelasticity_functional_control.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Meshopt
    {
      /// \cond internal
      /**
       * \brief Template wrapper for setting two parameters
       *
       * This is used to reduce the number of free template parameters of a Meshopt::HyperelasticityFunctional
       * class template from six to four by setting the (cell-) local Functional and the RefCellTrafo.
       */
      template
      <
        template<typename, typename, typename, typename, typename> class Functional_,
        typename CellFunctional_
      >
      class SetCellFunctional
      {
        public:
          template<typename DT_, typename IT_, typename Trafo_>
          using Functional = Functional_<DT_, IT_, Trafo_, CellFunctional_, FEAT::Meshopt::RumpfTrafo<Trafo_, typename Trafo_::CoordType>>;
      };
      /// \endcond

      /**
       * \brief Factory for MeshoptControl objects
       *
       * \tparam DT_
       * Floating point type for the solver.
       *
       * \tparam IT_
       * Index type for the solver.
       *
       * This class can construct MeshoptControl objects of different types at runtime based on PropertyMaps and
       * returns them through base class std::shared_ptrs.
       *
       * \note At this time, the only transformation available is Trafo::Standard.
       *
       */
      template<typename DT_, typename IT_>
      struct ControlFactory
      {
        /**
         * \brief Creates a DuDvFunctionalControl object
         *
         * \tparam DomCtrl_
         * The domain control type.
         *
         * \param[in] dom_ctrl
         * The DomainControl containing meshes on all levels, the atlas etc.
         *
         * \param[in] section_key
         * The name of the parameter section which configures this object.
         *
         * \param[in] meshopt_config
         * Meshopt configuration.
         *
         * \param[in] solver_config
         * Solver configuration.
         *
         * \returns
         * A BaseClass std::shared_ptr to the new object
         */
        template<typename DomCtrl_>
        static std::shared_ptr <Control::Meshopt::MeshoptControlBase<DomCtrl_>>
        create_dudv_control(DomCtrl_& dom_ctrl, const String& section_key, PropertyMap* meshopt_config, PropertyMap* solver_config)
        {
          std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_>> result(nullptr);

          // -1 causes the DuDvFunctionalControl to use max level of the underlying domain control
          int meshopt_lvl(-1);

          // Get Meshopt configuration section
          auto meshopt_section = meshopt_config->query_section(section_key);
          XASSERTM(meshopt_section != nullptr, "Application config is missing the mandatory MeshOptimizer section!");

          // Verify the type
          auto type_p = meshopt_section->query("type");
          XASSERTM(type_p.second, "MeshOptimizer section is missing the mandatory type!");
          XASSERTM(type_p.first == "DuDv", "Invalid type string!");

          // Get list of boundary conditions
          auto dirichlet_list_p = meshopt_section->query("dirichlet_boundaries");
          std::deque<String> dirichlet_list = dirichlet_list_p.first.split_by_whitespaces();

          // Get list of boundary conditions
          auto slip_list_p = meshopt_section->query("slip_boundaries");
          std::deque<String> slip_list = slip_list_p.first.split_by_whitespaces();

          // Get meshopt level (if any)
          auto meshopt_lvl_p = meshopt_section->query("meshopt_lvl");
          if(meshopt_lvl_p.second)
          {
            meshopt_lvl = std::stoi(meshopt_lvl_p.first);
          }

          auto config_section_p = meshopt_section->query("config_section");
          XASSERTM(config_section_p.second, "MeshOptimizer config section is missing config_section entry!");

          auto dudv_config_section = meshopt_config->query_section(config_section_p.first);
          if(dudv_config_section == nullptr)
          {
            XABORTM("config_section "+config_section_p.first+" not found!");
          }

          auto solver_p = dudv_config_section->query("solver_config");
          XASSERTM(solver_p.second, "DuDv config section is missing solver_config entry!");

          bool fixed_reference_domain(false);
          auto fixed_reference_domain_p = dudv_config_section->query("fixed_reference_domain");
          if(fixed_reference_domain_p.second)
          {
            fixed_reference_domain = (std::stoi(fixed_reference_domain_p.first) == 1);
          }

          typedef Control::Meshopt::DuDvFunctionalControl<DT_, IT_, DomCtrl_> DuDvCtrl;
          result = std::make_shared<DuDvCtrl>(
            dom_ctrl, meshopt_lvl,
            dirichlet_list, slip_list, solver_p.first, *solver_config, fixed_reference_domain);

          return result;
        } // create_dudv_control

        /**
         * \brief Creates a HyperelasticityFunctionalControl object
         *
         * \tparam DomCtrl_
         * The domain control type.
         *
         * \param[in] dom_ctrl
         * The DomainControl containing meshes on all levels, the atlas etc.
         *
         * \param[in] section_key
         * The name of the parameter section which configures this object.
         *
         * \param[in] meshopt_config
         * Meshopt configuration.
         *
         * \param[in] solver_config
         * Solver configuration.
         *
         * \returns
         * A BaseClass std::shared_ptr to the new object
         *
         * This is the first stage, where the (cell-)local functional is determined, configured, created and then
         * passed to the next stage.
         */
        template<typename DomCtrl_>
        static std::shared_ptr <Control::Meshopt::MeshoptControlBase<DomCtrl_>>
        create_hyperelasticity_control(DomCtrl_& dom_ctrl, const String& section_key, PropertyMap* meshopt_config,
        PropertyMap* solver_config)
        {
          typedef typename DomCtrl_::LevelType::TrafoType TrafoType;

          std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_>> result(nullptr);

          DT_ fac_norm(0);
          DT_ fac_det(0);
          DT_ fac_cof(0);
          DT_ fac_reg(0);

          // -1 causes the HyperelasticityFunctionalControl to use max level of the underlying domain control
          int meshopt_lvl(-1);

          //bool split_hypercubes(false);

          // Get Meshopt configuration section
          auto meshopt_section = meshopt_config->query_section(section_key);
          XASSERTM(meshopt_section != nullptr, "Application config is missing the mandatory MeshOptimizer section!");

          // Verify the type
          auto type_p = meshopt_section->query("type");
          XASSERTM(type_p.second, "MeshOptimizer section is missing the mandatory type!");
          XASSERTM(type_p.first == "Hyperelasticity", "Invalid type!");

          // Get list of boundary conditions
          auto dirichlet_list_p = meshopt_section->query("dirichlet_boundaries");
          std::deque<String> dirichlet_list = dirichlet_list_p.first.split_by_whitespaces();

          // Get list of boundary conditions
          auto slip_list_p = meshopt_section->query("slip_boundaries");
          std::deque<String> slip_list = slip_list_p.first.split_by_whitespaces();

          // Get meshopt level (if any)
          auto meshopt_lvl_p = meshopt_section->query("meshopt_lvl");
          if(meshopt_lvl_p.second)
          {
            meshopt_lvl = std::stoi(meshopt_lvl_p.first);
          }

          // Get the name of the MeshOptimizer configuration section
          auto config_section_p = meshopt_section->query("config_section");
          if(!config_section_p.second)
          {
            XABORTM("config_section "+config_section_p.first+" not found");
          }

          // Get the MeshOptimizer configuration section
          auto hyperelasticity_config_section = meshopt_config->query_section(config_section_p.first);

          // Get fac_norm
          auto fac_norm_p = hyperelasticity_config_section->query("fac_norm");
          if(!fac_norm_p.second)
          {
            XABORTM("config_section "+config_section_p.first+" is missing fac_norm entry!");
          }
          fac_norm = DT_(std::stod(fac_norm_p.first));

          // Get fac_det
          auto fac_det_p = hyperelasticity_config_section->query("fac_det");
          if(!fac_det_p.second)
          {
            XABORTM("config_section "+config_section_p.first+" is missing fac_det entry!");
          }
          fac_det = DT_(std::stod(fac_det_p.first));

          // Get fac_cof
          auto fac_cof_p = hyperelasticity_config_section->query("fac_cof");
          if(!fac_cof_p.second)
          {
            XABORTM("config_section "+config_section_p.first+" is missing fac_cof entry!");
          }
          fac_cof = DT_(std::stod(fac_cof_p.first));

          // Get fac_reg
          auto fac_reg_p = hyperelasticity_config_section->query("fac_reg");
          if(!fac_reg_p.second)
          {
            XABORTM("config_section "+config_section_p.first+" is missing fac_reg entry!");
          }
          fac_reg = DT_(std::stod(fac_reg_p.first));

          // Get the local functional
          auto cell_functional_p = hyperelasticity_config_section->query("cell_functional");
          if(!cell_functional_p.second)
          {
            XABORTM("config_section "+config_section_p.first+" is missing cell_functional entry!");
          }

          // Get the handling of the 1/det term
          int exponent_det(0);
          auto exponent_det_p = hyperelasticity_config_section->query("exponent_det");
          if(!exponent_det_p.second)
          {
            XABORTM("config_section "+config_section_p.first+" is missing exponent_det entry!");
          }
          exponent_det = std::stoi(exponent_det_p.first);

          // Get the local functional
          if(cell_functional_p.first == "RumpfFunctional")
          {
            typedef FEAT::Meshopt::RumpfFunctional<DT_, TrafoType> CellFunctionalType;
            std::shared_ptr<CellFunctionalType> my_functional = std::make_shared<CellFunctionalType>
              (fac_norm, fac_det, fac_cof, fac_reg, exponent_det);

            result = create_hyperelasticity_control_with_cell_functional
              (dom_ctrl, meshopt_lvl, hyperelasticity_config_section, meshopt_config, solver_config, my_functional,
              dirichlet_list, slip_list);
          }
          // This is disabled because the ***Unrolled classes produce huge object files, slow code and are a pain to compile.
          // If you need them for debugging purposes, use the code below
          //else if(cell_functional_p.first == "RumpfFunctionalUnrolled")
          //{
          //  typedef FEAT::Meshopt::RumpfFunctionalUnrolled<DT_, TrafoType> CellFunctionalType;
          //  std::shared_ptr<CellFunctionalType> my_functional = std::make_shared<CellFunctionalType>
          //    (fac_norm, fac_det, fac_cof, fac_reg, exponent_det);
          //  result = create_hyperelasticity_control_with_cell_functional(dom_ctrl, hyperelasticity_config_section,
          //  meshopt_config, solver_config, my_functional, dirichlet_list, slip_list);
          //}
          else
          {
            XABORTM("Unhandled cell_functional "+cell_functional_p.first);
          }

          return result;
        } // create_hyperelasticity_control

        /**
         * \brief Creates a HyperelasticityFunctionalControl object
         *
         * \tparam DomCtrl_
         * The domain control type.
         *
         * \tparam CellFunctional_
         * The cell-local functional's type
         *
         * \param[in] dom_ctrl
         * The DomainControl containing meshes on all levels, the atlas etc.
         *
         * \param[in] meshopt_lvl
         * The level (as in domain level) to optimize the mesh on. For all other levels, this solution gets
         * prolongated/restricted.
         *
         * \param[in] hyperelasticity_config_section
         * The PropertyMap(section) containing the configuration for the new object.
         *
         * \param[in] meshopt_config
         * Meshopt configuration PropertyMap.
         *
         * \param[in] solver_config
         * Solver configuration PropertyMap.
         *
         * \param[in] my_functional
         * The cell-local functional that was created at the first stage.
         *
         * \param[in] dirichlet_list
         * List of meshpart identifiers for Dirichlet boundary conditions
         *
         * \param[in] slip_list
         * List of meshpart identifiers for slip boundary conditions
         *
         * \returns
         * A BaseClass std::shared_ptr to the new object
         *
         * This is the second stage, where the (cell-) local functional is already finished and gets passed to the
         * HyperelasticityFunctional object's constructor.
         */
        template<typename DomCtrl_, typename CellFunctional_>
        static std::shared_ptr <Control::Meshopt::MeshoptControlBase<DomCtrl_>>
        create_hyperelasticity_control_with_cell_functional( DomCtrl_& dom_ctrl,
          const int meshopt_lvl,
          PropertyMap* hyperelasticity_config_section,
          PropertyMap* meshopt_config,
          PropertyMap* solver_config,
          std::shared_ptr<CellFunctional_> my_functional,
          const std::deque<String>& dirichlet_list, const std::deque<String>& slip_list)
        {
          typedef typename DomCtrl_::LevelType::TrafoType TrafoType;

          std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_>> result(nullptr);

          // Get the global mesh quality functional
          auto global_functional_p = hyperelasticity_config_section->query("global_functional");
          if(!global_functional_p.second)
          {
            XABORTM("Hyperelasticity config section is missing global_functional entry!");
          }

          // Get scale computation type, default is once_uniform
          FEAT::Meshopt::ScaleComputation scale_computation(FEAT::Meshopt::ScaleComputation::once_uniform);

          auto scale_computation_p = hyperelasticity_config_section->query("scale_computation");
          if(scale_computation_p.second)
          {
            scale_computation << scale_computation_p.first;
          }

          int align_mesh(0);
          auto align_mesh_p = hyperelasticity_config_section->query("align_mesh");
          if(align_mesh_p.second)
          {
            align_mesh = std::stoi(align_mesh_p.first);
          }

          // Get the solver config section
          auto solver_p = hyperelasticity_config_section->query("solver_config");
          if(!solver_p.second)
          {
            XABORTM("Hyperelasticity config section is missing solver entry!");
          }

          if(global_functional_p.first == "HyperelasticityFunctional")
          {
            typedef typename FEAT::Meshopt::
              HyperelasticityFunctional<DT_, IT_, TrafoType, CellFunctional_>::RefCellTrafo RefCellTrafo;

            std::shared_ptr<FEAT::Meshopt::MeshConcentrationFunctionBase<TrafoType, RefCellTrafo>>
              mesh_conc_func(nullptr);

            if( scale_computation == FEAT::Meshopt::ScaleComputation::current_concentration ||
                scale_computation == FEAT::Meshopt::ScaleComputation::iter_concentration ||
                align_mesh == 1)
                {
                  auto conc_func_section_p = hyperelasticity_config_section->query("conc_function");
                  XASSERTM(conc_func_section_p.second, "conc_function missing!");
                  mesh_conc_func = FEAT::Meshopt::MeshConcentrationFunctionFactory<TrafoType, RefCellTrafo>::
                    create(conc_func_section_p.first, meshopt_config);
                }

            result = std::make_shared<Control::Meshopt::HyperelasticityFunctionalControl
            <DT_, IT_, DomCtrl_,
            SetCellFunctional<FEAT::Meshopt::HyperelasticityFunctional, CellFunctional_>::template Functional>>
              (dom_ctrl, meshopt_lvl, dirichlet_list, slip_list, solver_p.first,
              *solver_config, my_functional, scale_computation, mesh_conc_func, DT_(align_mesh));
          }
          else
          {
            XABORTM("Unknown global_functional "+global_functional_p.first);
          }

          return result;
        } // create_hyperelasticity_control_with_functional

        /**
         * \brief Creates a MeshoptControlBase object according to a PropertyMap
         *
         * \tparam DomCtrl_
         * The type of domain control
         *
         * \param[in] dom_ctrl
         * The domain control containing the meshes etc. for the mesh optimization
         *
         * \param[in] section_key
         * The name of the configuration section for this object
         *
         * \param[in] meshopt_config
         * The PropertyMap containing the configuration referenced by section_key
         *
         * \param[in] solver_config
         * The PropertyMap containing the solver configuration
         *
         * \returns An std::shared_ptr<MeshoptControlBase> to the new object
         */
        template<typename DomCtrl_>
        static std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_>>
        create_meshopt_control(DomCtrl_& dom_ctrl, const String& section_key,
        PropertyMap* meshopt_config, PropertyMap* solver_config)
        {
          std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl_>> result(nullptr);

          String type("");

          auto meshopt_section = meshopt_config->query_section(section_key);
          if(meshopt_section == nullptr)
          {
            XABORTM("Meshopt config is missing the ["+section_key+"] section!");
          }
          else
          {
            // Get mandatory quality functional entry
            auto type_p = meshopt_section->query("type");
            if(!type_p.second)
            {
              XABORTM("MeshOptimizer section is missing the mandatory type!");
            }
            else
            {
              type = type_p.first;
            }
          }

          if(type == "DuDv")
          {
            result = create_dudv_control(dom_ctrl, section_key, meshopt_config, solver_config);
          }
          if(type == "Hyperelasticity")
          {
            result = create_hyperelasticity_control(dom_ctrl, section_key, meshopt_config, solver_config);
          }

          if(result == nullptr)
          {
            XABORTM("MeshOptimizer section has unhandled type "+type);
          }

          return result;
        }
      }; // struct ControlFactory

    } // namespace Meshopt
  } // namespace Control
} // namespace FEAT
