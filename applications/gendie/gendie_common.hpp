#pragma once

#include "applications/gendie/logger.hpp"
#include "kernel/util/property_map.hpp"
#include <applications/gendie/fbm_control_helper.hpp>
#include <kernel/runtime.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/atlas/cgal_surface_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/q1tbnp/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/assembly/stokes_fbm_assembler.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/checkpoint_control.hpp>
#include <control/statistics.hpp>

#include <applications/gendie/materials.hpp>
#include <applications/gendie/system_assembler.hpp>
#include <applications/gendie/defect_assembler.hpp>
#include <applications/gendie/steady_stokes_solver.hpp>
#include <applications/gendie/steady_flow_solver.hpp>

namespace Gendie
{
  using namespace FEAT;
  // define application dimension
  #define FEAT_GENDIE_APP_DIM 3
  static constexpr int dim = FEAT_GENDIE_APP_DIM;
  static_assert((dim == 2) || (dim == 3), "invalid dimension");


  /// our one and only index type
#ifdef FEAT_GENDIE_APP_INDEX32
  typedef std::uint32_t SystemIndexType;
  static constexpr int ix_size = 32;
  static const char* ix_typename = "uint32";
#else
  typedef std::uint64_t SystemIndexType;
  static constexpr int ix_size = 64;
  static const char* ix_typename = "uint64";
#endif

  // depending on whether FEAT_CCND_APP_QUADMATH is defined,
  // we use quad precision or double precision
#ifdef FEAT_CCND_APP_QUADMATH
#  define Q_(x) (x##Q)
  typedef __float128 SystemDataType;
  static constexpr int fp_num_digs = 35;
  static const char* fp_typename = "quadruple";
#else
#  define Q_(x) x
  typedef double SystemDataType;
  static constexpr int fp_num_digs = 15;
  static const char* fp_typename = "double";
#endif

  // our matrix types
  template<typename DT_, typename IT_> using LocalMatrixBlockA = LAFEM::SparseMatrixBCSR<DT_, IT_, dim, dim>;
  template<typename DT_, typename IT_> using LocalMatrixBlockB = LAFEM::SparseMatrixBCSR<DT_, IT_, dim, 1>;
  template<typename DT_, typename IT_> using LocalMatrixBlockD = LAFEM::SparseMatrixBCSR<DT_, IT_, 1, dim>;
  template<typename DT_, typename IT_> using LocalScalarMatrix = LAFEM::SparseMatrixCSR<DT_, IT_>;

  // our vector types
  template<typename DT_, typename IT_> using LocalVeloVector = LAFEM::DenseVectorBlocked<DT_, IT_, dim>;
  template<typename DT_, typename IT_> using LocalPresVector = LAFEM::DenseVector<DT_, IT_>;

  // define our mesh type and other geometry related classes
  typedef Shape::Hypercube<dim> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType, dim, SystemDataType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
  typedef Geometry::Atlas::ChartBase<MeshType> ChartBaseType;
  typedef Geometry::MeshAtlas<MeshType> MeshAtlasType;

  // define our trafo type: standard or isoparametric
#ifdef FEAT_GENDIE_APP_ISOPARAM
  static constexpr bool isoparam = true;
  typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
#else
  static constexpr bool isoparam = false;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
#endif

  // define FE space types
#if defined(FEAT_GENDIE_APP_Q1T_P0)
  typedef Space::CroRavRanTur::Element<TrafoType> SpaceVeloType;
  typedef Space::Discontinuous::ElementP0<TrafoType> SpacePresType;
#elif defined(FEAT_GENDIE_APP_Q1TBNP_P1DC)
  typedef Space::Q1TBNP::Element<TrafoType> SpaceVeloType;
  typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;
#else
  typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
  typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;
#endif

  typedef FEAT::Control::StokesFBMDomainLevelBase<SpaceVeloType, SpacePresType> DomainLevelBase;

  /**
   * \brief Navier-Stokes System Level base class
   */
  template<typename DT_, typename IT_>
  using SystemLevelBase = FEAT::Control::StokesFBMSystemLevelBase<dim, LocalMatrixBlockA<DT_, IT_>, LocalMatrixBlockB<DT_, IT_>, LocalMatrixBlockD<DT_, IT_>, LocalScalarMatrix<DT_, IT_>>;

  template<typename SystemLevel_, typename Domain_, typename InflowBounds_, typename ParamHolder_>
  class SteadyFlowSolverFactory
  {
  public:
    typedef std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>> SolverPointer;
    typedef typename SystemLevel_::DataType BaseDataType;
    typedef Domain_ DomainType;
    typedef std::vector<InflowBounds_> InflowBoundsType;
    typedef std::vector<Material<BaseDataType>> MaterialsType;

    enum solverType : std::uint8_t
    {
      alpine = 0,
      backtrace_newton = 1,
      pseudo_ts = 2,
      newton = 3,
      picard = 4,
      unknown = 5
    };

  protected:
    template<typename SolverDataType_>
    using SysBase = SystemLevelBase<SolverDataType_, typename SystemLevel_::IndexType>;

    template<typename SysBase_>
    using StS = MultigridVankaFlowSolver<SysBase_, true>;

    typedef CarreauSingleMatDefectFlowAssembler<Domain_, SystemLevel_> CaD;
    typedef CarreauSingleMatSystemFlowAssembler<Domain_> CaS;

    template<typename DT_>
    using LevelDeque = std::deque<std::shared_ptr<SysBase<DT_>>>;

    struct LevelDeques
    {
      LevelDeque<double> _solver_levels_double;
      LevelDeque<float> _solver_levels_float;

      template<typename DT_> LevelDeque<DT_>& get_level()
      {
        if constexpr(std::is_same_v<DT_, float>)
          return _solver_levels_float;
        if constexpr(std::is_same_v<DT_, double>)
          return _solver_levels_double;
      }

      void clear()
      {
        _solver_levels_double.clear();
        _solver_levels_float.clear();
      }
    };
    LevelDeques _level_deques;

    DomainType& _domain;
    const InflowBoundsType& _inflow_bounds;
    const MaterialsType& _materials;

    const ParamHolder_& _param_holder;

    std::shared_ptr<SystemLevel_> _system_level;
    std::shared_ptr<CarreauSingleMatDefectFlowAssembler<Domain_, SystemLevel_>> _defect_asm;
    std::shared_ptr<CarreauSingleMatSystemFlowAssembler<Domain_>> _system_asm;

    const PropertyMap* _config;
    bool _system_assembled;

    static constexpr char const * _cubature_transfer = "gauss-legendre:3";

    void _clear_levels()
    {
      _level_deques.clear();
    }


    template<typename SolverDataType_>
    LevelDeque<SolverDataType_>& _set_level(bool system_assembled = false)
    {
      constexpr bool system_is_solver_level = std::is_same_v<BaseDataType, SolverDataType_>;
      typedef SystemLevelBase<SolverDataType_, typename SystemLevel_::IndexType> SolverLevel;
      const String cubature_transfer(_cubature_transfer);
      // parse this from config
      const String& cubature_matrix_a = _param_holder.cubature_a;
      const String& cubature_matrix_m = _param_holder.cubature_m;
      const String& cubature_matrix_b = _param_holder.cubature_b;
      bool use_q2_fbm = _param_holder.assemble_q2_filter;
      bool use_coarse_fbm = true;
      std::size_t min_q2_fbm_level = ~std::size_t(0);
      BaseDataType mesh_unit_scale = _param_holder.mesh_unit_scale;

      //clear levels
      this->_clear_levels();
      //get our solver levels
      auto& solver_levels = _level_deques.template get_level<SolverDataType_>();

      const Index num_levels(this->_domain.size_physical());
      solver_levels.resize(num_levels);
      if constexpr(system_is_solver_level)
      {
        solver_levels.at(0) = this->_system_level;
        solver_levels.at(0)->clear_transfers();
      }
      else
      {
        solver_levels.at(0).reset(new SolverLevel());
        // and always convert from our system level
        solver_levels.at(0)->convert(*this->_system_level);
      }
      for (Index i(1); i < num_levels; ++i)
      {
        solver_levels.at(i).reset(new SolverLevel());
        solver_levels.at(i)->assemble_gates(this->_domain.at(i));
      }
      // assemble transfers
      for(Index i(0); (i < this->_domain.size_physical()) && ((i+1) < this->_domain.size_virtual()); ++i)
      {
        solver_levels.at(i)->assemble_coarse_muxers(this->_domain.at(i+1));
        if((i+1) < this->_domain.size_physical())
          solver_levels.at(i)->assemble_transfers_voxel(*solver_levels.at(i+1), this->_domain.at(i), this->_domain.at(i+1), cubature_transfer, true);
        else
          solver_levels.at(i)->assemble_transfers_voxel(this->_domain.at(i), this->_domain.at(i+1), cubature_transfer, true);
      }

      // extra logic for finest level for matrix assembly
      {
        // get domain and system levels
        auto& dom_lvl = *this->_domain.at(0);
        auto& sys_lvl = *solver_levels.at(0);

        // assemble matrix structures
        if(!system_assembled)
        {
          // if nothing is assembled, we have to do the deed now
          sys_lvl.assemble_velo_struct(dom_lvl.space_velo);
          sys_lvl.assemble_pres_struct(dom_lvl.space_pres);
          sys_lvl.assemble_velocity_laplace_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, SolverDataType_(this->_materials.front().get_viscosity_model()->get_data()[0] * 1E3/mesh_unit_scale), true, cubature_matrix_a);
          sys_lvl.assemble_velocity_mass_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, cubature_matrix_m);
          if constexpr(system_is_solver_level)
            sys_lvl.assemble_grad_div_matrices(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b);
          else
            sys_lvl.assemble_grad_div_matrices_high_prec(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b, BaseDataType(1));
        }
        else
        {
          if constexpr(!system_is_solver_level)
          {
            // in any other case, we already assembled the system, so we need to convert the missing members not cloned explicitly
            sys_lvl.velo_mass_matrix.convert(&sys_lvl.gate_velo, &sys_lvl.gate_velo , this->_system_level->velo_mass_matrix);
            sys_lvl.fbm_mask_velo = this->_system_level->fbm_mask_velo;
            sys_lvl.fbm_mask_pres = this->_system_level->fbm_mask_pres;
          }
        }

        // compile the system matrix
        sys_lvl.compile_system_matrix();

      }
      // assemble matrices
      for(Index i(1); i < num_levels; ++i)
      {
        // get domain and system levels
        auto& dom_lvl = *this->_domain.at(i);
        auto& sys_lvl = *solver_levels.at(i);

        sys_lvl.assemble_velo_struct(dom_lvl.space_velo);
        sys_lvl.assemble_pres_struct(dom_lvl.space_pres);
        sys_lvl.assemble_velocity_laplace_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, SolverDataType_(this->_materials.front().get_viscosity_model()->get_data()[0] * 1E+3/mesh_unit_scale), true, cubature_matrix_a);
        sys_lvl.assemble_velocity_mass_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, cubature_matrix_m);
        if constexpr(system_is_solver_level)
          sys_lvl.assemble_grad_div_matrices(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b);
        else
          sys_lvl.assemble_grad_div_matrices_high_prec(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b, BaseDataType(1));

        // compile the system matrix
        sys_lvl.compile_system_matrix();
      }

      // finally, we have to define our system filter
      // todo: export this to a component?
      {
        std::size_t start = system_assembled ? 1 : 0;
        if(system_assembled)
        {
          if constexpr(system_is_solver_level)
          {
            //nothing to do
          }
          else
          {
            // clone remaining filter
            solver_levels.front()->filter_interface_fbm.convert(this->_system_level->filter_interface_fbm);
            solver_levels.front()->compile_system_filter();
          }
        }

        // the names of the mesh parts on which to assemble
        std::deque<String> part_names = this->_domain.front()->get_mesh_node()->get_mesh_part_names(true);
        for(std::size_t i(start); i < solver_levels.size(); ++i)
        {
          // get our local velocity filters
          auto& filter_v_noflow = solver_levels.at(i)->get_local_velo_unit_filter_seq().find_or_add("noflow");
          auto& filter_v_inflow = solver_levels.at(i)->get_local_velo_unit_filter_seq().find_or_add("inflow");

          // create unit-filter assembler
          Assembly::UnitFilterAssembler<MeshType> /*unit_asm_inflow,*/ unit_asm_noflow;
          // loop over all boundary parts except for the right one, which is outflow
          for(const auto& name : part_names)
          {
            // skip non-boundary mesh-parts
            if((!name.starts_with("bnd:")))
              continue;

            // try to fetch the corresponding mesh part node
            auto* mesh_part_node = this->_domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
            XASSERT(mesh_part_node != nullptr);

            auto* mesh_part = mesh_part_node->get_mesh();
            String name_lower = name.lower();
            if (mesh_part != nullptr)
            {
              // we require the inflow to be of type "inflow_<type>_id"
              if((name_lower.replace_all(String("inflow"), String("")) > 0u))
              {
                auto name_deq = name.split_by_charset("_");
                if(name_deq.size() < 3)
                {
                  XABORTM("Can not parse string: " + name);
                }
                int id = -1;
                name_deq.back().parse(id);
                XASSERTM(id >= 0 && std::size_t(id) < this->_inflow_bounds.size(), "Id does not map to boundaries");
                Assembly::UnitFilterAssembler<MeshType> loc_unit_asm_inflow;

                // inflow
                loc_unit_asm_inflow.add_mesh_part(*mesh_part);
                if constexpr(system_is_solver_level)
                {
                  loc_unit_asm_inflow.assemble(filter_v_inflow, this->_domain.at(i)->space_velo, this->_inflow_bounds.at(std::size_t(id)).get_diri_inflow_function(this->_materials, mesh_unit_scale));
                }
                else
                {
                  // create temporary filter (of correct datatype)
                  typename SystemLevel_::LocalVeloUnitFilter fil_loc_v;
                  // directly assembly onto the inflow filter
                  loc_unit_asm_inflow.assemble(fil_loc_v, this->_domain.at(i)->space_velo, this->_inflow_bounds.at(std::size_t(id)).get_diri_inflow_function(this->_materials, mesh_unit_scale));
                  //and now convert
                  if(filter_v_inflow.size() == Index(0))
                    filter_v_inflow.convert(fil_loc_v);
                  else
                  {
                    if(fil_loc_v.used_elements() == 0u) // nothing to add here
                      continue;
                    const auto* indx = fil_loc_v.get_indices();
                    const auto* valx = fil_loc_v.get_values();
                    for(Index k = 0; k < fil_loc_v.used_elements(); ++k)
                    {
                      if constexpr(dim == 3)
                      {
                        FEAT::Tiny::Vector<SolverDataType_, dim> valn{SolverDataType_(valx[k][0]), SolverDataType_(valx[k][1]), SolverDataType_(valx[k][2])};
                        filter_v_inflow.add(indx[k], valn);
                      }
                    }
                  }
                }
              }
              else if((name != "bnd:n") && (!name.starts_with("bnd:out")) && ((name_lower.replace_all(String("outflow"), String("")) == 0u)))
              {
                // outflow
                unit_asm_noflow.add_mesh_part(*mesh_part);
              }
            }
          }
          {
            unit_asm_noflow.assemble(filter_v_noflow, this->_domain.at(i)->space_velo);
          }
          // assemble fbm parts
          {
            auto& fbm_asm = *this->_domain.at(i)->fbm_assembler;
            const bool q2_fbm = use_q2_fbm && ((unsigned short)(i) <= min_q2_fbm_level || use_coarse_fbm);
            solver_levels.at(i)->assemble_fbm_filters(fbm_asm, this->_domain.at(i)->space_velo, this->_domain.at(i)->space_pres, (i==0u), q2_fbm, false);
          }
          // compile system filter
          solver_levels.at(i)->compile_system_filter();
        }
      }

      return solver_levels;
      //todo: rescale interface filter of base system??

    }

    template<typename SolverDataType_>
    LevelDeque<SolverDataType_>& _get_level(bool system_assembled = false)
    {
      if(_level_deques.template get_level<SolverDataType_>().empty())
        return _set_level<SolverDataType_>(system_assembled);

      return _level_deques.template get_level<SolverDataType_>();
    }

    template<typename SolverDataType_>
    SolverPointer _create(solverType solver_type, const FEAT::PropertyMap* sol_config, const Gendie::Logger* logger = nullptr)
    {
      // constexpr bool system_is_solver_level = std::is_same_v<BaseDataType, SolverDataType_>;
      typedef SystemLevelBase<SolverDataType_, typename SystemLevel_::IndexType> SolverLevel;

      auto& solver_levels = _get_level<SolverDataType_>(this->_system_assembled);

      // system is assembled, now we need to create the components
      auto stokes_solver = std::make_shared<MultigridVankaFlowSolver<SolverLevel, true>>(solver_levels, this->_domain, sol_config, logger);
      auto scaling_factor = _defect_asm->calc_stiffness_factor();

      _defect_asm->set_stiffness_factor(scaling_factor);
      _system_asm->set_stiffness_factor(scaling_factor);

      // at this point, our base system is always assembled
      _system_assembled = true;

      switch(solver_type)
      {
        typedef StS<SysBase<SolverDataType_>> StSType;
        case solverType::alpine:
        {

          auto* sys_solver_ptr = new AlPiNeSteadyFlowSolver<StSType, CaD, CaS, decltype(this->_system_level->filter_sys)>(std::move(stokes_solver), this->_defect_asm, this->_system_asm, this->_system_level->filter_sys, logger);
          sys_solver_ptr->set_scaling_factor(scaling_factor);
          sys_solver_ptr->parse(sol_config);
          return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(sys_solver_ptr);
        }
        case solverType::backtrace_newton:
        {
          auto* sys_solver_ptr = new BackTraceNewtonSteadyFlowSolver<StSType, CaD, CaS, decltype(this->_system_level->filter_sys)>(std::move(stokes_solver), this->_defect_asm, this->_system_asm, this->_system_level->filter_sys, logger);
          sys_solver_ptr->set_scaling_factor(scaling_factor);
          sys_solver_ptr->parse(sol_config);
          return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(sys_solver_ptr);
        }
        case solverType::picard:
        {
          XABORTM("Picard Not implemented");
          return std::unique_ptr<NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(nullptr);
        }
        case solverType::newton:
        {
          XABORTM("Newton Not implemented");
          return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(nullptr);
        }
        case solverType::pseudo_ts:
        {
          auto* sys_solver_ptr = new PseudoUnsteadyFlowSolver<StSType, CaD, CaS, decltype(this->_system_level->filter_sys)>(std::move(stokes_solver), this->_defect_asm, this->_system_asm, this->_system_level->filter_sys, logger);
          sys_solver_ptr->set_scaling_factor(scaling_factor);
          sys_solver_ptr->parse(sol_config);
          return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(sys_solver_ptr);
        }
        case solverType::unknown:
        {
          XABORTM("Unkown solver Type");
          return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(nullptr);
        }
      }

      XABORTM("Option " + stringify(int(solver_type)) + String(" not known"));
      return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(nullptr);

    }




  public:

    explicit SteadyFlowSolverFactory(
      Domain_& domain, std::shared_ptr<SystemLevel_>& system_level, const std::vector<InflowBounds_>& inflow_bounds,
      const std::vector<Material<typename SystemLevel_::DataType>>& materials, const PropertyMap* config, bool system_assembled, const ParamHolder_& param_holder)
      :
       _level_deques(),
       _domain(domain),
       _inflow_bounds(inflow_bounds),
       _materials(materials),
       _param_holder(param_holder),
       _system_level(system_level),
       _defect_asm(nullptr),
       _system_asm(nullptr),
       _config(config),
       _system_assembled(system_assembled)
      {
        const String& cubature_matrix_a = param_holder.cubature_a;
        SystemDataType mesh_unit_scale = param_holder.mesh_unit_scale;
        _defect_asm = std::make_shared<CarreauSingleMatDefectFlowAssembler<Domain_, SystemLevel_>>(_domain, _system_level, _config, _materials.front(), _inflow_bounds.front().get_temperature(), mesh_unit_scale);
        _defect_asm->set_cubature(cubature_matrix_a);
        _system_asm = std::make_shared<CarreauSingleMatSystemFlowAssembler<Domain_>>(_domain, _config, _materials.front(), _inflow_bounds.front().get_temperature(), mesh_unit_scale);
        _system_asm->set_cubature(cubature_matrix_a);

      }



    SolverPointer create(solverType solver_type, Gendie::Logger* logger = nullptr)
    {
      // create the stokes solver
      String datatype = "double";
      if(this->_config)
        XASSERTM(this->_config->parse_entry("datatype", datatype, true), "Could not parse datatype");

      if(datatype.compare_no_case("double") == 0)
      {
        return _create<double>(solver_type, this->_config, logger);
      }
      else if(datatype.compare_no_case("float") == 0)
      {
        return _create<float>(solver_type, this->_config, logger);
      }
      return _create<SystemDataType>(solver_type, this->_config, logger);
    }

    template<typename DT_>
    SolverPointer create(solverType solver_type, DT_ DOXY(solver_dt), Gendie::Logger* logger = nullptr)
    {
      // create the stokes solver
      return _create<DT_>(solver_type, this->_config, logger);
    }

    static solverType parse_solver_type(const String& solver_string)
    {
      if(solver_string.compare_no_case("alpine") == 0)
      {
        return solverType::alpine;
      }
      else if(solver_string.compare_no_case("backtrace-newton") == 0)
      {
        return solverType::backtrace_newton;
      }
      else if(solver_string.compare_no_case("pseudo-ts") == 0)
      {
        return solverType::pseudo_ts;
      }
      else if(solver_string.compare_no_case("newton") == 0)
      {
        return solverType::newton;
      }
      else if(solver_string.compare_no_case("picard") == 0)
      {
        return solverType::picard;
      }
      return solverType::unknown;
    }

    /// for now we assume the create function is not called after this
    void reset()
    {
      this->_defect_asm.reset();
      this->_system_asm.reset();
      this->_level_deques.clear();
      this->_system_level.reset();
    }

  }; // class SteadyFlowFactory

}
