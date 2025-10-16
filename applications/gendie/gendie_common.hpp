#pragma once

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

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class DomainLevelBase :
    public Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType>
  {
  public:
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> BaseClass;

    /// the FBM assembler for this level
    std::unique_ptr<Assembly::StokesFBMAssembler<MeshType>> fbm_assembler;

    // inherit constructor
    using BaseClass::BaseClass;

    void create_fbm_assembler(const Dist::Comm& comm, const String& fbm_meshpart_names)
    {
      // get out mesh node
      const MeshNodeType& mesh_node = *this->get_mesh_node();

      // allocate assembler
      fbm_assembler.reset(new Assembly::StokesFBMAssembler<MeshType>(*mesh_node.get_mesh()));

      // loop over all FBM meshparts
      std::deque<String> fbm_name_deque = fbm_meshpart_names.split_by_whitespaces();
      for(const String& name : fbm_name_deque)
      {
        // find meshpart node
        const auto* part_node = mesh_node.find_mesh_part_node(name);
        XASSERTM(part_node != nullptr, "FBM meshpart node '" + name + "'not found");

        // get the meshpart if it exists
        const MeshPartType* fbm_meshpart = part_node->get_mesh();
        if(fbm_meshpart)
          fbm_assembler->add_fbm_meshpart(*fbm_meshpart);
      }
      // synchronize over comm
      fbm_assembler->sync(mesh_node, comm);

      // compile the assembler
      fbm_assembler->compile();
    }
  }; // class DomainLevel

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * \brief Navier-Stokes System Level base class
   */
  template<typename DT_, typename IT_>
  class SystemLevelBase :
    public Control::StokesBlockedCombinedSystemLevel<dim, DT_, IT_, LocalMatrixBlockA<DT_, IT_>, LocalMatrixBlockB<DT_, IT_>, LocalMatrixBlockD<DT_, IT_>, LocalScalarMatrix<DT_, IT_>>
  {
  public:
    /// out base class
    typedef Control::StokesBlockedCombinedSystemLevel<dim, DT_, IT_, LocalMatrixBlockA<DT_, IT_>, LocalMatrixBlockB<DT_, IT_>, LocalMatrixBlockD<DT_, IT_>, LocalScalarMatrix<DT_, IT_>> BaseClass;

    /// the filtered local system matrix for Vanka
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    /// the velocity mass matrix
    typename BaseClass::GlobalMatrixBlockA velo_mass_matrix;

    /// the local interface filter
    typename BaseClass::LocalVeloUnitFilter filter_interface_fbm;

    typedef DT_ DataType;
    typedef IT_ IndexType;
    typedef typename BaseClass::GlobalSystemVector GlobalSystemVector;
    typedef typename BaseClass::LocalVeloVector LocalVeloVector;

    /// FBM mask vectors of velocity and pressure
    std::vector<int> fbm_mask_velo, fbm_mask_pres;

    void assemble_velocity_laplace_matrix(Assembly::DomainAssembler<TrafoType>& dom_asm,
      const SpaceVeloType& space_velo, const DataType nu, bool defo, String cubature = "")
    {
      auto& loc_a = this->matrix_a.local();
      loc_a.format();
      Assembly::Common::LaplaceOperatorBlocked<dim> lapl_op;
      Assembly::Common::DuDvOperatorBlocked<dim> dudv_op;
      if(cubature.empty())
        cubature = "auto-degree:" + stringify(2*SpaceVeloType::local_degree+2);
      if(defo)
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_a, dudv_op, space_velo, cubature, nu);
      else
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_a, lapl_op, space_velo, cubature, nu);
    }

    void clear_velocity_laplace_matrix()
    {
      this->matrix_a.local() = BaseClass::LocalMatrixBlockA();
    }

    void assemble_velocity_mass_matrix(Assembly::DomainAssembler<TrafoType>& dom_asm, const SpaceVeloType& space_velo, String cubature = "")
    {
      this->velo_mass_matrix = this->matrix_a.clone(LAFEM::CloneMode::Weak);
      auto& loc_m = this->velo_mass_matrix.local();
      loc_m.format();

      if(cubature.empty())
        cubature = "auto-degree:" + stringify(2*SpaceVeloType::local_degree+2);
      Assembly::Common::IdentityOperatorBlocked<dim> id_op;
      Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_m, id_op, space_velo, cubature);
    }

    void compile_local_matrix()
    {
      // do we have to allocate the structures?
      if(this->local_matrix_sys.block_a().empty())
      {
        this->local_matrix_sys.block_a() = this->matrix_a.local().clone(LAFEM::CloneMode::Layout);
        this->local_matrix_sys.block_b() = this->matrix_b.local().clone(LAFEM::CloneMode::Layout);
        this->local_matrix_sys.block_d() = this->matrix_d.local().clone(LAFEM::CloneMode::Layout);
      }
      // copy local matrices (and sync 1 if necessary)
      this->matrix_a.convert_to_1(this->local_matrix_sys.block_a(), false);
      this->local_matrix_sys.block_b().copy(this->matrix_b.local());
      this->local_matrix_sys.block_d().copy(this->matrix_d.local());

      // apply velocity unit filters to A and B
      for(const auto& filter : this->get_local_velo_unit_filter_seq())
      {
        filter.second.filter_mat(this->local_matrix_sys.block_a());
        filter.second.filter_offdiag_row_mat(this->local_matrix_sys.block_b());
      }
    }

    void assemble_fbm_filters(Assembly::StokesFBMAssembler<MeshType>& fbm_asm, const SpaceVeloType& space_velo, const SpacePresType& space_pres, bool asm_mask, bool no_scale)
    {
      // assemble velocity and pressure unit filters
      fbm_asm.assemble_inside_filter(this->get_local_velo_unit_filter_seq().find_or_add("fbm"), space_velo);
      fbm_asm.assemble_inside_filter(this->get_local_pres_unit_filter(), space_pres);

      // assemble interface filter for velocity
      fbm_asm.assemble_interface_filter(filter_interface_fbm, space_velo, this->matrix_a, velo_mass_matrix, no_scale);
      //filter_interface_fbm = LocalVeloUnitFilter(space_velo.get_num_dofs());

      // assemble mask vectors on finest level
      if(asm_mask)
      {
        fbm_mask_velo.reserve(space_velo.get_num_dofs());
        for(int d(0); d <= dim; ++d)
        {
          for(auto k : fbm_asm.get_fbm_mask_vector(d))
            fbm_mask_velo.push_back(k);
        }
        fbm_mask_pres = fbm_asm.get_fbm_mask_vector(dim);
      }
    }

    void assemble_pressure_mean_filter(const SpacePresType& space_pres, bool enable_fbm)
    {
      typename BaseClass::GlobalPresVector  vec_glob_v(&this->gate_pres), vec_glob_w(&this->gate_pres);

      // fetch the local vectors
      LocalPresVector<DT_, IT_>& vec_loc_v = vec_glob_v.local();
      LocalPresVector<DT_, IT_>& vec_loc_w = vec_glob_w.local();

      // fetch the frequency vector of the pressure gate
      LocalPresVector<DT_, IT_>& vec_loc_f = this->gate_pres._freqs;

      // assemble the mean filter
      Assembly::MeanFilterAssembler::assemble(vec_loc_v, vec_loc_w, space_pres, "gauss-legendre:2");

      // synchronize the vectors
      vec_glob_v.sync_1();
      vec_glob_w.sync_0();

      // apply pressure unit filter if FBM is enabled
      if(enable_fbm)
      {
        const auto& fil_p = this->get_local_pres_unit_filter();
        fil_p.filter_cor(vec_loc_v);
        fil_p.filter_def(vec_loc_w);
      }

      // create mean filter
      this->get_local_pres_mean_filter() = LocalPresMeanFilter(vec_loc_v.clone(), vec_loc_w.clone(), vec_loc_f.clone(), this->gate_pres.get_comm());
    }

    void apply_fbm_filter_to_rhs(GlobalSystemVector& vec_rhs) const
    {
      this->apply_fbm_filter_to_rhs(vec_rhs.local().first());
    }

    void apply_fbm_filter_to_rhs(LocalVeloVector& vec_rhs_v) const
    {
      this->filter_interface_fbm.filter_def(vec_rhs_v);
    }

    void apply_fbm_filter_to_def(GlobalSystemVector& vec_def_v, const GlobalSystemVector& vec_sol_v, const DataType factor) const
    {
      this->apply_fbm_filter_to_def(vec_def_v.local().first(), vec_sol_v.local().first(), factor);
    }

    void apply_fbm_filter_to_def(LocalVeloVector& vec_def_v, const LocalVeloVector& vec_sol_v, const DataType factor) const
    {
      if(this->filter_interface_fbm.used_elements() == Index(0))
        return;

      auto* vdef = vec_def_v.elements();
      const auto* vsol = vec_sol_v.elements();
      const IndexType* fidx = this->filter_interface_fbm.get_indices();
      const auto* fval = this->filter_interface_fbm.get_values();
      const IndexType* row_ptr = this->velo_mass_matrix.local().row_ptr();
      const IndexType* col_idx = this->velo_mass_matrix.local().col_ind();
      const auto* mval = this->velo_mass_matrix.local().val();

      Index n = this->filter_interface_fbm.used_elements();
      for(Index i(0); i < n; ++i)
      {
        IndexType row = fidx[i];
        vdef[row] = DataType(0);
        for(IndexType j(row_ptr[row]); j < row_ptr[row+1]; ++j)
        {
          vdef[row] += mval[j] * vsol[col_idx[j]];
        }
        vdef[row] *= fval[i][0] * factor;
      }
    }


  }; // class SystemLevelBase


  template<typename SystemLevel_, typename Domain_, typename InflowBounds_, typename SolverDataType_, typename ParamHolder_>
  std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>> new_flow_solver_inner_type(
      Domain_& domain, std::shared_ptr<SystemLevel_>& system_level, const std::vector<InflowBounds_>& inflow_bounds,
      const std::vector<Material<typename SystemLevel_::DataType>>& materials, const PropertyMap* config, const Logger* logger, bool system_assembled, const ParamHolder_& param_holder)
  {
    constexpr bool system_is_solver_level = std::is_same_v<typename SystemLevel_::DataType, SolverDataType_>;
    typedef typename SystemLevel_::DataType BaseDataType;
    typedef SystemLevelBase<SolverDataType_, typename SystemLevel_::IndexType> SolverLevel;
    const String cubature_transfer("gauss-legendre:3");
    // parse this from config
    const String& cubature_matrix_a = param_holder.cubature_a;
    const String& cubature_matrix_m = param_holder.cubature_m;
    const String& cubature_matrix_b = param_holder.cubature_b;
    bool use_q2_fbm = param_holder.assemble_q2_filter;
    bool use_coarse_fbm = true;
    std::size_t min_q2_fbm_level = ~std::size_t(0);
    BaseDataType mesh_unit_scale = param_holder.mesh_unit_scale;

    //first of all, create solver_levels
    std::deque<std::shared_ptr<SolverLevel>> solver_levels;
    const Index num_levels(domain.size_physical());
    solver_levels.resize(num_levels);
    if constexpr(system_is_solver_level)
    {
      solver_levels.at(0) = system_level;
    }
    else
    {
      solver_levels.at(0).reset(new SolverLevel());
      // and always convert from our system level
      solver_levels.at(0)->convert(*system_level);
    }
    for (Index i(1); i < num_levels; ++i)
    {
      solver_levels.at(i).reset(new SolverLevel());
      solver_levels.at(i)->assemble_gates(domain.at(i));
    }
    // assemble transfers
    for(Index i(0); (i < domain.size_physical()) && ((i+1) < domain.size_virtual()); ++i)
    {
      solver_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
      if((i+1) < domain.size_physical())
        solver_levels.at(i)->assemble_transfers_voxel(*solver_levels.at(i+1), domain.at(i), domain.at(i+1), cubature_transfer, true);
      else
        solver_levels.at(i)->assemble_transfers_voxel(domain.at(i), domain.at(i+1), cubature_transfer, true);
    }

    // extra logic for finest level for matrix assembly
    {
      // get domain and system levels
      auto& dom_lvl = *domain.at(0);
      auto& sys_lvl = *solver_levels.at(0);

      // assemble matrix structures
      if(!system_assembled)
      {
        // if nothing is assembled, we have to do the deed now
        sys_lvl.assemble_velo_struct(dom_lvl.space_velo);
        sys_lvl.assemble_pres_struct(dom_lvl.space_pres);
        sys_lvl.assemble_velocity_laplace_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, SolverDataType_(materials.front().get_viscosity_model()->get_data()[0] * 1E3/mesh_unit_scale), true, cubature_matrix_a);
        sys_lvl.assemble_velocity_mass_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, cubature_matrix_m);
        if constexpr(system_is_solver_level)
          sys_lvl.assemble_grad_div_matrices(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b);
        else
          sys_lvl.assemble_grad_div_matrices_high_prec(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b, SystemDataType(1));
      }
      else
      {
        if constexpr(!system_is_solver_level)
        {
          // in any other case, we already assembled the system, so we need to convert the missing members not cloned explicitly
          sys_lvl.velo_mass_matrix.convert(&sys_lvl.gate_velo, &sys_lvl.gate_velo ,system_level->velo_mass_matrix);
          sys_lvl.fbm_mask_velo = system_level->fbm_mask_velo;
          sys_lvl.fbm_mask_pres = system_level->fbm_mask_pres;
        }
      }

      // compile the system matrix
      sys_lvl.compile_system_matrix();

    }
    // assemble matrices
    for(Index i(1); i < num_levels; ++i)
    {
      // get domain and system levels
      auto& dom_lvl = *domain.at(i);
      auto& sys_lvl = *solver_levels.at(i);

      sys_lvl.assemble_velo_struct(dom_lvl.space_velo);
      sys_lvl.assemble_pres_struct(dom_lvl.space_pres);
      sys_lvl.assemble_velocity_laplace_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, SolverDataType_(materials.front().get_viscosity_model()->get_data()[0] * 1E+3/mesh_unit_scale), true, cubature_matrix_a);
      sys_lvl.assemble_velocity_mass_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, cubature_matrix_m);
      if constexpr(system_is_solver_level)
        sys_lvl.assemble_grad_div_matrices(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b);
      else
        sys_lvl.assemble_grad_div_matrices_high_prec(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b, SystemDataType(1));

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
          solver_levels.front()->filter_interface_fbm.convert(system_level->filter_interface_fbm);
          solver_levels.front()->compile_system_filter();
        }
      }

      // the names of the mesh parts on which to assemble
      std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);
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
          auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
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
              XASSERTM(id >= 0 && std::size_t(id) < inflow_bounds.size(), "Id does not map to boundaries");
              Assembly::UnitFilterAssembler<MeshType> loc_unit_asm_inflow;

              // inflow
              loc_unit_asm_inflow.add_mesh_part(*mesh_part);
              if constexpr(system_is_solver_level)
              {
                loc_unit_asm_inflow.assemble(filter_v_inflow, domain.at(i)->space_velo, inflow_bounds.at(id).get_diri_inflow_function(materials, mesh_unit_scale));
              }
              else
              {
                // create temporary filter (of correct datatype)
                typename SystemLevel_::LocalVeloUnitFilter fil_loc_v;
                // directly assembly onto the inflow filter
                loc_unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, inflow_bounds.at(id).get_diri_inflow_function(materials, mesh_unit_scale));
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
                // comm.print("Assembled " + name);
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
          unit_asm_noflow.assemble(filter_v_noflow, domain.at(i)->space_velo);
        }
        // assemble fbm parts
        {
          auto& fbm_asm = *domain.at(i)->fbm_assembler;
          // auto& filter_fbm_p = solver_levels.at(i)->get_local_pres_unit_filter();
          // auto& filter_fbm_v = solver_levels.at(i)->get_local_velo_unit_filter_seq().find_or_add("fbm");
          auto& filter_fbm_int_v = solver_levels.at(i)->filter_interface_fbm;

          // assemble velocity unit filter
          if(use_q2_fbm && ((unsigned short)(i) <= min_q2_fbm_level || use_coarse_fbm))
          {
            // fbm_asm.assemble_inside_filter(filter_fbm_v, domain.at(i)->space_velo);
            // fbm_asm.assemble_inside_filter(filter_fbm_p, domain.at(i)->space_pres);
            // // if(correct_fbm_filter)
            // // {
            // //   //this probably does not work with interface filter...
            // //   correct_filters(*domain.at(i), filter_fbm_v, filter_fbm_p, *cgal_wrapper);
            // // }
            // fbm_asm.assemble_interface_filter(filter_fbm_int_v, domain.at(i)->space_velo, solver_levels.at(i)->matrix_a, solver_levels.at(i)->velo_mass_matrix, true);
            solver_levels.at(i)->assemble_fbm_filters(fbm_asm, domain.at(i)->space_velo, domain.at(i)->space_pres, (i==0u), false);
          }
          else
            filter_fbm_int_v = typename SolverLevel::LocalVeloUnitFilter(domain.at(i)->space_velo.get_num_dofs());

          // // assemble mask vectors on finest level
          // if(i == 0u)
          // {
          //   auto& mask_v = solver_levels.at(i)->fbm_mask_velo;
          //   mask_v.reserve(domain.at(i)->space_velo.get_num_dofs());
          //   for(int d(0); d <= dim; ++d)
          //   {
          //     for(auto k : fbm_asm.get_fbm_mask_vector(d))
          //       mask_v.push_back(k);
          //   }
          //   solver_levels.at(i)->fbm_mask_pres = fbm_asm.get_fbm_mask_vector(dim);
          // }

        }
        // compile system filter
        solver_levels.at(i)->compile_system_filter();
      }
    }

    // system is assembled, now we need to create the components
    auto stokes_solver = std::make_shared<MultigridVankaFBMFlowSolver<SolverLevel>>(solver_levels, domain, config->get_sub_section("solver-params"), logger);
    auto defect_asm = std::make_shared<CarreauSingleMatDefectFlowAssembler<Domain_, SystemLevel_>>(domain, system_level, config, materials.front(), inflow_bounds.front().get_temperature(), mesh_unit_scale);
    defect_asm->set_cubature(cubature_matrix_a);
    auto system_asm = std::make_shared<CarreauSingleMatSystemFlowAssembler<Domain_>>(domain, config, materials.front(), inflow_bounds.front().get_temperature(), mesh_unit_scale);
    system_asm->set_cubature(cubature_matrix_a);

    //rescale interface filter
    const SystemDataType scaling_factor = defect_asm->calc_stiffness_factor();
    {
      // our scaling factor for velocity <-> pressure
      for(auto& sys : stokes_solver->get_systems())
      {
        auto& fbm_int_filter = sys->filter_interface_fbm;
        // get prev value
        auto& _internal_vec = fbm_int_filter.get_filter_vector();
        if(_internal_vec.size() > 0u)
        {
          auto diag_a = sys->matrix_a.create_vector_l();
          auto diag_m = sys->velo_mass_matrix.create_vector_l();
          sys->matrix_a.extract_diag(diag_a, true);
          sys->velo_mass_matrix.extract_diag(diag_m, true);
          _internal_vec.format(diag_a.max_abs_element() / (diag_m.max_abs_element()*SolverDataType_(scaling_factor)));
        }
      }
      if constexpr(!system_is_solver_level)
      {
        auto& sys = defect_asm->system_level;

        auto& fbm_int_filter = sys->filter_interface_fbm;
        // get prev value
        auto& _internal_vec = fbm_int_filter.get_filter_vector();
        if(_internal_vec.size() > 0u)
        {
          auto diag_a = sys->matrix_a.create_vector_l();
          auto diag_m = sys->velo_mass_matrix.create_vector_l();
          sys->matrix_a.extract_diag(diag_a, true);
          sys->velo_mass_matrix.extract_diag(diag_m, true);
          _internal_vec.format(diag_a.max_abs_element() / (diag_m.max_abs_element()*scaling_factor));
        }
      }
      defect_asm->set_stiffness_factor(scaling_factor);
      system_asm->set_stiffness_factor(scaling_factor);
    }

    typedef MultigridVankaFBMFlowSolver<SolverLevel> StS;
    typedef CarreauSingleMatDefectFlowAssembler<Domain_, SystemLevel_> CaD;
    typedef CarreauSingleMatSystemFlowAssembler<Domain_> CaS;

    const auto* sol_config = config->get_sub_section("solver-params");
    if(sol_config && (sol_config->query("nl-solver", "pseudo-ts").compare_no_case("alpine") == 0))
    {
      auto* sys_solver_ptr = new AlPiNeSteadyFlowSolver<StS, CaD, CaS, decltype(system_level->filter_sys)>(stokes_solver, defect_asm, system_asm, system_level->filter_sys, logger);
      sys_solver_ptr->set_scaling_factor(scaling_factor);
      sys_solver_ptr->parse(sol_config);
      return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(sys_solver_ptr);
    }
    else if(sol_config && (sol_config->query("nl-solver", "pseudo-ts").compare_no_case("newton") == 0))
    {
      XABORTM("Not implemented");
      return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(nullptr);
    }
    else if(sol_config && (sol_config->query("nl-solver", "pseudo-ts").compare_no_case("picard") == 0))
    {
      XABORTM("Not implemented");
      return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(nullptr);
    }
    else if(sol_config && (sol_config->query("nl-solver", "pseudo-ts").compare_no_case("pseudo-ts") == 0))
    {
      auto* sys_solver_ptr = new PseudoUnsteadyFlowSolver<StS, CaD, CaS, decltype(system_level->filter_sys)>(stokes_solver, defect_asm, system_asm, system_level->filter_sys, logger);
      sys_solver_ptr->set_scaling_factor(scaling_factor);
      sys_solver_ptr->parse(sol_config);
      return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(sys_solver_ptr);
    }

    XABORTM("Option " + sol_config->query("nl-solver", "pseudo-ts") + String(" not known"));
    return std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>>(nullptr);
  }


  template<typename SystemLevel_, typename Domain_, typename InflowBounds_, typename ParamHolder_>
  std::unique_ptr<Gendie::NonlinearSteadyFlowSolverBase<typename SystemLevel_::GlobalSystemVector>> new_flow_solver(
      Domain_& domain, std::shared_ptr<SystemLevel_>& system_level, const std::vector<InflowBounds_>& inflow_bounds,
      const std::vector<Material<typename SystemLevel_::DataType>>& materials, const PropertyMap* config, const Logger* logger, bool system_assembled, const ParamHolder_& param_holder)
  {
    const auto* solver_config = config->get_sub_section("solver-params");
    auto datatype_entry = solver_config->query("datatype");
    if(datatype_entry.second)
    {
      if(datatype_entry.first.compare_no_case("double") == 0)
      {
        return new_flow_solver_inner_type<SystemLevel_, Domain_, InflowBounds_, double>(domain, system_level, inflow_bounds, materials, config, logger, system_assembled, param_holder);
      }
      else if(datatype_entry.first.compare_no_case("float") == 0)
      {
        return new_flow_solver_inner_type<SystemLevel_, Domain_, InflowBounds_, float>(domain, system_level, inflow_bounds, materials, config, logger, system_assembled, param_holder);
      }
      else if(datatype_entry.first.compare_no_case("same") == 0)
      {
        return new_flow_solver_inner_type<SystemLevel_, Domain_, InflowBounds_, typename SystemLevel_::DataType>(domain, system_level, inflow_bounds, materials, config, logger, system_assembled, param_holder);
      }
    }
    return new_flow_solver_inner_type<SystemLevel_, Domain_, InflowBounds_, typename SystemLevel_::DataType>(domain, system_level, inflow_bounds, materials, config, logger, system_assembled, param_holder);

  }
}
