#include <kernel/runtime.hpp>
#include <exception>
#include <iostream>

#include "gendie_common.hpp"


#include <applications/gendie/inflow_boundary.hpp>
#include <applications/gendie/materials.hpp>
#include <applications/gendie/logger.hpp>
#include <applications/gendie/system_assembler.hpp>
#include <applications/gendie/defect_assembler.hpp>
#include <applications/gendie/steady_stokes_solver.hpp>
#include <applications/gendie/steady_flow_solver.hpp>
#include <applications/gendie/cgal_meshpart_helper.hpp>
#include <applications/gendie/velo_analyser.hpp>
#include <kernel/assembly/discrete_projector.hpp>

#ifdef FEAT_HAVE_OMP
#include "omp.h"
#else
static inline int omp_get_max_threads()
{
  return 1;
}
#endif



namespace Gendie
{

  /// our domain level
  typedef DomainLevelBase DomainLevelBackground;
  typedef Control::Domain::VoxelDomainLevelWrapper<DomainLevelBackground> DomainLevel;

  template<typename DT_, typename IT_> using SystemLevel = SystemLevelBase<DT_, IT_>;

  /// our partidomaincontrol
  template<typename DomainLevel_>
  using PartitionControl = Control::Domain::VoxelDomainControl<DomainLevel_>;

  typedef Gendie::InflowBoundary<MeshType> InflowBound;

  static constexpr std::size_t pad_len = 30u;
  static constexpr char pad_char = '.';

  // class parsing and holding parameters used by our functions
  struct ParameterHolder
  {
    String surface_meshfile;
    String base_meshfile_names;
    String cubature_a = "gauss-legendre:5";
    String cubature_b = "gauss-legendre:3";
    String cubature_m = "gauss-legendre:3";
    String cubature_postproc = "gauss-legendre:5";
    String level_string;
    String voxel_map_file;
    String min_gap;
    String vtk_filename;
    String vtk_name;
    Tiny::Vector<SystemIndexType, dim> mesh_bb_num_cells;
    Tiny::Vector<SystemDataType, dim> mesh_bb_min, mesh_bb_max;
    SystemDataType resolution = SystemDataType(1);
    SystemDataType mesh_unit_scale = SystemDataType(1E+3);
    SystemDataType velo_scaling_factor = SystemDataType(1)/SystemDataType(1E+3);
    SystemDataType visc_scaling_factor = SystemDataType(1);
    SystemDataType pres_scaling_factor = SystemDataType(1)/SystemDataType(1E+5);




    bool use_base_mesh = false;
    bool use_voxel_map = false;
    bool meshpart_by_cgal_wrapper = false;
    bool assemble_q2_filter = false;
    bool want_vtk = false;
    bool refine_vtk = false;
    bool use_true_vtk_name = false;
    bool vtk_debug = false;
    bool measure_cgal_size = true;

    ParameterHolder() = default;

    ParameterHolder(const PropertyMap* config, const Logger* logger)
    {
      {
        if(const PropertyMap* problem_config = config->get_sub_section("problem_params"))
        {
          // parse off file
          std::pair<String, bool> off_string = problem_config->query("surface_mesh_filename");
          XASSERTM(off_string.second && !off_string.first.empty(), "Could not parse surface mesh file");
          //only first off file is parsed for now
          surface_meshfile = off_string.first.split_by_whitespaces().front();
        }
        else
        {
          XABORTM("Could not find problem_params section");
        }

      }

      // get the problem setup section
      if(auto* base_config = config->get_sub_section("problem-setup"))
      {
        {
          std::pair<String, bool> lvl_string_ = base_config->query("level");
          if(!lvl_string_.second)
          {
            logger->print("You must provide a level field inside the problem-setup section!", error);
            Runtime::abort();
          }
          level_string = lvl_string_.first;
          // TODO: somewhere at this point, we would need to set the base splitte if required

        }

        if(Gendie::check_for_config_option(base_config->query("use-base-mesh")))
        {
          resolution = SystemDataType(1);
          use_base_mesh = true;
          std::pair<String, bool> base_mesh_string = base_config->query("base-mesh");
          XASSERTM(base_mesh_string.second, "You must provide a mesh field inside the problem-setup section!");
          base_meshfile_names = base_mesh_string.first;

          XASSERTM(!base_meshfile_names.empty(), "Could not parse base mesh file name");
        }
        else
        {
          use_base_mesh = false;
          // parse base mesh boundingbox
          std::pair<String, bool> mesh_string = base_config->query("background-mesh");
          XASSERTM(mesh_string.second, "You must provide parameters for the background mesh!\n"
                            "Xmin Ymin [Zmin] Xmax Ymax [Zmax] n_x n_y [n_z]");
          std::deque<String> tmp_msh_string = mesh_string.first.split_by_whitespaces();
          if(tmp_msh_string.size() != Index(3*dim))
          {
            logger->print(String("Number of parameters for background-mesh does not fit to dimension!\n")
                      + String("Number of provided arguments: " + stringify(tmp_msh_string.size()) + " Number of expected arguments: "
                      + stringify(Index(3*dim))), error);
          }
          //first parse the min parameters
          for(int i = 0; i < dim; ++i)
          {
            if(!tmp_msh_string.at(i).parse(mesh_bb_min[i]))
            {
              logger->print("Could not parse " + stringify(i) + "th Parameter of BB min parameter!", error);
              Runtime::abort();
            }
          }
          for(int i = 0; i < dim; ++i)
          {
            if(!tmp_msh_string.at(i+dim).parse(mesh_bb_max[i]))
            {
              logger->print("Could not parse " + stringify(i) + "th Parameter of BB max parameter!", error);
              Runtime::abort();
            }
          }
          for(int i = 0; i < dim; ++i)
          {
            if(!tmp_msh_string.at(i+2*dim).parse(mesh_bb_num_cells[i]))
            {
              logger->print("Could not parse " + stringify(i) + "th Parameter of BB num cells parameter!", error);
              Runtime::abort();
            }
          }
          // calculate an approximate resolution
          resolution = Math::Limits<SystemDataType>::max();
          for(int i = 0; i < dim; ++i)
          {
            resolution = Math::min(resolution, Math::abs((mesh_bb_max[i]-mesh_bb_min[i])/ SystemDataType(mesh_bb_num_cells[i])));
          }

        }
        // refine resolution parameter depending on level
        resolution *= Math::pow(SystemDataType(2), -SystemDataType(std::stoi(level_string.split_by_whitespaces().front())));

        // create voxel map
        if(Gendie::check_for_config_option(base_config->query("voxel-map")))
        {
          use_voxel_map = true;
          voxel_map_file = base_config->query("voxel-map").first;
        }
        else
        {
          use_voxel_map = false;
          // parse resolution, if available
          const auto* spec_problem_def = config->get_sub_section("special-problem-definition");
          if(!spec_problem_def)
          {
            logger->print("No special problem definition section found", warning);
          }
          if(spec_problem_def && Gendie::check_for_config_option(spec_problem_def->query("min-gap")))
          {
            min_gap = spec_problem_def->query("min-gap").first;
            auto parse_ok = spec_problem_def->parse_entry("min-gap", resolution);
            XASSERTM(parse_ok, "Could not parse min-gap");
          }
        }
        if(const auto* spec_problem_def = config->get_sub_section("special-problem-definition"))
          assemble_q2_filter = Gendie::check_for_config_option(spec_problem_def->query("use-q2-fbm"));

        if(!base_config->parse_entry("mesh-unit-factor", mesh_unit_scale))
        {
          logger->print("Could not parse mesh-unit-factor", error);
          Runtime::abort();
        }

        velo_scaling_factor = SystemDataType(1)/mesh_unit_scale;
        visc_scaling_factor = SystemDataType(1E-3)*mesh_unit_scale;
        pres_scaling_factor = SystemDataType(1E-8)*mesh_unit_scale;

        meshpart_by_cgal_wrapper = Gendie::check_for_config_option(base_config->query("meshpart-by-cgal-wrapper"));
      }


      if(const auto* postproc_config = config->get_sub_section("postproc-params"))
      {
        measure_cgal_size = !Gendie::check_for_config_option(postproc_config->query("skip_cgal_meaure"));
        refine_vtk = Gendie::check_for_config_option(postproc_config->query("refine-vtk"));
        want_vtk = Gendie::check_for_config_option(postproc_config->query("vtk"));
        if(!want_vtk)
        {
          logger->print("Vtk output not requested!", warning);
          return;
        }
        std::pair<FEAT::String, bool> vtk_filename_query = postproc_config->query("vtk");
        if(want_vtk && vtk_filename_query.first.compare_no_case("true") != 0)
        {
          vtk_filename = vtk_filename_query.first;
        }
        use_true_vtk_name = Gendie::check_for_config_option(postproc_config->query("vtk-use-true-name"));

        vtk_debug = Gendie::check_for_config_option(postproc_config->query("vtk-debug-data"));
        // build VTK name
        vtk_name = vtk_filename;
        if(vtk_filename.empty())
        {
          vtk_name = "gendie_material";
          vtk_name += "-lvl" + level_string.split_by_whitespaces().front();
          vtk_name += "-n" + stringify(logger->comm.size());
        }
        else if(!use_true_vtk_name)
        {
          vtk_name += "_level_" + level_string.split_by_whitespaces().front();
        }
      }
      else
      {
        logger->print("No postproc-params section found!", warning);
      }

    }

    String format_string() const
    {
      String output;
      output += String("Surface Mesh File").pad_back(pad_len, pad_char) + ": " + surface_meshfile + "\n";
      output += String("Use Base Mesh").pad_back(pad_len, pad_char) + ": " + (use_base_mesh ? String("yes"):String("no")) + "\n";
      if(use_base_mesh)
      {
        output += String("Base Mesh Files").pad_back(pad_len, pad_char) + ": " +base_meshfile_names + "\n";
      }
      else
      {
        output += String("Background Mesh Num Cells").pad_back(pad_len, pad_char) + ": " + "[" + stringify(mesh_bb_num_cells[0]) + ", " + stringify(mesh_bb_num_cells[1]) + ", " + stringify(mesh_bb_num_cells[2]) + "]\n";
      }
      output += String("Background Mesh BB").pad_back(pad_len, pad_char) + ": " + "[" + stringify(mesh_bb_min[0]) + ", " + stringify(mesh_bb_max[0])
                      + "], [" + stringify(mesh_bb_min[1]) + ", " + stringify(mesh_bb_max[1]) + "], [" + stringify(mesh_bb_min[2]) + ", " + stringify(mesh_bb_max[2]) + "]\n";
      output += String("MinGap").pad_back(pad_len, pad_char) + ": " + min_gap + "\n";
      output += String("Use Voxelmap").pad_back(pad_len, pad_char) + ": " + (use_voxel_map ? String("yes"):String("no")) + "\n";
      if(use_voxel_map)
      {
        output += String("Voxelmap File").pad_back(pad_len, pad_char) + ": " + voxel_map_file + "\n";
      }
      else
      {
        output += String("Internal Voxel Map Resolution").pad_back(pad_len, pad_char) + ": " + stringify_fp_fix(resolution, 5) + "\n";
      }
      output += String("Use Voxel Map to Create Filters").pad_back(pad_len, pad_char) + ": " + (meshpart_by_cgal_wrapper ? String("no") : String("yes")) + "\n";
      output += String("Cubature Matrix a").pad_back(pad_len, pad_char) + ": " + cubature_a + "\n";
      output += String("Cubature Matrix b").pad_back(pad_len, pad_char) + ": " + cubature_b + "\n";
      output += String("Cubature Matrix m").pad_back(pad_len, pad_char) + ": " + cubature_m + "\n";
      output += String("Cubature Postproc").pad_back(pad_len, pad_char) + ": " + cubature_postproc + "\n";
      output += String("Use q2 fbm").pad_back(pad_len, pad_char) + ": " + (assemble_q2_filter ? String("yes"):String("no")) + "\n";
      output += String("Mesh Unit Scale").pad_back(pad_len, pad_char) + ": " + stringify_fp_sci(mesh_unit_scale, 2) + "\n";
      output += String("Write VTK").pad_back(pad_len, pad_char) + ": " + (want_vtk ? String("yes"):String("no")) + "\n";
      output += String("Vtk file path").pad_back(pad_len, pad_char) + ": " + vtk_filename + "\n";
      output += String("Vtk name").pad_back(pad_len, pad_char) + ": " + vtk_name + "\n";
      output += String("Refine VTK").pad_back(pad_len, pad_char) + ": " + (refine_vtk ? String("yes"):String("no")) + "\n";
      output += String("Debug VTK").pad_back(pad_len, pad_char) + ": " + (vtk_debug ? String("yes"):String("no")) + "\n";
      output += String("Velo Scaling Factor").pad_back(pad_len, pad_char) + ": " + stringify_fp_sci(velo_scaling_factor, 2) + "\n";
      output += String("Visco Scaling Factor").pad_back(pad_len, pad_char) + ": " + stringify_fp_sci(visc_scaling_factor, 2) + "\n";
      output += String("Pressure Scaling Factor").pad_back(pad_len, pad_char) + ": " + stringify_fp_sci(pres_scaling_factor, 2) + "\n";

      return output;

    }


  }; // class ParameterHolder

  template<typename ParamHolder_>
  void create_voxel_domain(PartitionControl<DomainLevel>& domain, const SimpleArgParser& /*args*/, const PropertyMap* config, ParamHolder_& param_holder, const Logger* logger)
  {
    logger->print("Going into create_voxel_domain", info);
    TimeStamp create_voxel_domain;
    {
      const auto* parti_config = config->get_sub_section("partition");
      bool success = true;
      if(parti_config)
        success = domain.parse_property_map(*parti_config);
      if(!success)
      {
        logger->print("ERROR: Parsing of partition segment went wrong!");
        Runtime::abort();
      }
    }

    // setup domain
    {
      domain.set_desired_levels(param_holder.level_string.split_by_whitespaces());

      // TODO: somewhere at this point, we would need to set the base splitte if required

    }

    if(param_holder.use_base_mesh)
    {
      std::deque<String> mesh_file_names = param_holder.base_meshfile_names.split_by_whitespaces();

      XASSERTM(!mesh_file_names.empty(), "Could not parse base mesh file name");

      // set base mesh for domain and set bounding boxes
      // Our mesh file reader
      Geometry::MeshFileReader mesh_reader;
      mesh_reader.add_mesh_files(domain.comm(), mesh_file_names);
      mesh_reader.read_root_markup();
      String mesh_file_type = mesh_reader.get_meshtype_string();
      XASSERTM(mesh_file_type == "conformal:hypercube:3:3", "Meshfile type has to be of conformal:hypercube:3:3");
      domain.set_base_mesh(mesh_reader.parse(domain.get_atlas()));
      param_holder.mesh_bb_min = domain.get_bounding_box_min();
      param_holder.mesh_bb_max = domain.get_bounding_box_max();
    }
    else
    {
      // so now we can create our background mesh
      domain.create_base_mesh_3d(param_holder.mesh_bb_num_cells[0], param_holder.mesh_bb_num_cells[1], param_holder.mesh_bb_num_cells[2], param_holder.mesh_bb_min[0],
                                param_holder.mesh_bb_max[0], param_holder.mesh_bb_min[1], param_holder.mesh_bb_max[1], param_holder.mesh_bb_min[2], param_holder.mesh_bb_max[2]);
    }

    // create slag mask
    {
      // create voxel map
      if(param_holder.use_voxel_map)
      {
        domain.read_voxel_map(param_holder.voxel_map_file);
      }
      else
      {
        domain.create_voxel_map_from_off(param_holder.surface_meshfile, false, param_holder.resolution);
      }
    }

    if(!param_holder.meshpart_by_cgal_wrapper)
    {
      domain.keep_voxel_map();
    }

    // now we can create our hirarchy
    domain.create_hierarchy();
    logger->print("Printing chosen parti info", info);

    if(domain.front_layer().comm().rank() == 0)
    {
      logger->print(domain.get_chosen_parti_info(), info);
      logger->flush_print();
    }

  }

  void parse_materials(std::vector<Material<SystemDataType>>& materials, const PropertyMap* mat_config)
  {
    XASSERTM(mat_config, "No material config list found");
    std::size_t i = 0u;
    while(true)
    {
      const auto* cur_mat = mat_config->query_section(stringify(i));
      if(cur_mat == nullptr)
      {
        XASSERTM(i > 0u, "Did not querry any materials!");
        break;
      }
      materials.emplace_back(cur_mat);
      ++i;
    }
  }

  void parse_inflow_boundaries(std::vector<InflowBound>& inflow_boundaries, const PropertyMap* inflow_list_config, const std::vector<Material<SystemDataType>>& /*materials*/,
                              const PartitionControl<DomainLevel>& domain)
  {
    XASSERTM(inflow_list_config, "Inflow config not found");
    // now, create our boundaries, for which we need our boundingbox
    Tiny::Vector<SystemDataType, dim> mesh_bb_min, mesh_bb_max;
    // get boundingboxes from our voxel domain
    mesh_bb_min = domain.get_bounding_box_min();
    mesh_bb_max = domain.get_bounding_box_max();

    for(auto it = inflow_list_config->begin_section(); it != inflow_list_config->end_section(); ++it)
    {
      int id = -1;
      it->first.parse(id);
      inflow_boundaries.emplace_back(it->second.get(), mesh_bb_min, mesh_bb_max, id);
    }
  }


  template<typename ParamHolder_>
  void create_boundaries(PartitionControl<DomainLevel>& domain, const std::vector<InflowBound>& inflow_boundaries, const ParamHolder_& param_holder, const Logger* logger)
  {
    // now, create our boundaries, for which we need our boundingbox
    Tiny::Vector<SystemDataType, dim> mesh_bb_min, mesh_bb_max;
    // get boundingboxes from our voxel domain
    mesh_bb_min = domain.get_bounding_box_min();
    mesh_bb_max = domain.get_bounding_box_max();

    const SystemDataType tol_adp = Math::pow(SystemDataType(2), -SystemDataType(domain.get_desired_level_max()));

    const SystemDataType z_max = mesh_bb_max[2]-SystemDataType(0.001) * tol_adp;
    logger->print(String("Zmax Outflow").pad_back(pad_len, pad_char) + ": " + stringify_fp_fix(z_max, 7));
    auto hit_r = [z_max](auto p) {return p[2] > z_max;};

    // create boundary meshparts on each level
    for(Index i(0); i < domain.size_physical(); ++i)
    {
      Geometry::RootMeshNode<MeshType>& mesh_node = *domain.at(i)->get_mesh_node();
      if(param_holder.use_base_mesh)
        mesh_node.remove_all_mesh_parts();
      // create mesh-part from chart by using a hit-test factory
      std::vector<std::pair<String, std::unique_ptr<MeshPartType>>> bnd_inflows;
      bnd_inflows.reserve(inflow_boundaries.size());
      for(const auto& b : inflow_boundaries)
      {
        //TODO: tolerance adaptive to mesh size?
        bnd_inflows.push_back(std::make_pair(b.get_bnd_name(), b.create_mesh_part(domain.at(i)->get_mesh(), SystemDataType(1E-2) * tol_adp)));
      }
      Geometry::HitTestFactory<decltype(hit_r), MeshType> factory_outflow(hit_r, domain.at(i)->get_mesh());

      // auto bnd_inflow = factory_inflow.make_unique();
      auto bnd_outflow = factory_outflow.make_unique();

      Geometry::GlobalMaskedBoundaryFactory<MeshType> boundary_factory(*mesh_node.get_mesh());
      // boundary_factory.add_mask_meshpart(*bnd_inflow);
      for(auto& [a,b] : bnd_inflows)
      {
        boundary_factory.add_mask_meshpart(*b);
      }

      boundary_factory.add_mask_meshpart(*bnd_outflow);
      for(const auto& v : mesh_node.get_halo_map())
        boundary_factory.add_halo(v.first,*v.second);
      boundary_factory.compile(domain.at(i).layer().comm());
      mesh_node.add_mesh_part("bnd:walls", boundary_factory.make_unique());
      for(auto& [a,b] : bnd_inflows)
      {
        mesh_node.add_mesh_part(a, std::move(b));
      }
      mesh_node.add_mesh_part("bnd:outflow", std::move(bnd_outflow));
    }
  }

  template<typename SystemLevel_, typename ParamHolder_>
  void  assemble_filters(SystemLevel_& system_level, const std::deque<String>& part_names, const DomainLevel& cur_dom, const std::vector<InflowBound>& inflow_boundaries, const std::vector<Material<SystemDataType>>& materials,
                            const ParamHolder_& param_holder, bool assemble_mask = true)
  {
    auto& filter_v_noflow = system_level.get_local_velo_unit_filter_seq().find_or_add("noflow");
    auto& filter_v_inflow = system_level.get_local_velo_unit_filter_seq().find_or_add("inflow");

    Assembly::UnitFilterAssembler<MeshType> unit_asm_noflow;

    // loop over all boundary parts except for the right one, which is outflow
    for(const auto& name : part_names)
    {
      // skip non-boundary mesh-parts
      if((!name.starts_with("bnd:")))
        continue;

      // try to fetch the corresponding mesh part node
      auto* mesh_part_node = cur_dom.get_mesh_node()->find_mesh_part_node(name);
      XASSERT(mesh_part_node != nullptr);

      auto* mesh_part = mesh_part_node->get_mesh();
      String name_lower = name.lower();
      if (mesh_part != nullptr)
      {
        // we require the inflow to be of type "inflow_<type>_id"
        // TODO: we should refactor this by saving the meshparts directly with the inflow boundary
        if((name_lower.replace_all(String("inflow"), String("")) > 0u))
        {
          auto name_deq = name.split_by_charset("_");
          if(name_deq.size() < 3)
          {
            XABORTM("Can not parse string: " + name);
          }
          int id = -1;
          name_deq.back().parse(id);
          XASSERTM(id >= 0 && std::size_t(id) < inflow_boundaries.size(), "Id does not map to boundaries");
          Assembly::UnitFilterAssembler<MeshType> loc_unit_asm_inflow;

          // inflow
          loc_unit_asm_inflow.add_mesh_part(*mesh_part);
          // directly assembly onto the inflow filter
          loc_unit_asm_inflow.assemble(filter_v_inflow, cur_dom.space_velo, inflow_boundaries.at(id).get_diri_inflow_function(materials, param_holder.mesh_unit_scale));
        }
        else if((name != "bnd:n") && (!name.starts_with("bnd:out")) && ((name_lower.replace_all(String("outflow"), String("")) == 0u)))
        {
          // outflow
          unit_asm_noflow.add_mesh_part(*mesh_part);
        }
      }
    }

    unit_asm_noflow.assemble(filter_v_noflow, cur_dom.space_velo);
    // assemble fbm filter
    {
      auto& fbm_asm = *cur_dom.fbm_assembler;
      // auto& filter_fbm_p = system_level.get_local_pres_unit_filter();
      // auto& filter_fbm_v = system_level.get_local_velo_unit_filter_seq().find_or_add("fbm");
      auto& filter_fbm_int_v = system_level.filter_interface_fbm;

      // assemble velocity unit filter
      if(param_holder.assemble_q2_filter)
      {
        // fbm_asm.assemble_inside_filter(filter_fbm_v, cur_dom.space_velo);
        // fbm_asm.assemble_inside_filter(filter_fbm_p, cur_dom.space_pres);
        // if(correct_fbm_filter)
        // {
        //   //this probably does not work with interface filter...
        //   correct_filters(*domain.at(i), filter_fbm_v, filter_fbm_p, *cgal_wrapper);
        // }
        // we assume matrix_a and velo matrix to be assembled at this point, but they can be freed after
        // fbm_asm.assemble_interface_filter(filter_fbm_int_v, cur_dom.space_velo, system_level.matrix_a, system_level.velo_mass_matrix, true);
        // fbm_asm.assemble_interface_filter(filter_fbm_int_v, cur_dom.space_velo, system_level.matrix_a, system_level.velo_mass_matrix, false);
        system_level.assemble_fbm_filters(fbm_asm, cur_dom.space_velo, cur_dom.space_pres, assemble_mask, false);
      }
      else
        filter_fbm_int_v = typename SystemLevel_::LocalVeloUnitFilter(cur_dom.space_velo.get_num_dofs());

      // assemble mask vectors on finest level
      // TODO: Only required for debugging purposes
      // if(assemble_mask)
      // {
      //   auto& mask_v = system_level.fbm_mask_velo;
      //   mask_v.reserve(cur_dom.space_velo.get_num_dofs());
      //   for(int d(0); d <= dim; ++d)
      //   {
      //     for(auto k : fbm_asm.get_fbm_mask_vector(d))
      //       mask_v.push_back(k);
      //   }
      //   system_level.fbm_mask_pres = fbm_asm.get_fbm_mask_vector(dim);
      // }

    }

    // compile system filter
    system_level.compile_system_filter();
  }

  template<typename SystemLevel_, typename DomainLevel_, typename VeloInfo_, typename ParamHolder_>
  void write_vtk(const typename SystemLevel_::GlobalSystemVector& vec_sol, const typename SystemLevel_::GlobalSystemVector& vec_rhs, const DomainLevel_& domain, SystemLevel_& system_level, const VeloInfo_& velo_info, const ParamHolder_& param_holder, const Logger* logger)
  {
    logger->print("Writing VTK output to '" + param_holder.vtk_name + ".pvtu'", info);

    std::unique_ptr<MeshType> refined_mesh;
    // get the mesh for the VTK output
    const MeshType* mesh = &domain.front()->get_mesh();
    if(param_holder.refine_vtk)
    {
      logger->print("Write out refined mesh", info);
      Geometry::StandardRefinery<MeshType> refinery(domain.front()->get_mesh());
      refined_mesh = refinery.make_unique();
      mesh = refined_mesh.get();
    }

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(*mesh);

    // write velocity
    exporter.add_vertex_vector("Velocity [m/s]", vec_sol.local().template at<0>(), param_holder.velo_scaling_factor);
    if(param_holder.vtk_debug)
    {
      auto vec_copyt = vec_sol.clone(LAFEM::CloneMode::Deep);
      vec_copyt.format(-1.5);
      system_level.filter_sys.filter_sol(vec_copyt);
      exporter.add_vertex_vector("Velocity testdebug", vec_copyt.local().template at<0>());
    }
    if(param_holder.vtk_debug)
      exporter.add_vertex_vector("Debug: rhs_v", vec_rhs.local().template at<0>());
    if(velo_info.vec_viscosity.size() == vec_sol.local().template at<0>().size())
    {
      exporter.add_vertex_scalar("Viscosity [Pa s]", velo_info.vec_viscosity.get_elements().front(), param_holder.visc_scaling_factor);
      exporter.add_vertex_scalar("NormShearRate [1/s]", velo_info.vec_norm_shear.get_elements().front());
    }


    if(!param_holder.surface_meshfile.empty())
    {
      // parse off file
      //only first off file is parsed for now
      const String& surface_mesh_file = param_holder.surface_meshfile;
      logger->print(String("CGAL mesh file for vtk").pad_back(pad_len, pad_char) + ": " + surface_mesh_file);
      //get filestream
      LAFEM::DenseVector<SystemDataType, Index> vec_dist(vec_sol.local().template at<0>().size());
      vec_dist.format();
      // create a new stream
      auto stream_cgal = std::make_shared<std::stringstream>();

      Geometry::CGALFileMode cgal_filemode;
      if(surface_mesh_file.ends_with(".off"))
        cgal_filemode = Geometry::CGALFileMode::fm_off;
      else if(surface_mesh_file.ends_with(".obj"))
        cgal_filemode = Geometry::CGALFileMode::fm_obj;
      else
        XABORTM("No valid file extension " + surface_mesh_file.split_by_charset(".").back());

      // read the stream
      DistFileIO::read_common(*stream_cgal, surface_mesh_file, domain.comm());

      //create cgal function
      Analytic::Distance::CGALDistanceFunction<SystemDataType> dist_func(*stream_cgal, cgal_filemode);

      //and now interpolate
      Assembly::Interpolator::project(vec_dist, dist_func, domain.front()->space_velo);
      // double metre_scaling_factor{1}; //for now hradcoded here...
      exporter.add_vertex_scalar("Shell", vec_dist.get_elements().front(), param_holder.velo_scaling_factor);
    }
    else
    {
      XABORTM("Could not find problem_params section");
    }

    // project pressure
    LAFEM::DenseVector<SystemDataType, Index> vtx_p, vtx_q;
    if(refined_mesh)
    {
      Assembly::DiscreteCellProjector::project_refined(vtx_p, vec_sol.local().template at<1>(), domain.front()->space_pres);
      Assembly::DiscreteCellProjector::project_refined(vtx_q, vec_rhs.local().template at<1>(), domain.front()->space_pres);
    }
    else
    {
      Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), domain.front()->space_pres, "gauss-legendre:2");
      Assembly::DiscreteCellProjector::project(vtx_q, vec_rhs.local().template at<1>(), domain.front()->space_pres, "gauss-legendre:2");
    }

    // write pressure
    exporter.add_cell_scalar("Pressure [bar]", vtx_p.elements(), param_holder.pres_scaling_factor);
    if(param_holder.vtk_debug)
      exporter.add_cell_scalar("Debug: rhs_p", vtx_q.elements());

    const bool use_q2_fbm = param_holder.assemble_q2_filter;

    // write FBM masks if
    if(refined_mesh)
    {
      if(param_holder.vtk_debug && use_q2_fbm)
        exporter.add_vertex_scalar("Debug: fbm_mask_v", system_level.fbm_mask_velo.data());

      if(param_holder.vtk_debug && use_q2_fbm)
      {
        const auto& filter_vec_velo = system_level.get_local_velo_unit_filter_seq().find_or_add("fbm").get_filter_vector();
        LAFEM::DenseVectorBlocked<SystemDataType, SystemIndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        for(Index k = 0; k < filter_vec_velo.used_elements(); ++k)
        {
          filter_vec(filter_vec_velo.indices()[k], filter_vec_velo.elements()[k]);
        }

        exporter.add_vertex_vector("Debug: fbm_filter_v", filter_vec);
      }

      if(param_holder.vtk_debug)
      {
        const auto& filter_vec_velo = system_level.get_local_velo_unit_filter_seq().find_or_add("noflow").get_filter_vector();
        LAFEM::DenseVectorBlocked<SystemDataType, SystemIndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        for(Index k = 0; k < filter_vec_velo.used_elements(); ++k)
        {
          filter_vec(filter_vec_velo.indices()[k], filter_vec_velo.elements()[k]);
        }

        exporter.add_vertex_vector("Debug: noflow_filter_v", filter_vec);
      }

      if(param_holder.vtk_debug && (system_level.filter_interface_fbm.get_filter_vector().used_elements() > 0))
      {
        const auto& filter_vec_velo = system_level.filter_interface_fbm.get_filter_vector();
        LAFEM::DenseVectorBlocked<SystemDataType, SystemIndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        for(Index k = 0; k < filter_vec_velo.used_elements(); ++k)
        {
          filter_vec(filter_vec_velo.indices()[k], filter_vec_velo.elements()[k]);
        }

        exporter.add_vertex_vector("Debug: interface_filter_v", filter_vec);
      }

      if(param_holder.vtk_debug)
      {
        const auto& filter_vec_velo = system_level.get_local_velo_unit_filter_seq().find_or_add("inflow").get_filter_vector();
        LAFEM::DenseVectorBlocked<SystemDataType, SystemIndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        for(Index k = 0; k < filter_vec_velo.used_elements(); ++k)
        {
          filter_vec(filter_vec_velo.indices()[k], filter_vec_velo.elements()[k]);
        }

        exporter.add_vertex_vector("Debug: inflow_filter_v", filter_vec);
      }

      const std::vector<int>& mask_p = domain.front()->fbm_assembler->get_fbm_mask_vector(dim);
      const int nc = 1 << dim;
      std::vector<int> mask_p_ref;
      mask_p_ref.reserve(mask_p.size() * std::size_t(nc));
      for(int i : mask_p)
        for(int k(0); k < nc; ++k)
          mask_p_ref.push_back(i);
      if(param_holder.vtk_debug)
        exporter.add_cell_scalar("Debug: fbm_mask_p", mask_p_ref.data());
    }
    else
    {
      if(param_holder.vtk_debug)
      {
        exporter.add_vertex_scalar("Debug: fbm_mask_v", domain.front()->fbm_assembler->get_fbm_mask_vector(0).data());
        exporter.add_cell_scalar("Debug: fbm_mask_p", domain.front()->fbm_assembler->get_fbm_mask_vector(dim).data());
      }
    }

    // finally, write the VTK file
    exporter.write(param_holder.vtk_name, domain.comm());

  }


  void main(int argc, char* argv[])
  {
    // we need a few timer
    FEAT::StopWatch watch_total, watch_create_domain, watch_create_system, watch_write_out, watch_create_filter, watch_create_solver, watch_analyze_sol;
    watch_total.start();
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());
    // setup our logging object
    print_level comm_out_level = print_level::info;
    Logger logger(comm, comm_out_level, print_level::off);
    logger.print("-------------------Welcome to FEAT3 Steady Gendie Simulator tool--------------------", info);

    // create arg parser
    SimpleArgParser args(argc, argv);
    args.support("config-file", "A list of config files in the ini file format.\n"
                                    "See the template ini file for details.");
    args.support("problem-file", "A list of additonal config files in the ini file format.");


    // our Propertymap
    std::unique_ptr<PropertyMap> config;

    // fill our propertymap
    if(auto parsed_file = args.query("config-file"))
    {
      // temp config
      PropertyMap config_in;
      XASSERTM(parsed_file->second.size() > std::size_t(0), "No config file given");
      logger.print(String("Config Files").pad_back(pad_len, pad_char) + ": " + std::accumulate(parsed_file->second.begin(), parsed_file->second.end(), String(""), [](const String& a, const String& b){return a + ", " + b;}));
      // read all files given and merge them together, last files taking prescidence
      for(auto& str : parsed_file->second)
      {
        PropertyMap temp_config;
        temp_config.read(comm, str);
        config_in.merge(temp_config, true);
      }
      if(auto parsed_problem_config = args.query("problem-file"))
      {
        logger.print(String("Problem Config file").pad_back(pad_len, pad_char) + ": " + std::accumulate(parsed_problem_config->second.begin(), parsed_problem_config->second.end(), String(""), [](const String& a, const String& b){return a + ", " + b;}));
        // also read inproblem config file and merge them in, last files taking prescidence
        for(auto& str : parsed_problem_config->second)
        {
          PropertyMap temp_config;
          temp_config.read(comm, str);
          config_in.merge(temp_config, true);
        }

      }
      //treeify our structure
      config = config_in.treeify_structures();
    }
    else
    {
      XABORTM("Did not provide a config file");
    }



    // start by constructing our domain and systemlevels
    logger.print("Number of Processes: " + stringify(comm.size()), info);
    logger.print("Number of OMP Threads: " + stringify(omp_get_max_threads()), info);
    logger.print("System DataType: " + String(fp_typename), info);
    logger.print("System IndexType: " + String(ix_typename), info);

    logger.print("Parsing Parameters...", info);
    // wrapper class for our parameters
    ParameterHolder param_holder(config.get(), &logger);

    if(param_holder.measure_cgal_size)
    {
      std::unique_ptr<Geometry::CGALWrapper<SystemDataType>> cgal = Gendie::create_cgal_wrapper(param_holder.surface_meshfile, comm);
      if(cgal)
      {
        std::size_t bytes_sum, bytes_min, bytes_max;
        bytes_sum = bytes_min = bytes_max = cgal->bytes();
        comm.allreduce(&bytes_min, &bytes_min, 1, FEAT::Dist::op_min);
        comm.allreduce(&bytes_max, &bytes_max, 1, FEAT::Dist::op_max);
        comm.allreduce(&bytes_sum, &bytes_sum, 1, FEAT::Dist::op_sum);
        logger.print(format_submemory_mm_t_kb("CGAL Size", bytes_sum, bytes_min, bytes_max, pad_len), info);
      }

    }
    /*----------------------------------------------------------------------------------------------*/
    /*------------------------- Domain and System Setup --------------------------------------------*/
    /*----------------------------------------------------------------------------------------------*/
    logger.print("Setting up domain", info);
    watch_create_domain.start();
    PartitionControl<DomainLevel> domain(comm, true);
    watch_create_domain.stop();

    // our materials
    std::vector<Material<SystemDataType>> materials;

    // our inflow boundaries
    std::vector<InflowBound> inflow_boundaries;

    watch_create_domain.start();
    //setup domain
    {
      create_voxel_domain(domain, args, config.get(), param_holder, &logger);


      if(const PropertyMap* problem_config = config->get_sub_section("problem_params"))
      {
        parse_materials(materials, problem_config->get_sub_section("list_materials"));
        parse_inflow_boundaries(inflow_boundaries, problem_config->get_sub_section("list_inflows"), materials, domain);
      }
      else
      {
        logger.print("Could not query problem_params section", error);
        Runtime::abort();
      }
      //print material and inflow information
      for(std::size_t i = 0; i < materials.size(); ++i)
      {
        logger.print(materials[i].format_string(int(i)), info);
      }
      for(const auto& in : inflow_boundaries)
      {
        logger.print(in.format_string(), info);
      }
      logger.print("Create Boundaries", info);
      // finally, create boundaries
      create_boundaries(domain, inflow_boundaries, param_holder, &logger);


      // now we create our fine resolved mesh part based on the surface file
      logger.print("Create FBM Meshpart and Assembler", info);
      {
        std::unique_ptr<Geometry::CGALWrapper<SystemDataType>> cgal_wrapper(nullptr);
        if(param_holder.meshpart_by_cgal_wrapper)
          cgal_wrapper = Gendie::create_cgal_wrapper(param_holder.surface_meshfile, comm);

        // create meshparts on all levels
        for(std::size_t i(0); i < domain.size_physical(); ++i)
        {

          // TODO: refactor with voxel map
          auto* ptr = param_holder.meshpart_by_cgal_wrapper ? Gendie::add_cgal_mesh_part(domain.at(i)->get_mesh_node(), cgal_wrapper.get(), "fbm")
                                                            : Gendie::add_voxel_mesh_part(domain.at(i)->get_mesh_node(), domain.get_voxel_map(), "fbm");
          if(ptr == nullptr)
          {
            logger.print("Meshpart is nullptr", error);
            Runtime::abort();
          }
          // // create FBM assembler
          domain.at(i)->create_fbm_assembler(domain.at(i).layer().comm(), "fbm");
        }
      }
      //and finally add trafo meshpart charts, only necessary if we use isoparametic trafos
      #ifdef FEAT_GENDIE_APP_ISOPARAM
      domain.add_trafo_mesh_part_charts();
      #endif

      watch_create_domain.stop();

      // at this point, our domain is compiled, so print chosen levels
      logger.print(String("Desired levels").pad_back(pad_len, pad_char) + ": " + domain.format_desired_levels(), info);
      logger.print(String("Chosen levels").pad_back(pad_len, pad_char) + ": " + domain.format_chosen_levels(), info);
      //print our config
      logger.print(param_holder.format_string(), info);
      //flush output
      logger.flush_print();

      watch_create_domain.stop();
      if(Gendie::check_for_config_option(config->query("postproc-params/mesh-hirarchy-vtk")))
      {
        const String& mesh_hirarchy_vtk = config->query("postproc-params/mesh-hirarchy-vtk").first;
        logger.print("Writing Mesh hirarchy output to '" + mesh_hirarchy_vtk + "_<LEVEL>.pvtu'", info);
        //todo: write out refinement hirarchy with mesh parts
        for(int i = 0; i < int(domain.size_physical()); ++i)
        {
          String cur_vtk_name = mesh_hirarchy_vtk + "_" + stringify(i).pad_front(2, '0');
          // get the mesh for the VTK output
          const MeshType* mesh = &domain.at(i)->get_mesh();
          const MeshNodeType* node = domain.at(i)->get_mesh_node();
          std::deque<String> part_names = node->get_mesh_part_names();

          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<MeshType> exporter(*mesh);

          std::vector<double> vtx_data(mesh->get_num_entities(0), 0.0);

          // loop over all mesh parts
          for(auto it = part_names.begin(); it != part_names.end(); ++it)
          {
            // get the mesh part
            const auto* part = node->find_mesh_part(*it);
            if(part == nullptr)
              continue;

            // get the vertex target set
            const auto& trg = part->template get_target_set<0>();

            // mark all indexes vertices
            for(Index k(0); k < trg.get_num_entities(); ++k)
              vtx_data[trg[k]] = 1.0;

            // add variable
            exporter.add_vertex_scalar(*it, vtx_data.data());

            // unmark vertices
            for(Index k(0); k < trg.get_num_entities(); ++k)
              vtx_data[trg[k]] = 0.0;
          }
          exporter.write(cur_vtk_name, domain.at(i).layer().comm());

        }
      }
    }


    logger.print("Compile Domain Assembler", info);
    {
      //create domain assembler for all levels
      const Index num_levels(domain.size_physical());
      for(Index i = 0; i < num_levels; ++i)
      {
        DomainLevel& dom_lvl = *domain.at(i);
        // in theory, we should also set a worker stretegy here, but is not necessary, due to omp...
        // if(num_worker_threads > 1)
        // {
        //   dom_lvl.domain_asm.set_max_worker_threads(num_worker_threads > 1 ? std::size_t(num_worker_threads) : std::size_t(0));
        //   dom_lvl.domain_asm.set_threading_strategy(Assembly::ThreadingStrategy::colored);
        // }
        dom_lvl.domain_asm.compile_all_elements();
      }

    }
    watch_create_domain.stop();

    logger.print("Compile Main Systemlevel", info);
    watch_create_system.start();
    //next up, create our defect system level, which does not require any mesh hirarchy, but at least the B matrices assembled
    std::shared_ptr<SystemLevel<SystemDataType, SystemIndexType>> system_level(new SystemLevel<SystemDataType, SystemIndexType>());
    system_level->assemble_gates(domain.front());
    system_level->assemble_pres_struct(domain.front()->space_pres);
    system_level->assemble_grad_div_matrices(domain.front()->domain_asm, domain.front()->space_velo, domain.front()->space_pres, param_holder.cubature_b);

    // TODO: mass matrix really required?
    system_level->assemble_velo_struct(domain.front()->space_velo);
    system_level->assemble_velocity_laplace_matrix(domain.front()->domain_asm, domain.front()->space_velo, materials.front().get_viscosity_model()->get_data()[0], true, param_holder.cubature_m);
    system_level->assemble_velocity_mass_matrix(domain.front()->domain_asm, domain.front()->space_velo, param_holder.cubature_m);

    //and compile the system matrix
    system_level->compile_system_matrix();

    watch_create_system.stop();

    // print system information
    {
      enum _SysID
      {
        elements = 0,
        dofs_v = 1,
        dofs_p = 2,
        num_entries = 3
      };

      Index dofs_total[_SysID::num_entries], dofs_min[_SysID::num_entries], dofs_max[_SysID::num_entries];

      {
        auto tv = system_level->gate_velo._freqs.clone(LAFEM::CloneMode::Deep);
        tv.format(1.0);
        dofs_total[_SysID::dofs_v] = Index(system_level->gate_velo.dot(tv, tv));
        dofs_min[_SysID::dofs_v] = dofs_max[_SysID::dofs_v] = Index(tv.template size<LAFEM::Perspective::pod>());
        dofs_total[_SysID::dofs_p] = dofs_min[_SysID::dofs_p] = dofs_max[_SysID::dofs_p] = system_level->gate_pres._freqs.size();
        dofs_total[_SysID::elements] = dofs_min[_SysID::elements] = dofs_max[_SysID::elements] = domain.front()->get_mesh().get_num_elements();
      }

      comm.allreduce(dofs_min, dofs_min, _SysID::num_entries, FEAT::Dist::op_min);
      comm.allreduce(dofs_max, dofs_max, _SysID::num_entries, FEAT::Dist::op_max);
      comm.allreduce(&dofs_total[_SysID::elements], &dofs_total[_SysID::elements], 1u, FEAT::Dist::op_sum);
      comm.allreduce(&dofs_total[_SysID::dofs_p], &dofs_total[_SysID::dofs_p], 1u, FEAT::Dist::op_sum);

      logger.print("\n----------------------------------------------------------------------\n", info);
      logger.print(format_index_mm("Element Count", dofs_total[_SysID::elements], dofs_min[_SysID::elements], dofs_max[_SysID::elements], pad_len), info);
      logger.print(format_index_mm("Velo Dofs", dofs_total[_SysID::dofs_v], dofs_min[_SysID::dofs_v], dofs_max[_SysID::dofs_v], pad_len), info);
      logger.print(format_index_mm("Pres Dofs", dofs_total[_SysID::dofs_p], dofs_min[_SysID::dofs_p], dofs_max[_SysID::dofs_p], pad_len), info);
      logger.print(format_index_mm("All Count", dofs_total[_SysID::dofs_v] + dofs_total[_SysID::dofs_p], dofs_min[_SysID::dofs_v] + dofs_min[_SysID::dofs_p], dofs_max[_SysID::dofs_v] + dofs_max[_SysID::dofs_p], pad_len), info);
    }

    logger.print("Compile Filter", info);
    // now we create our system defect filters
    watch_create_filter.start();
    assemble_filters(*system_level, domain.front()->get_mesh_node()->get_mesh_part_names(true), *domain.front(), inflow_boundaries, materials, param_holder, true);
    watch_create_filter.stop();
    // clear laplace matrix, since this will not be needed anymore
    // todo: unnessearry to save whole mass matrix for this
    // if(use_q2_fbm)
    // {
    //   system_level->clear_velocity_laplace_matrix();
    // }

    logger.print("Create Solver", info);
    // now we can construct our solver
    watch_create_solver.start();
    std::unique_ptr<NonlinearSteadyFlowSolverBase<SystemLevel<SystemDataType, SystemIndexType>::GlobalSystemVector>> flow_solver
        = Gendie::new_flow_solver(domain, system_level, inflow_boundaries, materials, config.get(), &logger,true, param_holder);
    watch_create_solver.stop();


    //print solver info
    logger.print(flow_solver->format_string(), info);

    logger.print("\n----------------------------------------------------------------------\n", info);
    logger.flush_print();

    // and now, init our solver
    flow_solver->init();

    logger.print("Initialize Solution", info);
    logger.print("\n--------------------------------------------------------------------------------------------------\n", info);
    logger.flush_print();
    //create solution vector and rhs
    typename SystemLevel<SystemDataType, SystemIndexType>::GlobalSystemVector vec_sol = system_level->create_global_vector_sys();
    typename SystemLevel<SystemDataType, SystemIndexType>::GlobalSystemVector vec_rhs = system_level->create_global_vector_sys();
    vec_sol.format();
    vec_rhs.format();
    // filter
    system_level->filter_sys.filter_sol(vec_sol);
    system_level->filter_sys.filter_rhs(vec_rhs);

    //init solution
    flow_solver->init_sol(vec_sol, vec_rhs);

    logger.print("Start Nonlinear Solver", info);
    logger.flush_print();
    //now solve
    FEAT::Solver::Status status = flow_solver->apply(vec_sol, vec_rhs);

    if(status != FEAT::Solver::Status::success)
    {
      logger.print("Solver did not converge with status code" + stringify(status), error);
    }

    logger.print("\n----------------------------------------------------------------------\n", info);
    // print timings
    logger.print(flow_solver->format_timings(), info);

    // and done
    flow_solver->done();

    double assembly_time_solver = flow_solver->get_assembly_time();
    double ls_time_solver = flow_solver->get_linear_solver_time();

    // reset solver to free memory
    flow_solver->reset();

    // and now, postprocessing, todo: refactor velo analyzer
    watch_analyze_sol.start();
    const auto mat_data = materials.front().get_viscosity_model()->get_data();
    auto velo_info = Gendie::CarreauVelocityAnalyser::compute(vec_sol.local().first(), domain.front()->space_velo, param_holder.cubature_postproc,
                                                materials.front().get_density_gram_per_unit(param_holder.mesh_unit_scale), inflow_boundaries.front().get_throughput(), mat_data[0] * param_holder.visc_scaling_factor, SystemDataType(0), mat_data[2], SystemDataType(2), mat_data[1], SystemDataType(40));
    {
      //create gate for communication
      // FEAT::Global::Gate<FEAT::LAFEM::DenseVector<SystemDataType,Index>, FEAT::LAFEM::VectorMirror<SystemDataType, Index>> tmp_gate;
      // FEAT::LAFEM::DenseVector<SystemDataType,Index> tmp_vec(system_level->gate_velo._freqs.size());
      // tmp_gate.convert(system_level->gate_velo, std::move(tmp_vec));
      // TODO: why not use scalar gate?
      velo_info.synchronize(comm, system_level->gate_scalar_velo);
    }
    watch_analyze_sol.stop();

    // print analysis
    logger.print(velo_info.format_string(), info);

    watch_write_out.start();
    if(param_holder.want_vtk)
      write_vtk(vec_sol, vec_rhs, domain, *system_level, velo_info, param_holder, &logger);
    watch_write_out.stop();


    // print memory usage
    {
      FEAT::MemoryUsage mi;
      std::size_t loc_peak_phys = mi.get_peak_physical();
      std::size_t loc_peak_virt = mi.get_peak_virtual();

      Index elemens, elements_total;
      elemens = elements_total = domain.front()->get_mesh().get_num_elements();

      std::size_t mem_per_cell_min, mem_per_cell_max;
      mem_per_cell_min = mem_per_cell_max = std::size_t(double(loc_peak_phys)/double(elemens));

      std::size_t peak_phys_min, peak_phys_max, peak_phys_total;
      std::size_t peak_virt_min, peak_virt_max, peak_virt_total;

      comm.allreduce(&loc_peak_phys, &peak_phys_min, 1u, FEAT::Dist::op_min);
      comm.allreduce(&loc_peak_phys, &peak_phys_max, 1u, FEAT::Dist::op_max);
      comm.allreduce(&loc_peak_phys, &peak_phys_total, 1u, FEAT::Dist::op_sum);
      comm.allreduce(&loc_peak_virt, &peak_virt_min, 1u, FEAT::Dist::op_min);
      comm.allreduce(&loc_peak_virt, &peak_virt_max, 1u, FEAT::Dist::op_max);
      comm.allreduce(&loc_peak_virt, &peak_virt_total, 1u, FEAT::Dist::op_sum);
      comm.allreduce(&mem_per_cell_min, &mem_per_cell_min, 1u, FEAT::Dist::op_min);
      comm.allreduce(&mem_per_cell_max, &mem_per_cell_max, 1u, FEAT::Dist::op_max);
      comm.allreduce(&elements_total, &elements_total, 1u, FEAT::Dist::op_sum);

      logger.print("\n----------------------------------------------------------------------\n", info);
      logger.print(format_submemory_mm_t("Peak Physical Memory", peak_phys_total, peak_phys_min, peak_phys_max, pad_len), info);
      logger.print(format_submemory_mm_t("Peak Virtual Memory", peak_virt_total, peak_virt_min, peak_virt_max, pad_len), info);
      logger.print(format_submemory_mm_t_kb("Peak Pyhsical Per Cell", std::size_t(double(peak_phys_total)/(double(elements_total))), mem_per_cell_min, mem_per_cell_max, pad_len), info);
    }

    logger.flush_print();

    watch_total.stop();

    {
      enum _TimeID
      {
        total = 0,
        create_domain = 1,
        create_system = 2,
        create_filter = 3,
        create_solver = 4,
        write_out = 5,
        analyze = 6,
        time_asm = 7,
        lin_sol = 8,
        num_entries = 9
      };
      FEAT::String s;
      double timings_max[num_entries], timings_min[num_entries], timings[num_entries];
      timings_min[_TimeID::total] = timings_max[_TimeID::total] = timings[_TimeID::total] = watch_total.elapsed();
      timings_min[_TimeID::create_domain] = timings_max[_TimeID::create_domain] = timings[_TimeID::create_domain] = watch_create_domain.elapsed();
      timings_min[_TimeID::create_system] = timings_max[_TimeID::create_system] = timings[_TimeID::create_system] = watch_create_system.elapsed();
      timings_min[_TimeID::create_filter] = timings_max[_TimeID::create_filter] = timings[_TimeID::create_filter] = watch_create_filter.elapsed();
      timings_min[_TimeID::create_solver] = timings_max[_TimeID::create_solver] = timings[_TimeID::create_solver] = watch_create_solver.elapsed();
      timings_min[_TimeID::write_out] = timings_max[_TimeID::write_out] = timings[_TimeID::write_out] = watch_write_out.elapsed();
      timings_min[_TimeID::analyze] = timings_max[_TimeID::analyze] = timings[_TimeID::analyze] = watch_analyze_sol.elapsed();
      timings_min[_TimeID::time_asm] = timings_max[_TimeID::time_asm] = timings[_TimeID::time_asm] = assembly_time_solver;
      timings_min[_TimeID::lin_sol] = timings_max[_TimeID::lin_sol] = timings[_TimeID::lin_sol] = ls_time_solver;
      // sync our timings
      comm.allreduce(timings_max, timings_max, _TimeID::num_entries, FEAT::Dist::op_max);
      comm.allreduce(timings_min, timings_min, _TimeID::num_entries, FEAT::Dist::op_min);


      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Overall Timings-------------------------------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");

      s += format_subtime_mm("Create Domain", timings[_TimeID::create_domain], timings[_TimeID::total], timings_min[_TimeID::create_domain], timings_max[_TimeID::create_domain], pad_len);
      s += format_subtime_mm("Create System", timings[_TimeID::create_system], timings[_TimeID::total], timings_min[_TimeID::create_system], timings_max[_TimeID::create_system], pad_len);
      s += format_subtime_mm("Create Filter", timings[_TimeID::create_filter], timings[_TimeID::total], timings_min[_TimeID::create_filter], timings_max[_TimeID::create_filter], pad_len);
      s += format_subtime_mm("Create Solver", timings[_TimeID::create_solver], timings[_TimeID::total], timings_min[_TimeID::create_solver], timings_max[_TimeID::create_solver], pad_len);
      s += format_subtime_mm("VTK Output", timings[_TimeID::write_out], timings[_TimeID::total], timings_min[_TimeID::write_out], timings_max[_TimeID::write_out], pad_len);
      s += format_subtime_mm("Solution Analyze", timings[_TimeID::analyze], timings[_TimeID::total], timings_min[_TimeID::analyze], timings_max[_TimeID::analyze], pad_len);
      s += format_subtime_mm("Total Time Assembly", timings[_TimeID::time_asm], timings[_TimeID::total], timings_min[_TimeID::time_asm], timings_max[_TimeID::time_asm], pad_len);
      s += format_subtime_mm("Total Time Linear Solver", timings[_TimeID::lin_sol], timings[_TimeID::total], timings_min[_TimeID::lin_sol], timings_max[_TimeID::lin_sol], pad_len);
      s += format_subtime_mm("Total Time", timings[_TimeID::total], timings[_TimeID::total], timings_min[_TimeID::total], timings_max[_TimeID::total], pad_len);

      logger.print(s, info);
    }

  }
} // namespace Gendie



int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  try
  {
    Gendie::main(argc, argv);
  }
  catch(std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return 0;
}
