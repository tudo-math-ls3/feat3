#include <kernel/runtime.hpp>
#include <exception>
#include <iostream>

#include "applications/gendie/parsing_helper.hpp"
#include "gendie_common.hpp"


#include <applications/gendie/inflow_boundary.hpp>
#include <applications/gendie/materials.hpp>
#include <applications/gendie/logger.hpp>
#include <applications/gendie/system_assembler.hpp>
#include <applications/gendie/defect_assembler.hpp>
#include <applications/gendie/steady_stokes_solver.hpp>
#include <applications/gendie/steady_flow_solver.hpp>
#include <applications/gendie/cgal_meshpart_helper.hpp>
#include <applications/gendie/shearrate_assembler.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/util/omp_util.hpp>

#include <control/domain/deslagged_domain_control.hpp>

namespace Gendie
{

  /// our domain level
  typedef DomainLevelBase DomainLevelBackground;
  typedef Control::Domain::DeslaggedDomainLevelWrapper<DomainLevelBackground> DomainLevel;
  typedef Control::Domain::CGALDeslaggingPredicate<MeshType, Real> DeslaggingPredicateType;

  template<typename DT_, typename IT_> using SystemLevel = SystemLevelBase<DT_, IT_>;

  /// our partidomaincontrol
  template<typename DomainLevel_>
  using PartitionControl = Control::Domain::DeslaggedDomainControl<DomainLevel_>;

  typedef Gendie::InflowBoundary<MeshType> InflowBound;

  static constexpr std::size_t pad_len = 30u;
  static constexpr char pad_char = '.';

  struct ParameterHolder : public ArgParser
  {
    struct MeshDef
    {
      Tiny::Vector<SystemIndexType, dim> num_cells;
      Tiny::Vector<SystemDataType, dim> min;
      Tiny::Vector<SystemDataType, dim> max;

      friend std::ostream& operator<<(std::ostream& s, const MeshDef& d)
      {
        s << "MeshDef{ num_cells: " << d.num_cells << ", min: " << d.min << ", max: " << d.max << " }";
        return s;
      }
    };

    Parameter<String> surface_meshfile = parameter(String(""))
                                                        .property("problem_params/surface_mesh_filename")
                                                        .required()
                                                        .validator([](const auto& v) { return !v.empty(); });
    Parameter<std::deque<String>> base_meshfile_names = parameter(std::deque<String>{})
                                                          .property("problem-setup/base-mesh");
    Parameter<String> cubature_a = parameter(String("gauss-legendre:5")).name("cubature_a");
    Parameter<String> cubature_b = parameter(String("gauss-legendre:3")).name("cubature_b");
    Parameter<String> cubature_m = parameter(String("gauss-legendre:3")).name("cubature_m");
    Parameter<String> cubature_postproc = parameter(String("gauss-legendre:5")).name("cubature_postproc");
    Parameter<String> level_string = parameter(String("")).property("problem-setup/level").required();
    Parameter<String> voxel_map_file = parameter(String("")).property("problem-setup/voxel-map");
    Parameter<Real> min_gap = parameter(Real(0)).property("special-problem-definition/min-gap");
    Parameter<String> vtk_filename = parameter(String("")).property("postproc-params/vtk");
    Parameter<MeshDef> mesh_def = parameter(MeshDef{})
                                    .property("problem-setup/background-mesh")
                                    .parser(
                                      [](const String& s)
                                      {
                                        MeshDef def;
                                        std::istringstream stream(s);

                                        for(int i(0); i < dim; i++)
                                        {
                                          stream >> def.min(i);
                                        }
                                        for(int i(0); i < dim; i++)
                                        {
                                          stream >> def.max(i);
                                        }
                                        for(int i(0); i < dim; i++)
                                        {
                                          stream >> def.num_cells(i);
                                        }

                                        return def;
                                      })
                                    .validator(
                                      [](const auto& d)
                                      {
                                        bool result = true;
                                        for(int i(0); i < dim; i++)
                                        {
                                          result &= d.min(i) < d.max(i);
                                        }
                                        return result;
                                      });

    Parameter<bool> use_base_mesh = parameter(false)
                                      .property("problem-setup/use-base-mesh")
                                      .needs_if(base_meshfile_names, [](const auto& b) { return b; })
                                      .needs_if(mesh_def, [](const auto& b) { return !b; });
    Parameter<SystemDataType> mesh_unit_scale = parameter(SystemDataType(1E+3)).property("problem-setup/mesh-unit-factor");//.name("Setup: Mesh Unit Factor");
    Parameter<SystemDataType> velo_scaling_factor = parameter(SystemDataType(1E-3)).property("postproc-params/velo-scale");//.name("Postproc: Velocity Scaling Factor");
    Parameter<SystemDataType> visc_scaling_factor = parameter(SystemDataType(1)).property("postproc-params/visc-scale");//.name("Postproc: Viscosity Scaling Factor");
    Parameter<SystemDataType> pres_scaling_factor = parameter(SystemDataType(1E-5)).property("postproc-params/pres-scale");//.name("Postproc: Pressure Scaling Factor");
    // Parameter<bool> meshpart_by_cgal_wrapper = parameter(false).property("problem-setup/meshpart-by-cgal-wrapper");
    Parameter<bool> assemble_q2_filter = parameter(false).property("special-problem-definition/use-q2-fbm");
    Parameter<bool> refine_vtk = parameter(false).property("postproc-params/refine-vtk");
    Parameter<bool> use_true_vtk_name = parameter(false).property("postproc-params/vtk-use-true-name");
    Parameter<bool> vtk_debug = parameter(false).property("postproc-params/vtk-debug-data");
    Parameter<bool> measure_cgal_size = parameter(true).property("postproc-params/skip-cgal-messure");
    Parameter<String> mesh_hierarchy_vtk = parameter(String("")).property("postproc-params/mesh-hierarchy-vtk");

    Parameter<bool> want_vtk = parameter(false).derives_from(
      vtk_filename,
      [](const String& filename) { return !Gendie::check_option_off(filename); });
    Parameter<String> vtk_name =
      parameter(String(""))
        .derives_from(
          vtk_filename,
          use_true_vtk_name,
          level_string,
          [](const String& filename, const bool& use_true_name, const String& ls)
          {
            if(filename.empty())
            {
              return String("gendie_material") + String("-lvl") + ls.split_by_whitespaces().front();
            }
            else if(!use_true_name)
            {
              return filename + "_level_" + ls.split_by_whitespaces().front();
            }
            else
            {
              return filename;
            }
          });
    Parameter<bool> cheap_weights = parameter(true).property("problem-setup/cheap-weights");
    Parameter<bool> always_flush_print = parameter(false).property("problem-setup/always-flush-print");
    Parameter<bool> exact_deslagging = parameter(false).property("problem-setup/exact-deslagging");

    Parameter<SystemDataType> resolution = parameter(0.0).derives_from(
      mesh_def,
      level_string,
      [](const MeshDef& md, const String& ls)
      {
        SystemDataType r = Math::Limits<SystemDataType>::max();
        for(int i = 0; i < dim; ++i)
        {
          r = Math::min(r, Math::abs((md.max[i] - md.min[i]) / SystemDataType(md.num_cells[i])));
        }
        r *= Math::pow(SystemDataType(2), -SystemDataType(std::stoi(ls.split_by_whitespaces().front())));
        return r;
      });
  };

  template<typename ParamHolder_>
  void create_voxel_domain(PartitionControl<DomainLevel>& domain, const SimpleArgParser& /*args*/, const PropertyMap* config, ParamHolder_& param_holder, const Logger* logger)
  {
    logger->print("Going into create_voxel_domain", PrintLevel::info);
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
      domain.set_desired_levels(param_holder.level_string->split_by_whitespaces());
    }

    logger->print("Creating mesh");

    if(*param_holder.use_base_mesh)
    {
      std::deque<String> mesh_file_names = *param_holder.base_meshfile_names;

      XASSERTM(!mesh_file_names.empty(), "Could not parse base mesh file name");

      // set base mesh for domain and set bounding boxes
      // Our mesh file reader
      Geometry::MeshFileReader mesh_reader;
      mesh_reader.add_mesh_files(domain.comm(), mesh_file_names);
      mesh_reader.read_root_markup();
      String mesh_file_type = mesh_reader.get_meshtype_string();
      XASSERTM(mesh_file_type == "conformal:hypercube:3:3", "Meshfile type has to be of conformal:hypercube:3:3");
      domain.set_base_mesh(mesh_reader.parse(domain.get_atlas()));
    }
    else
    {
      // so now we can create our background mesh
      domain.create_base_mesh_3d(param_holder.mesh_def->num_cells[0], param_holder.mesh_def->num_cells[1], param_holder.mesh_def->num_cells[2], param_holder.mesh_def->min[0],
                                param_holder.mesh_def->max[0], param_holder.mesh_def->min[1], param_holder.mesh_def->max[1], param_holder.mesh_def->min[2], param_holder.mesh_def->max[2]);
    }

    logger->print("Setting CGAL Wrapper");
    auto deslagging_predicate = std::make_unique<DeslaggingPredicateType>(domain.comm(), *param_holder.surface_meshfile, Geometry::CGALFileMode::fm_off);
    deslagging_predicate->get_cgal_wrapper().set_expensive_intersection(*param_holder.exact_deslagging);
    deslagging_predicate->cheap_weight_calculation(*param_holder.cheap_weights);

    domain.set_deslagging_predicate(std::move(deslagging_predicate));
    domain.keep_deslagging_predicate();

    // now we can create our hirarchy
    logger->print("Creating hirarchy");
    domain.create_hierarchy();

    logger->print("Printing chosen parti info", PrintLevel::info);

    if(domain.front_layer().comm().rank() == 0)
    {
      logger->print(domain.get_chosen_parti_info(), PrintLevel::info);
      logger->flush_print();
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
      if(*param_holder.use_base_mesh)
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
          loc_unit_asm_inflow.assemble(filter_v_inflow, cur_dom.space_velo, inflow_boundaries.at(std::size_t(id)).get_diri_inflow_function(materials, *param_holder.mesh_unit_scale));
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
      system_level.assemble_fbm_filters(fbm_asm, cur_dom.space_velo, cur_dom.space_pres, assemble_mask, *param_holder.assemble_q2_filter, false);
    }

    // compile system filter
    system_level.compile_system_filter();
  }

  template<typename SystemLevel_, typename DomainLevel_, typename ScalarVector_, typename ParamHolder_>
  void write_vtk(const typename SystemLevel_::GlobalSystemVector& vec_sol, const typename SystemLevel_::GlobalSystemVector& vec_rhs, const DomainLevel_& domain, SystemLevel_& system_level, const ScalarVector_& vec_shearrate, const ScalarVector_& vec_viscosity, const ParamHolder_& param_holder, const Logger* logger)
  {
    if(!*param_holder.want_vtk)
    {
      logger->print("Vtk output not requested!", PrintLevel::warning);
      return;
    }
    logger->print("Writing VTK output to '" + *param_holder.vtk_name + ".pvtu'", PrintLevel::info);

    std::unique_ptr<MeshType> refined_mesh;
    // get the mesh for the VTK output
    const MeshType* mesh = &domain.front()->get_mesh();
    if(*param_holder.refine_vtk)
    {
      logger->print("Write out refined mesh", PrintLevel::info);
      Geometry::StandardRefinery<MeshType> refinery(domain.front()->get_mesh());
      refined_mesh = refinery.make_unique();
      mesh = refined_mesh.get();
    }

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(*mesh);

    // write velocity
    exporter.add_vertex_vector("Velocity [m/s]", vec_sol.local().template at<0>(), *param_holder.velo_scaling_factor);
    if(*param_holder.vtk_debug)
    {
      auto vec_copyt = vec_sol.clone(LAFEM::CloneMode::Deep);
      vec_copyt.format(-1.5);
      system_level.filter_sys.filter_sol(vec_copyt);
      exporter.add_vertex_vector("Velocity testdebug", vec_copyt.local().template at<0>());
    }
    if(*param_holder.vtk_debug)
      exporter.add_vertex_vector("Debug: rhs_v", vec_rhs.local().template at<0>());
    if(vec_viscosity.local().size() == vec_sol.local().template at<0>().size())
    {
      exporter.add_vertex_scalar("Viscosity [Pa s]", vec_viscosity.local(), *param_holder.visc_scaling_factor);
      exporter.add_vertex_scalar("NormShearRate [1/s]", vec_shearrate.local());
    }


    if(!param_holder.surface_meshfile->empty())
    {
      // parse off file
      //only first off file is parsed for now
      const String& surface_mesh_file = *param_holder.surface_meshfile;
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
      exporter.add_vertex_scalar("Shell", vec_dist, *param_holder.velo_scaling_factor);
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
    exporter.add_cell_scalar("Pressure [bar]", vtx_p, *param_holder.pres_scaling_factor);
    if(*param_holder.vtk_debug)
      exporter.add_cell_scalar("Debug: rhs_p", vtx_q);

    const bool use_q2_fbm = *param_holder.assemble_q2_filter;

    // write FBM masks if
    if(refined_mesh)
    {
      if(*param_holder.vtk_debug && use_q2_fbm)
        exporter.add_vertex_scalar("Debug: fbm_mask_v", system_level.fbm_mask_velo.data());

      if(*param_holder.vtk_debug && use_q2_fbm)
      {
        const auto& filter_vec_velo = system_level.get_local_velo_unit_filter_seq().find_or_add("fbm");
        LAFEM::DenseVectorBlocked<SystemDataType, SystemIndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        filter_vec_velo.filter_sol(filter_vec);

        exporter.add_vertex_vector("Debug: fbm_filter_v", filter_vec);
      }

      if(*param_holder.vtk_debug)
      {
        const auto& filter_vec_velo = system_level.get_local_velo_unit_filter_seq().find_or_add("noflow");
        LAFEM::DenseVectorBlocked<SystemDataType, SystemIndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        filter_vec_velo.filter_sol(filter_vec);

        exporter.add_vertex_vector("Debug: noflow_filter_v", filter_vec);
      }

      if(*param_holder.vtk_debug && (system_level.filter_interface_fbm.get_filter_vector().num_nzes() > 0))
      {
        const auto& filter_vec_velo = system_level.filter_interface_fbm;
        LAFEM::DenseVectorBlocked<SystemDataType, SystemIndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        filter_vec_velo.filter_sol(filter_vec);

        exporter.add_vertex_vector("Debug: interface_filter_v", filter_vec);
      }

      if(*param_holder.vtk_debug)
      {
        const auto& filter_vec_velo = system_level.get_local_velo_unit_filter_seq().find_or_add("inflow");
        LAFEM::DenseVectorBlocked<SystemDataType, SystemIndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        filter_vec_velo.filter_sol(filter_vec);

        exporter.add_vertex_vector("Debug: inflow_filter_v", filter_vec);
      }

      const std::vector<int>& mask_p = domain.front()->fbm_assembler->get_fbm_mask_vector(dim);
      const int nc = 1 << dim;
      std::vector<int> mask_p_ref;
      mask_p_ref.reserve(mask_p.size() * std::size_t(nc));
      for(int i : mask_p)
        for(int k(0); k < nc; ++k)
          mask_p_ref.push_back(i);
      if(*param_holder.vtk_debug)
        exporter.add_cell_scalar("Debug: fbm_mask_p", mask_p_ref.data());
    }
    else
    {
      if(*param_holder.vtk_debug)
      {
        exporter.add_vertex_scalar("Debug: fbm_mask_v", domain.front()->fbm_assembler->get_fbm_mask_vector(0).data());
        exporter.add_cell_scalar("Debug: fbm_mask_p", domain.front()->fbm_assembler->get_fbm_mask_vector(dim).data());
      }
    }

    // finally, write the VTK file
    exporter.write(*param_holder.vtk_name, domain.comm());

  }


  void main(int argc, char* argv[])
  {
    // we need a few timer
    FEAT::StopWatch watch_total, watch_create_domain, watch_create_system, watch_write_out, watch_create_filter, watch_create_solver, watch_analyze_sol;
    watch_total.start();
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());
    // setup our logging object
    PrintLevel comm_out_level = PrintLevel::info;
    Logger logger(comm, comm_out_level, PrintLevel::off);
    logger.print("-------------------Welcome to FEAT3 Steady Gendie Simulator tool--------------------", PrintLevel::info);
    logger.flush_print();

    // create arg parser
    SimpleArgParser args(argc, argv);
    args.support("config-file", "A list of config files in the ini file format.\n"
                                    "See the template ini file for details.");
    args.support("problem-file", "A list of additonal config files in the ini file format.");
    args.support("log-level", "The log level to be used.");

    if(args.check("log-level") > 0)
    {
      logger.std_out_lvl = PrintLevel(args.parse_default("log-level", 2));
      logger.log_lvl = PrintLevel(args.parse_default("log-level", 2));
    }

    logger.print(FEAT::String("LogLevel: ") + stringify(logger.std_out_lvl));


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
    logger.print("Number of Processes: " + stringify(comm.size()), PrintLevel::info);
    logger.print("Number of OMP Threads: " + stringify(feat_omp_get_max_threads()), PrintLevel::info);
    logger.print("System DataType: " + String(fp_typename), PrintLevel::info);
    logger.print("System IndexType: " + String(ix_typename), PrintLevel::info);

    logger.print("Parsing Parameters...", PrintLevel::info);
    logger.flush_print();

    // wrapper class for our parameters
    ParameterHolder param_holder;
    if(!param_holder.parse(*config))
    {
      for(const String& s : param_holder.errors())
      {
        logger.print(s, PrintLevel::error);
      }
    }
    logger.print(param_holder.display(), PrintLevel::info);
    logger.always_flush = *param_holder.always_flush_print;

    /*----------------------------------------------------------------------------------------------*/
    /*------------------------- Domain and System Setup --------------------------------------------*/
    /*----------------------------------------------------------------------------------------------*/
    logger.print("Setting up domain", PrintLevel::info);
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

      if(*param_holder.measure_cgal_size)
      {
        std::size_t bytes_sum, bytes_min, bytes_max;
        bytes_sum = bytes_min = bytes_max = dynamic_cast<DeslaggingPredicateType&>(domain.get_deslagging_predicate()).get_cgal_wrapper().bytes();
        comm.allreduce(&bytes_min, &bytes_min, 1, FEAT::Dist::op_min);
        comm.allreduce(&bytes_max, &bytes_max, 1, FEAT::Dist::op_max);
        comm.allreduce(&bytes_sum, &bytes_sum, 1, FEAT::Dist::op_sum);
        logger.print(format_submemory_mm_t_kb("CGAL Size", bytes_sum, bytes_min, bytes_max, pad_len), PrintLevel::info);
      }

      logger.print(String("Reading materials and inflow boundaries"), PrintLevel::info);
      String material_info, inflow_info;
      if(const PropertyMap* problem_config = config->get_sub_section("problem_params"))
      {
        auto mat_parse = parse_materials(materials, problem_config->get_sub_section("list_materials"), material_info);
        if(mat_parse != MaterialParserError::success)
        {
          logger.print(String("Material parsing went wrong with: ") + parse_error(mat_parse), PrintLevel::error);
          Runtime::abort();
        }
        auto inflow_parse = parse_inflow_boundaries(inflow_boundaries, problem_config->get_sub_section("list_inflows"), materials, domain, inflow_info);
        if(inflow_parse != InflowParserError::success)
        {
          logger.print(String("Inflow parsing went wrong with") + parse_error(inflow_parse), PrintLevel::error);
          Runtime::abort();
        }
      }
      else
      {
        logger.print("Could not query problem_params section", PrintLevel::error);
        Runtime::abort();
      }
      logger.print(material_info, PrintLevel::info);
      logger.print(inflow_info, PrintLevel::info);
      logger.print("Create Boundaries", PrintLevel::info);
      // finally, create boundaries
      create_boundaries(domain, inflow_boundaries, param_holder, &logger);


      // now we create our fine resolved mesh part based on the surface file
      logger.print("Create FBM Meshpart and Assembler", PrintLevel::info);
      {
        // release here to free used memory of cgal afterwards
        std::unique_ptr<DeslaggingPredicateType> deslagging_predicate(dynamic_cast<DeslaggingPredicateType*>(domain.release_deslagging_predicate()));
        const auto& cgal_wrapper = deslagging_predicate->get_cgal_wrapper();

        // create meshparts on all levels
        for(std::size_t i(0); i < domain.size_physical(); ++i)
        {
          auto* ptr = Gendie::add_cgal_mesh_part(domain.at(i)->get_mesh_node(), &cgal_wrapper, "fbm");
          if(ptr == nullptr)
          {
            logger.print("Meshpart is nullptr", PrintLevel::error);
            Runtime::abort();
          }
          // create FBM assembler
          domain.at(i)->create_fbm_assembler(domain.at(i).layer().comm(), "fbm");
        }
      }
      //and finally add trafo meshpart charts, only necessary if we use isoparametic trafos
      #ifdef FEAT_GENDIE_APP_ISOPARAM
      domain.add_trafo_mesh_part_charts();
      #endif

      // at this point, our domain is compiled, so print chosen levels
      logger.print(String("Desired levels").pad_back(pad_len, pad_char) + ": " + domain.format_desired_levels(), PrintLevel::info);
      logger.print(String("Chosen levels").pad_back(pad_len, pad_char) + ": " + domain.format_chosen_levels(), PrintLevel::info);
      // print config
      // logger.print(param_holder.display());
      //flush output
      logger.flush_print();
    }

    watch_create_domain.stop();

    // TODO fix this as well
    // if(Gendie::check_for_config_option(config->query("postproc-params/mesh-hirarchy-vtk")))
    if(param_holder.mesh_hierarchy_vtk.set_by_user())
    {
      const String& mesh_hirarchy_vtk = *param_holder.mesh_hierarchy_vtk; // config->query("postproc-params/mesh-hirarchy-vtk").first;
      logger.print("Writing Mesh hirarchy output to '" + mesh_hirarchy_vtk + "_<LEVEL>.pvtu'", PrintLevel::info);
      //todo: write out refinement hirarchy with mesh parts
      for(std::size_t i = 0; i < std::size_t(domain.size_physical()); ++i)
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


    logger.print("Compile Domain Assembler", PrintLevel::info);
    watch_create_domain.start();
    {
      //create domain assembler for all levels
      const Index num_levels(domain.size_physical());
      for(Index i = 0; i < num_levels; ++i)
      {
        DomainLevel& dom_lvl = *domain.at(i);
        // in theory, we should also set a worker stretegy here, but is not necessary, due to omp...
        dom_lvl.domain_asm.set_max_worker_threads(std::size_t(0));
        // if(num_worker_threads > 1)
        // {
        //   dom_lvl.domain_asm.set_max_worker_threads(num_worker_threads > 1 ? std::size_t(num_worker_threads) : std::size_t(0));
        //   dom_lvl.domain_asm.set_threading_strategy(Assembly::ThreadingStrategy::colored);
        // }
        dom_lvl.domain_asm.compile_all_elements();
      }

    }
    watch_create_domain.stop();

    logger.print("Compile Main Systemlevel", PrintLevel::info);
    watch_create_system.start();
    //next up, create our defect system level, which does not require any mesh hirarchy, but at least the B matrices assembled
    std::shared_ptr<SystemLevel<SystemDataType, SystemIndexType>> system_level(new SystemLevel<SystemDataType, SystemIndexType>());
    system_level->assemble_gates(domain.front());
    system_level->assemble_grad_div_matrices(domain.front()->domain_asm, domain.front()->space_velo, domain.front()->space_pres, *param_holder.cubature_b);

    // TODO: mass matrix really required?
    system_level->assemble_velo_struct(domain.front()->space_velo);
    system_level->assemble_velocity_laplace_matrix(domain.front()->domain_asm, domain.front()->space_velo, materials.front().get_default_viscosity(*param_holder.mesh_unit_scale, SystemDataType(1)), true, *param_holder.cubature_m);
    system_level->assemble_velocity_mass_matrix(domain.front()->domain_asm, domain.front()->space_velo, *param_holder.cubature_m);

    //and compile the system matrix
    system_level->compile_system_matrix();

    watch_create_system.stop();

    // print system information
    {
      enum class _SysID : std::size_t
      {
        elements = 0,
        dofs_v,
        dofs_p,
        num_entries
      };

      Index dofs_total[to_underlying(_SysID::num_entries)], dofs_min[to_underlying(_SysID::num_entries)], dofs_max[to_underlying(_SysID::num_entries)];

      {
        auto tv = system_level->gate_velo._freqs.clone(LAFEM::CloneMode::Deep);
        tv.format(1.0);
        dofs_total[to_underlying(_SysID::dofs_v)] = Index(system_level->gate_velo.dot(tv, tv));
        dofs_min[to_underlying(_SysID::dofs_v)] = dofs_max[to_underlying(_SysID::dofs_v)] = Index(tv.size_raw());
        dofs_total[to_underlying(_SysID::dofs_p)] = dofs_min[to_underlying(_SysID::dofs_p)] = dofs_max[to_underlying(_SysID::dofs_p)] = system_level->gate_pres._freqs.size();
        dofs_total[to_underlying(_SysID::elements)] = dofs_min[to_underlying(_SysID::elements)] = dofs_max[to_underlying(_SysID::elements)] = domain.front()->get_mesh().get_num_elements();
      }

      comm.allreduce(dofs_min, dofs_min, to_underlying(_SysID::num_entries), FEAT::Dist::op_min);
      comm.allreduce(dofs_max, dofs_max, to_underlying(_SysID::num_entries), FEAT::Dist::op_max);
      comm.allreduce(&dofs_total[to_underlying(_SysID::elements)], &dofs_total[to_underlying(_SysID::elements)], 1u, FEAT::Dist::op_sum);
      comm.allreduce(&dofs_total[to_underlying(_SysID::dofs_p)], &dofs_total[to_underlying(_SysID::dofs_p)], 1u, FEAT::Dist::op_sum);

      logger.print("\n----------------------------------------------------------------------\n", PrintLevel::info);
      logger.print(format_index_mm("Element Count", dofs_total[to_underlying(_SysID::elements)], dofs_min[to_underlying(_SysID::elements)], dofs_max[to_underlying(_SysID::elements)], pad_len), PrintLevel::info);
      logger.print(format_index_mm("Velo Dofs", dofs_total[to_underlying(_SysID::dofs_v)], dofs_min[to_underlying(_SysID::dofs_v)], dofs_max[to_underlying(_SysID::dofs_v)], pad_len), PrintLevel::info);
      logger.print(format_index_mm("Pres Dofs", dofs_total[to_underlying(_SysID::dofs_p)], dofs_min[to_underlying(_SysID::dofs_p)], dofs_max[to_underlying(_SysID::dofs_p)], pad_len), PrintLevel::info);
      logger.print(format_index_mm("All Count", dofs_total[to_underlying(_SysID::dofs_v)] + dofs_total[to_underlying(_SysID::dofs_p)], dofs_min[to_underlying(_SysID::dofs_v)] + dofs_min[to_underlying(_SysID::dofs_p)], dofs_max[to_underlying(_SysID::dofs_v)] + dofs_max[to_underlying(_SysID::dofs_p)], pad_len), PrintLevel::info);
    }

    logger.print("Compile Filter", PrintLevel::info);
    // now we create our system defect filters
    watch_create_filter.start();
    assemble_filters(*system_level, domain.front()->get_mesh_node()->get_mesh_part_names(true), *domain.front(), inflow_boundaries, materials, param_holder, true);
    watch_create_filter.stop();

    logger.print("Create Solver", PrintLevel::info);
    // now we can construct our solver
    watch_create_solver.start();
    const PropertyMap* solver_config = config->query_section("solver-params");
    XASSERTM(solver_config, "No solver config section");
    Intern::SolverParser<SystemDataType> solver_parser;
    if(!solver_parser.parse(*solver_config))
    {
      for(const auto& err: solver_parser.errors())
      {
        logger.print(err, PrintLevel::error);
      }
      XABORTM("Could not parse solver section");
    }
    SteadyFlowSolverFactory flow_solver_factory(domain, system_level, inflow_boundaries, materials, solver_config, true, param_holder);
    auto solver_type = flow_solver_factory.parse_solver_type(*solver_parser.nl_solver);
    auto flow_solver = flow_solver_factory.create(solver_type, &logger);
    watch_create_solver.stop();

    if(materials.front().get_viscosity_model()->get_type() != MaterialType::pseudo_power_law)
    {
      watch_create_solver.start();
      flow_solver = flow_solver_factory.create(solver_type, &logger);
      watch_create_solver.stop();
    }
    else // we directly start with a beefed up ts solver
    {
      logger.print("Overwriting solver configuration for PowerLaw", PrintLevel::info);
      solver_type = decltype(solver_type)::pseudo_ts;
      auto ts_config = config->query_section("solver-params")->clone();
      watch_create_solver.start();
      Gendie::_adjust_solver_config(ts_config.get(), 64, 64);
      flow_solver = flow_solver_factory.create_overwrite_config(solver_type, ts_config.get(), &logger);
      watch_create_solver.stop();
    }

    //print solver info
    logger.print(flow_solver->format_string(), PrintLevel::info);
    logger.print(flow_solver->format_linear_solver(), PrintLevel::verbose);

    logger.print("\n----------------------------------------------------------------------\n", PrintLevel::info);
    logger.flush_print();

    // and now, init our solver
    flow_solver->init();

    logger.print("Initialize Solution", PrintLevel::info);
    logger.print("\n--------------------------------------------------------------------------------------------------\n", PrintLevel::info);
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

    logger.print("Start Nonlinear Solver", PrintLevel::info);
    logger.flush_print();
    //now solve
    NonLinearStatus status = flow_solver->apply(vec_sol, vec_rhs);

    if(status != NonLinearStatus::success)
    {
      logger.print("Solver did not converge with status code" + stringify(status), PrintLevel::warning);
      if(solver_type != decltype(solver_type)::pseudo_ts)
      {
        logger.print("Printing previous solver timings", PrintLevel::verbose);
        logger.print(flow_solver->format_timings(), PrintLevel::verbose);
        logger.print("\n----------------------------------------------------------------------\n", PrintLevel::info);
        logger.print("Retrying with pseudo timestepping", PrintLevel::info);
        logger.flush_print();
        flow_solver->done();
        flow_solver.reset();
        // always create ts solver with double solver precision in case we did not converge
        auto ts_config = config->query_section("solver-params")->clone();
        Gendie::_adjust_solver_config(ts_config.get(), 64, 64);
        flow_solver = flow_solver_factory.create_overwrite_config(decltype(solver_type)::pseudo_ts, ts_config.get(), &logger);
        logger.print(flow_solver->format_string(), PrintLevel::info);
        logger.print(flow_solver->format_linear_solver(), PrintLevel::verbose);
        flow_solver->init();
        // vec_sol.format();
        flow_solver->init_sol(vec_sol, vec_rhs);
        status = flow_solver->apply(vec_sol, vec_rhs);
        if(status != NonLinearStatus::success)
        {
          logger.print("Pseudo ts solver did not converge with status code" + stringify(status), PrintLevel::error);
        }
      }
    }

    logger.print("\n----------------------------------------------------------------------\n", PrintLevel::info);
    // print timings
    logger.print(flow_solver->format_timings(), PrintLevel::info);

    // and done
    flow_solver->done();

    double assembly_time_solver = flow_solver->get_assembly_time();
    double ls_time_solver = flow_solver->get_linear_solver_time();

    // reset solver to free memory
    flow_solver.reset();
    // also free solver factory
    flow_solver_factory.reset();

    FEAT::Global::Vector<FEAT::LAFEM::DenseVector<SystemDataType, SystemIndexType>, FEAT::LAFEM::VectorMirror<SystemDataType, SystemIndexType>> vec_shearrate{&system_level->gate_scalar_velo}, vec_viscosity{&system_level->gate_scalar_velo};
    // and now, postprocessing, todo: refactor velo analyzer
    watch_analyze_sol.start();
    if(*param_holder.want_vtk)
    {
      vec_shearrate.local() = FEAT::LAFEM::DenseVector<SystemDataType, SystemIndexType>(domain.front()->space_velo.get_num_dofs());
      vec_shearrate.format();
      {
        SystemLevel<SystemDataType, SystemIndexType>::GlobalVeloVector sol_velo(&system_level->gate_velo, vec_sol.local().first().clone(FEAT::LAFEM::CloneMode::Shallow));
        auto mass_vec = system_level->velo_mass_matrix.lump_rows(true);
        mass_vec.local().component_invert(mass_vec.local());
        Gendie::assemble_shear_vector(domain.front()->domain_asm, vec_shearrate.local(), sol_velo, mass_vec.local(), domain.front()->space_velo, *param_holder.cubature_postproc);
      }
      auto visco_func = materials.front().get_visc_func(*param_holder.mesh_unit_scale);
      vec_viscosity.local() = FEAT::LAFEM::DenseVector<SystemDataType, SystemIndexType>(domain.front()->space_velo.get_num_dofs());
      vec_viscosity.format();

      auto abiat_at = materials.front().get_abiat_aT(inflow_boundaries.front().get_temperature());

      {
        auto visco_val = vec_viscosity.local().elements_view_w();
        const auto shear_val = vec_shearrate.local().elements_view_r();

        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < vec_shearrate.local().size(); ++i)
        {
          visco_val[i] = visco_func(shear_val[i], abiat_at);
        }
      }

      // to be assure, average, so we have the same value everywhere
      vec_viscosity.sync_1();
      vec_shearrate.sync_1();
    }
    watch_analyze_sol.stop();

    watch_write_out.start();
    write_vtk(vec_sol, vec_rhs, domain, *system_level, vec_shearrate, vec_viscosity, param_holder, &logger);
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

      logger.print("\n----------------------------------------------------------------------\n", PrintLevel::info);
      logger.print(format_submemory_mm_t("Peak Physical Memory", peak_phys_total, peak_phys_min, peak_phys_max, pad_len), PrintLevel::info);
      logger.print(format_submemory_mm_t("Peak Virtual Memory", peak_virt_total, peak_virt_min, peak_virt_max, pad_len), PrintLevel::info);
      logger.print(format_submemory_mm_t_kb("Peak Pyhsical Per Cell", std::size_t(double(peak_phys_total)/(double(elements_total))), mem_per_cell_min, mem_per_cell_max, pad_len), PrintLevel::info);
    }

    logger.flush_print();

    watch_total.stop();

    {
      enum class _TimeID : std::size_t
      {
        total = 0,
        create_domain,
        create_system,
        create_filter,
        create_solver,
        write_out,
        analyze,
        time_asm,
        lin_sol,
        create_hirarch,
        num_entries
      };
      FEAT::String s;
      double timings_max[to_underlying(_TimeID::num_entries)], timings_min[to_underlying(_TimeID::num_entries)], timings[to_underlying(_TimeID::num_entries)];
      timings_min[to_underlying(_TimeID::total)] = timings_max[to_underlying(_TimeID::total)] = timings[to_underlying(_TimeID::total)] = watch_total.elapsed();
      timings_min[to_underlying(_TimeID::create_domain)] = timings_max[to_underlying(_TimeID::create_domain)] = timings[to_underlying(_TimeID::create_domain)] = watch_create_domain.elapsed();
      timings_min[to_underlying(_TimeID::create_system)] = timings_max[to_underlying(_TimeID::create_system)] = timings[to_underlying(_TimeID::create_system)] = watch_create_system.elapsed();
      timings_min[to_underlying(_TimeID::create_filter)] = timings_max[to_underlying(_TimeID::create_filter)] = timings[to_underlying(_TimeID::create_filter)] = watch_create_filter.elapsed();
      timings_min[to_underlying(_TimeID::create_solver)] = timings_max[to_underlying(_TimeID::create_solver)] = timings[to_underlying(_TimeID::create_solver)] = watch_create_solver.elapsed();
      timings_min[to_underlying(_TimeID::write_out)] = timings_max[to_underlying(_TimeID::write_out)] = timings[to_underlying(_TimeID::write_out)] = watch_write_out.elapsed();
      timings_min[to_underlying(_TimeID::analyze)] = timings_max[to_underlying(_TimeID::analyze)] = timings[to_underlying(_TimeID::analyze)] = watch_analyze_sol.elapsed();
      timings_min[to_underlying(_TimeID::time_asm)] = timings_max[to_underlying(_TimeID::time_asm)] = timings[to_underlying(_TimeID::time_asm)] = assembly_time_solver;
      timings_min[to_underlying(_TimeID::lin_sol)] = timings_max[to_underlying(_TimeID::lin_sol)] = timings[to_underlying(_TimeID::lin_sol)] = ls_time_solver;
      timings_min[to_underlying(_TimeID::create_hirarch)] = timings_max[to_underlying(_TimeID::create_hirarch)] = timings[to_underlying(_TimeID::create_hirarch)] = domain.get_watch_hierarchy().elapsed();
      // sync our timings
      comm.allreduce(timings_max, timings_max, to_underlying(_TimeID::num_entries), FEAT::Dist::op_max);
      comm.allreduce(timings_min, timings_min, to_underlying(_TimeID::num_entries), FEAT::Dist::op_min);


      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Overall Timings-------------------------------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");

      s += format_subtime_mm("Create Domain", timings[to_underlying(_TimeID::create_domain)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::create_domain)], timings_max[to_underlying(_TimeID::create_domain)], pad_len);
      s += format_subtime_mm("Create Domain Hirarch", timings[to_underlying(_TimeID::create_hirarch)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::create_hirarch)], timings_max[to_underlying(_TimeID::create_hirarch)], pad_len);
      s += format_subtime_mm("Create System", timings[to_underlying(_TimeID::create_system)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::create_system)], timings_max[to_underlying(_TimeID::create_system)], pad_len);
      s += format_subtime_mm("Create Filter", timings[to_underlying(_TimeID::create_filter)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::create_filter)], timings_max[to_underlying(_TimeID::create_filter)], pad_len);
      s += format_subtime_mm("Create Solver", timings[to_underlying(_TimeID::create_solver)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::create_solver)], timings_max[to_underlying(_TimeID::create_solver)], pad_len);
      s += format_subtime_mm("VTK Output", timings[to_underlying(_TimeID::write_out)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::write_out)], timings_max[to_underlying(_TimeID::write_out)], pad_len);
      s += format_subtime_mm("Solution Analyze", timings[to_underlying(_TimeID::analyze)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::analyze)], timings_max[to_underlying(_TimeID::analyze)], pad_len);
      s += format_subtime_mm("Total Time Assembly", timings[to_underlying(_TimeID::time_asm)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::time_asm)], timings_max[to_underlying(_TimeID::time_asm)], pad_len);
      s += format_subtime_mm("Total Time Linear Solver", timings[to_underlying(_TimeID::lin_sol)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::lin_sol)], timings_max[to_underlying(_TimeID::lin_sol)], pad_len);
      s += format_subtime_mm("Total Time", timings[to_underlying(_TimeID::total)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::total)], timings_max[to_underlying(_TimeID::total)], pad_len);

      logger.print(s, PrintLevel::info);
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
