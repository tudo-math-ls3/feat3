#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <control/domain/partitioner_domain_control.hpp>

#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_control_factory.hpp>

using namespace FEAT;

static void display_help();

template<typename Mem_, typename DT_, typename IT_, typename Mesh_>
struct MeshoptBoundaryApp
{
  /// The memory architecture. Although this looks freely chosable, it has to be Mem::Main for now because all the
  /// Hyperelasticity functionals are implemented for Mem::Main only
  typedef Mem_ MemType;
  /// The floating point type
  typedef DT_ DataType;
  /// The index type
  typedef IT_ IndexType;
  /// The type of mesh to use
  typedef Mesh_ MeshType;

  /// The only transformation available is the standard P1 or Q1 transformation
  typedef Trafo::Standard::Mapping<Mesh_> TrafoType;

#ifndef FEAST_HAVE_MPI
  // If we are in serial mode, there is no partitioning
  typedef Control::Domain::PartitionerDomainControl<Foundation::PExecutorNONE<DT_, IT_>, Mesh_> DomCtrl;
#else
#ifdef FEAT_HAVE_PARMETIS
  // If we have ParMETIS, we can use it for partitioning
  typedef Control::Domain::PartitionerDomainControl<Foundation::PExecutorParmetis<Foundation::ParmetisModePartKway>, Mesh_> DomCtrl;
#else
  // Otherwise we have to use the fallback partitioner
  typedef Control::Domain::PartitionerDomainControl<Foundation::PExecutorFallback<DT_, IT_>, Mesh_> DomCtrl;
#endif // FEAT_HAVE_PARMETIS
#endif // ! FEAST_HAVE_MPI

  /**
   * \brief Returns a descriptive string
   *
   * \returns The class name as string
   */
  static String name()
  {
    return "MeshoptBoundaryApp";
  }

  /**
   * \brief The routine that does the actual work
   */
  static int run(const String& meshopt_section_key, PropertyMap* meshopt_config, PropertyMap* solver_config,
  Geometry::MeshFileReader* mesh_file_reader, Geometry::MeshFileReader* chart_file_reader,
  int lvl_max, int lvl_min, const DataType delta_t, const DataType t_end,
  const bool write_vtk)
  {
    XASSERT(delta_t > DataType(0));
    XASSERT(t_end >= DataType(0));

    // This is the point we rotate our mesh around
    Tiny::Vector<DataType, MeshType::world_dim> midpoint(DataType(0.5));

    TimeStamp at;

    // Minimum number of cells we want to have in each patch
    Index part_min_elems(Util::Comm::size()*4);

    DomCtrl dom_ctrl(lvl_max, lvl_min, part_min_elems, mesh_file_reader, chart_file_reader);

    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl, TrafoType>> meshopt_ctrl(nullptr);
    meshopt_ctrl = Control::Meshopt::ControlFactory<Mem_, DT_, IT_, TrafoType>::create_meshopt_control(
      dom_ctrl, meshopt_section_key, meshopt_config, solver_config);

    // Print level information
    if(Util::Comm::rank() == 0)
    {
      std::cout << name() << " settings: " << std::endl;
      std::cout << "Timestep size: " << stringify_fp_fix(delta_t) << ", end time: " <<
        stringify_fp_fix(t_end) << std::endl;
      std::cout << "LVL-MAX: " <<
        dom_ctrl.get_levels().back()->get_level_index() << " [" << lvl_max << "]";
      std::cout << " LVL-MIN: " <<
        dom_ctrl.get_levels().front()->get_level_index() << " [" << lvl_min << "]" << std::endl;
    }

    String file_basename(name()+"_n"+stringify(Util::Comm::size()));

    // Get outer boundary MeshPart
    auto* outer_boundary = dom_ctrl.get_levels().back()->get_mesh_node()->find_mesh_part("outer");
    XASSERT(outer_boundary != nullptr);

    // Get vertex target set
    Geometry::TargetSet& boundary_set = outer_boundary->template get_target_set<0>();

    // A copy of the old vertex coordinates is kept here
    auto old_coords = meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep);

    // The mesh velocity is 1/delta_t*(coords_new - coords_old) and computed in each time step
    auto mesh_velocity = meshopt_ctrl->get_coords().clone();

    // Write initial vtk output
    if(write_vtk)
    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        String vtk_name = String(file_basename+"_pre_inital_lvl_"+stringify((*it)->get_level_index()));

        if(Util::Comm::rank() == 0)
          std::cout << "Writing " << vtk_name << std::endl;

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
        // Add everything from the MeshoptControl
        meshopt_ctrl->add_to_vtk_exporter(exporter, deque_position);
        exporter.write(vtk_name, int(Util::Comm::rank()), int(Util::Comm::size()));

        ++deque_position;
      }
    }

    // Copy the vertex coordinates to the buffer and get them via get_coords()
    meshopt_ctrl->mesh_to_buffer();
    auto new_coords = meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep);

    // Prepare the functional
    meshopt_ctrl->prepare(old_coords);
    // Optimise the mesh
    meshopt_ctrl->optimise();

    // Write vtk output
    if(write_vtk)
    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        String vtk_name = String(file_basename+"_post_inital_lvl_"+stringify((*it)->get_level_index()));

        if(Util::Comm::rank() == 0)
          std::cout << "Writing " << vtk_name << std::endl;

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
        // Add everything from the MeshoptControl
        meshopt_ctrl->add_to_vtk_exporter(exporter, deque_position);
        exporter.write(vtk_name, int(Util::Comm::rank()), int(Util::Comm::size()));

        ++deque_position;
      }
    }

    // Initial time
    DataType time(0);
    // Counter for timesteps
    Index n(0);

    while(time < t_end)
    {
      n++;
      time+= delta_t;

      if(Util::Comm::rank() == 0)
        std::cout << "Timestep " << n << " t = " << stringify_fp_fix(time) <<std::endl;

      // Save old vertex coordinates
      meshopt_ctrl->mesh_to_buffer();
      old_coords.clone(meshopt_ctrl->get_coords());

      // Get coords for modification
      auto& coords = (meshopt_ctrl->get_coords());
      auto& coords_loc = *coords;

      // Update the boundary
      for(Index i(0); i < boundary_set.get_num_entities(); ++i)
      {
        Index j = boundary_set[i];
        Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp0(coords_loc(j));
        Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp1(tmp0);

        tmp1(0) += delta_t * (DataType(1)*(tmp0(0) - midpoint(0))
            + Math::pow(DataType(1)*(tmp0(1) - midpoint(1)),DataType(3) ) );
        tmp1(1) -= delta_t * (DataType(1)*(tmp1(0) - midpoint(0))
            + Math::pow(DataType(1)*(tmp0(0) - midpoint(0)),DataType(3) ) );

        coords_loc(j, tmp1);
      }

      meshopt_ctrl->prepare(coords);
      meshopt_ctrl->optimise();

      // Compute mesh velocity
      mesh_velocity.axpy(meshopt_ctrl->get_coords(), old_coords, DataType(-1));
      mesh_velocity.scale(mesh_velocity, DataType(1)/delta_t);

      // Compute maximum of the mesh velocity
      DataType max_mesh_velocity(0);
      for(IT_ i(0); i < (*mesh_velocity).size(); ++i)
        max_mesh_velocity = Math::max(max_mesh_velocity, (*mesh_velocity)(i).norm_euclid());

      if(Util::Comm::rank() == 0)
        std::cout << "max. mesh velocity: " << stringify_fp_sci(max_mesh_velocity) << std::endl;

      if(write_vtk)
      {
        String vtk_name(file_basename+"_post_"+stringify(n));

        if(Util::Comm::rank() == 0)
          std::cout << "Writing " << vtk_name << std::endl;

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_ctrl.get_levels().back()->get_mesh());
        // Add mesh velocity
        exporter.add_vertex_vector("mesh_velocity", *mesh_velocity);
        // Add everything from the MeshoptControl
        meshopt_ctrl->add_to_vtk_exporter(exporter, int(dom_ctrl.get_levels().size())-1);
        // Write the file
        exporter.write(vtk_name, int(Util::Comm::rank()), int(Util::Comm::size()));
      }

    } // time loop

    if(Util::Comm::rank() == 0)
    {
      TimeStamp bt;
      std::cout << "Elapsed time: " << bt.elapsed(at) << std::endl;
    }

    return 0;

  }
}; // struct MeshSmootherApp

int main(int argc, char* argv[])
{
  // Even though this *looks* configureable, it is not: All HyperelasticityFunctionals are implemented for Mem::Main
  // only
  typedef Mem::Main MemType;
  // Floating point type
  typedef double DataType;
  // Index type
  typedef Index IndexType;

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> H2M2D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, 3, Real> S2M3D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, Real> S3M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, 1, Real> H1M1D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, 2, Real> H1M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, 3, Real> H1M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, 3, Real> H2M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, Real> H3M3D;

  int rank(0);
  int nprocs(0);

  // initialise
  FEAT::Runtime::initialise(argc, argv, rank, nprocs);
#ifdef FEAT_HAVE_MPI
  if (rank == 0)
  {
    std::cout << "NUM-PROCS: " << nprocs << std::endl;
  }
#endif

  // Mininum refinement level, parsed from the application config file
  int lvl_min(-1);
  // Maximum refinement level, parsed from the application config file
  int lvl_max(-1);
  // Timestep size, parsed from the application config file
  DataType delta_t(0);
  // End time, parsed from the application config file
  DataType t_end(0);
  // Filename to read the mesh from, parsed from the application config file
  String mesh_filename("");
  // Filename to read the seperate chart from (if any), parsed from the application config file
  String chart_filename("");
  // String containing the mesh type, read from the header of the mesh file
  String mesh_type("");
  // Do we want to write vtk files. Read from the command line arguments
  bool write_vtk(false);

  // Streams for synchronising information read from files
  std::stringstream synchstream_mesh;
  std::stringstream synchstream_chart;
  std::stringstream synchstream_app_config;
  std::stringstream synchstream_meshopt_config;
  std::stringstream synchstream_solver_config;

  // Create a parser for command line arguments.
  SimpleArgParser args(argc, argv);
  args.support("application_config");
  args.support("vtk");
  args.support("help");

  if( args.check("help") > -1 || args.num_args()==1)
  {
    display_help();
  }

  // Get unsupported command line arguments
  std::deque<std::pair<int,String> > unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
  }

  // Check if we want to write vtk files
  if(args.check("vtk") >= 0 )
  {
    write_vtk = true;
  }


  // Read the application config file on rank 0
  if(Util::Comm::rank() == 0)
  {
    // Input application configuration file name, required
    String application_config_filename("");
    // Check and parse --application_config
    if(args.check("application_config") != 1 )
    {
      std::cout << "You need to specify a application configuration file with --application_config.";
      throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --application_config");
    }
    else
    {
      args.parse("application_config", application_config_filename);
      std::cout << "Reading application configuration from file " << application_config_filename << std::endl;
      std::ifstream ifs(application_config_filename);
      if(!ifs.good())
        throw InternalError(__func__, __FILE__, __LINE__, "config file "+application_config_filename+" not found!");
      else
        synchstream_app_config << ifs.rdbuf();
    }
  }

#ifdef FEAT_HAVE_MPI
  // If we are in parallel mode, we need to synchronise the stream
  Util::Comm::synch_stringstream(synchstream_app_config);
#endif

  // Parse the application config from the (synchronised) stream
  PropertyMap* application_config = new PropertyMap;
  application_config->parse(synchstream_app_config, true);

  // Get the application settings section
  auto app_settings_section = application_config->query_section("ApplicationSettings");
  if(app_settings_section == nullptr)
    throw InternalError(__func__,__FILE__,__LINE__,
    "Application config is missing the mandatory ApplicationSettings section!");

  // We read the files only on rank 0. After reading, we synchronise the streams like above.
  if(Util::Comm::rank() == 0)
  {
    // Read the mesh file to stream
    auto mesh_filename_p = app_settings_section->query("mesh_file");
    if(!mesh_filename_p.second)
    {
      throw InternalError(__func__,__FILE__,__LINE__,
      "ApplicationSettings section is missing the mandatory mesh_file entry!");
    }
    else
    {
      mesh_filename = mesh_filename_p.first;
      std::ifstream ifs(mesh_filename);
      if(!ifs.good())
        throw InternalError(__func__, __FILE__, __LINE__, "mesh file "+mesh_filename+" not found!");
      else
      {
        std::cout << "Reading mesh from file " << mesh_filename << std::endl;
        synchstream_mesh << ifs.rdbuf();
      }
    }

    // Read the chart file to stream
    auto chart_filename_p = app_settings_section->query("chart_file");
    if(chart_filename_p.second)
    {
      chart_filename = chart_filename_p.first;
      std::ifstream ifs(chart_filename);
      if(!ifs.good())
        throw InternalError(__func__, __FILE__, __LINE__, "chart file "+chart_filename+" not found!");
      else
      {
        std::cout << "Reading charts from file " << chart_filename << std::endl;
        synchstream_chart << ifs.rdbuf();
      }
    }

    // Read configuration for mesh optimisation to stream
    auto meshopt_config_filename_p = app_settings_section->query("meshopt_config_file");
    if(!meshopt_config_filename_p.second)
    {
      throw InternalError(__func__,__FILE__,__LINE__,
      "ApplicationConfig section is missing the mandatory meshopt_config_file entry!");
    }
    else
    {
      std::ifstream ifs(meshopt_config_filename_p.first);
      if(!ifs.good())
        throw InternalError(__func__, __FILE__, __LINE__,
        "config file "+meshopt_config_filename_p.first+" not found!");
      else
      {
        std::cout << "Reading mesh optimisation config from file " <<meshopt_config_filename_p.first << std::endl;
        synchstream_meshopt_config << ifs.rdbuf();
      }
    }

    // Read solver configuration to stream
    auto solver_config_filename_p = app_settings_section->query("solver_config_file");
    if(!solver_config_filename_p.second)
    {
      throw InternalError(__func__,__FILE__,__LINE__,
      "ApplicationConfig section is missing the mandatory solver_config_file entry!");
    }
    else
    {
      std::ifstream ifs(solver_config_filename_p.first);
      if(!ifs.good())
        throw InternalError(__func__, __FILE__, __LINE__,
        "config file "+solver_config_filename_p.first+" not found!");
      else
      {
        std::cout << "Reading solver config from file " << solver_config_filename_p.first << std::endl;
        synchstream_solver_config << ifs.rdbuf();
      }
    }
  } // Util::Comm::rank() == 0

#ifdef FEAT_HAVE_MPI
  // Synchronise all those streams in parallel mode
  Util::Comm::synch_stringstream(synchstream_mesh);
  Util::Comm::synch_stringstream(synchstream_chart);
  Util::Comm::synch_stringstream(synchstream_meshopt_config);
  Util::Comm::synch_stringstream(synchstream_solver_config);
#endif

  // Create a MeshFileReader and parse the mesh stream
  Geometry::MeshFileReader* mesh_file_reader(new Geometry::MeshFileReader(synchstream_mesh));
  mesh_file_reader->read_root_markup();

  // Create a MeshFileReader and parse the chart stream
  Geometry::MeshFileReader* chart_file_reader(nullptr);
  if(!synchstream_chart.str().empty())
    chart_file_reader = new Geometry::MeshFileReader(synchstream_chart);

  // Create PropertyMaps and parse the configuration streams
  PropertyMap* meshopt_config = new PropertyMap;
  meshopt_config->parse(synchstream_meshopt_config, true);

  PropertyMap* solver_config = new PropertyMap;
  solver_config->parse(synchstream_solver_config, true);

  // Get the coarse mesh and finest mesh levels from the application settings
  auto lvl_min_p = app_settings_section->query("lvl_min");
  if(!lvl_min_p.second)
    lvl_min = 0;
  else
    lvl_min = std::stoi(lvl_min_p.first);

  auto lvl_max_p = app_settings_section->query("lvl_max");
  if(!lvl_max_p.second)
    lvl_max = 0;
  else
    lvl_max = std::stoi(lvl_max_p.first);

  // Get timestep size
  auto delta_t_p = app_settings_section->query("delta_t");
  if(!delta_t_p.second)
    throw InternalError(__func__,__FILE__,__LINE__,
    "ApplicationConfig section is missing the mandatory delta_t entry!");
  else
    delta_t = std::stod(delta_t_p.first);

  // Get end time
  auto t_end_p = app_settings_section->query("t_end");
  if(!delta_t_p.second)
    throw InternalError(__func__,__FILE__,__LINE__,
    "ApplicationConfig section is missing the mandatory t_end entry!");
  else
    t_end = std::stod(t_end_p.first);

  // Get the mesh optimiser key from the application settings
  auto meshoptimiser_key_p = app_settings_section->query("mesh_optimiser");
  if(!meshoptimiser_key_p.second)
    throw InternalError(__func__,__FILE__,__LINE__,
    "ApplicationConfig section is missing the mandatory meshoptimiser entry!");

  int ret(1);

  // Get the mesh type sting from the parsed mesh so we know with which template parameter to call the application
  mesh_type = mesh_file_reader->get_meshtype_string();

  // Call the appropriate class' run() function
  if(mesh_type == "conformal:hypercube:2:2")
  {
    ret = MeshoptBoundaryApp<MemType, DataType, IndexType, H2M2D>::run(
      meshoptimiser_key_p.first, meshopt_config, solver_config, mesh_file_reader, chart_file_reader,
      lvl_max, lvl_min, delta_t, t_end, write_vtk);
  }

  if(mesh_type == "conformal:simplex:2:2")
  {
    ret = MeshoptBoundaryApp<MemType, DataType, IndexType, S2M2D>::run(
      meshoptimiser_key_p.first, meshopt_config, solver_config, mesh_file_reader, chart_file_reader,
      lvl_max, lvl_min, delta_t, t_end, write_vtk);
  }

  delete mesh_file_reader;
  if(chart_file_reader != nullptr)
    delete chart_file_reader;

  FEAT::Runtime::finalise();
  return ret;
}

static void display_help()
{
  if(Util::Comm::rank() == 0)
  {
    std::cout << "meshopt_boundary-app: Moving the boundary of a mesh and computing an extension into the interiour"
    << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << " --application_config: Path to the application configuration file" << std::endl;
      std::cout << "Optional arguments:" << std::endl;
    std::cout << " --vtk: If this is set, vtk files are written" << std::endl;
    std::cout << " --help: Displays this text" << std::endl;
  }
}
