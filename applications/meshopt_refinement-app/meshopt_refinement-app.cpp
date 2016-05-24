#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_quality_heuristic.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <control/domain/partitioner_domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_control_factory.hpp>

using namespace FEAST;

static void display_help();

#ifndef SERIAL
static void synch_stringstream(std::stringstream& iss, Util::Communicator comm = Util::Communicator(MPI_COMM_WORLD))
{
  Index me(Util::Comm::rank(comm));
  Index size;
  std::string str;

  //bcast size
  if(me == 0)
  {
    str = (iss.str());
    size = str.length();
  }
  // synchronize length
  Util::Comm::bcast(&size, 1, 0, comm);

  //allocate
  char* buf = new char[size + 1];

  //fill
  if(me == 0) //master
  {
    std::strcpy(buf, str.c_str());
  }

  //bcast data
  Util::Comm::bcast(buf, size, 0, comm);

  //convert
  if(me != 0)
  {
    std::string res_str(buf, size);
    iss << res_str;
  }

  delete[] buf;
}
#endif


template<typename Mem_, typename DT_, typename IT_, typename Mesh_>
struct MeshRefinementOptimiserApp
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

#ifdef SERIAL
  // If we are in serial mode, there is no partitioning
  typedef Control::Domain::PartitionerDomainControl<Foundation::PExecutorNONE<DT_, IT_>, Mesh_> DomCtrl;
#else
#ifdef FEAST_HAVE_PARMETIS
  // If we have ParMETIS, we can use it for partitioning
  typedef Control::Domain::PartitionerDomainControl
  <
    Foundation::PExecutorParmetis<Foundation::ParmetisModePartKway>,
    Mesh_
  >
  DomCtrl;
#else
  // Otherwise we have to use the fallback partitioner
  typedef Control::Domain::PartitionerDomainControl<Foundation::PExecutorFallback<DT_, IT_>, Mesh_> DomCtrl;
#endif // FEAST_HAVE_PARMETIS
#endif // SERIAL

  /**
   * \brief Returns a descriptive string
   *
   * \returns The class name as string
   */
  static String name()
  {
    return "MeshoptRefinementApp";
  }

  static int run(const String& meshopt_section_key, PropertyMap* meshopt_config, PropertyMap* solver_config,
  Geometry::MeshFileReader* mesh_file_reader, Geometry::MeshFileReader* chart_file_reader,
  int lvl_max, int lvl_min, const bool write_vtk)
  {
    // Base for vtk filenames, we attach some more information to that
    String file_basename(name()+"_n"+stringify(Util::Comm::size()));

    // Minimum number of cells we want to have in each patch
    Index part_min_elems(Util::Comm::size()*4);

    // Create the DomainControl
    DomCtrl dom_ctrl(lvl_max, lvl_min, part_min_elems, mesh_file_reader, chart_file_reader, Geometry::AdaptMode::none);
    // Print level information
    if(Util::Comm::rank() == 0)
    {
      std::cout << name() << " settings: " << std::endl;
      std::cout << "LVL-MAX: " <<
        dom_ctrl.get_levels().back()->get_level_index() << " [" << lvl_max << "]";
      std::cout << " LVL-MIN: " <<
        dom_ctrl.get_levels().front()->get_level_index() << " [" << lvl_min << "]" << std::endl;
      std::cout << "Cells: " << dom_ctrl.get_levels().back()->get_mesh().get_num_entities(MeshType::shape_dim) <<
        " vertices: " << dom_ctrl.get_levels().back()->get_mesh().get_num_entities(0) <<
        " DoF: " << MeshType::world_dim*dom_ctrl.get_levels().back()->get_mesh().get_num_entities(0) << std::endl;

    }

    // Create the MeshoptControl
    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl, TrafoType>> meshopt_ctrl(nullptr);
    meshopt_ctrl = Control::Meshopt::ControlFactory<Mem_, DT_, IT_, TrafoType>::create_meshopt_control(
      dom_ctrl, meshopt_section_key, meshopt_config, solver_config);

    // Adapt the finest level
    dom_ctrl.get_levels().back()->get_mesh_node()->adapt();

    meshopt_ctrl->mesh_to_buffer();
    auto new_coords = meshopt_ctrl->get_coords().clone();

    meshopt_ctrl->prepare(new_coords);

    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());

        // Write initial vtk output
        if(write_vtk)
        {
          String vtk_name = String(file_basename+"_pre_lvl_"+stringify(lvl_index));
          if(Util::Comm::rank() == 0)
            std::cout << "Writing " << vtk_name << std::endl;

          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
          meshopt_ctrl->add_to_vtk_exporter(exporter, deque_position);
          exporter.write(vtk_name, int(Util::Comm::rank()), int(Utill::Comm::size()));
        }

        auto quality = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::compute(
          (*it)->get_mesh().template get_index_set<MeshType::shape_dim, 0>(), (*it)->get_mesh().get_vertex_set());

        auto angle = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::angle(
          (*it)->get_mesh().template get_index_set<MeshType::shape_dim, 0>(), (*it)->get_mesh().get_vertex_set());

        auto min_quality(quality);
        auto min_angle(angle);
#ifdef FEAST_HAVE_MPI
        Util::Comm::allreduce(&quality, Index(1), &min_quality, MPI_MIN);
        Util::Comm::allreduce(&angle, Index(1), &min_angle, MPI_MIN);
#endif
        if(Util::Comm::rank() == 0)
          std::cout << "Pre: Level " << lvl_index << ": Quality indicator " << " " <<
            stringify_fp_sci(min_quality) << ", minimum angle " << stringify_fp_fix(min_angle) << std::endl;

        ++deque_position;
      }
    }

    TimeStamp pre_opt;
    meshopt_ctrl->optimise();
    TimeStamp post_opt;

    if(Util::Comm::rank() == 0)
      std::cout << "Solve time: " << post_opt.elapsed(pre_opt) << std::endl;

    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());

        if(write_vtk)
        {
          String vtk_name = String(file_basename+"_post_lvl_"+stringify(lvl_index));

          if(Util::Comm::rank() == 0)
            std::cout << "Writing " << vtk_name << std::endl;

          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
          meshopt_ctrl->add_to_vtk_exporter(exporter, deque_position);
          exporter.write(vtk_name, int(Comm::rank()), int(Comm::size()));
        }

        auto quality = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::compute(
          (*it)->get_mesh().template get_index_set<MeshType::shape_dim, 0>(), (*it)->get_mesh().get_vertex_set());

        auto angle = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::angle(
          (*it)->get_mesh().template get_index_set<MeshType::shape_dim, 0>(), (*it)->get_mesh().get_vertex_set());

        auto min_quality(quality);
        auto min_angle(angle);
#ifdef FEAST_HAVE_MPI
        Util::Comm::allreduce(&quality, Index(1), &min_quality, MPI_MIN);
        Util::Comm::allreduce(&angle, Index(1), &min_angle, MPI_MIN);
#endif
        if(Util::Comm::rank() == 0)
          std::cout << "Post: Level " << lvl_index << ": Quality indicator " << " " <<
            stringify_fp_sci(min_quality) << ", minimum angle " << stringify_fp_fix(min_angle) << std::endl;

        ++deque_position;
      }
    }

    return 0;

  }
}; // struct MeshOptimiserRefinementApp

int main(int argc, char* argv[])
{
  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> H2M2D;
  // These mesh types do exist, but no MeshQuality functional is implemented for them
  //typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, 3, Real> S2M3D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, Real> S3M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, 1, Real> H1M1D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, 2, Real> H1M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, 3, Real> H1M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, 3, Real> H2M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, Real> H3M3D;

  // Measure total execution time
  TimeStamp ts_start;

  int rank(0);
  int nprocs(0);

  // initialise
  FEAST::Runtime::initialise(argc, argv, rank, nprocs);
#ifndef SERIAL
  if (rank == 0)
  {
    std::cout << "NUM-PROCS: " << nprocs << std::endl;
  }
#endif

  // Mininum refinement level, parsed from the application config file
  int lvl_min(-1);
  // Maximum refinement level, parsed from the application config file
  int lvl_max(-1);
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

#ifndef SERIAL
  // If we are in parallel mode, we need to synchronise the stream
  synch_stringstream(synchstream_app_config);
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
    String mesh_filename("");
    String chart_filename("");
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

#ifndef SERIAL
  // Synchronise all those streams in parallel mode
  synch_stringstream(synchstream_mesh);
  synch_stringstream(synchstream_chart);
  synch_stringstream(synchstream_meshopt_config);
  synch_stringstream(synchstream_solver_config);
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

  // Get the mesh optimiser key from the application settings
  auto mesh_optimiser_key_p = app_settings_section->query("mesh_optimiser");
  if(!mesh_optimiser_key_p.second)
    throw InternalError(__func__,__FILE__,__LINE__,
    "ApplicationConfig section is missing the mandatory mesh_optimiser entry!");

  // Get the mesh type sting from the parsed mesh so we know with which template parameter to call the application
  const String mesh_type(mesh_file_reader->get_meshtype_string());

  int ret(1);

  // Call the run() method of the appropriate wrapper class
  if(mesh_type == "conformal:hypercube:2:2")
  {
    ret = MeshRefinementOptimiserApp<MemType, DataType, IndexType, H2M2D>::run(
      mesh_optimiser_key_p.first, meshopt_config, solver_config, mesh_file_reader, chart_file_reader,
      lvl_max, lvl_min, write_vtk);
  }
  else if(mesh_type == "conformal:simplex:2:2")
  {
    ret = MeshRefinementOptimiserApp<MemType, DataType, IndexType, S2M2D>::run(
      mesh_optimiser_key_p.first, meshopt_config, solver_config, mesh_file_reader, chart_file_reader,
      lvl_max, lvl_min, write_vtk);
  }
  else
    throw InternalError(__func__,__FILE__,__LINE__,"Unknown mesh type string "+mesh_type);

  delete mesh_file_reader;
  if(chart_file_reader != nullptr)
    delete chart_file_reader;

  // Measure total execution time
  TimeStamp ts_end;

  if(Util::Comm::rank() == 0)
    std::cout << "Total time: " << ts_end.elapsed(ts_start) << std::endl;

  FEAST::Runtime::finalise();
  return ret;
}

static void display_help()
{
  if(Util::Comm::rank() == 0)
  {
    std::cout << "meshopt_refinement-app: This refines a mesh without boundary adaption, then just adapts the finest mesh and uses a mesh optimiser on this"
    << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << " --application_config: Path to application configuration file" << std::endl;
      std::cout << "Optional arguments:" << std::endl;
    std::cout << " --vtk: If this is set, vtk files are written" << std::endl;
    std::cout << " --help: Displays this text" << std::endl;
  }
}
