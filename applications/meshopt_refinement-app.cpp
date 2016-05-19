#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <control/domain/partitioner_domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_control_factory.hpp>

using namespace FEAST;

static void display_help();

#ifndef SERIAL
static void synch_stringstream(std::stringstream& iss, Communicator comm = Communicator(MPI_COMM_WORLD))
{
  Index me(Comm::rank(comm));
  Index size;
  std::string str;

  //bcast size
  if(me == 0)
  {
    str = (iss.str());
    size = str.length();
  }
  // synchronize length
  Comm::bcast(&size, 1, 0, comm);

  //allocate
  char* buf = new char[size + 1];

  //fill
  if(me == 0) //master
  {
    std::strcpy(buf, str.c_str());
  }

  //bcast data
  Comm::bcast(buf, size, 0, comm);

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
    String file_basename(name()+"_n"+stringify(Comm::size()));

    // Minimum number of cells we want to have in each patch
    Index part_min_elems(Comm::size()*4);

    // Create the DomainControl
    DomCtrl dom_ctrl(lvl_max, lvl_min, part_min_elems, mesh_file_reader, chart_file_reader, Geometry::AdaptMode::none);

    // Create the MeshoptControl
    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl, TrafoType>> meshopt_ctrl(nullptr);
    meshopt_ctrl = Control::Meshopt::ControlFactory<Mem_, DT_, IT_, TrafoType>::create_meshopt_control(
      dom_ctrl, meshopt_section_key, meshopt_config, solver_config);

    // Adapt the finest level
    dom_ctrl.get_levels().back()->get_mesh_node()->adapt();

    meshopt_ctrl->mesh_to_buffer();
    auto new_coords = meshopt_ctrl->get_coords().clone();

    meshopt_ctrl->prepare(new_coords);

    // Write initial vtk output
    if(write_vtk)
    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        String vtk_name = String(file_basename+"_pre_lvl_"+stringify((*it)->get_level_index()));

        if(Comm::rank() == 0)
          std::cout << "Writing " << vtk_name << std::endl;

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
        meshopt_ctrl->add_to_vtk_exporter(exporter, deque_position);
        exporter.write(vtk_name, int(Comm::rank()), int(Comm::size()));
        ++deque_position;
      }
    }


    TimeStamp pre_opt;
    meshopt_ctrl->optimise();
    TimeStamp post_opt;

    if(Comm::rank() == 0)
      std::cout << "Solve time: " << post_opt.elapsed(pre_opt) << std::endl;

    // Write initial vtk output
    if(write_vtk)
    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());
        String vtk_name = String(file_basename+"_post_lvl_"+stringify(lvl_index));

        if(Comm::rank() == 0)
          std::cout << "Writing " << vtk_name << std::endl;

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
        meshopt_ctrl->add_to_vtk_exporter(exporter, lvl_index);
        exporter.write(vtk_name, int(Comm::rank()), int(Comm::size()));
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
  //typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, 3, Real> S2M3D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, Real> S3M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, 1, Real> H1M1D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, 2, Real> H1M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, 3, Real> H1M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> H2M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, 3, Real> H2M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, Real> H3M3D;

  // Measure total execution time
  TimeStamp ts_start;

  int rank(0);
  int nprocs(0);

  int lvl_min(-1);
  int lvl_max(-1);
  String mesh_optimiser_key("");

  bool write_vtk(false);

  // initialise
  FEAST::Runtime::initialise(argc, argv, rank, nprocs);
#ifndef SERIAL
  if (rank == 0)
  {
    std::cout << "NUM-PROCS: " << nprocs << std::endl;
  }
#endif

  // Creata a parser for command line arguments.
  SimpleArgParser args(argc, argv);
  args.support("help");
  args.support("level");
  args.support("meshfile");
  args.support("meshopt_config");
  args.support("mesh_optimiser");
  args.support("solver_config");
  args.support("vtk");

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

  args.parse("level", lvl_max, lvl_min);
  args.parse("mesh_optimiser", mesh_optimiser_key);

  // Check if we want to write vtk files
  if(args.check("vtk") >= 0 )
  {
    write_vtk = true;
  }

  std::stringstream synchstream_mesh;
  std::stringstream synchstream_chart;
  std::stringstream synchstream_meshopt_config;
  std::stringstream synchstream_solver_config;

  // Only read files on rank 0
  if(Comm::rank() == 0)
  {
    String mesh_filename("");
    String chart_filename("");
    String meshopt_config_filename("");
    String solver_config_filename("");

    // Read the mesh file to stream
    if(args.check("meshfile") != 1 )
    {
      throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --meshfile");
    }
    else
    {
      args.parse("meshfile",mesh_filename);
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
    if(args.check("chartfile") == 1 )
    {
      args.parse("chartfile",chart_filename);
      std::ifstream ifs(chart_filename);
      if(!ifs.good())
        throw InternalError(__func__, __FILE__, __LINE__, "chart file "+chart_filename+" not found!");
      else
      {
        std::cout << "Reading charts from file " << chart_filename << std::endl;
        synchstream_chart << ifs.rdbuf();
      }
    }
    // Get meshopt_config
    if(args.check("meshopt_config") != 1 )
    {
      std::cout << "You need to specify a mesh optimiser configuration file with --meshopt_config.";
      throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --meshopt_config");
    }
    else
    {
      args.parse("meshopt_config", meshopt_config_filename);
      std::cout << "Reading application configuration from file " << meshopt_config_filename << std::endl;

      std::ifstream ifs(meshopt_config_filename);
      if(!ifs.good())
        throw InternalError(__func__, __FILE__, __LINE__,
        "config file "+meshopt_config_filename+" not found!");
      else
        synchstream_meshopt_config << ifs.rdbuf();
    }

    if(args.check("solver_config") != 1 )
    {
      std::cout << "You need to specify a solver configuration file with --solver_config.";
      throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --solver_config");
    }
    else
    {
      args.parse("solver_config", solver_config_filename);
      std::cout << "Reading application configuration from file " << solver_config_filename << std::endl;

      std::ifstream ifs(solver_config_filename);
      if(!ifs.good())
        throw InternalError(__func__, __FILE__, __LINE__,
        "config file "+solver_config_filename+" not found!");
      else
        synchstream_solver_config << ifs.rdbuf();
    }
  }

#ifndef SERIAL
  // Synchronise all those streams in parallel mode
  synch_stringstream(synchstream_mesh);
  synch_stringstream(synchstream_chart);
  synch_stringstream(synchstream_meshopt_config);
  synch_stringstream(synchstream_solver_config);
#endif

  // Create PropertyMaps and parse the configuration streams
  PropertyMap* meshopt_config = new PropertyMap;
  meshopt_config->parse(synchstream_meshopt_config, true);

  PropertyMap* solver_config = new PropertyMap;
  solver_config->parse(synchstream_solver_config, true);

  // Create a MeshFileReader and parse the mesh stream
  Geometry::MeshFileReader* mesh_file_reader(new Geometry::MeshFileReader(synchstream_mesh));
  mesh_file_reader->read_root_markup();

  // Create a MeshFileReader and parse the chart stream
  Geometry::MeshFileReader* chart_file_reader(nullptr);
  if(!synchstream_chart.str().empty())
    chart_file_reader = new Geometry::MeshFileReader(synchstream_chart);

  // Get the mesh type sting from the parsed mesh so we know with which template parameter to call the application
  const String mesh_type(mesh_file_reader->get_meshtype_string());

  int ret(1);

  // Call the run() method of the appropriate wrapper class
  if(mesh_type == "conformal:hypercube:2:2")
  {
    ret = MeshRefinementOptimiserApp<MemType, DataType, IndexType, H2M2D>::run(
      mesh_optimiser_key, meshopt_config, solver_config, mesh_file_reader, chart_file_reader,
      lvl_max, lvl_min, write_vtk);
  }
  else if(mesh_type == "conformal:simplex:2:2")
  {
    ret = MeshRefinementOptimiserApp<MemType, DataType, IndexType, S2M2D>::run(
      mesh_optimiser_key, meshopt_config, solver_config, mesh_file_reader, chart_file_reader,
      lvl_max, lvl_min, write_vtk);
  }
  else
    throw InternalError(__func__,__FILE__,__LINE__,"Unknown mesh type string "+mesh_type);

  delete mesh_file_reader;
  if(chart_file_reader != nullptr)
    delete chart_file_reader;

  // Measure total execution time
  TimeStamp ts_end;

  if(Comm::rank() == 0)
    std::cout << "Total time: " << ts_end.elapsed(ts_start) << std::endl;

  FEAST::Runtime::finalise();
  return ret;
}

static void display_help()
{
  if(Comm::rank() == 0)
  {
    std::cout << "meshopt_refinement-app: This refines a mesh without boundary adaption, then just adapts the finest mesh and uses a mesh optimiser on this"
    << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << " --meshfile: Path to the mesh file" << std::endl;
    std::cout << " --mesh_optimiser: Which mesh quality functional to use" << std::endl;
    std::cout << " --meshopt_config: Path to the mesh optimiser configuration file" << std::endl;
    std::cout << " --solver_config: Path to the solver configuration file" << std::endl;
    std::cout << " --level [LVL-MAX, LVL-MIN]: Maximum refinement level and coarses level" << std::endl;
      std::cout << "Optional arguments:" << std::endl;
    std::cout << " --vtk: If this is set, vtk files are written" << std::endl;
    std::cout << " --help: Displays this text" << std::endl;
  }
}
