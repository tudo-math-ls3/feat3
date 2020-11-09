// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/geometry/mesh_quality_heuristic.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_control_factory.hpp>

using namespace FEAT;

static void display_help(const Dist::Comm&);
static void read_test_application_config(std::stringstream&, const int);
static void read_test_meshopt_config(std::stringstream&, const int);
static void read_test_solver_config(std::stringstream&, const int);
static void read_test_mesh_file_names(std::deque<String>&, const int);

template<typename Mem_, typename DT_, typename IT_, typename Mesh_>
struct MeshoptRefinementApp
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
  /// The shape type of the mesh's cells
  typedef typename Mesh_::ShapeType ShapeType;

  /// Type for points in the mesh
  typedef Tiny::Vector<DataType, MeshType::world_dim> WorldPoint;

  /// The only transformation available is the standard P1 or Q1 transformation
  typedef Trafo::Standard::Mapping<Mesh_> TrafoType;
  /// FE space for the transformation. The mesh optimization problem is solved on this
  typedef typename Meshopt::Intern::TrafoFE<TrafoType>::Space TrafoFESpace;

  /// The domain level, including trafo and FE space
  typedef Control::Domain::SimpleDomainLevel<Mesh_, TrafoType, TrafoFESpace> DomLvl;
  /// Domain Control Type
  typedef Control::Domain::PartiDomainControl<DomLvl> DomCtrl;

  /**
   * \brief Returns a descriptive string
   *
   * \returns The class name as string
   */
  static String name()
  {
    return "MeshoptRefinementApp";
  }

  /**
   * \brief The routine that does the actual work
   */
  static int run(const SimpleArgParser& args, Dist::Comm& comm, PropertyMap& application_config,
    PropertyMap& meshopt_config, PropertyMap& solver_config, Geometry::MeshFileReader& mesh_file_reader)
  {

    static constexpr int pad_width = 30;

    int ret(0);

    // Mininum refinement level, parsed from the application config file
    int lvl_min(-1);
    // Maximum refinement level, parsed from the application config file
    int lvl_max(-1);
    // Do we want to write vtk files. Read from the command line arguments
    bool write_vtk(false);
    // Do we want to write the optimized mesh to xml. Read from the command line arguments
    bool write_xml(false);
    // Is the application running as a test? Read from the command line arguments
    int test_number(0);

    // Check if we want to write vtk files
    if(args.check("vtk") >= 0 )
    {
      write_vtk = true;
    }

    // Check if we want to write xml
    if(args.check("xml") >= 0 )
    {
      write_xml = true;
    }

    // Check if we are to perform test 1 or test 2, if any
    if( args.check("test") >=0 )
    {
      comm.print("Running in test mode, all other command line arguments and configuration files are ignored.");

      if(args.check("test") > 1)
      {
        XABORTM("Too many options for --test");
      }

      args.parse("test",test_number);
      if(test_number != 1 && test_number != 2)
      {
        XABORTM("Encountered unhandled test number "+stringify(test_number));
      }
    }

    // Get the application settings section
    auto app_settings_section = application_config.query_section("ApplicationSettings");
    XASSERTM(app_settings_section != nullptr,
    "Application config is missing the mandatory ApplicationSettings section!");

    // Get the mesh optimizer key from the application settings
    auto meshoptimizer_key_p = app_settings_section->query("mesh_optimizer");
    XASSERTM(meshoptimizer_key_p.second,
    "ApplicationConfig section is missing the mandatory meshoptimizer entry!");

    // Get the application settings section
    auto domain_control_settings_section = application_config.query_section("DomainControlSettings");
    XASSERTM(domain_control_settings_section != nullptr,
    "DomainControl config is missing the mandatory DomainControlSettings section!");

    // Get the coarse mesh and finest mesh levels from the application settings
    auto lvl_min_p = domain_control_settings_section->query("lvl_min");
    if(lvl_min_p.second)
    {
      lvl_min = std::stoi(lvl_min_p.first);
    }
    else
    {
      lvl_min = 0;
    }

    auto lvl_max_p = domain_control_settings_section->query("lvl_max");
    if(lvl_max_p.second)
    {
      lvl_max = std::stoi(lvl_max_p.first);
    }
    else
    {
      lvl_max = lvl_min;
    }

    // Get the mode for adapting the mesh upon refinement
    Geometry::AdaptMode refinement_adapt_mode(Geometry::AdaptMode::none);
    auto adapt_mode_p = domain_control_settings_section->query("refinement_adapt_mode");
    if(adapt_mode_p.second)
    {
      refinement_adapt_mode << adapt_mode_p.first;
    }

    // Get the mode for adapting the finest mesh
    Geometry::AdaptMode finest_adapt_mode(Geometry::AdaptMode::none);
    adapt_mode_p = domain_control_settings_section->query("finest_adapt_mode");
    if(adapt_mode_p.second)
    {
      finest_adapt_mode << adapt_mode_p.first;
    }

    TimeStamp at;

    // Create domain control
    DomCtrl dom_ctrl(comm, false);
    dom_ctrl.set_adapt_mode(refinement_adapt_mode);
    dom_ctrl.parse_property_map(*domain_control_settings_section);
    dom_ctrl.set_desired_levels(lvl_max, lvl_min);
    dom_ctrl.create(mesh_file_reader);

    // Mesh on the finest level, mainly for computing quality indicators
    const auto& finest_mesh = dom_ctrl.front()->get_mesh();

    // Print level information
    comm.print(name()+" settings:");
    comm.print("LVL-MAX "+stringify(dom_ctrl.max_level_index())
        +" [" +stringify(lvl_max) + "] "
        +"LVL-MIN "+stringify(dom_ctrl.min_level_index())+" [" +stringify(lvl_min) + "]");
    dom_ctrl.print();

    // Create MeshoptControl
    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl>> meshopt_ctrl(nullptr);
    meshopt_ctrl = Control::Meshopt::ControlFactory<Mem_, DT_, IT_>::create_meshopt_control(
      dom_ctrl, meshoptimizer_key_p.first, &meshopt_config, &solver_config);

    String file_basename(name()+"_n"+stringify(comm.size()));

    // Adapt the finest level
    dom_ctrl.set_adapt_mode(finest_adapt_mode);
    dom_ctrl.front()->get_mesh_node()->adapt();

    // Save new coordinates. We need them for calling prepare() to set the initial guess
    meshopt_ctrl->mesh_to_buffer();
    auto new_coords(meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep));

    // Now call prepare
    meshopt_ctrl->prepare(new_coords);

    // For the tests these have to have function global scope
    DT_ qi_min(0);
    DT_ qi_mean(0);
    DataType* qi_cellwise(new DataType[finest_mesh.get_num_entities(MeshType::shape_dim)]);

    DT_ edge_angle(0);
    DataType* edge_angle_cellwise(new DataType[finest_mesh.get_num_entities(MeshType::shape_dim)]);

    DT_ cell_size_defect(0);

    // Write initial vtk output
    if(write_vtk)
    {
      for(size_t l(0); l < dom_ctrl.size_physical(); ++l)
      {
        auto& dom_lvl = dom_ctrl.at(l);

        int lvl_index(dom_lvl->get_level_index());

        String vtk_name = String(file_basename+"_pre_initial_lvl_"+stringify(lvl_index));
        comm.print("Writing "+vtk_name);

        // Compute mesh quality on this level
        meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise, lvl_index);

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_lvl->get_mesh());

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise);
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise);

        meshopt_ctrl->add_to_vtk_exporter(exporter, lvl_index);

        exporter.write(vtk_name, comm.rank(), comm.size());
      }
    }

    // Compute and print quality indicators on the finest level only
    {
      DT_ lambda_min(Math::huge<DT_>());
      DT_ lambda_max(0);
      DT_ vol(0);
      DT_ vol_min(Math::huge<DT_>());
      DT_ vol_max(0);

      cell_size_defect = meshopt_ctrl->compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max, vol);

      // If we did not compute this for the vtk output, we have to do it here
      if(!write_vtk)
      {
        meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);
      }

      String msg("");
      comm.print(msg);

      msg = String("Initial total volume").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(vol);
      comm.print(msg);

      msg = String("Initial QI min/mean").pad_back(pad_width,' ') + String(": ") + stringify_fp_sci(qi_min) + String(" / ") + stringify_fp_sci(qi_mean);
      comm.print(msg);

      msg = String("Initial worst edge angle").pad_back(pad_width, ' ' ) + String(": ") + stringify_fp_fix(edge_angle);
      comm.print(msg);

      msg = String("Initial cell size defect").pad_back(pad_width, ' ' ) + String(": ") + stringify_fp_sci(cell_size_defect);
      comm.print(msg);

      msg = String("Initial lambda min/max").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(lambda_min) + String(" / ") + stringify_fp_sci(lambda_max) ;
      comm.print(msg);

      msg = String("Initial vol fraction min/max").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(vol_min) + " / " + stringify_fp_sci(vol_max);
      comm.print(msg);

      msg = String("");
      comm.print(msg);

    }

    // Check for the hard coded settings for test mode
    if(test_number == 1)
    {
      if( Math::abs(edge_angle - DT_(45)) > Math::sqrt(Math::eps<DataType>()))
      {
        comm.print("FAILED: Initial worst angle should be = "+stringify_fp_fix(45)+ " but is "
            +stringify_fp_fix(edge_angle));
        ret++;
      }
    }
    else if(test_number == 2)
    {
      if( Math::abs(edge_angle - DT_(26.103429982165846)) > Math::sqrt(Math::eps<DataType>()))
      {
        comm.print("FAILED: Initial worst angle should be >= "+stringify_fp_fix(26.103429982165846)+ " but is "
            +stringify_fp_fix(edge_angle));
        ret++;
      }
    }

    // Optimize the mesh
    meshopt_ctrl->optimize();

    // Write output again
    if(write_vtk)
    {
      for(size_t l(0); l < dom_ctrl.size_physical(); ++l)
      {
        auto& dom_lvl = dom_ctrl.at(l);

        int lvl_index(dom_lvl->get_level_index());

        String vtk_name = String(file_basename+"_post_initial_lvl_"+stringify(lvl_index));
        comm.print("Writing "+vtk_name);

        // Compute mesh quality on this level
        meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise, lvl_index);

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_lvl->get_mesh());

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise);
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise);

        meshopt_ctrl->add_to_vtk_exporter(exporter, lvl_index);

        exporter.write(vtk_name, comm.rank(), comm.size());
      }
    }

    // Compute and print quality indicators on the finest level only
    {
      DT_ lambda_min(Math::huge<DT_>());
      DT_ lambda_max(0);
      DT_ vol(0);
      DT_ vol_min(Math::huge<DT_>());
      DT_ vol_max(0);

      cell_size_defect = meshopt_ctrl->compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max, vol);

      // If we did not compute this for the vtk output, we have to do it here
      if(!write_vtk)
      {
        meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);
      }

      String msg("");
      comm.print(msg);

      msg = String("Final total volume").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(vol);
      comm.print(msg);

      msg = String("Final QI min/mean").pad_back(pad_width,' ') + String(": ") + stringify_fp_sci(qi_min) + String(" / ") + stringify_fp_sci(qi_mean);
      comm.print(msg);

      msg = String("Final worst edge angle").pad_back(pad_width, ' ' ) + String(": ") + stringify_fp_fix(edge_angle);
      comm.print(msg);

      msg = String("Final cell size defect").pad_back(pad_width, ' ' ) + String(": ") + stringify_fp_sci(cell_size_defect);
      comm.print(msg);

      msg = String("Final lambda min/max").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(lambda_min) + String(" / ") + stringify_fp_sci(lambda_max) ;
      comm.print(msg);

      msg = String("Final vol fraction min/max").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(vol_min) + " / " + stringify_fp_sci(vol_max);
      comm.print(msg);

      msg = String("");
      comm.print(msg);

    }

    if(write_xml)
    {

      if(comm.size() != 1)
      {
        XABORTM("Writing the mesh to xml is only possible with 1 process!");
      }


      for(size_t l(0); l < dom_ctrl.size_physical(); ++l)
      {
        auto& dom_lvl = dom_ctrl.at(l);

        int lvl_index(dom_lvl->get_level_index());
        String xml_name(file_basename+"_"+stringify(lvl_index)+".xml");

        comm.print("Writing "+xml_name);

        std::ofstream ofs(xml_name);

        ofs << std::setprecision(16);

        if(!ofs.is_open() || !ofs.good())
        {
          throw FileError(String("Failed to open output file: ") + xml_name);
        }

        Geometry::MeshFileWriter mesh_writer(ofs, false);
        mesh_writer.write_mesh(dom_lvl->get_mesh());

      }
    }

    // Check for the hard coded settings for test mode
    if(test_number == 1)
    {
      if(edge_angle < DT_(55.1))
      {
        comm.print("FAILED: Post Initial worst angle should be >= "+stringify_fp_fix(55.1)+
            " but is "+stringify_fp_fix(edge_angle));
        ret++;
      }
      if(qi_min < DT_(4.12e-1))
      {
        comm.print("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(4.12e-1)+
            " but is "+stringify_fp_fix(qi_min)+"\n");
        ret++;
      }
      if(cell_size_defect > DT_(2.6e-1))
      {
        comm.print("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(2.6e-1)+
            " but is "+stringify_fp_fix(cell_size_defect)+"\n");
        ret++;
      }
    }
    else if(test_number == 2)
    {
      if(edge_angle < DT_(22))
      {
        comm.print("FAILED: Post Initial worst angle should be >= "+stringify_fp_fix(22)+
            " but is "+stringify_fp_fix(edge_angle)+"\n");
        ret++;
      }
      if(qi_min < DT_(6.4e-1))
      {
        comm.print("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(6.4e-1)+
            " but is "+stringify_fp_fix(qi_min)+"\n");
        ret++;
      }
      if(cell_size_defect > DT_(1.8e-1))
      {
        comm.print("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(1.8e-1)+
            " but is "+stringify_fp_fix(cell_size_defect)+"\n");
        ret++;
      }
    }

    comm.print(meshopt_ctrl->info());

    // Print success or not
    if(ret == 0)
    {
      comm.print("\nFinished successfully!");
    }
    else
    {
      String msg("\nFAILED: "+stringify(ret) + " check");
      if(ret > 1)
      {
        msg+="s";
      }
      comm.print(msg);
    }

    TimeStamp bt;
    comm.print("Elapsed time: "+stringify(bt.elapsed(at)));

    delete[] qi_cellwise;
    delete[] edge_angle_cellwise;

    return ret;

  }
}; // struct MeshoptRefinementApp

int run_app(int argc, char* argv[])
{
  // Even though this *looks* configurable, it is not: All HyperelasticityFunctionals are implemented for Mem::Main
  // only
  typedef Mem::Main MemType;
  // Floating point type
  typedef double DataType;
  // Index type
  typedef Index IndexType;

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, Real> S2M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, Real> H1M1D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, Real> H1M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, Real> H1M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, Real> H2M3D;

  // create world communicator
  Dist::Comm comm(Dist::Comm::world());

  // print number of processes
  comm.print("Number of Processes: " + stringify(comm.size()));

  // Filenames to read the mesh from, parsed from the application config file
  std::deque<String> mesh_files;
  String mesh_path;
  // String containing the mesh type, read from the header of the mesh file
  String mesh_type("");
  // Is the application running as a test? Read from the command line arguments
  int test_number(0);

  // Streams for synchronising information read from files
  std::stringstream synchstream_app_config;
  std::stringstream synchstream_meshopt_config;
  std::stringstream synchstream_solver_config;

  // Create a parser for command line arguments.
  SimpleArgParser args(argc, argv);
  args.support("application_config");
  args.support("mesh-path");
  args.support("help");
  args.support("test");
  args.support("vtk");
  args.support("xml");

  if( args.check("help") > -1 || args.num_args()==1)
  {
    display_help(comm);
  }

  // Get unsupported command line arguments
  std::deque<std::pair<int,String> > unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
  }

  if( args.check("test") >=0 )
  {
    if(args.check("test") > 1)
    {
      XABORTM("Too many options for --test");
    }

    args.parse("test",test_number);
    if(test_number != 1 && test_number != 2)
    {
      XABORTM("Encountered unhandled test number "+stringify(test_number));
    }
  }

  // parse mesh directory
  if(args.check("mesh-path") >= 0)
  {
    if(args.parse("mesh-path", mesh_path) != 1)
    {
      XABORTM("Failed to parse mesh directory");
    }
  }

  // create a mesh file reader
  Geometry::MeshFileReader mesh_file_reader;

  // Application settings, has to be created here because it gets filled differently according to test
  PropertyMap application_config;
  PropertyMap meshopt_config;
  PropertyMap solver_config;

  // If we are not in test mode, parse command line arguments, read files, synchronize streams
  if(test_number == 0)
  {
    // Read the application config file, required
    String application_config_filename("");
    // Check and parse --application_config
    if(args.check("application_config") != 1 )
    {
      comm.print("You need to specify a application configuration file with --application_config.");
      XABORTM("Invalid option for --application_config");
    }
    else
    {
      args.parse("application_config", application_config_filename);
      comm.print("Reading application configuration from file "+application_config_filename);

      DistFileIO::read_common(synchstream_app_config, application_config_filename);
    }

    // Parse the application config from the (synchronized) stream
    application_config.read(synchstream_app_config, true);

    // Get the application settings section
    auto app_settings_section = application_config.query_section("ApplicationSettings");
    XASSERTM(app_settings_section != nullptr,
    "Application config is missing the mandatory ApplicationSettings section!");

    auto mesh_files_p = app_settings_section->query("mesh_files");
    mesh_files = mesh_files_p.first.split_by_whitespaces();

    // Read configuration for mesh optimization to stream
    auto meshopt_config_filename_p = app_settings_section->query("meshopt_config_file");

    XASSERTM(meshopt_config_filename_p.second,
    "ApplicationConfig section is missing the mandatory meshopt_config_file entry!");

    comm.print("Reading mesh optimization config from file "+meshopt_config_filename_p.first);
    DistFileIO::read_common(synchstream_meshopt_config, meshopt_config_filename_p.first);
    meshopt_config.read(synchstream_meshopt_config, true);

    // Read solver configuration to stream
    auto solver_config_filename_p = app_settings_section->query("solver_config_file");

    XASSERTM(solver_config_filename_p.second,
    "ApplicationConfig section is missing the mandatory solver_config_file entry!");

    comm.print("Reading solver config from file "+solver_config_filename_p.first);
    DistFileIO::read_common(synchstream_solver_config, solver_config_filename_p.first);
    solver_config.read(synchstream_solver_config, true);
  }
  // If we are in test mode, all streams are filled by the hard coded stuff below
  else
  {
    read_test_application_config(synchstream_app_config, test_number);
    application_config.read(synchstream_app_config, true);

    read_test_meshopt_config(synchstream_meshopt_config, test_number);
    meshopt_config.read(synchstream_meshopt_config, true);

    read_test_solver_config(synchstream_solver_config, test_number);
    solver_config.read(synchstream_solver_config, true);

    read_test_mesh_file_names(mesh_files, test_number);
  }

  // Now we have all configurations in the corresponding streams and know the mesh file names

  // Read all mesh files
  mesh_file_reader.add_mesh_files(comm, mesh_files, mesh_path);

  int ret(1);

  // Get the mesh type sting from the parsed mesh so we know with which template parameter to call the application
  mesh_file_reader.read_root_markup();
  mesh_type = mesh_file_reader.get_meshtype_string();

  // Call the appropriate class' run() function
  if(mesh_type == "conformal:hypercube:2:2")
  {
    ret = MeshoptRefinementApp<MemType, DataType, IndexType, H2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  //else if(mesh_type == "conformal:hypercube:3:3")
  //{
  //  ret = MeshoptRefinementApp<MemType, DataType, IndexType, H3M3D>::run(
  //    args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  //}
  else if(mesh_type == "conformal:simplex:2:2")
  {
    ret = MeshoptRefinementApp<MemType, DataType, IndexType, S2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  //else if(mesh_type == "conformal:simplex:3:3")
  //{
  //  ret = MeshoptRefinementApp<MemType, DataType, IndexType, S3M3D>::run(
  //    args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  //}
  else
  {
    XABORTM("Unhandled mesh type "+mesh_type);
  }

  return ret;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialize(argc, argv);
  int ret = run_app(argc, argv);
  FEAT::Runtime::finalize();
  return ret;
}

static void read_test_application_config(std::stringstream& iss, const int test_number)
{
  if(test_number == 1)
  {
    iss << "[ApplicationSettings]" << std::endl;
    iss << "mesh_optimizer = HyperelasticityDefault" << std::endl;
    iss << "solver_config_file = ./solver_config.ini" << std::endl;

    iss << "[DomainControlSettings]" << std::endl;
    iss << "refinement_adapt_mode = none" << std::endl;
    iss << "finest_adapt_mode = chart" << std::endl;
    iss << "lvl_min = 1" << std::endl;
    iss << "lvl_max = 3" << std::endl;
  }
  else if(test_number == 2)
  {
    iss << "[ApplicationSettings]" << std::endl;
    iss << "mesh_optimizer = HyperelasticityDefault" << std::endl;
    iss << "solver_config_file = ./solver_config.ini" << std::endl;

    iss << "[DomainControlSettings]" << std::endl;
    iss << "refinement_adapt_mode = none" << std::endl;
    iss << "finest_adapt_mode = chart" << std::endl;
    iss << "lvl_min = 1" << std::endl;
    iss << "lvl_max = 3" << std::endl;
  }
  else
  {
    XABORTM("Unknown test number: "+stringify(test_number));
  }
}

static void read_test_meshopt_config(std::stringstream& iss, const int test_number)
{
  if(test_number == 1)
  {
    iss << "[HyperElasticityDefault]" << std::endl;
    iss << "type = Hyperelasticity" << std::endl;
    iss << "config_section = HyperelasticityDefaultParameters" << std::endl;
    iss << "dirichlet_boundaries = bnd:o" << std::endl;

    iss << "[HyperelasticityDefaultParameters]" << std::endl;
    iss << "global_functional = HyperelasticityFunctional" << std::endl;
    iss << "cell_functional = RumpfFunctional" << std::endl;
    iss << "solver_config = NLCG" << std::endl;
    iss << "fac_norm = 1.0" << std::endl;
    iss << "fac_det = 1.0" << std::endl;
    iss << "fac_cof = 0.0" << std::endl;
    iss << "fac_reg = 4e-8" << std::endl;
    iss << "exponent_det = 1" << std::endl;
    iss << "scale_computation = once_uniform" << std::endl;
  }
  else if(test_number == 2)
  {
    iss << "[HyperElasticityDefault]" << std::endl;
    iss << "type = Hyperelasticity" << std::endl;
    iss << "config_section = HyperelasticityDefaultParameters" << std::endl;
    iss << "slip_boundaries = bnd:o" << std::endl;

    iss << "[HyperelasticityDefaultParameters]" << std::endl;
    iss << "global_functional = HyperelasticityFunctional" << std::endl;
    iss << "cell_functional = RumpfFunctional" << std::endl;
    iss << "solver_config = NLCG" << std::endl;
    iss << "fac_norm = 1.0" << std::endl;
    iss << "fac_det = 1.0" << std::endl;
    iss << "fac_cof = 0.0" << std::endl;
    iss << "fac_reg = 2e-8" << std::endl;
    iss << "exponent_det = 2" << std::endl;
    iss << "scale_computation = current_uniform" << std::endl;
  }
  else
  {
    XABORTM("Unknown test number "+stringify(test_number));
  }
}

static void read_test_solver_config(std::stringstream& iss, const int test_number)
{
  if(test_number == 1)
  {
    iss << "[NLCG]" << std::endl;
    iss << "type = NLCG" << std::endl;
    iss << "precon = DuDvPrecon" << std::endl;
    iss << "plot_mode = all" << std::endl;
    iss << "tol_rel = 1e-8" << std::endl;
    iss << "max_iter = 1000" << std::endl;
    iss << "linesearch = MQCLinesearch" << std::endl;
    iss << "direction_update = DYHSHybrid" << std::endl;
    iss << "keep_iterates = 0" << std::endl;

    iss << "[MQCLinesearch]" << std::endl;
    iss << "type = MQCLinesearch" << std::endl;
    iss << "plot_mode = none" << std::endl;
    iss << "max_iter = 20" << std::endl;
    iss << "tol_decrease = 1e-3" << std::endl;
    iss << "tol_curvature = 0.3" << std::endl;
    iss << "keep_iterates = 0" << std::endl;

    iss << "[DuDvPrecon]" << std::endl;
    iss << "type = DuDvPrecon" << std::endl;
    iss << "dirichlet_boundaries = bnd:o" << std::endl;
    iss << "fixed_reference_domain = 1" << std::endl;
    iss << "linear_solver = PCG-MG" << std::endl;

    iss << "[PCG-JAC]" << std::endl;
    iss << "type = pcg" << std::endl;
    iss << "max_iter = 10" << std::endl;
    iss << "tol_rel = 1e-8" << std::endl;
    iss << "precon = jac" << std::endl;

    iss << "[PCG-MG]" << std::endl;
    iss << "type = pcg" << std::endl;
    iss << "max_iter = 2" << std::endl;
    iss << "tol_rel = 1e-8" << std::endl;
    iss << "plot_mode = summary" << std::endl;
    iss << "precon = MG1" << std::endl;

    iss << "[cg]" << std::endl;
    iss << "type = pcg" << std::endl;
    iss << "max_iter = 4" << std::endl;
    iss << "min_iter = 4" << std::endl;

    iss << "[rich]" << std::endl;
    iss << "type = richardson" << std::endl;
    iss << "max_iter = 4" << std::endl;
    iss << "min_iter = 4" << std::endl;
    iss << "precon = jac" << std::endl;

    iss << "[jac]" << std::endl;
    iss << "type = jacobi" << std::endl;
    iss << "omega = 0.5" << std::endl;

    iss << "[MG1]" << std::endl;
    iss << "type = mg" << std::endl;
    iss << "hierarchy = s:rich-c:pcg" << std::endl;
    iss << "lvl_min = -1" << std::endl;
    iss << "lvl_max = 0" << std::endl;
    iss << "cycle = w" << std::endl;

    iss << "[s:rich-c:pcg]" << std::endl;
    iss << "type = hierarchy" << std::endl;
    iss << "smoother = rich" << std::endl;
    iss << "coarse = PCG-JAC" << std::endl;

  }
  else if(test_number == 2)
  {
    iss << "[NLCG]" << std::endl;
    iss << "type = NLCG" << std::endl;
    iss << "precon = DuDvPrecon" << std::endl;
    iss << "plot_mode = all" << std::endl;
    iss << "tol_rel = 1e-8" << std::endl;
    iss << "max_iter = 100" << std::endl;
    iss << "linesearch = MQCLinesearch" << std::endl;
    iss << "direction_update = DYHSHybrid" << std::endl;
    iss << "keep_iterates = 0" << std::endl;

    iss << "[DuDvPrecon]" << std::endl;
    iss << "type = DuDvPrecon" << std::endl;
    iss << "slip_boundaries = bnd:o" << std::endl;
    iss << "fixed_reference_domain = 1" << std::endl;
    iss << "linear_solver = PCG-MG" << std::endl;

    iss << "[PCG-JAC]" << std::endl;
    iss << "type = pcg" << std::endl;
    iss << "max_iter = 100" << std::endl;
    iss << "tol_rel = 1e-8" << std::endl;
    iss << "precon = jac" << std::endl;

    iss << "[PCG-MG]" << std::endl;
    iss << "type = pcg" << std::endl;
    iss << "max_iter = 10" << std::endl;
    iss << "tol_rel = 1e-3" << std::endl;
    iss << "plot_mode = summary" << std::endl;
    iss << "precon = MG1" << std::endl;

    iss << "[MQCLinesearch]" << std::endl;
    iss << "type = MQCLinesearch" << std::endl;
    iss << "plot_mode = none" << std::endl;
    iss << "max_iter = 20" << std::endl;
    iss << "tol_decrease = 1e-3" << std::endl;
    iss << "tol_curvature = 0.3" << std::endl;
    iss << "keep_iterates = 0" << std::endl;

    iss << "[rich]" << std::endl;
    iss << "type = richardson" << std::endl;
    iss << "max_iter = 4" << std::endl;
    iss << "min_iter = 4" << std::endl;
    iss << "precon = jac" << std::endl;

    iss << "[cg]" << std::endl;
    iss << "type = pcg" << std::endl;
    iss << "max_iter = 4" << std::endl;
    iss << "min_iter = 4" << std::endl;

    iss << "[jac]" << std::endl;
    iss << "type = jacobi" << std::endl;
    iss << "omega = 0.7" << std::endl;

    iss << "[MG1]" << std::endl;
    iss << "type = mg" << std::endl;
    iss << "hierarchy = s:cg-c:pcg" << std::endl;
    iss << "lvl_min = -1" << std::endl;
    iss << "lvl_max = 0" << std::endl;
    iss << "cycle = v" << std::endl;

    iss << "[s:cg-c:pcg]" << std::endl;
    iss << "type = hierarchy" << std::endl;
    iss << "smoother = cg" << std::endl;
    iss << "coarse = PCG-JAC" << std::endl;

    iss << "[s:rich-c:pcg]" << std::endl;
    iss << "type = hierarchy" << std::endl;
    iss << "smoother = rich" << std::endl;
    iss << "coarse = PCG-JAC" << std::endl;

  }
  else
  {
    XABORTM("Unknown test number "+stringify(test_number));
  }
}

static void read_test_mesh_file_names(std::deque<String>& mesh_files, const int test_number)
{
  if(test_number == 1)
  {
    mesh_files.push_back("unit-circle-quad.xml");
  }
  else if(test_number == 2)
  {
    mesh_files.push_back("unit-circle-tria.xml");
  }
  else
  {
    XABORTM("Encountered unhandled test "+stringify(test_number));
  }
}

static void display_help(const Dist::Comm& comm)
{
  if(comm.rank() == 0)
  {
    std::cout << "meshopt_refinement-app: This refines a mesh without boundary adaption, then just adapts the" <<
      " finest mesh and uses a mesh optimizer on this" << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << " --application_config: Path to the application configuration file" << std::endl;
    std::cout << " --mesh-path: Path to the mesh directory" << std::endl;
    std::cout << "Optional arguments:" << std::endl;
    std::cout << " --test: Run as a test. Ignores configuration files and uses hard coded settings." << std::endl;
    std::cout << " --test [1 or 2]: Run as a test. Ignores configuration files and uses hard coded settings. " <<
      "Test 1 is quadrilateral cells, test 2 is triangular cells" << std::endl;
    std::cout << " --vtk: If this is set, vtk files are written" << std::endl;
    std::cout << " --help: Displays this text" << std::endl;
  }
}
