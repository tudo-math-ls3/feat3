#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_quality_heuristic.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_control_factory.hpp>

using namespace FEAT;

static void display_help(const Dist::Comm&);
static void read_test_application_config(std::stringstream&);
static void read_test_meshopt_config(std::stringstream&, const int);
static void read_test_solver_config(std::stringstream&);
static void read_test_mesh_file_names(std::deque<String>&, const int);

template<typename Mem_, typename DT_, typename IT_, typename Mesh_>
struct MeshoptRAdaptApp
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
  /// FE space for the transformation. The mesh optimisation problem is solved on this
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
    return "MeshoptRAdaptApp";
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
    // Timestep size, parsed from the application config file
    DataType delta_t(0);
    // End time, parsed from the application config file
    DataType t_end(0);
    // Do we want to write vtk files. Read from the command line arguments
    bool write_vtk(false);
    // If write_vtk is set, we write out every vtk_freq time steps
    Index vtk_freq(1);
    // Is the application running as a test? Read from the command line arguments
    int test_number(0);

    // Check if we want to write vtk files and at what frequency
    if(args.check("vtk") >= 0 )
    {
      write_vtk = true;

      if(args.check("vtk") > 1)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "Too many options for --vtk");
      }

      args.parse("vtk",vtk_freq);
    }

    // Check if we are to perform test 1 or test 2, if any
    if( args.check("test") >=0 )
    {
      if(args.check("test") > 1)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "Too many options for --test");
      }

      args.parse("test",test_number);
      if(test_number != 1 && test_number != 2)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "Encountered unhandled test number "+stringify(test_number));
      }
    }

    // Get the application settings section
    auto app_settings_section = application_config.query_section("ApplicationSettings");
    XASSERTM(app_settings_section != nullptr,
    "Application config is missing the mandatory ApplicationSettings section!\n");

    // Get timestep size
    auto delta_t_p = app_settings_section->query("delta_t");
    XASSERTM(delta_t_p.second, "ApplicationConfig section is missing the mandatory delta_t entry!\n");
    delta_t = std::stod(delta_t_p.first);
    XASSERT(delta_t > DataType(0));

    // Get end time
    auto t_end_p = app_settings_section->query("t_end");
    XASSERTM(delta_t_p.second, "ApplicationConfig section is missing the mandatory t_end entry!");
    t_end = std::stod(t_end_p.first);
    XASSERT(t_end >= DataType(0));

    // Get the mesh optimiser key from the application settings
    auto meshoptimiser_key_p = app_settings_section->query("mesh_optimiser");
    XASSERTM(meshoptimiser_key_p.second,
    "ApplicationConfig section is missing the mandatory meshoptimiser entry!");

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

    TimeStamp at;

    // Create domain control
    DomCtrl dom_ctrl(comm);
    dom_ctrl.read_mesh(mesh_file_reader);
    dom_ctrl.parse_property_map(domain_control_settings_section);
    dom_ctrl.create_partition();
    dom_ctrl.create_hierarchy(lvl_max, lvl_min);

    // Mesh on the finest level, mainly for computing quality indicators
    const auto& finest_mesh = dom_ctrl.front()->get_mesh();

    // Print level information
    comm.print(name()+" settings:");
    comm.print("LVL-MAX "+stringify(dom_ctrl.max_level_index())
        +" [" +stringify(lvl_max) + "] "
        +"LVL-MIN "+stringify(dom_ctrl.min_level_index())+" [" +stringify(lvl_min) + "]");
    comm.print("Timestep size: "+stringify_fp_fix(delta_t)+", end time: "+ stringify_fp_fix(t_end));
    dom_ctrl.print();

    // Create MeshoptControl
    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl>> meshopt_ctrl(nullptr);
    meshopt_ctrl = Control::Meshopt::ControlFactory<Mem_, DT_, IT_>::create_meshopt_control(
      dom_ctrl, meshoptimiser_key_p.first, &meshopt_config, &solver_config);

    String file_basename(name()+"_n"+stringify(comm.size()));

    // Copy the vertex coordinates to the buffer and get them via get_coords()
    meshopt_ctrl->mesh_to_buffer();
    // A copy of the old vertex coordinates is kept here
    auto old_coords = meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep);
    auto new_coords = meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep);

    // Prepare the functional
    meshopt_ctrl->prepare(old_coords);

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
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise, lvl_index);

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
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);
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
      if( Math::abs(edge_angle - DT_(90)) > Math::sqrt(Math::eps<DataType>()))
      {
        comm.print("FAILED: Initial worst angle should be = "+stringify_fp_fix(90)+ " but is "
            + stringify_fp_fix(edge_angle));
        ++ret;
      }
    }
    else if(test_number == 2)
    {
      if( Math::abs(edge_angle - DT_(45)) > Math::sqrt(Math::eps<DataType>()))
      {
        comm.print("FAILED: Initial worst angle should be >= "+stringify_fp_fix(45)+ " but is "
            + stringify_fp_fix(edge_angle));
        ++ret;
      }
    }

    // Optimise the mesh
    meshopt_ctrl->optimise();

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
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise, lvl_index);

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

      dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);

      String msg("");
      comm.print(msg);

      msg = String("Optimised total volume").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(vol);
      comm.print(msg);

      msg = String("Optimised QI min/mean").pad_back(pad_width,' ') + String(": ") + stringify_fp_sci(qi_min) + String(" / ") + stringify_fp_sci(qi_mean);
      comm.print(msg);

      msg = String("Optimised worst edge angle").pad_back(pad_width, ' ' ) + String(": ") + stringify_fp_fix(edge_angle);
      comm.print(msg);

      msg = String("Optimised cell size defect").pad_back(pad_width, ' ' ) + String(": ") + stringify_fp_sci(cell_size_defect);
      comm.print(msg);

      msg = String("Optimised lambda min/max").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(lambda_min) + String(" / ") + stringify_fp_sci(lambda_max) ;
      comm.print(msg);

      msg = String("Optimised vol fraction min/max").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(vol_min) + " / " + stringify_fp_sci(vol_max);
      comm.print(msg);

      msg = String("");
      comm.print(msg);

    }

    // Check for the hard coded settings for test mode
    if(test_number == 1)
    {
      if(edge_angle < DT_(50))
      {
        comm.print("FAILED: Post Initial worst angle should be >= "+stringify_fp_fix(50)+
            " but is "+stringify_fp_fix(edge_angle));
        ++ret;
      }
      if(qi_min < DT_(3.4e-1))
      {
        comm.print("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(3.4e-1)+
            " but is "+stringify_fp_fix(qi_min));
        ++ret;
      }
      if(cell_size_defect > DT_(4.6e-2))
      {
        comm.print("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(4.6e-2)+
            " but is "+stringify_fp_fix(cell_size_defect));
        ++ret;
      }
    }
    else if(test_number == 2)
    {
      if(edge_angle < DT_(24))
      {
        comm.print("FAILED: Post Initial worst angle should be >= "+stringify_fp_fix(24)+
            " but is "+stringify_fp_fix(edge_angle));
        ++ret;
      }
      if(qi_min < DT_(5.8e-1))
      {
        comm.print("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(5.8e-1)+
            " but is "+stringify_fp_fix(qi_min));
        ++ret;
      }
      if(cell_size_defect > DT_(8.4e-2))
      {
        comm.print("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(8.4e-2)+
            " but is "+stringify_fp_fix(cell_size_defect));
        ++ret;
      }
    }

    // Initial time
    DataType time(0);
    // Counter for timesteps
    Index n(0);

    // The mesh velocity is 1/delta_t*(coords_new - coords_old) and computed in each time step
    auto mesh_velocity = meshopt_ctrl->get_coords().clone();

    // This is the centre reference point
    WorldPoint midpoint(DataType(0));
    midpoint(0) = DataType(0.25) *(DataType(2) + Math::cos(time));
    midpoint(1) = DataType(0.25) *(DataType(2) + Math::sin(DataType(3)*time));

    WorldPoint rotation_centre(DataType(0.5));
    DataType rotation_speed(DataType(2)*Math::pi<DataType>());
    WorldPoint rotation_angles(DataType(0));
    rotation_angles(0) = rotation_speed*delta_t;

    WorldPoint dir(delta_t/DataType(2));

    while(time < t_end)
    {
      n++;
      time+= delta_t;

      bool abort(false);

      WorldPoint old_midpoint(midpoint);

      midpoint(0) = DataType(0.25) *(DataType(2) + Math::cos(time));
      midpoint(1) = DataType(0.25) *(DataType(2) + Math::sin(DataType(3)*time));

      comm.print("Timestep "+stringify(n)+": t = "+stringify_fp_fix(time)+", midpoint = "+stringify(midpoint));

      // Save old vertex coordinates
      meshopt_ctrl->mesh_to_buffer();
      old_coords.copy(meshopt_ctrl->get_coords());

      for(auto& it:dom_ctrl.get_atlas().get_mesh_chart_map())
      {
        if(it.first.find("moving_") != String::npos)
        {
          comm.print(it.first+" by "+stringify(midpoint-old_midpoint));
          it.second->move_by(midpoint-old_midpoint);
        }

        if(it.first.find("pos_merging_") != String::npos)
        {
          comm.print(it.first+" by "+stringify(dir));
          it.second->move_by(dir);
        }

        if(it.first.find("neg_merging_") != String::npos)
        {
          comm.print(it.first+" by "+stringify(DataType(-1)*dir));
          it.second->move_by(DataType(-1)*dir);
        }

        if(it.first.find("rotating_") != String::npos)
        {
          comm.print(it.first+" around "+stringify(rotation_centre)+" by "+stringify_fp_fix(rotation_angles(0)));
          it.second->rotate(rotation_centre, rotation_angles);
        }
      }

      // Get coords for modification
      new_coords.copy(meshopt_ctrl->get_coords());

      // Now prepare the functional
      meshopt_ctrl->prepare(new_coords);

      meshopt_ctrl->optimise();

      // Compute mesh velocity
      {
        mesh_velocity.axpy(old_coords, meshopt_ctrl->get_coords(), DataType(-1));
        mesh_velocity.scale(mesh_velocity, DataType(1)/delta_t);

        // Compute maximum of the mesh velocity
        DataType max_mesh_velocity(0);
        for(IT_ i(0); i < (*mesh_velocity).size(); ++i)
        {
          max_mesh_velocity = Math::max(max_mesh_velocity, (*mesh_velocity)(i).norm_euclid());
        }

        comm.print("");
        String msg = String("max. mesh velocity").pad_back(pad_width, ' ') + ": "
          + stringify_fp_sci(max_mesh_velocity);
        comm.print(msg);
      }

      // Write output on the finest level only
      if(write_vtk && ( (n%vtk_freq == 0) || abort ))
      {
        // Compute mesh quality on the finest
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);

        String vtk_name = String(file_basename+"_post_"+stringify(n));
        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_ctrl.front()->get_mesh());

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise);
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise);

        meshopt_ctrl->add_to_vtk_exporter(exporter, dom_ctrl.max_level_index());

        exporter.write(vtk_name, comm.rank(), comm.size());
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
        if(! (write_vtk && ( (n%vtk_freq == 0) || abort ) ) )
        {
          dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);
        }

        String msg;
        msg = String("Post total volume").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(vol);
        comm.print(msg);

        msg = String("Post QI min/mean").pad_back(pad_width,' ') + String(": ") + stringify_fp_sci(qi_min) + String(" / ") + stringify_fp_sci(qi_mean);
        comm.print(msg);

        msg = String("Post worst edge angle").pad_back(pad_width, ' ' ) + String(": ") + stringify_fp_fix(edge_angle);
        comm.print(msg);

        msg = String("Post cell size defect").pad_back(pad_width, ' ' ) + String(": ") + stringify_fp_sci(cell_size_defect);
        comm.print(msg);

        msg = String("Post lambda min/max").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(lambda_min) + String(" / ") + stringify_fp_sci(lambda_max) ;
        comm.print(msg);

        msg = String("Post vol fraction min/max").pad_back(pad_width, ' ') + String(": ") + stringify_fp_sci(vol_min) + " / " + stringify_fp_sci(vol_max);
        comm.print(msg);

        msg = String("");
        comm.print(msg);

      }

      if(edge_angle < DT_(1))
      {
        comm.print("FAILED: Mesh deteriorated, stopping");
        ++ret;
        abort = true;
      }

      // compress all statistics from the current timestep for further analysis after the solution is finished
      FEAT::Statistics::compress_solver_expressions();

      if(abort)
      {
        break;
      }

    } // time loop

    meshopt_ctrl->print();

    // Write final vtk output
    if(write_vtk)
    {
      for(size_t l(0); l < dom_ctrl.size_physical(); ++l)
      {
        auto& dom_lvl = dom_ctrl.at(l);

        int lvl_index(dom_lvl->get_level_index());

        String vtk_name = String(file_basename+"_final_lvl_"+stringify(lvl_index));
        comm.print("Writing "+vtk_name);

        // Compute mesh quality on this level
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise, lvl_index);

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_lvl->get_mesh());

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise);
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise);

        meshopt_ctrl->add_to_vtk_exporter(exporter, lvl_index);

        exporter.write(vtk_name, comm.rank(), comm.size());
      }
    }


    // Check for the hard coded settings for test mode
    if(test_number == 1)
    {
      if(edge_angle < DT_(50))
      {
        comm.print("FAILED: Final worst angle should be >= "+stringify_fp_fix(50)+
            " but is "+stringify_fp_fix(edge_angle));
        ++ret;
      }
      if(qi_min < DT_(3.45e-1))
      {
        comm.print("FAILED: Final worst shape quality should be >= "+stringify_fp_fix(3.45e-1)+
            " but is "+stringify_fp_fix(qi_min));
        ++ret;
      }
      if(cell_size_defect > DT_(4.6e-2))
      {
        comm.print("FAILED: Final cell size distribution defect should be < "+stringify_fp_fix(4.6e-2)+
            " but is "+stringify_fp_fix(cell_size_defect));
        ++ret;
      }
    }
    else if(test_number == 2)
    {
      if(edge_angle < DT_(28))
      {
        comm.print("FAILED: Final worst angle should be >= "+stringify_fp_fix(28)+
            " but is "+stringify_fp_fix(edge_angle));
        ++ret;
      }
      if(qi_min < DT_(6.5e-1))
      {
        comm.print("FAILED: Final worst shape quality should be >= "+stringify_fp_fix(6.5e-1)+
            " but is "+stringify_fp_fix(qi_min));
        ++ret;
      }
      if(cell_size_defect > DT_(6.6e-2))
      {
        comm.print("FAILED: Final cell size distribution defect should be < "+stringify_fp_fix(6.6e-2)+
            " but is "+stringify_fp_fix(cell_size_defect));
        ++ret;
      }
    }

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
    comm.print("Elapsed time: "+ stringify(bt.elapsed(at)));

    delete[] qi_cellwise;
    delete[] edge_angle_cellwise;

    return ret;

  }
}; // struct MeshoptRAdaptApp

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
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> H2M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, Real> H3M3D;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> S2M2D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, Real> S3M3D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, 3, Real> S2M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, 1, Real> H1M1D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, 2, Real> H1M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, 3, Real> H1M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, 3, Real> H2M3D;

  // create world communicator
  Dist::Comm comm(Dist::Comm::world());
  comm.print("NUM-PROCS: "+stringify(comm.size()));

  // Filenames to read the mesh from, parsed from the application config file
  std::deque<String> mesh_files;
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
      throw InternalError(__func__, __FILE__, __LINE__, "Too many options for --test");
    }

    args.parse("test",test_number);
    if(test_number != 1 && test_number != 2)
    {
      throw InternalError(__func__, __FILE__, __LINE__, "Encountered unhandled test number "+stringify(test_number));
    }
  }

  // create a mesh file reader
  Geometry::MeshFileReader mesh_file_reader;

  // Application settings, has to be created here because it gets filled differently according to test
  PropertyMap application_config;
  PropertyMap meshopt_config;
  PropertyMap solver_config;

  // If we are not in test mode, parse command line arguments, read files, synchronise streams
  if(test_number == 0)
  {
    // Read the application config file, required
    String application_config_filename("");
    // Check and parse --application_config
    if(args.check("application_config") != 1 )
    {
      comm.print("You need to specify a application configuration file with --application_config.");
      throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --application_config");
    }
    else
    {
      args.parse("application_config", application_config_filename);
      comm.print("Reading application configuration from file "+application_config_filename);

      DistFileIO::read_common(synchstream_app_config, application_config_filename);
    }

    // Parse the application config from the (synchronised) stream
    application_config.parse(synchstream_app_config, true);

    // Get the application settings section
    auto app_settings_section = application_config.query_section("ApplicationSettings");
    XASSERTM(app_settings_section != nullptr,
    "Application config is missing the mandatory ApplicationSettings section!");

    auto mesh_files_p = app_settings_section->query("mesh_files");
    mesh_files_p.first.split_by_charset(mesh_files, " ");

    // Read configuration for mesh optimisation to stream
    auto meshopt_config_filename_p = app_settings_section->query("meshopt_config_file");

    XASSERTM(meshopt_config_filename_p.second,
    "ApplicationConfig section is missing the mandatory meshopt_config_file entry!");

    comm.print("Reading mesh optimisation config from file "+meshopt_config_filename_p.first);
    DistFileIO::read_common(synchstream_meshopt_config, meshopt_config_filename_p.first);
    meshopt_config.parse(synchstream_meshopt_config, true);

    // Read solver configuration to stream
    auto solver_config_filename_p = app_settings_section->query("solver_config_file");

    XASSERTM(solver_config_filename_p.second,
    "ApplicationConfig section is missing the mandatory solver_config_file entry!");

    comm.print("Reading solver config from file "+solver_config_filename_p.first);
    DistFileIO::read_common(synchstream_solver_config, solver_config_filename_p.first);
    solver_config.parse(synchstream_solver_config, true);
  }
  // If we are in test mode, all streams are filled by the hard coded stuff below
  else
  {
    read_test_application_config(synchstream_app_config);
    application_config.parse(synchstream_app_config, true);

    read_test_meshopt_config(synchstream_meshopt_config, test_number);
    meshopt_config.parse(synchstream_meshopt_config, true);

    read_test_solver_config(synchstream_solver_config);
    solver_config.parse(synchstream_solver_config, true);

    read_test_mesh_file_names(mesh_files, test_number);
  }

  // Now we have all configurations in the corresponding streams and know the mesh file names

  // Read all mesh files
  std::deque<std::stringstream> mesh_streams(mesh_files.size());
  for(std::size_t i(0); i < mesh_files.size(); ++i)
  {
    // Read the stream
    comm.print("Reading mesh file "+mesh_files.at(i));
    DistFileIO::read_common(mesh_streams.at(i), mesh_files.at(i));

    // Add to mesh reader
    mesh_file_reader.add_stream(mesh_streams.at(i));
  }

  int ret(1);

  // Get the mesh type sting from the parsed mesh so we know with which template parameter to call the application
  mesh_file_reader.read_root_markup();
  mesh_type = mesh_file_reader.get_meshtype_string();

  // Call the appropriate class' run() function
  if(mesh_type == "conformal:hypercube:2:2")
  {
    ret = MeshoptRAdaptApp<MemType, DataType, IndexType, H2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  //else if(mesh_type == "conformal:hypercube:3:3")
  //{
  //  ret = MeshoptRAdaptApp<MemType, DataType, IndexType, H3M3D>::run(
  //    args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  //}
  else if(mesh_type == "conformal:simplex:2:2")
  {
    ret = MeshoptRAdaptApp<MemType, DataType, IndexType, S2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  //else if(mesh_type == "conformal:simplex:3:3")
  //{
  //  ret = MeshoptRAdaptApp<MemType, DataType, IndexType, S3M3D>::run(
  //    args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  //}
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"Unhandled mesh type "+mesh_type);
  }

  return ret;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialise(argc, argv);
  int ret = run_app(argc, argv);
  FEAT::Runtime::finalise();
  return ret;
}

static void display_help(const Dist::Comm& comm)
{
  if(comm.rank() == 0)
  {
    std::cout << "meshopt_r_adapt-app: Chart distance based r-adaptivity and surface alignment"
    << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << " --application_config: Path to the application configuration file" << std::endl;
    std::cout << "Optional arguments:" << std::endl;
    std::cout << " --test [1 or 2]: Run as a test. Ignores configuration files and uses hard coded settings. " <<
      "Test 1 is r-adaptivity, test 2 is surface alignment" << std::endl;
    std::cout << " --vtk [freq]: If this is set, vtk files are written every freq time steps. freq defaults to 1" <<
      std::endl;
    std::cout << " --help: Displays this text" << std::endl;
  }
}

static void read_test_application_config(std::stringstream& iss)
{
  iss << "[ApplicationSettings]" << std::endl;
  //iss << "mesh_file = ./unit-square-tria.xml" << std::endl;
  //iss << "chart_file = ./moving_circle_chart.xml" << std::endl;
  //iss << "meshopt_config_file = ./meshopt_config.ini" << std::endl;
  iss << "mesh_optimiser = HyperelasticityDefault" << std::endl;
  iss << "solver_config_file = ./solver_config.ini" << std::endl;
  iss << "delta_t = 1e-2" << std::endl;
  iss << "t_end = 2e-2" << std::endl;
  iss << "[DomainControlSettings]" << std::endl;
  iss << "parti-type = fallback parmetis" << std::endl;
  iss << "parti-rank-elems = 4" << std::endl;
  iss << "lvl_min = 3" << std::endl;
  iss << "lvl_max = 3" << std::endl;

}

static void read_test_meshopt_config(std::stringstream& iss, const int test)
{
  if(test == 1)
  {
    iss << "[HyperElasticityDefault]" << std::endl;
    iss << "type = Hyperelasticity" << std::endl;
    iss << "config_section = HyperelasticityDefaultParameters" << std::endl;
    iss << "dirichlet_boundaries = bnd:b bnd:t bnd:l bnd:r" << std::endl;

    iss << "[HyperelasticityDefaultParameters]" << std::endl;
    iss << "global_functional = HyperelasticityFunctional" << std::endl;
    iss << "cell_functional = RumpfFunctional" << std::endl;
    iss << "solver_config = NLCG" << std::endl;
    iss << "fac_norm = 1e-2" << std::endl;
    iss << "fac_det = 1.0" << std::endl;
    iss << "fac_cof = 0.0" << std::endl;
    iss << "fac_reg = 5e-8" << std::endl;
    iss << "exponent_det = 1" << std::endl;
    iss << "scale_computation = iter_concentration" << std::endl;
    iss << "conc_function = OuterDist" << std::endl;

    iss << "[OuterDist]" << std::endl;
    iss << "type = ChartDistance" << std::endl;
    iss << "chart_list = moving_circle" << std::endl;
    iss << "operation = min" << std::endl;
    iss << "function_type = PowOfDist" << std::endl;
    iss << "minval = 2e-2" << std::endl;
    iss << "exponent = 0.5" << std::endl;
    iss << "use_derivative = 1" << std::endl;
  }
  else if(test == 2)
  {
    iss << "[HyperElasticityDefault]" << std::endl;
    iss << "type = Hyperelasticity" << std::endl;
    iss << "config_section = HyperelasticityDefaultParameters" << std::endl;
    iss << "dirichlet_boundaries = bnd:b bnd:t bnd:l bnd:r" << std::endl;

    iss << "[HyperelasticityDefaultParameters]" << std::endl;
    iss << "global_functional = HyperelasticityFunctional" << std::endl;
    iss << "cell_functional = RumpfFunctional" << std::endl;
    iss << "solver_config = QPenalty" << std::endl;
    iss << "fac_norm = 1.0" << std::endl;
    iss << "fac_det = 1.0" << std::endl;
    iss << "fac_cof = 0.0" << std::endl;
    iss << "fac_reg = 5e-8" << std::endl;
    iss << "exponent_det = 1" << std::endl;
    iss << "scale_computation = once_uniform" << std::endl;
    iss << "conc_function = OuterDist" << std::endl;
    iss << "align_mesh = 1" << std::endl;

    iss << "[OuterDist]" << std::endl;
    iss << "type = ChartDistance" << std::endl;
    iss << "chart_list = moving_circle" << std::endl;
    iss << "operation = min" << std::endl;
    iss << "function_type = default" << std::endl;

  }
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"Unknown test: "+stringify(test));
  }
}

static void read_test_solver_config(std::stringstream& iss)
{
  iss << "[NLCG]" << std::endl;
  iss << "type = NLCG" << std::endl;
  iss << "precon = none" << std::endl;
  iss << "plot = iter" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "max_iter = 500" << std::endl;
  iss << "linesearch = MQCLinesearch" << std::endl;
  iss << "direction_update = DYHSHybrid" << std::endl;
  iss << "keep_iterates = 0" << std::endl;

  iss << "[QPenalty]" << std::endl;
  iss << "type = QPenalty" << std::endl;
  iss << "max_iter = 10" << std::endl;
  iss << "tol_rel = 1e5" << std::endl;
  iss << "tol_abs = 1e-8" << std::endl;
  iss << "initial_penalty_param = 1e4" << std::endl;
  iss << "plot = iter" << std::endl;
  iss << "inner_solver = NLCG" << std::endl;

  iss << "[MQCLinesearch]" << std::endl;
  iss << "type = MQCLinesearch" << std::endl;
  iss << "plot = none" << std::endl;
  iss << "max_iter = 20" << std::endl;
  iss << "tol_decrease = 1e-3" << std::endl;
  iss << "tol_curvature = 0.3" << std::endl;
  iss << "keep_iterates = 0" << std::endl;
}

static void read_test_mesh_file_names(std::deque<String>& mesh_files, const int test_number)
{
  String chart_filename(FEAT_SRC_DIR);
  chart_filename +="/applications/meshopt_r_adapt-app/moving-circle-chart.xml";
  mesh_files.push_back(chart_filename);

  String mesh_filename(FEAT_SRC_DIR);
  if(test_number == 1)
  {
    mesh_filename +="/data/meshes/unit-square-quad.xml";
  }
  else if(test_number == 2)
  {
    mesh_filename +="/data/meshes/unit-square-tria.xml";
  }
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"Encountered unhandled test "+stringify(test_number));
  }

  mesh_files.push_back(mesh_filename);
}
