#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_quality_heuristic.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/mpi_cout.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_control_factory.hpp>

using namespace FEAT;

static void display_help(const Dist::Comm&);
static void read_test_application_config(std::stringstream&, const int);
static void read_test_meshopt_config(std::stringstream&, const int);
static void read_test_solver_config(std::stringstream&);
static void read_test_mesh_file_names(std::deque<String>&, const int);

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
  /// The shape type of the mesh's cells
  typedef typename Mesh_::ShapeType ShapeType;

  /// The only transformation available is the standard P1 or Q1 transformation
  typedef Trafo::Standard::Mapping<Mesh_> TrafoType;

  /// Type for points in the mesh
  typedef Tiny::Vector<DataType, MeshType::world_dim> WorldPoint;

  /// Domain Control Type
  typedef Control::Domain::PartiDomainControl<MeshType> DomCtrl;

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
  static int run(const SimpleArgParser& args, Dist::Comm& comm, PropertyMap* application_config,
    PropertyMap* meshopt_config, PropertyMap* solver_config, Geometry::MeshFileReader& mesh_file_reader)
  {

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
        throw InternalError(__func__, __FILE__, __LINE__, "Too many options for --vtk");

      args.parse("vtk",vtk_freq);
    }

    // Check if we are to perform test 1 or test 2, if any
    if( args.check("test") >=0 )
    {
      Util::mpi_cout("Running in test mode, all other command line arguments and configuration files are ignored.\n");

      if(args.check("test") > 1)
        throw InternalError(__func__, __FILE__, __LINE__, "Too many options for --test");

      args.parse("test",test_number);
      if(test_number != 1 && test_number != 2)
        throw InternalError(__func__, __FILE__, __LINE__, "Encountered unhandled test number "+stringify(test_number));
    }

    // Get the application settings section
    auto app_settings_section = application_config->query_section("ApplicationSettings");
    XASSERTM(app_settings_section != nullptr,
    "Application config is missing the mandatory ApplicationSettings section!");

    // Get timestep size
    auto delta_t_p = app_settings_section->query("delta_t");
    XASSERTM(delta_t_p.second, "ApplicationConfig section is missing the mandatory delta_t entry!");
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

    // Get the centre of rotation, default to (0,0)^T
    auto midpoint_p = app_settings_section->query("midpoint");
    std::deque<String> midpoint_deque;
    if(midpoint_p.second)
    {
      midpoint_p.first.split_by_charset(midpoint_deque," ");
    }

    Tiny::Vector<DataType,2> midpoint(DataType(0));
    if(midpoint_deque.size() > size_t(0))
    {
      XASSERTM(midpoint_deque.size() == size_t(MeshType::world_dim),"midpoint has invalid number of components!");
      for(size_t i(0); i < midpoint_deque.size(); ++i)
      {
        midpoint(int(i)) = DataType(std::stod(midpoint_deque.at(i)));
      }
    }

    // Get the application settings section
    auto domain_control_settings_section = application_config->query_section("DomainControlSettings");
    XASSERTM(domain_control_settings_section != nullptr,
    "DomainControl config is missing the mandatory DomainControlSettings section!");

    // Get the coarse mesh and finest mesh levels from the application settings
    auto lvl_min_p = domain_control_settings_section->query("lvl_min");
    if(!lvl_min_p.second)
      lvl_min = 0;
    else
      lvl_min = std::stoi(lvl_min_p.first);

    auto lvl_max_p = domain_control_settings_section->query("lvl_max");
    if(!lvl_max_p.second)
      lvl_max = lvl_min;
    else
      lvl_max = std::stoi(lvl_max_p.first);

    TimeStamp at;

    DomCtrl dom_ctrl(comm);
    dom_ctrl.read_mesh(mesh_file_reader);
    dom_ctrl.parse_property_map(domain_control_settings_section);
    dom_ctrl.create_partition();
    dom_ctrl.create_hierarchy(lvl_max, lvl_min);

    // Mesh on the finest level, mainly for computing quality indicators
    const auto& finest_mesh = dom_ctrl.get_levels().back()->get_mesh();

    Index ncells(dom_ctrl.get_levels().back()->get_mesh().get_num_entities(MeshType::shape_dim));
    comm.allreduce(&ncells, &ncells, std::size_t(1), Dist::op_sum);

    // Print level information
    if(comm.rank() == 0)
    {
      std::cout << name() << "settings:" << std::endl;
      std::cout << "Timestep size: " << stringify_fp_fix(delta_t) << ", end time: " <<
        stringify_fp_fix(t_end) << std::endl;
      std::cout << "LVL-MAX: " <<
        dom_ctrl.get_levels().back()->get_level_index() << " [" << lvl_max << "]";
      std::cout << " LVL-MIN: " <<
        dom_ctrl.get_levels().front()->get_level_index() << " [" << lvl_min << "]" << std::endl;
      std::cout << "Cells: " << ncells << std::endl;
    }

    // Create MeshoptControl
    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl, TrafoType>> meshopt_ctrl(nullptr);
    meshopt_ctrl = Control::Meshopt::ControlFactory<Mem_, DT_, IT_, TrafoType>::create_meshopt_control(
      dom_ctrl, meshoptimiser_key_p.first, meshopt_config, solver_config);

    String file_basename(name()+"_n"+stringify(comm.size()));

    // Copy the vertex coordinates to the buffer and get them via get_coords()
    meshopt_ctrl->mesh_to_buffer();
    // A copy of the old vertex coordinates is kept here
    auto old_coords = meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep);
    auto new_coords = meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep);

    // Prepare the functional
    meshopt_ctrl->prepare(old_coords);

    // Write initial vtk output
    if(write_vtk)
    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());

        String vtk_name = String(file_basename+"_pre_lvl_"+stringify(lvl_index));
        Util::mpi_cout("Writing "+vtk_name+"\n");

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
        meshopt_ctrl->add_to_vtk_exporter(exporter, deque_position);
        exporter.write(vtk_name, comm.rank(), comm.size());

        ++deque_position;
      }

    }

    // For the tests these have to have function global scope
    DT_ qual_min(0);
    DT_ qual_sum(0);
    DT_ worst_angle(0);
    DT_ cell_size_defect(0);

    // Compute quality indicators
    {
      Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::compute(qual_min, qual_sum,
      finest_mesh.template get_index_set<MeshType::shape_dim, 0>(), finest_mesh.get_vertex_set());

      worst_angle = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::angle(
        finest_mesh.template get_index_set<MeshType::shape_dim, 0>(), finest_mesh.get_vertex_set());

      comm.allreduce(&qual_min, &qual_min, std::size_t(1), Dist::op_min);
      comm.allreduce(&qual_sum, &qual_sum, std::size_t(1), Dist::op_sum);
      comm.allreduce(&worst_angle, &worst_angle, std::size_t(1), Dist::op_min);

      DataType qual_avg(qual_sum/DataType(ncells));

      DT_ lambda_min;
      DT_ lambda_max;
      DT_ vol_min;
      DT_ vol_max;
      cell_size_defect = meshopt_ctrl->compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max);

      if(comm.rank() == 0)
      {
        std::cout << "Pre initial quality indicator: " << stringify_fp_sci(qual_min) <<
          " / " << stringify_fp_sci(qual_avg) << " worst angle: " << stringify_fp_fix(worst_angle) << std::endl;
        std::cout << "Pre initial cell size defect: " << stringify_fp_sci(cell_size_defect) <<
          " lambda: " << stringify_fp_sci(lambda_min) << " " << stringify_fp_sci(lambda_max) <<
          " vol: " << stringify_fp_sci(vol_min) << " " << stringify_fp_sci(vol_max) << std::endl;
      }
    }

    // Check for the hard coded settings for test mode
    if(test_number == 1)
    {
      if( Math::abs(worst_angle - DT_(45)) > Math::sqrt(Math::eps<DataType>()))
      {
        Util::mpi_cout("FAILED:");
        throw InternalError(__func__,__FILE__,__LINE__,
        "Initial worst angle should be = "+stringify_fp_fix(45)+ " but is "+stringify_fp_fix(worst_angle)+"\n");
      }
      else if(test_number == 2)
      {
        if( Math::abs(worst_angle - DT_(45)) > Math::sqrt(Math::eps<DataType>()))
        {
          Util::mpi_cout("FAILED:");
          throw InternalError(__func__,__FILE__,__LINE__,
          "Initial worst angle should be >= "+stringify_fp_fix(45)+ " but is "+stringify_fp_fix(worst_angle)+"\n");
        }
      }
    }

    // Optimise the mesh
    meshopt_ctrl->optimise();

    // Write output again
    if(write_vtk)
    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());

        String vtk_name = String(file_basename+"_post_lvl_"+stringify(lvl_index));
        Util::mpi_cout("Writing "+vtk_name+"\n");

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
        meshopt_ctrl->add_to_vtk_exporter(exporter, deque_position);
        exporter.write(vtk_name, comm.rank(), comm.size());

        ++deque_position;
      }

    }

    // Compute quality indicators
    {
      Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::compute(qual_min, qual_sum,
      finest_mesh.template get_index_set<MeshType::shape_dim, 0>(), finest_mesh.get_vertex_set());

      worst_angle = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::angle(
        finest_mesh.template get_index_set<MeshType::shape_dim, 0>(), finest_mesh.get_vertex_set());

      comm.allreduce(&qual_min, &qual_min, std::size_t(1), Dist::op_min);
      comm.allreduce(&qual_sum, &qual_sum, std::size_t(1), Dist::op_sum);
      comm.allreduce(&worst_angle, &worst_angle, std::size_t(1), Dist::op_min);

      DataType qual_avg(qual_sum/DataType(ncells));

      DT_ lambda_min;
      DT_ lambda_max;
      DT_ vol_min;
      DT_ vol_max;
      cell_size_defect = meshopt_ctrl->compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max);

      if(comm.rank() == 0)
      {
        std::cout << "Post initial quality indicator: " << stringify_fp_sci(qual_min) <<
          " / " << stringify_fp_sci(qual_avg ) << " worst angle: " << stringify_fp_fix(worst_angle) << std::endl;
        std::cout << "Post initial cell size defect: " << stringify_fp_sci(cell_size_defect) <<
          " lambda: " << stringify_fp_sci(lambda_min) << " " << stringify_fp_sci(lambda_max) <<
          " vol: " << stringify_fp_sci(vol_min) << " " << stringify_fp_sci(vol_max) << std::endl;
      }
    }

    // Check for the hard coded settings for test mode
    if(test_number == 1)
    {
      if(worst_angle < DT_(44))
      {
        Util::mpi_cout("FAILED: Post Initial worst angle should be >= "+stringify_fp_fix(44)+
            " but is "+stringify_fp_fix(worst_angle)+"\n");
        return 1;
      }
      if(qual_min < DT_(2.4e-1))
      {
        Util::mpi_cout("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(2.4e-1)+
            " but is "+stringify_fp_fix(qual_min)+"\n");
        return 1;
      }
      if(cell_size_defect > DT_(2.2))
      {
        Util::mpi_cout("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(2.2)+
            " but is "+stringify_fp_fix(cell_size_defect)+"\n");
        return 1;
      }
    }
    else if(test_number == 2)
    {
      if(worst_angle < DT_(22))
      {
        Util::mpi_cout("FAILED: Post Initial worst angle should be >= "+stringify_fp_fix(22)+
            " but is "+stringify_fp_fix(worst_angle)+"\n");
        return 1;
      }
      if(qual_min < DT_(6.4e-1))
      {
        Util::mpi_cout("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(6.4e-1)+
            " but is "+stringify_fp_fix(qual_min)+"\n");
        return 1;
      }
      if(cell_size_defect > DT_(1.2e-1))
      {
        Util::mpi_cout("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(1.2e-1)+
            " but is "+stringify_fp_fix(cell_size_defect)+"\n");
        return 1;
      }
    }

    // Initial time
    DataType time(0);
    // Counter for timesteps
    Index n(0);

    // The mesh velocity is 1/delta_t*(coords_new - coords_old) and computed in each time step
    auto mesh_velocity = meshopt_ctrl->get_coords().clone();

    int return_value(0);

    bool* todo_boundary(new bool[(*old_coords).size()]);
    for(Index i(0); i < (*old_coords).size(); ++i)
      todo_boundary[i] = true;

    // Get a list of the dirichlet boundaries for modifying the boundary coordinates
    std::deque<String> dirichlet_boundaries(meshopt_ctrl->get_dirichlet_boundaries());

    while(time < t_end)
    {
      n++;
      time+= delta_t;

      if(comm.rank() == 0)
        std::cout << "Timestep " << n << " t = " << stringify_fp_fix(time) <<std::endl;

      // Clear statistics data so it does not eat us alive
      FEAT::Statistics::reset_solver_statistics();

      bool abort(false);

      // Save old vertex coordinates
      meshopt_ctrl->mesh_to_buffer();
      old_coords.copy(meshopt_ctrl->get_coords());

      // Get coords for modification
      new_coords.copy(meshopt_ctrl->get_coords());
      auto& coords_loc = *new_coords;

      for(Index i(0); i < (*old_coords).size(); ++i)
        todo_boundary[i] = true;

      // Update the boundary
      for(const auto& it:dirichlet_boundaries)
      {
        // Get outer boundary MeshPart
        auto* boundary_meshpart = dom_ctrl.get_levels().back()->get_mesh_node()->find_mesh_part(it);

        // If the meshpart is not there, our patch does not lie on that boundary, which is fine
        if(boundary_meshpart == nullptr)
          continue;

        // Get vertex target set
        const Geometry::TargetSet& boundary_set = boundary_meshpart->template get_target_set<0>();

        for(Index i(0); i < boundary_set.get_num_entities(); ++i)
        {
          Index j = boundary_set[i];
          if(todo_boundary[j])
          {
            todo_boundary[j] = false;
            Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp0(coords_loc(j));
            Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp1(tmp0);

            tmp1(0) += delta_t * (DataType(1)*(tmp0(0) - midpoint(0))
                + Math::pow(DataType(1)*(tmp0(1) - midpoint(1)),DataType(3) ) );
            tmp1(1) -= delta_t * (DataType(1)*(tmp1(0) - midpoint(1))
                + Math::pow(DataType(1)*(tmp0(0) - midpoint(0)),DataType(3) ) );

            coords_loc(j, tmp1);
          }
        }
      }

      meshopt_ctrl->prepare(new_coords);
      meshopt_ctrl->optimise();

      // Compute mesh velocity
      mesh_velocity.axpy(old_coords, meshopt_ctrl->get_coords(), DataType(-1));
      mesh_velocity.scale(mesh_velocity, DataType(1)/delta_t);

      // Compute maximum of the mesh velocity
      DataType max_mesh_velocity(0);
      for(IT_ i(0); i < (*mesh_velocity).size(); ++i)
        max_mesh_velocity = Math::max(max_mesh_velocity, (*mesh_velocity)(i).norm_euclid());

      if(comm.rank() == 0)
        std::cout << "max. mesh velocity: " << stringify_fp_sci(max_mesh_velocity) << std::endl;

      // Compute quality indicators
      {
        Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::compute(qual_min, qual_sum,
        finest_mesh.template get_index_set<MeshType::shape_dim, 0>(), finest_mesh.get_vertex_set());

        worst_angle = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::angle(
          finest_mesh.template get_index_set<MeshType::shape_dim, 0>(), finest_mesh.get_vertex_set());

        comm.allreduce(&qual_min, &qual_min, std::size_t(1), Dist::op_min);
        comm.allreduce(&qual_sum, &qual_sum, std::size_t(1), Dist::op_sum);
        comm.allreduce(&worst_angle, &worst_angle, std::size_t(1), Dist::op_min);

        DataType qual_avg(qual_sum/DataType(ncells));

        DT_ lambda_min;
        DT_ lambda_max;
        DT_ vol_min;
        DT_ vol_max;
        cell_size_defect = meshopt_ctrl->compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max);

        if(comm.rank() == 0)
        {
          std::cout << "Quality indicator: " << stringify_fp_sci(qual_min) <<
            " / " << stringify_fp_sci(qual_avg ) << " worst angle: " << stringify_fp_fix(worst_angle) << std::endl;
          std::cout << "Cell size defect: " << stringify_fp_sci(cell_size_defect) <<
            " lambda: " << stringify_fp_sci(lambda_min) << " " << stringify_fp_sci(lambda_max) <<
            " vol: " << stringify_fp_sci(vol_min) << " " << stringify_fp_sci(vol_max) << std::endl;
        }
      }

      if(worst_angle < DT_(1))
      {
        Util::mpi_cout("Mesh deteriorated, stopping.\n");
        return_value = 1;
        abort = true;
      }

      if(write_vtk && ( (n%vtk_freq == 0) || abort ))
      {
        String vtk_name(file_basename+"_post_"+stringify(n));

        Util::mpi_cout("Writing "+vtk_name+"\n");

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_ctrl.get_levels().back()->get_mesh());
        // Add mesh velocity
        exporter.add_vertex_vector("mesh_velocity", *mesh_velocity);
        // Add everything from the MeshoptControl
        meshopt_ctrl->add_to_vtk_exporter(exporter, int(dom_ctrl.get_levels().size())-1);
        // Write the file
        exporter.write(vtk_name, comm.rank(), comm.size());
      }

      if(abort)
        break;

    } // time loop

    Util::mpi_cout("Finished!\n");
    meshopt_ctrl->print();

    // Write final vtk output
    if(write_vtk)
    {
      int deque_position(0);
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());

        String vtk_name = String(file_basename+"_final_lvl_"+stringify(lvl_index));
        Util::mpi_cout("Writing "+vtk_name+"\n");

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));
        meshopt_ctrl->add_to_vtk_exporter(exporter, deque_position);
        exporter.write(vtk_name, comm.rank(), comm.size());

        ++deque_position;
      }

    }

    // Check for the hard coded settings for test mode
    if(test_number == 1)
    {
      if(worst_angle < DT_(43))
      {
        Util::mpi_cout("FAILED: Final worst angle should be >= "+stringify_fp_fix(43)+
            " but is "+stringify_fp_fix(worst_angle)+"\n");
        return 1;
      }
      if(qual_min < DT_(2.3e-1))
      {
        Util::mpi_cout("FAILED: Final worst shape quality should be >= "+stringify_fp_fix(2.3e-1)+
            " but is "+stringify_fp_fix(qual_min)+"\n");
        return 1;
      }
      if(cell_size_defect > DT_(2.3))
      {
        Util::mpi_cout("FAILED: Final cell size distribution defect should be < "+stringify_fp_fix(2.3)+
            " but is "+stringify_fp_fix(cell_size_defect)+"\n");
        return 1;
      }
    }
    else if(test_number == 2)
    {
      if(worst_angle < DT_(22))
      {
        Util::mpi_cout("FAILED: Final worst angle should be >= "+stringify_fp_fix(28)+
            " but is "+stringify_fp_fix(worst_angle)+"\n");
        return 1;
      }
      if(qual_min < DT_(6.4e-1))
      {
        Util::mpi_cout("FAILED: Final worst shape quality should be >= "+stringify_fp_fix(6.4e-1)+
            " but is "+stringify_fp_fix(qual_min)+"\n");
        return 1;
      }
      if(cell_size_defect > DT_(1.2))
      {
        Util::mpi_cout("FAILED: Final cell size distribution defect should be < "+stringify_fp_fix(1.2)+
            " but is "+stringify_fp_fix(cell_size_defect)+"\n");
        return 1;
      }
    }

    if(comm.rank() == 0)
    {
      TimeStamp bt;
      std::cout << "Elapsed time: " << bt.elapsed(at) << std::endl;
    }

    return return_value;

  }
}; // struct MeshoptBoundaryApp

int run_app(int argc, char* argv[])
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

  // create world communicator
  Dist::Comm comm(Dist::Comm::world());

  // initialise
#ifdef FEAT_HAVE_MPI
  if (comm.rank() == 0)
  {
    std::cout << "NUM-PROCS: " << comm.size() << std::endl;
  }
#endif

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

  if( args.check("help") > -1 || args.num_args()==1)
    display_help(comm);

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
    Util::mpi_cout("Running in test mode, all other command line arguments and configuration files are ignored.\n");

    if(args.check("test") > 1)
      throw InternalError(__func__, __FILE__, __LINE__, "Too many options for --test");

    args.parse("test",test_number);
    if(test_number != 1 && test_number != 2)
      throw InternalError(__func__, __FILE__, __LINE__, "Encountered unhandled test number "+stringify(test_number));
  }

  // Application settings, has to be created here because it gets filled differently according to test
  PropertyMap* application_config = new PropertyMap;

  // create a mesh file reader
  Geometry::MeshFileReader mesh_file_reader;

  // If we are not in test mode, parse command line arguments, read files, synchronise streams
  if(test_number == 0)
  {
    // Read the application config file on rank 0
    if(comm.rank() == 0)
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
          throw FileNotFound(application_config_filename);

        synchstream_app_config << ifs.rdbuf();
      }
    }

    // If we are in parallel mode, we need to synchronise the stream
    comm.bcast_stringstream(synchstream_app_config);

    // Parse the application config from the (synchronised) stream
    application_config->parse(synchstream_app_config, true);

    // Get the application settings section
    auto app_settings_section = application_config->query_section("ApplicationSettings");
    XASSERTM(app_settings_section != nullptr,
    "Application config is missing the mandatory ApplicationSettings section!");

    auto mesh_files_p = app_settings_section->query("mesh_files");
    mesh_files_p.first.split_by_charset(mesh_files, " ");

    // We read the files only on rank 0. After reading, we synchronise the streams like above.
    if(comm.rank() == 0)
    {
      // Read configuration for mesh optimisation to stream
      auto meshopt_config_filename_p = app_settings_section->query("meshopt_config_file");
      XASSERTM(meshopt_config_filename_p.second,
      "ApplicationConfig section is missing the mandatory meshopt_config_file entry!");
      {
        std::ifstream ifs(meshopt_config_filename_p.first);
        if(!ifs.good())
          throw FileNotFound(meshopt_config_filename_p.first);

        std::cout << "Reading mesh optimisation config from file " <<meshopt_config_filename_p.first << std::endl;
        synchstream_meshopt_config << ifs.rdbuf();
      }

      // Read solver configuration to stream
      auto solver_config_filename_p = app_settings_section->query("solver_config_file");
      XASSERTM(solver_config_filename_p.second,
      "ApplicationConfig section is missing the mandatory solver_config_file entry!");
      {
        std::ifstream ifs(solver_config_filename_p.first);
        if(!ifs.good())
          throw FileNotFound(solver_config_filename_p.first);
        else
        {
          std::cout << "Reading solver config from file " << solver_config_filename_p.first << std::endl;
          synchstream_solver_config << ifs.rdbuf();
        }
      }
    } // comm.rank() == 0

    // Synchronise all those streams in parallel mode
    comm.bcast_stringstream(synchstream_meshopt_config);
    comm.bcast_stringstream(synchstream_solver_config);
  }
  // If we are in test mode, all streams are filled by the hard coded stuff below
  else
  {
    read_test_application_config(synchstream_app_config, test_number);
    // Parse the application config from the (synchronised) stream
    application_config->parse(synchstream_app_config, true);

    read_test_meshopt_config(synchstream_meshopt_config, test_number);
    read_test_solver_config(synchstream_solver_config);

    read_test_mesh_file_names(mesh_files, test_number);
  }
  // Now we have all configurations in the corresponding streams and know the mesh file names

  // Create PropertyMaps and parse the configuration streams
  PropertyMap* meshopt_config = new PropertyMap;
  meshopt_config->parse(synchstream_meshopt_config, true);

  PropertyMap* solver_config = new PropertyMap;
  solver_config->parse(synchstream_solver_config, true);

  std::deque<std::stringstream> mesh_streams(mesh_files.size());

  // read all files
  for(std::size_t i(0); i < mesh_files.size(); ++i)
  {
    // read the stream
    DistFileIO::read_common(mesh_streams.at(i), mesh_files.at(i));

    // add to mesh reader
    mesh_file_reader.add_stream(mesh_streams.at(i));
  }

  int ret(1);

  // Get the mesh type sting from the parsed mesh so we know with which template parameter to call the application
  mesh_file_reader.read_root_markup();
  mesh_type = mesh_file_reader.get_meshtype_string();

  // Call the appropriate class' run() function
  if(mesh_type == "conformal:hypercube:2:2")
  {
    ret = MeshoptBoundaryApp<MemType, DataType, IndexType, H2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  else if(mesh_type == "conformal:simplex:2:2")
  {
    ret = MeshoptBoundaryApp<MemType, DataType, IndexType, S2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  else
    throw InternalError(__func__,__FILE__,__LINE__,"Unhandled mesh type "+mesh_type);

  delete application_config;
  delete meshopt_config;
  delete solver_config;

  return ret;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialise(argc, argv);
  int ret = run_app(argc, argv);
  FEAT::Runtime::finalise();
  return ret;
}

static void read_test_application_config(std::stringstream& iss, const int test_number)
{
  if(test_number == 1)
  {
    iss << "[ApplicationSettings]" << std::endl;
    iss << "mesh_optimiser = DuDvDefault" << std::endl;
    iss << "solver_config_file = ./solver_config.ini" << std::endl;
    iss << "delta_t = 1e-2" << std::endl;
    iss << "t_end = 2e-2" << std::endl;
    iss << "midpoint = 0.0 0.0" << std::endl;

    iss << "[DomainControlSettings]" << std::endl;
    iss << "parti-type = fallback parmetis" << std::endl;
    iss << "parti-rank-elems = 4" << std::endl;
    iss << "lvl_min = 1" << std::endl;
    iss << "lvl_max = 5" << std::endl;
  }
  else if(test_number == 2)
  {
    iss << "[ApplicationSettings]" << std::endl;
    iss << "mesh_optimiser = HyperelasticityDefault" << std::endl;
    iss << "solver_config_file = ./solver_config.ini" << std::endl;
    iss << "delta_t = 1e-2" << std::endl;
    iss << "t_end = 2e-2" << std::endl;
    iss << "midpoint = 0.0 0.0" << std::endl;

    iss << "[DomainControlSettings]" << std::endl;
    iss << "parti-type = fallback parmetis" << std::endl;
    iss << "parti-rank-elems = 4" << std::endl;
    iss << "lvl_min = 1" << std::endl;
    iss << "lvl_max = 3" << std::endl;
  }
  else
    throw InternalError(__func__,__FILE__,__LINE__,"Unknown test number: "+stringify(test_number));
}

static void read_test_meshopt_config(std::stringstream& iss, const int test_number)
{
  if(test_number == 1)
  {
    iss << "[DuDvDefault]" << std::endl;
    iss << "type = DuDv" << std::endl;
    iss << "config_section = DuDvDefaultParameters" << std::endl;
    iss << "fixed_reference_domain = 1" << std::endl;
    iss << "dirichlet_boundaries = outer" << std::endl;

    iss << "[DuDvDefaultParameters]" << std::endl;
    iss << "solver_config = PCG-MGV" << std::endl;
  }
  else if(test_number == 2)
  {
    iss << "[HyperElasticityDefault]" << std::endl;
    iss << "type = Hyperelasticity" << std::endl;
    iss << "config_section = HyperelasticityDefaultParameters" << std::endl;
    iss << "dirichlet_boundaries = bottom top left right" << std::endl;

    iss << "[HyperelasticityDefaultParameters]" << std::endl;
    iss << "global_functional = HyperelasticityFunctional" << std::endl;
    iss << "local_functional = RumpfFunctional" << std::endl;
    iss << "solver_config = NLCG" << std::endl;
    iss << "fac_norm = 1.0" << std::endl;
    iss << "fac_det = 1.0" << std::endl;
    iss << "fac_cof = 0.0" << std::endl;
    iss << "fac_reg = 1e-8" << std::endl;
    iss << "exponent_det = 2" << std::endl;
    iss << "scale_computation = once_uniform" << std::endl;
  }
  else
    throw InternalError(__func__,__FILE__,__LINE__,"Unknown test number "+stringify(test_number));
}

static void read_test_solver_config(std::stringstream& iss)
{
  iss << "[NLCG]" << std::endl;
  iss << "type = NLCG" << std::endl;
  iss << "precon = DuDvPrecon" << std::endl;
  iss << "plot = 1" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "max_iter = 1000" << std::endl;
  iss << "linesearch = StrongWolfeLinesearch" << std::endl;
  iss << "direction_update = DYHSHybrid" << std::endl;
  iss << "keep_iterates = 0" << std::endl;

  iss << "[DuDvPrecon]" << std::endl;
  iss << "type = DuDvPrecon" << std::endl;
  iss << "dirichlet_boundaries = bottom top left right" << std::endl;
  iss << "fixed_reference_domain = 1" << std::endl;
  iss << "linear_solver = PCG-MGV" << std::endl;

  iss << "[PCG-MGV]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "max_iter = 20" << std::endl;
  iss << "min_stag_iter = 1" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "plot = 1" << std::endl;
  iss << "precon = mgv" << std::endl;

  iss << "[strongwolfelinesearch]" << std::endl;
  iss << "type = StrongWolfeLinesearch" << std::endl;
  iss << "plot = 0" << std::endl;
  iss << "max_iter = 20" << std::endl;
  iss << "tol_decrease = 1e-3" << std::endl;
  iss << "tol_curvature = 0.3" << std::endl;
  iss << "keep_iterates = 0" << std::endl;

  iss << "[rich]" << std::endl;
  iss << "type = richardson" << std::endl;
  iss << "max_iter = 4" << std::endl;
  iss << "min_iter = 4" << std::endl;
  iss << "precon = jac" << std::endl;

  iss << "[jac]" << std::endl;
  iss << "type = jac" << std::endl;
  iss << "omega = 0.5" << std::endl;

  iss << "[mgv]" << std::endl;
  iss << "type = mgv" << std::endl;
  iss << "smoother = rich" << std::endl;
  iss << "coarse = pcg" << std::endl;

  iss << "[pcg]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "max_iter = 50" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "precon = jac" << std::endl;
}

static void read_test_mesh_file_names(std::deque<String>& mesh_files, const int test_number)
{
  String mesh_filename(FEAT_SRC_DIR);
  if(test_number == 1)
    mesh_filename +="/data/meshes/unit-circle-quad.xml";
  else if(test_number == 2)
    mesh_filename +="/data/meshes/unit-axe-tria.xml";
  else
    throw InternalError(__func__,__FILE__,__LINE__,"Encountered unhandled test "+stringify(test_number));

  mesh_files.push_back(mesh_filename);
}

static void display_help(const Dist::Comm& comm)
{
  if(comm.rank() == 0)
  {
    std::cout << "meshopt_boundary-app: Moving the boundary of a mesh and computing an extension into the interiour"
    << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << " --application_config: Path to the application configuration file" << std::endl;
    std::cout << "Optional arguments:" << std::endl;
    std::cout << " --test: Run as a test. Ignores configuration files and uses hard coded settings." << std::endl;
    std::cout << " --test [1 or 2]: Run as a test. Ignores configuration files and uses hard coded settings. " <<
      "Test 1 is quadrilateral cells, test 2 is triangular cells" << std::endl;
    std::cout << " --vtk <FREQ>: If this is set, vtk files are written every <FREQ> time steps." << std::endl;
    std::cout << " --help: Displays this text" << std::endl;
  }
}
