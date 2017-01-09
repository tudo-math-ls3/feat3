#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/assembly/fe_interpolator.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_quality_heuristic.hpp>
#include <kernel/geometry/mesh_extruder.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_control_factory.hpp>
#include <control/stokes_blocked.hpp>

using namespace FEAT;

static void display_help(const Dist::Comm&);
static void read_test_application_config(std::stringstream&, const int);
static void read_test_meshopt_config(std::stringstream&, const int);
static void read_test_solver_config(std::stringstream&, const int);
static void read_test_mesh_file_names(std::deque<String>&, const int);

template<typename Mesh_>
struct MeshExtrudeHelper
{
  typedef Mesh_ MeshType;
  typedef Mesh_ ExtrudedMeshType;
  typedef typename MeshType::CoordType CoordType;

  Geometry::RootMeshNode<MeshType>* extruded_mesh_node;

  explicit MeshExtrudeHelper(Geometry::RootMeshNode<MeshType>* DOXY(rmn), Index DOXY(slices), CoordType DOXY(z_min), CoordType DOXY(z_max), const String& DOXY(z_min_part_name), const String& DOXY(z_max_part_name)) :
    extruded_mesh_node(nullptr)
    {
    }

  MeshExtrudeHelper(const MeshExtrudeHelper&) = delete;

  ~MeshExtrudeHelper()
  {
  }

  void extrude_vertex_set(const typename MeshType::VertexSetType& DOXY(vtx))
  {
  }

};

template<typename Coord_>
struct MeshExtrudeHelper<Geometry::ConformalMesh<Shape::Hypercube<2>,2,2,Coord_>>
{
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>,2,2,Coord_> MeshType;
  typedef Coord_ CoordType;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>,3,3,Coord_> ExtrudedMeshType;
  typedef Geometry::MeshAtlas<ExtrudedMeshType> ExtrudedAtlasType;

  Geometry::MeshExtruder<MeshType> mesh_extruder;
  ExtrudedAtlasType* extruded_atlas;
  Geometry::RootMeshNode<ExtrudedMeshType>* extruded_mesh_node;

  explicit MeshExtrudeHelper(Geometry::RootMeshNode<MeshType>* rmn, Index slices, CoordType z_min, CoordType z_max, const String& z_min_part_name, const String& z_max_part_name) :
    mesh_extruder(slices, z_min, z_max, z_min_part_name, z_max_part_name),
    extruded_atlas(new ExtrudedAtlasType),
    extruded_mesh_node(new Geometry::RootMeshNode<ExtrudedMeshType>(nullptr, extruded_atlas))
  {
    mesh_extruder.extrude_atlas(*extruded_atlas, *(rmn->get_atlas()));
    mesh_extruder.extrude_root_node(*extruded_mesh_node, *rmn, extruded_atlas);
  }

  MeshExtrudeHelper(const MeshExtrudeHelper&) = delete;

  ~MeshExtrudeHelper()
  {
    if(extruded_atlas != nullptr)
    {
      delete extruded_atlas;
    }
    if(extruded_mesh_node != nullptr)
    {
      delete extruded_mesh_node;
    }
  }

  void extrude_vertex_set(const typename MeshType::VertexSetType& vtx)
  {
    mesh_extruder.extrude_vertex_set(extruded_mesh_node->get_mesh()->get_vertex_set(), vtx);
  }
};

static inline void dump_time(const Dist::Comm& comm, String s, double t, double total)
{
  comm.print(s.pad_back(30, '.') + ": " + stringify_fp_fix(t, 3, 10)
      + " (" + stringify_fp_fix(100.0*t/total,3,7) + "%)");
}


/**
 * \brief Navier-Stokes System Level class
 *
 * This extends the StokesBlockedSystemLevel by the corresponding filters for
 * the velocity and pressure sub-systems.
 */
template
<
  int dim_,
  typename MemType_ = Mem::Main,
  typename DataType_ = Real,
  typename IndexType_ = Index,
  typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
  typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
  typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
  typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
>
class NavierStokesBlockedSystemLevel :
  public Control::StokesBlockedUnitVeloMeanPresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
{
  public:
    typedef Control::StokesBlockedUnitVeloMeanPresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;

}; // class NavierStokesBlockedSystemLevel

template<typename Mem_, typename DT_, typename IT_, typename Mesh_>
struct NavierStokesScrewsApp
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
  /// FE space for the velocity
  typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
  /// FE space for the pressure
  typedef Space::Lagrange1::Element<TrafoType> SpacePresType;
  /// FE space for the transformation in the mesh optimisation problem
  typedef Space::Lagrange1::Element<TrafoType> SpaceTrafoType;

  /// The domain level with FE spaces for velocity and pressure
  typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
  /// Domain Control Type
  typedef Control::Domain::PartiDomainControl<DomainLevelType> DomCtrl;

  /// Type of the extruded mesh, which is Simplex<2> for Simplex<2> meshes (no extrusion) and Hypercube<3> for
  /// Hypercube<2> meshes
  typedef typename MeshExtrudeHelper<MeshType>::ExtrudedMeshType ExtrudedMeshType;

  /// This is how far the inner screw's centre deviates from the outer screw's
  static constexpr DataType excentricity_inner = DataType(0.2833);

  /**
   * \brief Returns a descriptive string
   *
   * \returns The class name as string
   */
  static String name()
  {
    return "NavierStokesScrewsApp";
  }

  /**
   * \brief The routine that does the actual work
   */
  static int run(const SimpleArgParser& args, Dist::Comm& comm, PropertyMap* application_config,
  PropertyMap* meshopt_config, PropertyMap* solver_config, Geometry::MeshFileReader& mesh_file_reader)
  {

    // Create a time-stamp
    TimeStamp stamp_start;

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

    DataType reynolds(1);
    bool use_deformation(false);
    Index nl_its_max(0);
    Index uzawa_its_max(0);

    // Need some pi for all the angles
    DataType pi(Math::pi<DataType>());


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
        throw InternalError(__func__, __FILE__, __LINE__,
        "Encountered unhandled test number "+stringify(test_number));
      }
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

    // Get reynolds number
    auto reynolds_p = app_settings_section->query("reynolds");
    XASSERTM(reynolds_p.second, "ApplicationConfig section is missing the mandatory reynolds entry!");
    reynolds = std::stod(reynolds_p.first);
    XASSERT(reynolds > DataType(0));

    // Get the number of nonlinear iterations
    auto nl_its_max_p = app_settings_section->query("nl_its_max");
    if(nl_its_max_p.second)
    {
      nl_its_max = Index(std::stoul(nl_its_max_p.first));
    }

    // Get the number of Uzawa iterations
    auto uzawa_its_max_p = app_settings_section->query("uzawa_its_max");
    if(uzawa_its_max_p.second)
    {
      uzawa_its_max = Index(std::stoul(uzawa_its_max_p.first));
    }

    // Get the application settings section
    auto domain_control_settings_section = application_config->query_section("DomainControlSettings");
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

    MeshExtrudeHelper<MeshType> extruder(dom_ctrl.front()->get_mesh_node(),
    Index(10*(lvl_max+1)), DataType(0), DataType(1), "bnd:b", "bnd:t");

    // Get outer boundary MeshPart. Can be nullptr if this process' patch does not lie on that boundary
    auto* outer_boundary_part = dom_ctrl.front()->get_mesh_node()->find_mesh_part("bnd:o");
    Geometry::TargetSet* outer_indices(nullptr);
    if(outer_boundary_part != nullptr)
    {
      outer_indices = &(outer_boundary_part->template get_target_set<0>());
    }

    // This is the centre point of the rotation of the outer screw
    WorldPoint centre_outer(DataType(0));
    DT_ angular_velocity_outer(DT_(2)*pi);
    auto* outer_chart = dom_ctrl.get_atlas().find_mesh_chart("screw:o");

    // Get inner boundary MeshPart. Can be nullptr if this process' patch does not lie on that boundary
    auto* inner_boundary = dom_ctrl.front()->get_mesh_node()->find_mesh_part("bnd:i");
    Geometry::TargetSet* inner_indices(nullptr);
    if(inner_boundary != nullptr)
    {
      inner_indices = &(inner_boundary->template get_target_set<0>());
    }

    // This is the centre point of the rotation of the inner screw
    WorldPoint centre_inner(DataType(0));
    centre_inner.v[0] = -excentricity_inner;
    DT_ angular_velocity_inner(angular_velocity_outer*DT_(7)/DT_(6));
    auto* inner_chart = dom_ctrl.get_atlas().find_mesh_chart("screw:i");

    // Create MeshoptControl
    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl>> meshopt_ctrl(nullptr);
    meshopt_ctrl = Control::Meshopt::ControlFactory<Mem_, DT_, IT_>::create_meshopt_control(
      dom_ctrl, meshoptimiser_key_p.first, meshopt_config, solver_config);

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

    //// Check for the hard coded settings for test mode
    //if( (test_number == 1) && ( edge_angle < DT_(5.69)) )
    //{
    //  comm.print("FAILED: Initial worst edge angle should be >= "+stringify_fp_fix(5.69)
    //      + " but is "+stringify_fp_fix(edge_angle));
    //  ++ret;
    //}
    //else if( (test_number == 2) && ( edge_angle < DT_(2.83)) )
    //{
    //  comm.print("FAILED: Initial worst edge angle should be >= "+stringify_fp_fix(2.83)
    //      + " but is "+stringify_fp_fix(edge_angle));
    //}

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

      // If we did not compute this for the vtk output, we have to do it here
      if(!write_vtk)
      {
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);
      }

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
    //if(test_number == 1)
    //{
    //  if(edge_angle < DT_(5.69))
    //  {
    //    comm.print("FAILED: Post Initial worst edge angle should be >= "+stringify_fp_fix(5.69)+
    //        " but is "+stringify_fp_fix(edge_angle));
    //    ++ret;
    //  }
    //  if(qi_min < DT_(2.5e-2))
    //  {
    //    comm.print("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(2.5e-2)+
    //        " but is "+stringify_fp_fix(qi_min));
    //    ++ret;
    //  }
    //  if(cell_size_defect > DT_(2.5))
    //  {
    //    comm.print("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(2.5)+
    //        " but is "+stringify_fp_fix(cell_size_defect));
    //    ++ret;
    //  }
    //}
    //else if(test_number == 2)
    //{
    //  if(edge_angle < DT_(3.2))
    //  {
    //    comm.print("FAILED: Post Initial worst edge angle should be >= "+stringify_fp_fix(3.2)+
    //        " but is "+stringify_fp_fix(edge_angle));
    //    ++ret;
    //  }
    //  if(qi_min < DT_(2.5e-1))
    //  {
    //    comm.print("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(2.5e-1)+
    //        " but is "+stringify_fp_fix(qi_min));
    //    ++ret;
    //  }
    //  if(cell_size_defect > DT_(3.8e-1))
    //  {
    //    comm.print("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(3.8e-1)+
    //        " but is "+stringify_fp_fix(cell_size_defect));
    //    ++ret;
    //  }
    //}

    // define our velocity and pressure system levels
    typedef NavierStokesBlockedSystemLevel<ShapeType::dimension, MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = dom_ctrl.size_physical();

    // create a batch of stop-watches
    StopWatch watch_total, watch_asm_rhs, watch_asm_mat, watch_calc_def,
    watch_sol_init, watch_solver_a, watch_solver_s, watch_vtk, watch_meshopt;

    watch_total.start();

    // create stokes and system levels
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("auto-degree:7");

    /* ***************************************************************************************** */

    comm.print("Assembling gates...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gates(dom_ctrl.at(i));
    }

    /* ***************************************************************************************** */

    comm.print("Assembling transfers...");

    for (Index i(0); (i+1) < dom_ctrl.size_virtual(); ++i)
    {
      system_levels.at(i)->assemble_coarse_muxers(dom_ctrl.at(i+1));
      system_levels.at(i)->assemble_transfers(dom_ctrl.at(i), dom_ctrl.at(i+1), cubature);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling basic matrices...");

    for(Index i(0); i < num_levels; ++i)
    {
      // assemble velocity matrix structure
      system_levels.at(i)->assemble_velo_struct(dom_ctrl.at(i)->space_velo);
      // assemble pressure matrix structure
      system_levels.at(i)->assemble_pres_struct(dom_ctrl.at(i)->space_pres);
      // assemble pressure laplace matrix
      system_levels.at(i)->matrix_s.local().format();
      Assembly::Common::LaplaceOperator laplace_op;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(system_levels.at(i)->matrix_s.local(),
      laplace_op, dom_ctrl.at(i)->space_pres, cubature);
    }

    // assemble B/D matrices on finest level
    system_levels.front()->assemble_grad_div_matrices(dom_ctrl.front()->space_velo, dom_ctrl.front()->space_pres, cubature);

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    //Analytic::Common::XYPlaneRotation<DataType, MeshType::world_dim> rotation_inner(angular_velocity_inner, centre_inner);
    //Analytic::Common::XYPlaneRotation<DataType, MeshType::world_dim> rotation_outer(angular_velocity_outer, centre_outer);
    Tiny::Vector<DataType, 2> zeros_y(0);
    zeros_y(0) = -DataType(1);
    zeros_y(1) = DataType(1);
    DataType amplitude(1);
    Analytic::Common::YZPlaneParabolic<DataType, MeshType::world_dim> rotation_outer(amplitude, zeros_y);

    std::vector<Assembly::UnitFilterAssembler<MeshType>> unit_asm_velo_i(num_levels), unit_asm_velo_o(num_levels);

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local system filters
      typename SystemLevelType::LocalVeloFilter& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // Add inner boundary to assembler
      {
        auto* mesh_part_node = dom_ctrl.at(i)->get_mesh_node()->find_mesh_part_node("bnd:i");

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        if(mesh_part_node != nullptr)
        {
          auto* mesh_part = mesh_part_node->get_mesh();

          if(mesh_part != nullptr)
          {
            unit_asm_velo_i[i].add_mesh_part(*mesh_part);
          }
        }
      }

      // Add outer boundary to assembler
      {
        auto* mesh_part_node = dom_ctrl.at(i)->get_mesh_node()->find_mesh_part_node("bnd:o");

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        if(mesh_part_node != nullptr)
        {
          auto* mesh_part = mesh_part_node->get_mesh();

          if(mesh_part != nullptr)
          {
            unit_asm_velo_o[i].add_mesh_part(*mesh_part);
          }
        }
      }

      // assemble the velocity filter
      //unit_asm_velo_i[i].assemble(fil_loc_v, dom_ctrl.at(i)->space_velo, rotation_inner);
      unit_asm_velo_o[i].assemble(fil_loc_v, dom_ctrl.at(i)->space_velo, rotation_outer);

      // assembly of the pressure filter is done in the system level
      system_levels.at(i)->assemble_pressure_mean_filter(dom_ctrl.at(i)->space_pres, cubature);

    } // all levels

    /* ***************************************************************************************** */

    // get our vector types
    typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *dom_ctrl.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // get our fine-level matrices
    typename SystemLevelType::GlobalMatrixBlockA& matrix_a = the_system_level.matrix_a;
    typename SystemLevelType::GlobalMatrixBlockB& matrix_b = the_system_level.matrix_b;
    typename SystemLevelType::GlobalMatrixBlockD& matrix_d = the_system_level.matrix_d;
    typename SystemLevelType::GlobalSchurMatrix& matrix_s = the_system_level.matrix_s;

    // get out fine-level filters
    typename SystemLevelType::GlobalVeloFilter& filter_v = the_system_level.filter_velo;
    typename SystemLevelType::GlobalPresFilter& filter_p = the_system_level.filter_pres;

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("Setting up Velocity Multigrid...");

    Solver::MatrixStock<typename SystemLevelType::GlobalMatrixBlockA, typename SystemLevelType::GlobalVeloFilter,
    typename SystemLevelType::GlobalVeloTransfer> matrix_stock_velo;
    for (auto & system_level: system_levels)
    {
      matrix_stock_velo.systems.push_back(system_level->matrix_a.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_velo.gates_row.push_back(&system_level->gate_velo);
      matrix_stock_velo.gates_col.push_back(&system_level->gate_velo);
      matrix_stock_velo.filters.push_back(system_level->filter_velo.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_velo.muxers.push_back(&system_level->coarse_muxer_velo);
      matrix_stock_velo.transfers.push_back(system_level->transfer_velo.clone(LAFEM::CloneMode::Shallow));
    }

    auto tsolver_a = Control::SolverFactory::create_scalar_solver(matrix_stock_velo, solver_config, "linsolver_a");
    Solver::PreconditionedIterativeSolver<typename decltype(tsolver_a)::element_type::VectorType>* solver_a =
      (Solver::PreconditionedIterativeSolver<typename decltype(tsolver_a)::element_type::VectorType>*) &(*tsolver_a);
    matrix_stock_velo.hierarchy_init_symbolic();
    solver_a->init_symbolic();


    /* ***************************************************************************************** */

    comm.print("Setting up Pressure Multigrid...");

    Solver::MatrixStock<typename SystemLevelType::GlobalSchurMatrix, typename SystemLevelType::GlobalPresFilter,
    typename SystemLevelType::GlobalPresTransfer> matrix_stock_pres;
    for (auto & system_level: system_levels)
    {
      matrix_stock_pres.systems.push_back(system_level->matrix_s.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_pres.gates_row.push_back(&system_level->gate_pres);
      matrix_stock_pres.gates_col.push_back(&system_level->gate_pres);
      matrix_stock_pres.filters.push_back(system_level->filter_pres.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_pres.muxers.push_back(&system_level->coarse_muxer_pres);
      matrix_stock_pres.transfers.push_back(system_level->transfer_pres.clone(LAFEM::CloneMode::Shallow));
    }

    auto tsolver_s = Control::SolverFactory::create_scalar_solver(matrix_stock_pres, solver_config, "linsolver_s");
    Solver::PreconditionedIterativeSolver<typename decltype(tsolver_s)::element_type::VectorType>* solver_s =
      (Solver::PreconditionedIterativeSolver<typename decltype(tsolver_s)::element_type::VectorType>*) &(*tsolver_s);
    matrix_stock_pres.hierarchy_init_symbolic();
    solver_s->init_symbolic();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("\n");

    // create RHS and SOL vectors
    GlobalVeloVector vec_sol_v = matrix_a.create_vector_l();
    GlobalPresVector vec_sol_p = matrix_s.create_vector_l();
    GlobalVeloVector vec_rhs_v = matrix_a.create_vector_l();
    GlobalPresVector vec_rhs_p = matrix_s.create_vector_l();

    // create defect and correction vectors
    GlobalVeloVector vec_def_v = matrix_a.create_vector_l();
    GlobalPresVector vec_def_p = matrix_s.create_vector_l();
    GlobalVeloVector vec_cor_v = matrix_a.create_vector_l();
    GlobalPresVector vec_cor_p = matrix_s.create_vector_l();
    // create convection vector
    GlobalVeloVector vec_conv = matrix_a.create_vector_l();

    // format the vectors
    vec_sol_v.format();
    vec_sol_p.format();
    vec_rhs_v.format();
    vec_rhs_p.format();

    // apply filter onto solution vector
    filter_v.filter_sol(vec_sol_v);

    // create solution backup vectors; these store vec_sol_v/p of the last two time-steps
    GlobalVeloVector vec_sol_v_1 = vec_sol_v.clone();
    GlobalVeloVector vec_sol_v_2 = vec_sol_v.clone();
    GlobalPresVector vec_sol_p_1 = vec_sol_p.clone();

    // set up a burgers assembler for the RHS vector
    Assembly::BurgersAssembler<DataType, IndexType, 2> burgers_rhs;
    burgers_rhs.deformation = use_deformation;
    burgers_rhs.nu = -delta_t/reynolds;
    burgers_rhs.beta = -delta_t;
    burgers_rhs.theta = DataType(1);

    // set up a burgers assembler for the velocity matrix
    Assembly::BurgersAssembler<DataType, IndexType, 2> burgers_mat;
    burgers_mat.deformation = use_deformation;
    burgers_mat.nu = delta_t/reynolds;
    burgers_mat.beta = delta_t;
    burgers_mat.theta = DataType(1);

    // Initial time
    DataType time(0);
    // Counter for timesteps
    Index time_step(0);

    // The mesh velocity is 1/delta_t*(coords_new - coords_old) and computed in each time step
    auto mesh_velocity = meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep);
    mesh_velocity.format();

    GlobalVeloVector mesh_velocity_q2 = matrix_a.create_vector_l();

    // This is the absolute turning angle of the screws
    DataType alpha(0);
    WorldPoint rotation_angles(0);

    while(time < t_end)
    {
      ++time_step;
      time += delta_t;

      // Clear statistics data so it does not eat us alive
      FEAT::Statistics::reset_solver_statistics();

      DataType alpha_old = alpha;
      alpha = angular_velocity_outer*time;
      DataType delta_alpha = alpha - alpha_old;

      bool failure(false);

      comm.print("Timestep "+stringify(time_step)+": t = "+stringify_fp_fix(time)+", angle = "
          +stringify_fp_fix(alpha/(DataType(2)*pi)*DataType(360)) + " degrees");

      watch_meshopt.start();
      // Save old vertex coordinates
      meshopt_ctrl->mesh_to_buffer();
      old_coords.copy(meshopt_ctrl->get_coords());
      new_coords.copy(meshopt_ctrl->get_coords());
      auto& coords_loc(new_coords.local());

      // Rotate the charts
      rotation_angles(0) = delta_t*angular_velocity_outer;
      rotation_angles(1) = DataType(0);

      if(outer_chart != nullptr)
      {
        outer_chart->transform(centre_outer, rotation_angles, centre_outer);
      }

      rotation_angles(0) = delta_t*angular_velocity_inner;
      if(inner_chart != nullptr)
      {
        inner_chart->transform(centre_inner, rotation_angles, centre_inner);
      }

      //// Update boundary of the inner screw
      //// This is the 2x2 matrix representing the turning by the angle delta_alpha of the inner screw
      //if(inner_indices != nullptr)
      //{
      //  Tiny::Matrix<DataType, 2, 2> rot(DataType(0));

      //  rot.set_rotation_2d(delta_t*angular_velocity_inner);

      //  WorldPoint tmp(DataType(0));
      //  WorldPoint tmp2(DataType(0));

      //  for(Index i(0); i < inner_indices->get_num_entities(); ++i)
      //  {
      //    // Index of boundary vertex i in the mesh
      //    Index j(inner_indices->operator[](i));
      //    // Translate the point to the centre of rotation
      //    tmp = coords_loc(j) - centre_inner;
      //    // Rotate
      //    tmp2.set_mat_vec_mult(rot, tmp);
      //    // Translate the point by the new centre of rotation
      //    coords_loc(j, centre_inner + tmp2);
      //  }
      //}

      //// The outer screw has 7 teeth as opposed to the inner screw with 6, and it rotates at 6/7 of the speed
      //if(outer_indices != nullptr)
      //{
      //  Tiny::Matrix<DataType, 2, 2> rot(DataType(0));
      //  rot.set_rotation_2d(delta_t*angular_velocity_outer);

      //  WorldPoint tmp(DataType(0));
      //  WorldPoint tmp2(DataType(0));

      //  // The outer screw rotates centrically, so centre_outer remains the same at all times
      //  for(Index i(0); i < outer_indices->get_num_entities(); ++i)
      //  {
      //    // Index of boundary vertex i in the mesh
      //    Index j(outer_indices->operator[](i));
      //    tmp = coords_loc(j) - centre_outer;

      //    tmp2.set_mat_vec_mult(rot, tmp);

      //    coords_loc(j, centre_outer+tmp2);
      //  }
      //}

      // Now prepare the functional
      meshopt_ctrl->prepare(new_coords);

      meshopt_ctrl->optimise();
      watch_meshopt.stop();

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
        String msg = String("max. mesh velocity").pad_back(pad_width, ' ') + ": " + stringify_fp_sci(max_mesh_velocity);
        comm.print(msg);
      }

      // Now the flow problem

      // Assemble filters on all levels
      for(Index i(0); i < num_levels; ++i)
      {
        // get our local system filters
        typename SystemLevelType::LocalVeloFilter& fil_loc_v = system_levels.at(i)->filter_velo.local();

        // assemble the velocity filter
        //unit_asm_velo_i.at(i).assemble(fil_loc_v, dom_ctrl.at(i)->space_velo, rotation_inner);
        unit_asm_velo_o.at(i).assemble(fil_loc_v, dom_ctrl.at(i)->space_velo, rotation_outer);

        // assembly of the pressure filter is done in the system level
        system_levels.at(i)->assemble_pressure_mean_filter(dom_ctrl.at(i)->space_pres, cubature);

      } // all levels

      // apply filter onto solution vector
      filter_v.filter_sol(vec_sol_v);
      filter_p.filter_sol(vec_sol_p);

      // assemble B/D matrices on finest level
      system_levels.front()->assemble_grad_div_matrices(dom_ctrl.front()->space_velo, dom_ctrl.front()->space_pres, cubature);

      matrix_stock_pres.hierarchy_init_numeric();
      solver_s->init_numeric();

      // assemble RHS vector
      watch_asm_rhs.start();
      vec_rhs_v.format();
      vec_rhs_p.format();
      burgers_rhs.assemble(the_domain_level.space_velo, cubature, vec_sol_v.local(), nullptr, &vec_rhs_v.local());
      vec_rhs_v.sync_0();


      // subtract pressure (?)
      //matrix_b.apply(vec_rhs_v, vec_sol_p, vec_rhs_v, -DataType(1));
      watch_asm_rhs.stop();

      // apply RHS filter
      filter_v.filter_rhs(vec_rhs_v);
      filter_p.filter_rhs(vec_rhs_p);

      // non-linear loop
      for(Index nl_it(0); nl_it < nl_its_max; ++nl_it)
      {
        // Phase 1: compute convection vector
        // extrapolate previous time-step solution in first NL step
        if((time_step > Index(2)) && (nl_it == Index(0)))
        {
          // linear extrapolation of solution in time
          vec_conv.scale(vec_sol_v_1, DataType(2));
          vec_conv.axpy(vec_sol_v_2, vec_conv, -DataType(1));
        }
        else
        {
          // constant extrapolation of solution in time
          vec_conv.copy(vec_sol_v);
        }

        // Subtract ALE velocity
        Assembly::FEInterpolator<SpaceVeloType, SpaceTrafoType>::interpolate(mesh_velocity_q2.local(), mesh_velocity.local(), the_domain_level.space_velo, the_domain_level.space_pres);
        vec_conv.axpy(mesh_velocity_q2, vec_conv, -DataType(1));

        // Phase 2: loop over all levels and assemble the burgers matrices
        watch_asm_mat.start();
        //if(cfg.multigrid_a)
        //{
        // assemble burgers matrices on all levels
        for(std::size_t i(0); i < system_levels.size(); ++i)
        {
          auto& loc_mat_a = system_levels.at(i)->matrix_a.local();
          typename GlobalVeloVector::LocalVectorType vec_cv(vec_conv.local(), loc_mat_a.rows(), IndexType(0));
          loc_mat_a.format();
          burgers_mat.assemble(dom_ctrl.at(i)->space_velo, cubature, vec_cv, &loc_mat_a, nullptr);
        }

        //}
        //else
        //{
        //  // assemble burgers matrices on finest level
        //  the_system_level.matrix_a.local().format();
        //  burgers_mat.assemble(the_domain_level.space_velo, cubature, vec_conv.local(),
        //    &the_system_level.matrix_a.local(), nullptr);
        //}
        watch_asm_mat.stop();

        // Phase 4: initialise linear solvers
        watch_sol_init.start();
        matrix_stock_velo.hierarchy_init_numeric();
        solver_a->init_numeric();
        watch_sol_init.stop();

        // Phase 3: compute non-linear defects
        watch_calc_def.start();
        matrix_a.apply(vec_def_v, vec_sol_v, vec_rhs_v, -DataType(1));
        matrix_b.apply(vec_def_v, vec_sol_p, vec_def_v, -DataType(1));
        matrix_d.apply(vec_def_p, vec_sol_v, vec_rhs_p, -DataType(1));
        filter_v.filter_def(vec_def_v);
        filter_p.filter_def(vec_def_p);

        // compute defect norms
        const DataType def_nl1_v = vec_def_v.norm2();
        const DataType def_nl1_p = vec_def_p.norm2();
        watch_calc_def.stop();

        // console output, part 1
        if(comm.rank() == 0)
        {
          String line("NVS-NL intermediate:");
          line += stringify(nl_it).pad_front(4);
          line += " : ";
          line += stringify_fp_sci(def_nl1_v) + " ";
          line += stringify_fp_sci(def_nl1_p);
          std::cout << line << std::endl;
        }

        // linear solver iterations counts
        Index iter_v(0), iter_p(0);

        // Phase 5: linear DPM loop
        // Note: we need to perform one additional velocity solve,
        // so the break condition of the loop is inside...
        for(Index dpm_step(0); ; ++dpm_step)
        {
          // solve velocity system
          FEAT::Statistics::expression_target = "solver_a";
          watch_solver_a.start();
          Solver::Status status_a = solver_a->apply(vec_cor_v, vec_def_v);
          watch_solver_a.stop();
          if(!Solver::status_success(status_a))
          {
            comm.print("\n\nERROR: velocity solver broke down!\n");
            failure = true;
            break;
          }
          iter_v += solver_a->get_num_iter();

          // update velocity solution
          vec_sol_v.axpy(vec_cor_v, vec_sol_v);

          // are we done yet?
          if(dpm_step >= uzawa_its_max)
            break;

          // update pressure defect
          watch_calc_def.start();
          matrix_d.apply(vec_def_p, vec_sol_v, vec_rhs_p, -DataType(1));
          DT_ def_div_v = vec_def_p.norm2();
          filter_p.filter_def(vec_def_p);
          DT_ def_div_v_filtered = vec_def_p.norm2();

          watch_calc_def.stop();

          // solve pressure system
          FEAT::Statistics::expression_target = "solver_s";
          watch_solver_s.start();
          Solver::Status status_s = solver_s->apply(vec_cor_p, vec_def_p);
          watch_solver_s.stop();
          if(!Solver::status_success(status_s))
          {
            comm.print("\n\nERROR: pressure solver broke down!\n");
            failure = true;
            break;
          }
          iter_p += solver_s->get_num_iter();

          // update pressure solution
          filter_p.filter_cor(vec_cor_p);
          vec_sol_p.axpy(vec_cor_p, vec_sol_p, -DataType(1));

          // compute new defect
          watch_calc_def.start();
          matrix_a.apply(vec_def_v, vec_sol_v, vec_rhs_v, -DataType(1));
          matrix_b.apply(vec_def_v, vec_sol_p, vec_def_v, -DataType(1));
          filter_v.filter_def(vec_def_v);
          filter_p.filter_def(vec_def_p);
          watch_calc_def.stop();

          // console output, part 2
          if(comm.rank() == 0)
          {
            String line("DPM : ");
            line += stringify(dpm_step).pad_front(4);
            line += ": div(v) "+ stringify_fp_sci(def_div_v);
            line += ": filterd "+ stringify_fp_sci(def_div_v_filtered);
            std::cout << line << std::endl;
          }
        } // inner Uzawa loop

        // Phase 6: release linear solvers
        solver_a->done_numeric();
        //if(cfg.multigrid_a)
        //  solver_a->get_hierarchy()->done_numeric();
        matrix_stock_velo.hierarchy_done_numeric();

        // epic fail?
        if(failure)
        {
          break;
        }

        // Phase 7: compute final defect and norms (only for console output)
        watch_calc_def.start();
        matrix_a.apply(vec_def_v, vec_sol_v, vec_rhs_v, -DataType(1));
        matrix_b.apply(vec_def_v, vec_sol_p, vec_def_v, -DataType(1));
        matrix_d.apply(vec_def_p, vec_sol_v, vec_rhs_p, -DataType(1));
        filter_v.filter_def(vec_def_v);
        filter_p.filter_def(vec_def_p);

        const DataType def_nl2_v = vec_def_v.norm2();
        const DataType def_nl2_p = vec_def_p.norm2();
        watch_calc_def.stop();

        // console output, part 2
        if(comm.rank() == 0)
        {
          String line("NVS-NL projected:");
          line += stringify(uzawa_its_max).pad_front(4);
          line += ": v "+ stringify_fp_sci(def_nl2_v);
          line += "(" + stringify(iter_v).pad_front(4)+")";
          line += ": p " + stringify_fp_sci(def_nl2_p);
          line += "(" + stringify(iter_p).pad_front(4)+")";
          std::cout << line << std::endl;
        }
      } // non-linear loop

      solver_s->done_numeric();
      matrix_stock_pres.hierarchy_done_numeric();

      // epic fail?
      if(failure)
        break;

      // Write output on the finest level only
      if(write_vtk && ( (time_step%vtk_freq == 0) || failure))
      {

        watch_vtk.start();

        // Compute mesh quality on the finest
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);

        String vtk_name = String(file_basename+"_post_"+stringify(time_step));
        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_ctrl.front()->get_mesh());

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise);
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise);

        meshopt_ctrl->add_to_vtk_exporter(exporter, -1);

        exporter.add_vertex_vector("mesh velocity", mesh_velocity.local());

        exporter.add_vertex_vector("v", vec_sol_v.local());
        exporter.add_vertex_scalar("def_p", vec_def_p.local().elements());

        // compute and write time-derivatives
        GlobalVeloVector vec_der_v = vec_sol_v.clone();
        GlobalPresVector vec_der_p = vec_sol_p.clone();

        // Scale the pressure with the missing factor 1/delta_t. We abuse vec_der_p for that
        vec_der_p.scale(vec_der_p, DataType(1) / delta_t);
        exporter.add_vertex_scalar("p", vec_der_p.local().elements());

        vec_der_v.axpy(vec_sol_v_1, vec_der_v, -DataType(1));
        vec_der_p.axpy(vec_sol_p_1, vec_der_p, -DataType(1));

        vec_der_v.scale(vec_der_v, DataType(1) / delta_t);
        vec_der_p.scale(vec_der_p, DataType(1) / delta_t);

        exporter.add_vertex_vector("v_dt", vec_der_v.local());
        exporter.add_vertex_scalar("p_dt", vec_der_p.local().elements());

        exporter.write(vtk_name, comm);

        if(extruder.extruded_mesh_node != nullptr)
        {
          //extruder.extrude_vertex_set(finest_mesh.get_vertex_set());
          //// Create a VTK exporter for our mesh
          //String extruded_vtk_name = String(file_basename+"_post_extruded_"+stringify(time_step));

          //comm.print("Writing "+extruded_vtk_name);

          //Geometry::ExportVTK<ExtrudedMeshType> extruded_exporter(*(extruder.extruded_mesh_node->get_mesh()));
          //extruded_exporter.write(extruded_vtk_name, comm);
        }

        watch_vtk.stop();
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
        if(! (write_vtk && ( (time_step%vtk_freq == 0) || failure) ) )
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
        failure = true;
      }

      if(failure)
      {
        break;
      }

      // finally, update our solution vector backups
      vec_sol_v_2.copy(vec_sol_v_1);
      vec_sol_v_1.copy(vec_sol_v);
      vec_sol_p_1.copy(vec_sol_p);

      // continue with next time-step

    } // time loop

    solver_a->done_symbolic();
    solver_s->done_symbolic();
    matrix_stock_velo.hierarchy_done_symbolic();
    matrix_stock_pres.hierarchy_done_symbolic();

    meshopt_ctrl->print();

    // Write final vtk output
    if(write_vtk)
    {
      watch_vtk.start();
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

        exporter.add_vertex_vector("v", vec_sol_v.local());

        // compute and write time-derivatives
        GlobalVeloVector vec_der_v = vec_sol_v.clone();
        GlobalPresVector vec_der_p = vec_sol_p.clone();

        // Scale the pressure with the missing factor 1/delta_t. We abuse vec_der_p for that
        vec_der_p.scale(vec_der_p, DataType(1) / delta_t);
        exporter.add_vertex_scalar("p", vec_der_p.local().elements());

        vec_der_v.axpy(vec_sol_v_1, vec_der_v, -DataType(1));
        vec_der_p.axpy(vec_sol_p_1, vec_der_p, -DataType(1));

        vec_der_v.scale(vec_der_v, DataType(1) / delta_t);
        vec_der_p.scale(vec_der_p, DataType(1) / delta_t);

        exporter.add_vertex_vector("v_dt", vec_der_v.local());
        exporter.add_vertex_scalar("p_dt", vec_der_p.local().elements());

        exporter.write(vtk_name, comm.rank(), comm.size());
      }
      watch_vtk.stop();
    }

    // Check for the hard coded settings for test mode
    if(test_number == 1)
    {
      if(edge_angle < DT_(5.53))
      {
        comm.print("FAILED: Final worst edge angle should be >= "+stringify_fp_fix(5.53)+
            " but is "+stringify_fp_fix(edge_angle)+"\n");
        ++ret;
      }
      if(qi_min < DT_(2.51e-2))
      {
        comm.print("FAILED: Final worst shape quality should be <= "+stringify_fp_fix(2.51e-2)+
            " but is "+stringify_fp_fix(qi_min)+"\n");
        ++ret;
      }
      if(cell_size_defect > DT_(2.5))
      {
        comm.print("FAILED: Final cell size distribution defect should be < "+stringify_fp_fix(2.5)+
            " but is "+stringify_fp_fix(cell_size_defect)+"\n");
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

    watch_total.stop();

    double t_total = watch_total.elapsed();
    double t_asm_mat = watch_asm_mat.elapsed();
    double t_asm_rhs = watch_asm_rhs.elapsed();
    double t_calc_def = watch_calc_def.elapsed();
    double t_sol_init = watch_sol_init.elapsed();
    double t_solver_a = watch_solver_a.elapsed();
    double t_solver_s = watch_solver_s.elapsed();
    double t_vtk = watch_vtk.elapsed();
    double t_meshopt = watch_meshopt.elapsed();
    double t_sum = t_asm_mat + t_asm_rhs + t_calc_def + t_sol_init + t_solver_a + t_solver_s + t_vtk + t_meshopt;

    // write timings
    if(comm.rank() == 0)
    {
      comm.print("");
      dump_time(comm, "Total Solver Time", t_total, t_total);
      dump_time(comm, "Matrix Assembly Time", t_asm_mat, t_total);
      dump_time(comm, "Vector Assembly Time", t_asm_rhs, t_total);
      dump_time(comm, "Defect-Calc Time", t_calc_def, t_total);
      dump_time(comm, "Solver-A Init Time", t_sol_init, t_total);
      dump_time(comm, "Solver-A Time", t_solver_a, t_total);
      dump_time(comm, "Solver-S Time", t_solver_s, t_total);
      dump_time(comm, "VTK-Write Time", t_vtk, t_total);
      dump_time(comm, "Mesh Optimisation Time", t_meshopt, t_total);
      dump_time(comm, "Other Time", t_total-t_sum, t_total);
    }

    TimeStamp bt;
    comm.print("Elapsed time: "+ stringify(bt.elapsed(at)));

    delete[] qi_cellwise;
    delete[] edge_angle_cellwise;

    return ret;

  }
}; // struct NavierStokesScrewsApp

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
  //typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> S2M2D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, 3, Real> S2M3D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, Real> S3M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, 1, Real> H1M1D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, 2, Real> H1M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, 3, Real> H1M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, 3, Real> H2M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, Real> H3M3D;

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
    {
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
    }
  }

  if( args.check("test") >=0 )
  {
    comm.print("Running in test mode, all other command line arguments and configuration files are ignored.");

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
        {
          throw FileNotFound(application_config_filename);
        }

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
        {
          throw FileNotFound(meshopt_config_filename_p.first);
        }

        std::cout << "Reading mesh optimisation config from file " <<meshopt_config_filename_p.first << std::endl;
        synchstream_meshopt_config << ifs.rdbuf();
      }

      // Read solver configuration to stream
      auto solver_config_filename_p = app_settings_section->query("solver_config_file");
      XASSERTM(solver_config_filename_p.second,
      "ApplicationConfig section is missing the mandatory solver_config_file entry!");
      {
        std::ifstream ifs(solver_config_filename_p.first);
        if(ifs.good())
        {
          std::cout << "Reading solver config from file " << solver_config_filename_p.first << std::endl;
          synchstream_solver_config << ifs.rdbuf();
        }
        else
        {
          throw FileNotFound(solver_config_filename_p.first);
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
    read_test_solver_config(synchstream_solver_config, test_number);

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

    comm.print("Reading mesh file "+mesh_files.at(i));
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
    ret = NavierStokesScrewsApp<MemType, DataType, IndexType, H2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"This supports only conformal:hypercube:2:2, but got "+mesh_type);
  }

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
    iss << "delta_t = 1e-3" << std::endl;
    iss << "t_end = 1e-3" << std::endl;
    iss << "midpoint = 0.0 0.0" << std::endl;

    iss << "[DomainControlSettings]" << std::endl;
    iss << "parti-type = fallback parmetis" << std::endl;
    iss << "parti-rank-elems = 4" << std::endl;
    iss << "lvl_min = 0" << std::endl;
    iss << "lvl_max = 2" << std::endl;
  }
  else if(test_number == 2)
  {
    iss << "[ApplicationSettings]" << std::endl;
    iss << "mesh_optimiser = HyperelasticityDefault" << std::endl;
    iss << "solver_config_file = ./solver_config.ini" << std::endl;
    iss << "delta_t = 1e-4" << std::endl;
    iss << "t_end = 0" << std::endl;
    iss << "midpoint = 0.0 0.0" << std::endl;

    iss << "[DomainControlSettings]" << std::endl;
    iss << "parti-type = fallback parmetis" << std::endl;
    iss << "parti-rank-elems = 4" << std::endl;
    iss << "lvl_min = 0" << std::endl;
    iss << "lvl_max = 0" << std::endl;
  }
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"Unknown test number: "+stringify(test_number));
  }
}

static void read_test_meshopt_config(std::stringstream& iss, const int test_number)
{
  if(test_number == 1)
  {
    iss << "[DuDvDefault]" << std::endl;
    iss << "type = DuDv" << std::endl;
    iss << "config_section = DuDvDefaultParameters" << std::endl;
    iss << "fixed_reference_domain = 1" << std::endl;
    iss << "dirichlet_boundaries = bnd:i bnd:o" << std::endl;

    iss << "[DuDvDefaultParameters]" << std::endl;
    iss << "solver_config = PCG-MG" << std::endl;
  }
  else if(test_number == 2)
  {
    iss << "[HyperElasticityDefault]" << std::endl;
    iss << "type = Hyperelasticity" << std::endl;
    iss << "config_section = HyperelasticityDefaultParameters" << std::endl;
    iss << "slip_boundaries = bnd:i bnd:o" << std::endl;

    iss << "[HyperelasticityDefaultParameters]" << std::endl;
    iss << "global_functional = HyperelasticityFunctional" << std::endl;
    iss << "local_functional = RumpfFunctional" << std::endl;
    iss << "solver_config = NLCG" << std::endl;
    iss << "fac_norm = 1e-1" << std::endl;
    iss << "fac_det = 1.0" << std::endl;
    iss << "fac_cof = 0.0" << std::endl;
    iss << "fac_reg = 1e-8" << std::endl;
    iss << "exponent_det = 2" << std::endl;
    iss << "scale_computation = current_concentration" << std::endl;
    iss << "conc_function = GapWidth" << std::endl;

    iss << "[GapWidth]" << std::endl;
    iss << "type = ChartDistance" << std::endl;
    iss << "operation = add" << std::endl;
    iss << "chart_list = screw:i screw:o" << std::endl;
    iss << "function_type = default" << std::endl;

  }
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"Unknown test number "+stringify(test_number));
  }
}

static void read_test_solver_config(std::stringstream& iss, const int test_number)
{
  if(test_number == 1)
  {
    iss << "[PCG-JAC]" << std::endl;
    iss << "type = pcg" << std::endl;
    iss << "max_iter = 1000" << std::endl;
    iss << "tol_rel = 1e-8" << std::endl;
    iss << "precon = jac" << std::endl;

    iss << "[PCG-MG]" << std::endl;
    iss << "type = pcg" << std::endl;
    iss << "max_iter = 75" << std::endl;
    iss << "tol_rel = 1e-8" << std::endl;
    iss << "precon = MG1" << std::endl;
    iss << "plot = 1" << std::endl;

    iss << "[rich]" << std::endl;
    iss << "type = richardson" << std::endl;
    iss << "max_iter = 4" << std::endl;
    iss << "min_iter = 4" << std::endl;
    iss << "precon = jac" << std::endl;

    iss << "[jac]" << std::endl;
    iss << "type = jac" << std::endl;
    iss << "omega = 0.5" << std::endl;

    iss << "[MG1]" << std::endl;
    iss << "type = mg" << std::endl;
    iss << "hierarchy = s:rich-c:pcg" << std::endl;
    iss << "lvl_min = 0" << std::endl;
    iss << "lvl_max = -1" << std::endl;
    iss << "cycle = v" << std::endl;

    iss << "[s:rich-c:pcg]" << std::endl;
    iss << "smoother = rich" << std::endl;
    iss << "coarse = PCG-JAC" << std::endl;
  }
  else if (test_number == 2)
  {
    iss << "[NLCG]" << std::endl;
    iss << "type = NLCG" << std::endl;
    iss << "precon = none" << std::endl;
    iss << "plot = 1" << std::endl;
    iss << "tol_rel = 1e-8" << std::endl;
    iss << "max_iter = 1000" << std::endl;
    iss << "linesearch = MQCLinesearch" << std::endl;
    iss << "direction_update = DYHSHybrid" << std::endl;
    iss << "keep_iterates = 0" << std::endl;

    iss << "[MQCLinesearch]" << std::endl;
    iss << "type = MQCLinesearch" << std::endl;
    iss << "plot = 0" << std::endl;
    iss << "max_iter = 20" << std::endl;
    iss << "tol_decrease = 1e-3" << std::endl;
    iss << "tol_curvature = 0.3" << std::endl;
    iss << "keep_iterates = 0" << std::endl;
  }
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"Encountered unhandled test "+stringify(test_number));
  }
}

static void read_test_mesh_file_names(std::deque<String>& mesh_files, const int test_number)
{
  String mesh_filename(FEAT_SRC_DIR);
  String chart_filename(mesh_filename+"/data/meshes/screws_2d_chart_bezier_24_28.xml");
  mesh_files.push_back(chart_filename);

  if(test_number == 1)
  {
    mesh_filename +="/data/meshes/screws_2d_mesh_quad_360_1.xml";
  }
  else if(test_number == 2)
  {
    mesh_filename +="/data/meshes/screws_2d_mesh_tria_360_1.xml";
  }
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"Encountered unhandled test "+stringify(test_number));
  }

  mesh_files.push_back(mesh_filename);
}

static void display_help(const Dist::Comm& comm)
{
  if(comm.rank() == 0)
  {
    std::cout << "Navier Stokes Screws Application"
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
