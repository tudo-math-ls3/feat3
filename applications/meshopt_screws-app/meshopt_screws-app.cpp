#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

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

#include <control/domain/parti_domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_control_factory.hpp>

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

template<typename Mem_, typename DT_, typename IT_, typename Mesh_>
struct MeshoptScrewsApp
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

  /// Type of the extruded mesh, which is Simplex<2> for Simplex<2> meshes (no extrusion) and Hypercube<3> for
  /// Hypercube<2> meshes
  typedef typename MeshExtrudeHelper<MeshType>::ExtrudedMeshType ExtrudedMeshType;

  /// Type for points in the mesh
  typedef Tiny::Vector<DataType, MeshType::world_dim> WorldPoint;

  /// Domain Control Type
  typedef Control::Domain::PartiDomainControl<MeshType> DomCtrl;

  /// This is how far the inner screw's centre deviates from the outer screw's
  static constexpr DataType excentricity_inner = DataType(0.2833);

  /**
   * \brief Returns a descriptive string
   *
   * \returns The class name as string
   */
  static String name()
  {
    return "MeshoptScrewsApp";
  }

  /**
   * \brief The routine that does the actual work
   */
  static int run(const SimpleArgParser& args, Dist::Comm& comm, PropertyMap* application_config,
    PropertyMap* meshopt_config, PropertyMap* solver_config, Geometry::MeshFileReader& mesh_file_reader)
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
    const auto& finest_mesh = dom_ctrl.get_levels().back()->get_mesh();

    // Print level information
    comm.print(name()+" settings:");
    comm.print("LVL-MAX "+stringify(dom_ctrl.get_levels().back()->get_level_index())
        +" [" +stringify(lvl_max) + "] "
        +"LVL-MIN "+stringify(dom_ctrl.get_levels().front()->get_level_index())+" [" +stringify(lvl_min) + "]");
    comm.print("Timestep size: "+stringify_fp_fix(delta_t)+", end time: "+ stringify_fp_fix(t_end));
    dom_ctrl.print();

    MeshExtrudeHelper<MeshType> extruder(dom_ctrl.get_levels().back()->get_mesh_node(),
    Index(10*(lvl_max+1)), DataType(0), DataType(1), "bnd:b", "bnd:t");

    // Get inner boundary MeshPart. Can be nullptr if this process' patch does not lie on that boundary
    auto* inner_boundary = dom_ctrl.get_levels().back()->get_mesh_node()->find_mesh_part("bnd:i");
    Geometry::TargetSet* inner_indices(nullptr);
    if(inner_boundary != nullptr)
    {
      inner_indices = &(inner_boundary->template get_target_set<0>());
    }

    // This is the centre point of the rotation of the inner screw
    WorldPoint centre_inner(DataType(0));
    centre_inner.v[0] = -excentricity_inner;
    auto* inner_chart = dom_ctrl.get_atlas()->find_mesh_chart("screw:i");

    // Get outer boundary MeshPart. Can be nullptr if this process' patch does not lie on that boundary
    auto* outer_boundary_part = dom_ctrl.get_levels().back()->get_mesh_node()->find_mesh_part("bnd:o");
    Geometry::TargetSet* outer_indices(nullptr);
    if(outer_boundary_part != nullptr)
    {
      outer_indices = &(outer_boundary_part->template get_target_set<0>());
    }

    // This is the centre point of the rotation of the outer screw
    WorldPoint centre_outer(DataType(0));
    auto* outer_chart = dom_ctrl.get_atlas()->find_mesh_chart("screw:o");

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
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());

        String vtk_name = String(file_basename+"_pre_initial_lvl_"+stringify(lvl_index));
        comm.print("Writing "+vtk_name);

        // Compute mesh quality on this level
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise, lvl_index);
        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));

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
    if( (test_number == 1) && ( edge_angle < DT_(5.69)) )
    {
      comm.print("FAILED: Initial worst edge angle should be >= "+stringify_fp_fix(5.69)
          + " but is "+stringify_fp_fix(edge_angle));
      ++ret;
    }
    else if( (test_number == 2) && ( edge_angle < DT_(2.83)) )
    {
      comm.print("FAILED: Initial worst edge angle should be >= "+stringify_fp_fix(2.83)
          + " but is "+stringify_fp_fix(edge_angle));
    }

    // Optimise the mesh
    meshopt_ctrl->optimise();

    // Write output again
    if(write_vtk)
    {
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());

        String vtk_name = String(file_basename+"_post_initial_lvl_"+stringify(lvl_index));
        comm.print("Writing "+vtk_name);

        // Compute mesh quality on this level
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise, lvl_index);

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));

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
    if(test_number == 1)
    {
      if(edge_angle < DT_(5.69))
      {
        comm.print("FAILED: Post Initial worst edge angle should be >= "+stringify_fp_fix(5.69)+
            " but is "+stringify_fp_fix(edge_angle));
        ++ret;
      }
      if(qi_min < DT_(2.5e-2))
      {
        comm.print("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(2.5e-2)+
            " but is "+stringify_fp_fix(qi_min));
        ++ret;
      }
      if(cell_size_defect > DT_(2.5))
      {
        comm.print("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(2.5)+
            " but is "+stringify_fp_fix(cell_size_defect));
        ++ret;
      }
    }
    else if(test_number == 2)
    {
      if(edge_angle < DT_(3.2))
      {
        comm.print("FAILED: Post Initial worst edge angle should be >= "+stringify_fp_fix(3.2)+
            " but is "+stringify_fp_fix(edge_angle));
        ++ret;
      }
      if(qi_min < DT_(2.5e-1))
      {
        comm.print("FAILED: Post Initial worst shape quality should be >= "+stringify_fp_fix(2.5e-1)+
            " but is "+stringify_fp_fix(qi_min));
        ++ret;
      }
      if(cell_size_defect > DT_(3.8e-1))
      {
        comm.print("FAILED: Post Initial cell size distribution defect should be <= "+stringify_fp_fix(3.8e-1)+
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

    // This is the absolute turning angle of the screws
    DataType alpha(0);
    // Need some pi for all the angles
    DataType pi(Math::pi<DataType>());

    WorldPoint rotation_angles(0);

    while(time < t_end)
    {
      n++;
      time+= delta_t;

      // Clear statistics data so it does not eat us alive
      FEAT::Statistics::reset_solver_statistics();

      DataType alpha_old = alpha;
      alpha = -DataType(2)*pi*time;
      DataType delta_alpha = alpha - alpha_old;

      bool abort(false);

      comm.print("Timestep "+stringify(n)+": t = "+stringify_fp_fix(time)+", angle = "
          +stringify_fp_fix(alpha/(DataType(2)*pi)*DataType(360)) + " degrees");

      // Save old vertex coordinates
      meshopt_ctrl->mesh_to_buffer();
      old_coords.copy(meshopt_ctrl->get_coords());
      new_coords.copy(meshopt_ctrl->get_coords());
      auto& coords_loc(*new_coords);

      // Rotate the charts
      rotation_angles(0) = delta_alpha;
      rotation_angles(1) = DataType(0);

      outer_chart->rotate(centre_outer, rotation_angles);

      rotation_angles(0) *= DataType(7)/DataType(6);
      inner_chart->rotate(centre_inner, rotation_angles);

      // Update boundary of the inner screw
      // This is the 2x2 matrix representing the turning by the angle delta_alpha of the inner screw
      if(inner_indices != nullptr)
      {
        Tiny::Matrix<DataType, 2, 2> rot(DataType(0));

        rot(0,0) = Math::cos(delta_alpha*DataType(7)/DataType(6));
        rot(0,1) = - Math::sin(delta_alpha*DataType(7)/DataType(6));
        rot(1,0) = -rot(0,1);
        rot(1,1) = rot(0,0);

        WorldPoint tmp(DataType(0));
        WorldPoint tmp2(DataType(0));

        for(Index i(0); i < inner_indices->get_num_entities(); ++i)
        {
          // Index of boundary vertex i in the mesh
          Index j(inner_indices->operator[](i));
          // Translate the point to the centre of rotation
          tmp = coords_loc(j) - centre_inner;
          // Rotate
          tmp2.set_vec_mat_mult(tmp, rot);
          // Translate the point by the new centre of rotation
          coords_loc(j, centre_inner + tmp2);
        }
      }

      // The outer screw has 7 teeth as opposed to the inner screw with 6, and it rotates at 6/7 of the speed
      if(outer_indices != nullptr)
      {
        Tiny::Matrix<DataType, 2, 2> rot(DataType(0));
        rot(0,0) = Math::cos(delta_alpha);
        rot(0,1) = - Math::sin(delta_alpha);
        rot(1,0) = -rot(0,1);
        rot(1,1) = rot(0,0);

        WorldPoint tmp(DataType(0));
        WorldPoint tmp2(DataType(0));

        // The outer screw rotates centrically, so centre_outer remains the same at all times
        for(Index i(0); i < outer_indices->get_num_entities(); ++i)
        {
          // Index of boundary vertex i in the mesh
          Index j(outer_indices->operator[](i));
          tmp = coords_loc(j) - centre_outer;

          tmp2.set_vec_mat_mult(tmp, rot);

          coords_loc(j, centre_outer+tmp2);
        }
      }

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
        String msg = String("max. mesh velocity").pad_back(pad_width, ' ') + ": " + stringify_fp_sci(max_mesh_velocity);
        comm.print(msg);
      }

      // Write output on the finest level only
      if(write_vtk && ( (n%vtk_freq == 0) || abort ))
      {
        // Compute mesh quality on the finest
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise);

        String vtk_name = String(file_basename+"_post_"+stringify(n));
        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_ctrl.get_levels().back()->get_mesh());

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise);
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise);

        meshopt_ctrl->add_to_vtk_exporter(exporter, int(dom_ctrl.get_levels().size())-1);

        exporter.write(vtk_name, comm);

        if(extruder.extruded_mesh_node != nullptr)
        {
          extruder.extrude_vertex_set(finest_mesh.get_vertex_set());
          // Create a VTK exporter for our mesh
          String extruded_vtk_name = String(file_basename+"_post_extruded_"+stringify(n));

          comm.print("Writing "+extruded_vtk_name);

          Geometry::ExportVTK<ExtrudedMeshType> extruded_exporter(*(extruder.extruded_mesh_node->get_mesh()));
          extruded_exporter.write(extruded_vtk_name, comm);
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

      if(abort)
      {
        break;
      }

    } // time loop

    meshopt_ctrl->print();

    // Write final vtk output
    if(write_vtk)
    {
      for(auto it = dom_ctrl.get_levels().begin(); it !=  dom_ctrl.get_levels().end(); ++it)
      {
        int lvl_index((*it)->get_level_index());

        String vtk_name = String(file_basename+"_final_lvl_"+stringify(lvl_index));
        comm.print("Writing "+vtk_name);

        // Compute mesh quality on this level
        dom_ctrl.compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise, qi_cellwise, lvl_index);

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(((*it)->get_mesh()));

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise);
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise);

        meshopt_ctrl->add_to_vtk_exporter(exporter, lvl_index);
        exporter.write(vtk_name, comm.rank(), comm.size());
      }
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


    TimeStamp bt;
    comm.print("Elapsed time: "+ stringify(bt.elapsed(at)));

    delete[] qi_cellwise;
    delete[] edge_angle_cellwise;

    return ret;

  }
}; // struct MeshoptScrewsApp

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
    ret = MeshoptScrewsApp<MemType, DataType, IndexType, H2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  else if(mesh_type == "conformal:simplex:2:2")
  {
    ret = MeshoptScrewsApp<MemType, DataType, IndexType, S2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  else
  {
    throw InternalError(__func__,__FILE__,__LINE__,"Unhandled mesh type "+mesh_type);
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
    iss << "linesearch = StrongWolfeLinesearch" << std::endl;
    iss << "direction_update = DYHSHybrid" << std::endl;
    iss << "keep_iterates = 0" << std::endl;

    iss << "[strongwolfelinesearch]" << std::endl;
    iss << "type = StrongWolfeLinesearch" << std::endl;
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
