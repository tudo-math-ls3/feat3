// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/assembly/analytic_projector.hpp>
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/fe_interpolator.hpp>
#include <kernel/assembly/grad_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/geometry/common_factories.hpp>
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
#include <control/time/nvs_bdf_q.hpp>

// For validating the time discretisation with analytic solutions
//#define ANALYTIC_SOLUTION
using namespace FEAT;

static void display_help(const Dist::Comm&);
static void read_test_application_config(std::stringstream&);
static void read_test_meshopt_config(std::stringstream&);
static void read_test_solver_config(std::stringstream&);
static void read_test_mesh_file_names(std::deque<String>&);

template<typename Mesh_, bool extrude>
struct MeshExtrudeHelper
{
  typedef Mesh_ MeshType;
  typedef Mesh_ ExtrudedMeshType;
  typedef typename MeshType::CoordType CoordType;
  typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
  typedef Geometry::RootMeshNode<ExtrudedMeshType> ExtrudedMeshNodeType;
  typedef Geometry::MeshAtlas<ExtrudedMeshType> ExtrudedAtlasType;

  static void create(
    std::shared_ptr<ExtrudedMeshNodeType>& extruded_mesh_node,
    const std::shared_ptr<MeshNodeType> rmn,
    Index DOXY(slices), CoordType DOXY(z_min), CoordType DOXY(z_max), const String& DOXY(z_min_part_name), const String& DOXY(z_max_part_name))
    {
      extruded_mesh_node = rmn;
    }

  static void clear(ExtrudedMeshNodeType* DOXY(extruded_mesh_node))
  {
  }
};

template<typename Coord_>
struct MeshExtrudeHelper<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Coord_>, true>
{
  typedef Coord_ CoordType;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Coord_> MeshType;
  typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Coord_> ExtrudedMeshType;
  typedef Geometry::RootMeshNode<ExtrudedMeshType> ExtrudedMeshNodeType;
  typedef Geometry::MeshAtlas<ExtrudedMeshType> ExtrudedAtlasType;

  static void create(
    std::shared_ptr<ExtrudedMeshNodeType>& extruded_mesh_node,
    const std::shared_ptr<MeshNodeType> rmn,
    Index slices, CoordType z_min, CoordType z_max, const String& z_min_part_name, const String& z_max_part_name)
    {
      XASSERT(extruded_mesh_node == nullptr);

      ExtrudedAtlasType* extruded_atlas(new ExtrudedAtlasType);
      extruded_mesh_node = std::make_shared<ExtrudedMeshNodeType>(nullptr, extruded_atlas);

      Geometry::MeshExtruder<MeshType> mesh_extruder(slices, z_min, z_max, z_min_part_name, z_max_part_name);

      mesh_extruder.extrude_atlas(*extruded_atlas, *(rmn->get_atlas()));
      mesh_extruder.extrude_root_node(*extruded_mesh_node, *rmn, extruded_atlas);
    }

  static void clear(ExtrudedMeshNodeType* extruded_mesh_node)
  {
    delete extruded_mesh_node->get_atlas();
  }
};

static inline void dump_time(const Dist::Comm& comm, String s, double t, double total)
{
  comm.print(s.pad_back(30, '.') + ": " + stringify_fp_fix(t, 3, 12)
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
  public Control::StokesBlockedSlipUnitVeloMeanPresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
{
  public:
    typedef Control::StokesBlockedSlipUnitVeloMeanPresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;
}; // class NavierStokesBlockedSystemLevel
//class NavierStokesBlockedSystemLevel :
//  public Control::StokesBlockedUnitVeloMeanPresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
//{
//  public:
//    typedef Control::StokesBlockedUnitVeloMeanPresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;
//
//}; // class NavierStokesBlockedSystemLevel

template<typename DomainLevel_, typename DomCtrl2d_, bool extrude>
class ExtrudedPartiDomainControl :
  public Control::Domain::DomainControl<DomainLevel_>
{
  public:
    /// Our base class
    typedef Control::Domain::DomainControl<DomainLevel_> BaseClass;
    typedef DomCtrl2d_ DomCtrl2d;
    typedef typename DomCtrl2d_::MeshType Mesh2d;
    /// our domain level type
    typedef DomainLevel_ LevelType;
    /// our domain layer type
    typedef typename BaseClass::LayerType LayerType;
    /// our mesh type
    typedef typename BaseClass::MeshType MeshType;
    /// our atlas type
    typedef typename BaseClass::AtlasType AtlasType;
    /// our root mesh node type
    typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

    typedef typename MeshType::CoordType CoordType;

    typedef MeshExtrudeHelper<Mesh2d, extrude> MeshExtruder;

    std::vector<Index> _vtx_map;
    const DomCtrl2d& _dom_ctrl_2d;

  public:
    /// default constructor
    explicit ExtrudedPartiDomainControl(const Dist::Comm& comm_, const DomCtrl2d& dom_ctrl_2d,
    Index slices, CoordType z_min, CoordType z_max, const String& z_min_part_name, const String& z_max_part_name ) :
      BaseClass(comm_),
      _vtx_map(),
      _dom_ctrl_2d(dom_ctrl_2d)
      {
        //XASSERTM(!this->_have_hierarchy, "domain control already has a hierarchy!");

        this->_layers.push_back(std::make_shared<LayerType>(comm_.comm_dup(), 0));
        this->_layer_levels.resize(std::size_t(1));

        auto coarse_mesh_node_2d = _dom_ctrl_2d.back()->get_mesh_node_ptr();
        int lvl_min(_dom_ctrl_2d.min_level_index());
        int lvl_max(_dom_ctrl_2d.max_level_index());

        // Create coarse mesh node for this domain control
        std::shared_ptr<MeshNodeType> mesh_node(nullptr);
        MeshExtruder::create(mesh_node, coarse_mesh_node_2d, slices, z_min, z_max, z_min_part_name, z_max_part_name);
        // Copy the neighbours information from the lower dimensional domain control
        this->_layers.front()->set_neighbour_ranks(_dom_ctrl_2d.front().layer().get_neighbour_ranks());

        // Add coarse mesh node to layer_levels
        auto& laylevs = this->_layer_levels.front();
        laylevs.push_front(std::make_shared<LevelType>(lvl_min, mesh_node));

        // Refine up to desired maximum level
        for(int lvl(lvl_min); lvl < lvl_max; ++lvl)
        {
          std::shared_ptr<MeshNodeType> coarse_node = mesh_node;
          mesh_node = std::shared_ptr<MeshNodeType>(coarse_node->refine(/*this->_adapt_mode*/));

          laylevs.push_front(std::make_shared<LevelType>(lvl+1, mesh_node));
        }

        // Finally, compile the virtual levels
        this->compile_virtual_levels();

        // Build vertex mapping by brute force
        const auto& finest_mesh_2d = _dom_ctrl_2d.front()->get_mesh();
        const auto& finest_extruded_mesh = laylevs.front()->get_mesh();

        const auto& vtx_2d = finest_mesh_2d.get_vertex_set();
        const auto& extruded_vtx = finest_extruded_mesh.get_vertex_set();

        Index num_verts(finest_mesh_2d.get_num_entities(0));
        Index extruded_num_verts(finest_extruded_mesh.get_num_entities(0));

        _vtx_map.resize(extruded_num_verts);

        const CoordType tol(Math::eps<CoordType>());

        for(Index i_extruded(0); i_extruded < extruded_num_verts; ++i_extruded)
        {
          bool found(false);
          Tiny::Vector<CoordType, 2> tmp(0);
          tmp(0) = extruded_vtx[i_extruded](0);
          tmp(1) = extruded_vtx[i_extruded](1);
          CoordType min_dist(Math::huge<CoordType>());
          //Index best_index(0);
          for(Index j(0); j < num_verts; ++j)
          {
            CoordType my_dist((tmp - vtx_2d[j]).norm_euclid_sqr());
            if(my_dist < min_dist)
            {
              //best_index = j;
              min_dist = my_dist;
            }
            if(my_dist <= tol )
            {
              found = true;
              _vtx_map.at(i_extruded) = j;
              break;
            }
          }
          if(!found)
          {
            XABORTM("Could not find 2d point for "+stringify(extruded_vtx[i_extruded]));
          }
        }

      }

    /// virtual destructor
    virtual ~ExtrudedPartiDomainControl()
    {
      MeshExtruder::clear(this->front()->get_mesh_node());
    }

    /**
     * \brief Extrudes the 2d vertex set to the 3d mesh
     */
    void extrude_vertex_sets()
    {

      for(Index l(0); l < this->size_physical(); ++l)
      {
        const auto& vtx = _dom_ctrl_2d.at(l)->get_mesh().get_vertex_set();
        auto& extruded_vtx = this->at(l)->get_mesh().get_vertex_set();

        for(Index i(0); i < extruded_vtx.get_num_vertices(); ++i)
        {
          Index j(_vtx_map[i]);
          for(int d(0); d < Mesh2d::world_dim; ++d)
          {
            extruded_vtx[i][d] = vtx[j](d);
          }
        }
      }
    }

    /**
     * \brief Extrudes a 2d vertex vector to 3d
     */
    template<typename DT_, typename IT_>
    void extrude_vertex_vector(LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, MeshType::world_dim>& v,
    const LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, Mesh2d::world_dim>& v_in)
    {
      Tiny::Vector<DT_, MeshType::world_dim> tmp(DT_(0));

      for(Index i(0); i < v.size(); ++i)
      {
        Index j(_vtx_map[i]);
        for(int d(0); d < Mesh2d::world_dim; ++d)
        {
          tmp(d) = v_in(j)(d);
        }
        v(i, tmp);
      }
    }
};

template<typename Mem_, typename DT_, typename IT_, typename Mesh_>
struct NavierStokesScrewsApp
{
  static constexpr bool extrude = false;

  /// The memory architecture. Although this looks freely chosable, it has to be Mem::Main for now because all the
  /// Hyperelasticity functionals are implemented for Mem::Main only
  typedef Mem_ MemType;
  /// The floating point type
  typedef DT_ DataType;
  /// The index type
  typedef IT_ IndexType;
  /// The type of mesh to use
  typedef Mesh_ MeshType;

  /// Type for points in the mesh
  typedef Tiny::Vector<DataType, MeshType::world_dim> WorldPoint;

  /// Type of the extruded mesh, which is Simplex<2> for Simplex<2> meshes (no extrusion) and Hypercube<3> for
  /// Hypercube<2> meshes
  typedef MeshExtrudeHelper<Mesh_, extrude> MeshExtruder;

  /// The only transformation available is the standard P1 or Q1 transformation
  typedef Trafo::Standard::Mapping<Mesh_> TrafoType;
  /// FE space for the transformation. The mesh optimisation problem is solved on this
  typedef typename Meshopt::Intern::TrafoFE<TrafoType>::Space TrafoFESpace;

  /// The domain level, including trafo and FE space
  typedef Control::Domain::SimpleDomainLevel<Mesh_, TrafoType, TrafoFESpace> DomainLevelType;

  typedef Control::Domain::PartiDomainControl<DomainLevelType> DomCtrl;

  typedef typename MeshExtruder::ExtrudedMeshType ExtrudedMeshType;
  /// Type for points in the mesh
  typedef Tiny::Vector<DataType, ExtrudedMeshType::world_dim> ExtrudedWorldPoint;
  /// The shape type of the mesh's cells
  typedef typename ExtrudedMeshType::ShapeType ExtrudedShapeType;

  /// The only transformation available is the standard P1 or Q1 transformation
  typedef Trafo::Standard::Mapping<ExtrudedMeshType> ExtrudedTrafoType;
  /// FE space for the velocity
  typedef Space::Lagrange2::Element<ExtrudedTrafoType> SpaceVeloType;
  /// FE space for the pressure
  typedef Space::Lagrange1::Element<ExtrudedTrafoType> SpacePresType;
  /// For the extruded mesh velocity, we will need this
  typedef typename Meshopt::Intern::TrafoFE<ExtrudedTrafoType>::Space ExtrudedTrafoFESpace;

  /// The domain level with FE spaces for velocity and pressure
  typedef Control::Domain::StokesDomainLevel
  <
    ExtrudedMeshType,
    ExtrudedTrafoType,
    SpaceVeloType, SpacePresType
  > ExtrudedDomainLevelType;
  /// Domain Control Type
  typedef ExtrudedPartiDomainControl<ExtrudedDomainLevelType, DomCtrl, extrude> ExtrudedDomCtrl;


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
  static int run(const SimpleArgParser& args, Dist::Comm& comm, PropertyMap& application_config,
  PropertyMap& meshopt_config, PropertyMap& solver_config, Geometry::MeshFileReader& mesh_file_reader)
  {

    // Create a time-stamp
    TimeStamp stamp_start;

    // create a batch of stop-watches
    StopWatch watch_total, watch_asm_rhs, watch_asm_mat, watch_asm_fil,
      watch_sol_init, watch_solver_a, watch_solver_s, watch_solver_m_p, watch_vtk,
      watch_meshopt, watch_quality, watch_meshopt_preproc;

    watch_total.start();

    static constexpr int pad_width = 30;

    int failed_checks(0);

    // Mininum refinement level, parsed from the application config file
    int lvl_min(-1);
    // Maximum refinement level, parsed from the application config file
    int lvl_max(-1);
    DataType z_min(0);
    DataType z_max(1);
    Index slices(1);
    // End time, parsed from the application config file
    DataType t_end(0);
    // Do we want to write vtk files. Read from the command line arguments
    bool write_vtk(false);
    // If write_vtk is set, we write out every vtk_freq time steps
    Index vtk_freq(1);
    // Do we want to write xml files. Read from the command line arguments
    bool write_xml(false);
    // If write_xml is set, we write out every xml_freq time steps
    Index xml_freq(1);
    // Is the application running as a test? Read from the command line arguments
    bool test(false);

    // Solve the flow problem?
    bool solve_flow(false);
    // Optimise the mesh in every time step?
    bool solve_mesh_optimisation(false);

    // The Reynolds number for the flow problem
    DataType reynolds(1);
    // Use the deformation tensor-based bilinear form for the viscous term in the Navier-Stokes equations?
    bool use_deformation(false);

    // Need some pi for all the angles
    DataType pi(Math::pi<DataType>());

    // Check if we want to write vtk files and at what frequency
    if(args.check("vtk") >= 0 )
    {
      write_vtk = true;

      if(args.check("vtk") > 1)
      {
        XABORTM("Too many options for --vtk");
      }

      args.parse("vtk",vtk_freq);
      XASSERT(vtk_freq > Index(0));
    }

    // Check if we want to write xml files and at what frequency
    if(args.check("xml") >= 0 )
    {
      write_xml = true;

      if(args.check("xml") > 1)
      {
        XABORTM("Too many options for --xml");
      }

      args.parse("xml",xml_freq);
      XASSERT(xml_freq > Index(0));
    }

    // Check if we are running in test mode
    if( args.check("test") >= 0 )
    {
      test = true;
    }

    // Get the application settings section
    auto app_settings_section = application_config.query_section("ApplicationSettings");
    XASSERTM(app_settings_section != nullptr,
    "Application config is missing the mandatory ApplicationSettings section!");

    // Get the mesh optimiser key from the application settings
    auto meshoptimiser_key_p = app_settings_section->query("mesh_optimiser");
    XASSERTM(meshoptimiser_key_p.second,
    "ApplicationConfig section is missing the mandatory meshoptimiser entry!");

    // End time
    auto t_end_p = app_settings_section->query("t_end");
    XASSERTM(t_end_p.second, "Application config section is missing the mandatory t_end entry!");
    t_end = DataType(std::stod(t_end_p.first));

    // Solve the mesh_optimisation problem?
    auto solve_mesh_optimisation_p = app_settings_section->query("solve_mesh_optimisation");
    XASSERTM(solve_mesh_optimisation_p.second,
    "Application config section is missing the mandatory solve_mesh_optimisation entry!");
    solve_mesh_optimisation = (std::stoi(solve_mesh_optimisation_p.first) == 1);

    // Solve the flow problem?
    auto solve_flow_p = app_settings_section->query("solve_flow");
    XASSERTM(solve_flow_p.second, "Application config section is missing the mandatory solve_flow entry!");
    solve_flow = (std::stoi(solve_flow_p.first) == 1);

    if(solve_flow)
    {
      // Get reynolds number
      auto reynolds_p = app_settings_section->query("reynolds");
      XASSERTM(reynolds_p.second, "ApplicationConfig section is missing the mandatory reynolds entry!");
      reynolds = std::stod(reynolds_p.first);
      XASSERT(reynolds > DataType(0));

      // Use the D(u) : D(v) bilinear form instead of the grad(u):grad(v) bilinear form?
      auto use_deformation_p = app_settings_section->query("use_deformation");
      XASSERTM(use_deformation_p.second,
      "ApplicationSettings section is missing the mandatory use_deformation entry!");
      XASSERTM(std::stoi(use_deformation_p.first) == 0 || std::stoi(use_deformation_p.first) == 1,
      "use_deformation has to be set to 0 or 1");
      use_deformation = std::stoi(use_deformation_p.first) == 1;
    }

    // Get the domain control settings section
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

    auto z_min_p = domain_control_settings_section->query("z_min");
    if(z_min_p.second)
    {
      z_min = DataType(std::stod(z_min_p.first));
    }

    auto z_max_p = domain_control_settings_section->query("z_max");
    if(z_max_p.second)
    {
      z_max = DataType(std::stod(z_max_p.first));
    }

    auto slices_p = domain_control_settings_section->query("slices");
    if(slices_p.second)
    {
      slices = Index(std::stoul(slices_p.first));
    }

    // Get the time discretisation settings section
    auto time_disc_settings = application_config.query_section("TimeDiscretisation");
    XASSERTM(time_disc_settings!= nullptr,
    "Application config is missing the mandatory TimeDiscretisation section!");

    // Create BDF time discretisation coefficients in startup mode
    Control::Time::NvsBdfQ<DataType> time_disc(time_disc_settings, reynolds, use_deformation, true);

    // Create domain control
    DomCtrl dom_ctrl(comm, false);
    dom_ctrl.parse_property_map(*domain_control_settings_section);
    dom_ctrl.set_desired_levels(lvl_max, lvl_min);
    dom_ctrl.create(mesh_file_reader);

    ExtrudedDomCtrl extruded_dom_ctrl(comm, dom_ctrl, slices, z_min, z_max, "bnd:b", "bnd:t");

    // Mesh on the finest level, mainly for computing quality indicators
    const auto& finest_extruded_mesh = extruded_dom_ctrl.front()->get_mesh();

    // Print configuration
    {
      String msg("");

      comm.print(name()+" settings:");

      msg = String("Maximum Level").pad_back(pad_width, '.') + String(": "+ stringify(dom_ctrl.max_level_index()))
        +String(" [" +stringify(lvl_max) + "]");
      comm.print(msg);
      msg = String("Minimum Level").pad_back(pad_width, '.') + String(": "+ stringify(dom_ctrl.min_level_index()))
        +String(" [" +stringify(lvl_min) + "]");
      comm.print(msg);

      msg = String("delta_t").pad_back(pad_width, '.') + String(": "+ stringify_fp_fix(time_disc.delta_t()));
      comm.print(msg);
      msg = String("t_end").pad_back(pad_width, '.') + String(": "+ stringify_fp_fix(t_end));
      comm.print(msg);

      if(solve_flow)
      {
        msg = String("Reynolds").pad_back(pad_width, '.') + String(": "+ stringify_fp_fix(reynolds));
        comm.print(msg);

        msg = String("ALE term").pad_back(pad_width, '.') + String(": "+stringify(time_disc.ale_handling));
        comm.print(msg);
        msg = String("Convection term").pad_back(pad_width, '.') + String(": "+stringify(time_disc.conv_handling));
        comm.print(msg);
        msg = String("Viscous term").pad_back(pad_width, '.') + String(": "+stringify(time_disc.visc_handling));
        comm.print(msg);
      }

      dom_ctrl.print();
    }

    // Get outer boundary MeshPart. Can be nullptr if this process' patch does not lie on that boundary
    auto* outer_boundary_part = dom_ctrl.front()->get_mesh_node()->find_mesh_part("bnd:o");
    Geometry::TargetSet* outer_indices(nullptr);
    if(outer_boundary_part != nullptr)
    {
      outer_indices = &(outer_boundary_part->template get_target_set<0>());
    }

    // This is the centre point of the rotation of the outer screw
    WorldPoint centre_outer(DataType(0));

    ExtrudedWorldPoint extruded_centre_outer(DataType(0));
    for(int d(0); d < MeshType::world_dim; ++d)
    {
      extruded_centre_outer(d) = centre_outer(d);
    }

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
    centre_inner(0) = -excentricity_inner;

    ExtrudedWorldPoint extruded_centre_inner(DataType(0));
    for(int d(0); d < MeshType::world_dim; ++d)
    {
      extruded_centre_inner(d) = centre_inner(d);
    }

    DT_ angular_velocity_inner(angular_velocity_outer*DT_(7)/DT_(6));
    auto* inner_chart = dom_ctrl.get_atlas().find_mesh_chart("screw:i");

    // Create MeshoptControl
    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl>> meshopt_ctrl(nullptr);

    meshopt_ctrl = Control::Meshopt::ControlFactory<Mem_, DT_, IT_>::create_meshopt_control(
      dom_ctrl, meshoptimiser_key_p.first, &meshopt_config, &solver_config);

    std::shared_ptr<Control::Meshopt::MeshoptControlBase<DomCtrl>> meshopt_preproc(nullptr);
    // Get mesh optimiser settings section
    auto meshopt_settings_section = meshopt_config.query_section(meshoptimiser_key_p.first);
    if(meshopt_settings_section != nullptr)
    {
      auto preproc_p = meshopt_settings_section->query("preprocessor_config");
      if(preproc_p.second)
      {
        meshopt_preproc = Control::Meshopt::ControlFactory<Mem_, DT_, IT_>::create_meshopt_control(
          dom_ctrl, preproc_p.first, &meshopt_config, &solver_config);
        XASSERT(meshopt_preproc != nullptr);
        comm.print("Using "+preproc_p.first+" for meshopt preprocessing");

      }
    }

    comm.print(String("\n") + meshopt_ctrl->info() + String("\n"));

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
    std::vector<DataType> qi_cellwise(finest_extruded_mesh.get_num_entities(ExtrudedMeshType::shape_dim));

    DT_ edge_angle(0);
    std::vector<DataType> edge_angle_cellwise(finest_extruded_mesh.get_num_entities(ExtrudedMeshType::shape_dim));

    DT_ cell_size_defect(0);

    // Write initial vtk output
    if(write_vtk)
    {
      watch_vtk.start();
      for(size_t l(0); l < dom_ctrl.size_physical(); ++l)
      {
        auto& dom_lvl = dom_ctrl.at(l);

        int lvl_index(dom_lvl->get_level_index());

        String vtk_name = String("mesh_pre_initial_n"+stringify(comm.size())+"_lvl_"+stringify(lvl_index));
        comm.print("Writing "+vtk_name+".vtk");

        // Compute mesh quality on this level
        meshopt_ctrl->compute_mesh_quality(
          edge_angle, qi_min, qi_mean, edge_angle_cellwise.data(), qi_cellwise.data(), lvl_index);

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_lvl->get_mesh());

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise.data());
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise.data());

        meshopt_ctrl->add_to_vtk_exporter(exporter, lvl_index);

        exporter.write(vtk_name, comm.rank(), comm.size());
      }
      watch_vtk.stop();
    }

    // Compute and print quality indicators on the finest level only
    {
      watch_quality.start();
      DT_ lambda_min(Math::huge<DT_>());
      DT_ lambda_max(0);
      DT_ vol(0);
      DT_ vol_min(Math::huge<DT_>());
      DT_ vol_max(0);

      cell_size_defect = meshopt_ctrl->compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max, vol);

      meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise.data(), qi_cellwise.data());

      String msg("");

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

      watch_quality.stop();
    }

    // Check for the hard coded settings for test mode
    if(test)
    {
      if( edge_angle < DT_(28))
      {
        comm.print("FAILED: Initial worst edge angle should be >= "+stringify_fp_fix(28)
            + " but is "+stringify_fp_fix(edge_angle));
        ++failed_checks;
      }
      if( qi_min < DT_(3e-2))
      {
        comm.print("FAILED: Initial minimal quality indicator should be >= "+stringify_fp_fix(3e-2)
            + " but is "+stringify_fp_fix(edge_angle));
        ++failed_checks;
      }
      if( cell_size_defect > DT_(6.8e-1))
      {
        comm.print("FAILED: Initial cell size defect should be <= "+stringify_fp_fix(6.8e-1)
            + " but is "+stringify_fp_fix(edge_angle));
        ++failed_checks;
      }
    }

    if(solve_mesh_optimisation)
    {
      // Optimise the mesh
      FEAT::Statistics::expression_target = "meshopt_preproc";
      watch_meshopt.start();
      meshopt_ctrl->optimise();
      extruded_dom_ctrl.extrude_vertex_sets();
      watch_meshopt.stop();

      // Compute and print quality indicators on the finest level only

      watch_quality.start();
      DT_ lambda_min(Math::huge<DT_>());
      DT_ lambda_max(0);
      DT_ vol(0);
      DT_ vol_min(Math::huge<DT_>());
      DT_ vol_max(0);

      cell_size_defect = meshopt_ctrl->compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max, vol);

      meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise.data(), qi_cellwise.data());

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
      watch_quality.stop();


      // Check for the hard coded settings for test mode
      if(test)
      {
        if( edge_angle < DT_(28))
        {
          comm.print("FAILED: Optimised worst edge angle should be >= "+stringify_fp_fix(28)
              + " but is "+stringify_fp_fix(edge_angle));
          ++failed_checks;
        }
        if( qi_min < DT_(3e-2))
        {
          comm.print("FAILED: Optimised minimal quality indicator should be >= "+stringify_fp_fix(3e-2)
              + " but is "+stringify_fp_fix(edge_angle));
          ++failed_checks;
        }
        if( cell_size_defect > DT_(6.8e-1))
        {
          comm.print("FAILED: Optimised cell size defect should be <= "+stringify_fp_fix(6.8e-1)
              + " but is "+stringify_fp_fix(edge_angle));
          ++failed_checks;
        }
      }
    }

    // define our velocity and pressure system levels
    typedef NavierStokesBlockedSystemLevel
    <
      ExtrudedShapeType::dimension, MemType, DataType, IndexType
    > SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = dom_ctrl.size_physical();

    // create stokes and system levels
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("auto-degree:7");

    watch_asm_mat.start();
    /* ***************************************************************************************** */

    comm.print("Assembling gates...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gates(extruded_dom_ctrl.at(i));
    }

    /* ***************************************************************************************** */

    comm.print("Assembling transfers...");

    for (Index i(0); (i < num_levels) && ((i+1) < extruded_dom_ctrl.size_virtual()); ++i)
    {
      system_levels.at(i)->assemble_coarse_muxers(extruded_dom_ctrl.at(i+1));
      system_levels.at(i)->assemble_transfers(extruded_dom_ctrl.at(i), extruded_dom_ctrl.at(i+1), cubature);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling basic matrices...");

    //// Create assembler for the Robin BCs for the pressure Poisson problem - only needed in the presence of
    //// do-nothing boundaries
    //std::vector<std::shared_ptr<Assembly::TraceAssembler<ExtrudedTrafoType>>> robin_asm_p(num_levels);
    //for(Index i(0); i < num_levels; ++i)
    //{
    //  robin_asm_p.at(i) =
    //    std::make_shared<Assembly::TraceAssembler<ExtrudedTrafoType>>(extruded_dom_ctrl.at(i)->trafo);

    //  auto* mpp = extruded_dom_ctrl.at(i)->get_mesh_node()->find_mesh_part("bnd:b");
    //  if(mpp != nullptr)
    //  {
    //    robin_asm_p.at(i)->add_mesh_part(*mpp);
    //  }

    //  mpp = extruded_dom_ctrl.at(i)->get_mesh_node()->find_mesh_part("bnd:t");
    //  if(mpp != nullptr)
    //  {
    //    robin_asm_p.at(i)->add_mesh_part(*mpp);
    //  }
    //}

    for(Index i(0); i < num_levels; ++i)
    {
      // assemble velocity matrix structure
      system_levels.at(i)->assemble_velo_struct(extruded_dom_ctrl.at(i)->space_velo);
      // assemble pressure matrix structure
      system_levels.at(i)->assemble_pres_struct(extruded_dom_ctrl.at(i)->space_pres);
      // assemble pressure laplace matrix
      system_levels.at(i)->matrix_s.local().format();
      Assembly::Common::LaplaceOperator laplace_op;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(system_levels.at(i)->matrix_s.local(),
      laplace_op, extruded_dom_ctrl.at(i)->space_pres, cubature);

      //// Add Robin boundary terms to the lhs
      //Assembly::Common::IdentityOperator mass_op;
      //robin_asm_p.at(i)->assemble_operator_matrix1(
      //  system_levels.at(i)->matrix_s.local(), mass_op, extruded_dom_ctrl.at(i)->space_pres, cubature, DT_(0e3));
    }

    // assemble B/D matrices on finest level
    system_levels.front()->assemble_grad_div_matrices(extruded_dom_ctrl.front()->space_velo, extruded_dom_ctrl.front()->space_pres, cubature);

    // Pressure mass matrix on the finest level
    typename SystemLevelType::GlobalSchurMatrix matrix_m_p = system_levels.front()->matrix_s.clone(LAFEM::CloneMode::Layout);
    matrix_m_p.local().format();
    Assembly::Common::IdentityOperator mass_p_op;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_m_p.local(), mass_p_op, extruded_dom_ctrl.front()->space_pres, cubature);
    typedef Global::Filter<LAFEM::NoneFilter<MemType, DataType, IndexType>, typename SystemLevelType::ScalarMirror> MassPFilter;

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

#ifndef ANALYTIC_SOLUTION
    Analytic::Common::XYPlaneRotation<DataType, ExtrudedMeshType::world_dim> rotation_inner(angular_velocity_inner, extruded_centre_inner);
    Analytic::Common::XYPlaneRotation<DataType, ExtrudedMeshType::world_dim> rotation_outer(angular_velocity_outer, extruded_centre_outer);
#else
    //WorldPoint zeros_y(0);
    //zeros_y(0) = -DataType(1);
    //zeros_y(1) = DataType(1);
    //DataType amplitude(1);
    //Analytic::Common::YZPlaneParabolic<DataType, MeshType::world_dim> rotation_outer(amplitude, zeros_y);

    //typedef Analytic::Common::SinYT0<DataType, MeshType::world_dim> AnalyticSolV;
    //typedef Analytic::Common::ConstantFunction<MeshType::world_dim, DataType> AnalyticSolP;
    //AnalyticSolV velo_sol;
    //AnalyticSolP pres_sol(DataType(0));
    //Analytic::Common::SinYT0StokesRhs<DataType, MeshType::world_dim> velo_rhs(reynolds);

    typedef Analytic::Common::GuermondStokesSol<DataType, ExtrudedMeshType::world_dim> AnalyticSolV;
    typedef Analytic::Common::GuermondStokesSolPressure<DataType, ExtrudedMeshType::world_dim> AnalyticSolP;
    AnalyticSolV velo_sol;
    AnalyticSolP pres_sol;
    Analytic::Common::GuermondStokesSolRhs<DataType, ExtrudedMeshType::world_dim> velo_rhs(reynolds);
#endif

    std::vector<Assembly::UnitFilterAssembler<ExtrudedMeshType>>
      unit_asm_velo_i(num_levels), unit_asm_velo_o(num_levels);

    std::vector<std::shared_ptr<Assembly::SlipFilterAssembler<ExtrudedMeshType>>>
      slip_asm_velo_b(num_levels), slip_asm_velo_t(num_levels);

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local system filters
      typename SystemLevelType::LocalVeloUnitFilter& unit_fil_loc_v = system_levels.at(i)->filter_velo.local().template at<1>();
      typename SystemLevelType::LocalVeloSlipFilter& slip_fil_loc_v = system_levels.at(i)->filter_velo.local().template at<0>();

      // Add inner boundary to assembler
      {
        auto* mesh_part_node = extruded_dom_ctrl.at(i)->get_mesh_node()->find_mesh_part_node("bnd:i");

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
        auto* mesh_part_node = extruded_dom_ctrl.at(i)->get_mesh_node()->find_mesh_part_node("bnd:o");

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

      slip_asm_velo_b.at(i) = std::make_shared<Assembly::SlipFilterAssembler<ExtrudedMeshType>>(extruded_dom_ctrl.at(i)->get_mesh());

      // Add bottom boundary to assembler
      {
        auto* mesh_part_node = extruded_dom_ctrl.at(i)->get_mesh_node()->find_mesh_part_node("bnd:b");

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        if(mesh_part_node != nullptr)
        {
          auto* mesh_part = mesh_part_node->get_mesh();

          if(mesh_part != nullptr)
          {
            slip_asm_velo_b.at(i)->add_mesh_part(*mesh_part);
          }
        }
      }

      slip_asm_velo_t.at(i) = std::make_shared<Assembly::SlipFilterAssembler<ExtrudedMeshType>>(extruded_dom_ctrl.at(i)->get_mesh());

      // Add bottom boundary to assembler
      {
        auto* mesh_part_node = extruded_dom_ctrl.at(i)->get_mesh_node()->find_mesh_part_node("bnd:t");

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        if(mesh_part_node != nullptr)
        {
          auto* mesh_part = mesh_part_node->get_mesh();

          if(mesh_part != nullptr)
          {
            slip_asm_velo_t.at(i)->add_mesh_part(*mesh_part);
          }
        }
      }

      // assemble the velocity filter
#ifndef ANALYTIC_SOLUTION
      unit_asm_velo_i[i].assemble(unit_fil_loc_v, extruded_dom_ctrl.at(i)->space_velo, rotation_inner);
      unit_asm_velo_o[i].assemble(unit_fil_loc_v, extruded_dom_ctrl.at(i)->space_velo, rotation_outer);
#else
      unit_asm_velo_o[i].assemble(unit_fil_loc_v, extruded_dom_ctrl.at(i)->space_velo, velo_sol);
#endif
      slip_asm_velo_b.at(i)->assemble(slip_fil_loc_v, extruded_dom_ctrl.at(i)->space_velo);
      slip_asm_velo_t.at(i)->assemble(slip_fil_loc_v, extruded_dom_ctrl.at(i)->space_velo);

      // assembly of the pressure filter is done in the system level
      system_levels.at(i)->assemble_global_filters(extruded_dom_ctrl.at(i)->space_pres, cubature);

    } // all levels

    watch_asm_mat.stop();
    /* ***************************************************************************************** */

    // get our vector types
    typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    ExtrudedDomainLevelType& the_domain_level = *extruded_dom_ctrl.front();
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

    watch_sol_init.start();
    Solver::MatrixStock<
      typename SystemLevelType::GlobalMatrixBlockA,
      typename SystemLevelType::GlobalVeloFilter,
      typename SystemLevelType::GlobalVeloTransfer> matrix_stock_velo(extruded_dom_ctrl.size_virtual());
    for (auto & system_level: system_levels)
    {
      matrix_stock_velo.systems.push_back(system_level->matrix_a.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_velo.gates_row.push_back(&system_level->gate_velo);
      matrix_stock_velo.gates_col.push_back(&system_level->gate_velo);
      matrix_stock_velo.filters.push_back(system_level->filter_velo.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_velo.muxers.push_back(&system_level->coarse_muxer_velo);
      matrix_stock_velo.transfers.push_back(system_level->transfer_velo.clone(LAFEM::CloneMode::Shallow));
    }

    auto tsolver_a = Solver::SolverFactory::create_scalar_solver(matrix_stock_velo, &solver_config, "linsolver_a");
    Solver::PreconditionedIterativeSolver<typename decltype(tsolver_a)::element_type::VectorType>* solver_a =
      (Solver::PreconditionedIterativeSolver<typename decltype(tsolver_a)::element_type::VectorType>*) &(*tsolver_a);
    matrix_stock_velo.hierarchy_init_symbolic();
    solver_a->init_symbolic();
    solver_a->set_plot_name(solver_a->name()+" (V-Bur)");

    /* ***************************************************************************************** */

    comm.print("Setting up Pressure Multigrid...");

    Solver::MatrixStock<
      typename SystemLevelType::GlobalSchurMatrix,
      typename SystemLevelType::GlobalPresFilter,
      typename SystemLevelType::GlobalPresTransfer> matrix_stock_pres(extruded_dom_ctrl.size_virtual());
    for (auto & system_level: system_levels)
    {
      matrix_stock_pres.systems.push_back(system_level->matrix_s.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_pres.gates_row.push_back(&system_level->gate_pres);
      matrix_stock_pres.gates_col.push_back(&system_level->gate_pres);
      matrix_stock_pres.filters.push_back(system_level->filter_pres.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_pres.muxers.push_back(&system_level->coarse_muxer_pres);
      matrix_stock_pres.transfers.push_back(system_level->transfer_pres.clone(LAFEM::CloneMode::Shallow));
    }

    auto tsolver_s = Solver::SolverFactory::create_scalar_solver(matrix_stock_pres, &solver_config, "linsolver_s");
    Solver::PreconditionedIterativeSolver<typename decltype(tsolver_s)::element_type::VectorType>* solver_s =
      (Solver::PreconditionedIterativeSolver<typename decltype(tsolver_s)::element_type::VectorType>*) &(*tsolver_s);
    matrix_stock_pres.hierarchy_init_symbolic();
    solver_s->init_symbolic();
    solver_s->set_plot_name(solver_s->name()+" (P-Lap)");

    Solver::MatrixStock<
      typename SystemLevelType::GlobalSchurMatrix,
      MassPFilter,
      typename SystemLevelType::GlobalPresTransfer> ms_mass_p(std::size_t(1));
    MassPFilter mass_p_filter;
    ms_mass_p.systems.push_back(matrix_m_p.clone(LAFEM::CloneMode::Shallow));
    ms_mass_p.gates_row.push_back(&the_system_level.gate_pres);
    ms_mass_p.gates_col.push_back(&the_system_level.gate_pres);
    ms_mass_p.filters.push_back(mass_p_filter.clone(LAFEM::CloneMode::Shallow));

    auto tsolver_m_p = Solver::SolverFactory::create_scalar_solver(ms_mass_p, &solver_config, "solver_m_p");
    Solver::PreconditionedIterativeSolver<typename decltype(tsolver_m_p)::element_type::VectorType>* solver_m_p =
      (Solver::PreconditionedIterativeSolver<typename decltype(tsolver_m_p)::element_type::VectorType>*) &(*tsolver_m_p);
    solver_m_p->init_symbolic();
    solver_m_p->set_plot_name(solver_s->name()+" (P-Mass)");
    watch_sol_init.stop();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("\n");
    for(Index i(0); i < num_levels; ++i)
    {
      Index ndof_v(system_levels.at(i)->matrix_a.columns());
      Index ndof_p(system_levels.at(i)->matrix_s.columns());
      String msg("Level "+stringify(extruded_dom_ctrl.at(i)->get_level_index())+" DoF: "+stringify(ndof_v)+ " / " + stringify(ndof_p));
      comm.print(msg);
    }

    // Create solution vectors
    std::vector<GlobalVeloVector> vec_sol_v(3);
    for(size_t i(0); i < vec_sol_v.size(); ++i)
    {
      vec_sol_v.at(i) = std::move(matrix_a.create_vector_r());
      vec_sol_v.at(i).format();
    }

    // This will hold D v because it gets used several times if the rotational version of the splitting is used
    GlobalPresVector vec_D_v = matrix_d.create_vector_l();
    vec_D_v.format();

    // create convection vector
    GlobalVeloVector vec_conv = matrix_a.create_vector_l();
    // The mesh velocity is 1/delta_t*(coords_new - coords_old) and computed in each time step
    auto mesh_velocity = meshopt_ctrl->get_coords().clone(LAFEM::CloneMode::Deep);
    mesh_velocity.format();

    GlobalVeloVector extruded_mesh_velocity(vec_sol_v.at(0).get_gate(), the_domain_level.get_mesh().get_num_entities(0));
    extruded_mesh_velocity.format();

    // The mesh velocity interpolated to the velocity space is needed for the convection vector if the mesh moves
    GlobalVeloVector mesh_velocity_q2 = matrix_a.create_vector_l();

    // We need at Max(2, p_extrapolation_steps+1) pressure vectors
    std::vector<GlobalPresVector> vec_sol_p(Math::max(time_disc.get_max_p_extrapolation_steps(),Index(1))+Index(1));
    for(size_t i(0); i < vec_sol_p.size(); ++i)
    {
      vec_sol_p.at(i) = std::move(matrix_s.create_vector_r());
      vec_sol_p.at(i).format();
    }

    GlobalPresVector vec_extrapolated_p = matrix_s.create_vector_l();
    vec_extrapolated_p.format();

    // We need several solution vectors for the auxillary Poisson problem
    std::vector<GlobalPresVector> vec_sol_phi(time_disc.get_max_num_steps()+Index(1));
    for(size_t i(0); i < vec_sol_phi.size(); ++i)
    {
      vec_sol_phi.at(i) = std::move(matrix_s.create_vector_r());
      vec_sol_phi.at(i).format();
    }

    // RHS vectors
    GlobalVeloVector vec_rhs_v = matrix_a.create_vector_l();
    GlobalPresVector vec_rhs_p = matrix_s.create_vector_l();
    GlobalPresVector vec_rhs_phi = matrix_s.create_vector_l();
    vec_rhs_v.format();
    vec_rhs_p.format();
    vec_rhs_phi.format();

    // Vectors for the analytical solution
    GlobalVeloVector vec_sol_v_analytic = matrix_a.create_vector_l();
    GlobalPresVector vec_sol_p_analytic = matrix_s.create_vector_l();

    // Vectors for time derivatives
    GlobalVeloVector vec_der_v = vec_sol_v.at(0).clone();
    GlobalPresVector vec_der_p = vec_sol_p.at(0).clone();

    // Initial time
    DataType time(0);
    // Counter for timesteps
    Index time_step(0);

    // set up a burgers assembler for the velocity matrix
    Assembly::BurgersAssembler<DataType, IndexType, ExtrudedMeshType::world_dim> burgers_lhs;
    burgers_lhs.deformation = use_deformation;
    burgers_lhs.theta = time_disc.coeff_lhs_v.at(0);
    burgers_lhs.beta = time_disc.coeff_lhs_v.at(1);
    burgers_lhs.nu = time_disc.coeff_lhs_v.at(2);

    // Set up a burgers assembler for the RHS vector that uses the last velocity vector (BDF1 and BDF2)
    Assembly::BurgersAssembler<DataType, IndexType, ExtrudedMeshType::world_dim> burgers_rhs_1;
    burgers_rhs_1.deformation = use_deformation;
    burgers_rhs_1.theta = time_disc.coeff_rhs_v_1.at(0);
    burgers_rhs_1.beta = time_disc.coeff_rhs_v_1.at(1);
    burgers_rhs_1.nu = time_disc.coeff_rhs_v_1.at(2);

    // Set up a burgers assembler for the RHS vector that uses the previous velocity vector (BDF2)
    Assembly::BurgersAssembler<DataType, IndexType, ExtrudedMeshType::world_dim> burgers_rhs_2;
    burgers_rhs_2.deformation = use_deformation;
    burgers_rhs_2.theta = time_disc.coeff_rhs_v_2;
    burgers_rhs_2.beta = DataType(0);
    burgers_rhs_2.nu = DataType(0);

    // This is the absolute turning angle of the screws
    DataType alpha(0);
    WorldPoint rotation_angles(0);

#ifdef ANALYTIC_SOLUTION
    // Hackish way to start with the analytic solution at something other than t=0
    DataType start_time(0);

    velo_sol.set_time(start_time);
    pres_sol.set_time(start_time);
    velo_rhs.set_time(start_time);

    Assembly::Interpolator::project(vec_sol_v.at(0).local(), velo_sol, the_domain_level.space_velo);
    // For starting with an exact pressure or phi
    //Assembly::Interpolator::project(vec_sol_p.at(0).local(), pres_sol, the_domain_level.space_pres);
    //Assembly::Interpolator::project(vec_sol_phi.at(0).local(), pres_sol, the_domain_level.space_pres);

    DataType err_v_L2_max(0);
    DataType err_v_H1_max(0);
    DataType err_p_L2_max(0);

    DataType err_v_l2_L2(0);
    DataType err_v_l2_H1(0);
    DataType err_p_l2_L2(0);
#endif

    // Write output again
    if(write_vtk)
    {
      watch_vtk.start();
      for(size_t l(0); l < dom_ctrl.size_physical(); ++l)
      {
        auto& dom_lvl = dom_ctrl.at(l);

        int lvl_index(dom_lvl->get_level_index());

        String vtk_name = String("mesh_post_initial_n"+stringify(comm.size())+"_lvl_"+stringify(lvl_index));
        comm.print("Writing "+vtk_name+".vtk");

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(dom_lvl->get_mesh());

        exporter.add_cell_scalar("Worst angle", edge_angle_cellwise.data());
        exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise.data());

        meshopt_ctrl->add_to_vtk_exporter(exporter, lvl_index);

        exporter.write(vtk_name, comm.rank(), comm.size());
      }

      if(solve_flow)
      {
        for(size_t l(0); l < extruded_dom_ctrl.size_physical(); ++l)
        {
          auto& dom_lvl = extruded_dom_ctrl.at(l);

          int lvl_index(dom_lvl->get_level_index());

          String vtk_name = String("flow_initial_n"+stringify(comm.size())+"_lvl_"+stringify(lvl_index));
          comm.print("Writing "+vtk_name+".vtk");

          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<ExtrudedMeshType> exporter(dom_lvl->get_mesh());

          //exporter.add_vertex_scalar("extrude_idx", extruded_dom_ctrl._vtx_map.data());
          exporter.add_vertex_vector("v", vec_sol_v.at(0).local());
          exporter.add_vertex_scalar("p", vec_sol_p.at(0).local().elements());
          //exporter.add_vertex_vector("nu", system_levels.at(l)->filter_velo.local().template at<0>().get_filter_vector());
#ifdef ANALYTIC_SOLUTION
          Assembly::AnalyticVertexProjector::project(vec_sol_v_analytic.local(), velo_sol, the_domain_level.trafo);
          exporter.add_vertex_vector("v_analytic", vec_sol_v_analytic.local());

          Assembly::AnalyticVertexProjector::project(vec_sol_v_analytic.local(), velo_rhs, the_domain_level.trafo);
          exporter.add_vertex_vector("v_rhs_analytic", vec_sol_v_analytic.local());

          Assembly::AnalyticVertexProjector::project(vec_sol_p_analytic.local(), pres_sol, the_domain_level.trafo);
          exporter.add_vertex_scalar("p_analytic", vec_sol_p_analytic.local().elements());
#endif

          exporter.write(vtk_name, comm.rank(), comm.size());
        }
      }

      watch_vtk.stop();
    }

    // Write xml
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
        String xml_name(file_basename+"_post_initial_"+stringify(lvl_index));

        comm.print("Writing "+xml_name+".xml");

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

    while(time < t_end)
    {
      ++time_step;
      time += time_disc.delta_t();

      //DataType alpha_old = alpha;
      alpha = angular_velocity_outer*time;
      //DataType delta_alpha = alpha - alpha_old;

      bool failure(false);

      comm.print("");
      comm.print(String(80u, '#'));
      comm.print("Timestep "+stringify(time_step)+": t = "+stringify_fp_fix(time)+", angle = "
          +stringify_fp_fix(alpha/(DataType(2)*pi)*DataType(360)) + " degrees\n");

      if(solve_mesh_optimisation)
      {
        watch_meshopt.start();
        // Save old vertex coordinates
        meshopt_ctrl->mesh_to_buffer();
        old_coords.copy(meshopt_ctrl->get_coords());
        new_coords.copy(meshopt_ctrl->get_coords());
        auto& coords_loc(new_coords.local());

        // Rotate the charts
        rotation_angles(0) = time_disc.delta_t()*angular_velocity_outer;
        rotation_angles(1) = DataType(0);

        if(outer_chart != nullptr)
        {
          outer_chart->transform(centre_outer, rotation_angles, centre_outer);
        }

        rotation_angles(0) = time_disc.delta_t()*angular_velocity_inner;
        if(inner_chart != nullptr)
        {
          inner_chart->transform(centre_inner, rotation_angles, centre_inner);
        }

        // Update boundary of the inner screw
        // This is the 2x2 matrix representing the turning by the angle delta_alpha of the inner screw
        if(inner_indices != nullptr)
        {
          Tiny::Matrix<DataType, 2, 2> rot(DataType(0));

          rot.set_rotation_2d(time_disc.delta_t()*angular_velocity_inner);

          WorldPoint tmp(DataType(0));
          WorldPoint tmp2(DataType(0));

          for(Index i(0); i < inner_indices->get_num_entities(); ++i)
          {
            // Index of boundary vertex i in the mesh
            Index j(inner_indices->operator[](i));
            // Translate the point to the centre of rotation
            tmp = coords_loc(j) - centre_inner;
            // Rotate
            tmp2.set_mat_vec_mult(rot, tmp);
            // Translate the point by the new centre of rotation
            coords_loc(j, centre_inner + tmp2);
          }
        }

        // The outer screw has 7 teeth as opposed to the inner screw with 6, and it rotates at 6/7 of the speed
        if(outer_indices != nullptr)
        {
          Tiny::Matrix<DataType, 2, 2> rot(DataType(0));
          rot.set_rotation_2d(time_disc.delta_t()*angular_velocity_outer);

          WorldPoint tmp(DataType(0));
          WorldPoint tmp2(DataType(0));

          // The outer screw rotates centrically, so centre_outer remains the same at all times
          for(Index i(0); i < outer_indices->get_num_entities(); ++i)
          {
            // Index of boundary vertex i in the mesh
            Index j(outer_indices->operator[](i));
            tmp = coords_loc(j) - centre_outer;

            tmp2.set_mat_vec_mult(rot, tmp);

            coords_loc(j, centre_outer+tmp2);
          }
        }

        watch_meshopt.stop();
        if(meshopt_preproc!=nullptr)
        {
          watch_meshopt_preproc.start();
          comm.print("Meshopt preprocessor:");
          FEAT::Statistics::expression_target = "meshopt_preproc";
          meshopt_preproc->prepare(new_coords);
          meshopt_preproc->optimise();
          comm.print("");
          new_coords.copy(meshopt_preproc->get_coords());
          watch_meshopt_preproc.stop();

          //if(write_vtk && ( (time_step%vtk_freq == 0) || failure))
          //{
          //  watch_vtk.start();

          //  // Compute mesh quality on the finest
          //  meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise.data(), qi_cellwise.data());

          //  String vtk_name = String("mesh_pre_n"+stringify(comm.size())+"_"+stringify(time_step));
          //  // Create a VTK exporter for our mesh
          //  Geometry::ExportVTK<MeshType> exporter(dom_ctrl.front()->get_mesh());

          //  exporter.add_cell_scalar("Worst angle", edge_angle_cellwise.data());
          //  exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise.data());

          //  meshopt_ctrl->add_to_vtk_exporter(exporter, dom_ctrl.front()->get_level_index());

          //  mesh_velocity.axpy(old_coords, new_coords, DataType(-1));
          //  mesh_velocity.scale(mesh_velocity, DataType(1)/time_disc.delta_t());

          //  exporter.add_vertex_vector("mesh velocity", mesh_velocity.local());

          //  exporter.write(vtk_name, comm);

          //  watch_vtk.stop();
          //}
        }
        watch_meshopt.start();

        comm.print("Mesh optimisation:");
        // Now prepare the functional
        FEAT::Statistics::expression_target = "meshopt";
        meshopt_ctrl->prepare(new_coords);
        meshopt_ctrl->optimise();
        new_coords.copy(meshopt_ctrl->get_coords());
        extruded_dom_ctrl.extrude_vertex_sets();
        comm.print("");
        watch_meshopt.stop();

        // Compute mesh velocity
        {
          mesh_velocity.axpy(old_coords, new_coords, DataType(-1));
          mesh_velocity.scale(mesh_velocity, DataType(1)/time_disc.delta_t());

          // Compute maximum of the mesh velocity
          DataType max_mesh_velocity(0);

          for(IT_ i(0); i < mesh_velocity.local().size(); ++i)
          {
            max_mesh_velocity = Math::max(max_mesh_velocity, (mesh_velocity.local())(i).norm_euclid());
          }

          if(mesh_velocity.get_gate() != nullptr)
          {
            max_mesh_velocity = mesh_velocity.get_gate()->max(max_mesh_velocity);
          }

          String msg = String("max. mesh velocity").pad_back(pad_width, ' ') + ": " + stringify_fp_sci(max_mesh_velocity)+"\n";
          comm.print(msg);

          extruded_dom_ctrl.extrude_vertex_vector(extruded_mesh_velocity.local(), mesh_velocity.local());
        }
      } // solve_mesh_optimisation

      // Now the flow problem
      if(solve_flow)
      {
        comm.print("Solving flow problem:");
#ifdef ANALYTIC_SOLUTION
        velo_sol.set_time(time+start_time);
        pres_sol.set_time(time+start_time);
        velo_rhs.set_time(time+start_time);
#endif

        if((time_step == time_disc.get_max_num_steps())
            && (time_disc.get_num_steps() != time_disc.get_max_num_steps()))
            {
              comm.print("Switching time discretisation to BDF"+stringify(time_disc.get_max_num_steps()));
              time_disc.finish_startup();

              // Update the coefficients in the Burgers matrix assembler
              burgers_lhs.theta = time_disc.coeff_lhs_v.at(0);
              burgers_lhs.beta = time_disc.coeff_lhs_v.at(1);
              burgers_lhs.nu = time_disc.coeff_lhs_v.at(2);

              burgers_rhs_1.theta = time_disc.coeff_rhs_v_1.at(0);
              burgers_rhs_1.beta = time_disc.coeff_rhs_v_1.at(1);
              burgers_rhs_1.nu = time_disc.coeff_rhs_v_1.at(2);

              burgers_rhs_2.theta = time_disc.coeff_rhs_v_2;
              burgers_rhs_2.beta = DataType(0);
              burgers_rhs_2.nu = DataType(0);

            }

        // Save previous velocity solutions
        for(Index step(time_disc.get_num_steps()); step > Index(0); --step)
        {
          vec_sol_v.at(step).copy(vec_sol_v.at(step-Index(1)));
        }

        // Save previous phi solutions
        for(Index step(time_disc.get_num_steps()); step > Index(0); --step)
        {
          vec_sol_phi.at(step).copy(vec_sol_phi.at(step-Index(1)));
        }

        // Save previous pressure solutions
        for(Index step(time_disc.get_p_extrapolation_steps()); step > Index(0); --step)
        {
          vec_sol_p.at(step).copy(vec_sol_p.at(step-Index(1)));
        }

        // Extrapolate previous pressure solutions
        vec_extrapolated_p.format();
        for(Index step(0); step < time_disc.get_p_extrapolation_steps(); ++step)
        {
          vec_extrapolated_p.axpy(
            vec_sol_p.at(step+Index(1)), vec_extrapolated_p, time_disc.coeff_extrapolation_p.at(step));
        }

        // Assemble filters on all levels
        watch_asm_fil.start();
        for(Index i(0); i < num_levels; ++i)
        {
          // get our local system filters
          typename SystemLevelType::LocalVeloUnitFilter& fil_loc_v = system_levels.at(i)->filter_velo.local().template at<1>();
          //typename SystemLevelType::LocalVeloFilter& fil_loc_v = system_levels.at(i)->filter_velo.local();

          // assemble the velocity filter
#ifndef ANALYTIC_SOLUTION
          unit_asm_velo_i.at(i).assemble(fil_loc_v, extruded_dom_ctrl.at(i)->space_velo, rotation_inner);
          unit_asm_velo_o.at(i).assemble(fil_loc_v, extruded_dom_ctrl.at(i)->space_velo, rotation_outer);
#else
          unit_asm_velo_o.at(i).assemble(fil_loc_v, extruded_dom_ctrl.at(i)->space_velo, velo_sol);
#endif

          // assembly of the pressure filter is done in the system level
          system_levels.at(i)->assemble_global_filters(extruded_dom_ctrl.at(i)->space_pres, cubature);

        } // all levels
        watch_asm_fil.stop();

        // apply filter onto solution vector
        filter_v.filter_sol(vec_sol_v.at(0));
        filter_p.filter_sol(vec_sol_p.at(0));

        // Does the old velocity have to be filtered with the new BCs? - I think not!
        //for (Index step(0); step< time_disc.get_num_steps(); ++step)
        //{
        //  filter_v.filter_sol(vec_sol_v.at(step+Index(1)));
        //}

        // assemble B/D matrices on finest level
        watch_asm_mat.start();
        system_levels.front()->assemble_grad_div_matrices(
          extruded_dom_ctrl.front()->space_velo, extruded_dom_ctrl.front()->space_pres, cubature);

        // Assemble pressure mass matrix on the finest level
        Assembly::BilinearOperatorAssembler::assemble_matrix1(
          matrix_m_p.local(), mass_p_op, extruded_dom_ctrl.front()->space_pres, cubature);

        // assemble pressure Laplace matrices on all levels
        for(std::size_t i(0); i < system_levels.size(); ++i)
        {
          auto& loc_mat_s = system_levels.at(i)->matrix_s.local();
          loc_mat_s.format();

          Assembly::Common::LaplaceOperator laplace_op;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(
            loc_mat_s, laplace_op, extruded_dom_ctrl.at(i)->space_pres, cubature);

          //// Add Robin boundary terms to the lhs
          //Assembly::Common::IdentityOperator mass_op;
          //robin_asm_p.at(i)->assemble_operator_matrix1(
          //  system_levels.at(i)->matrix_s.local(), mass_op, extruded_dom_ctrl.at(i)->space_pres, cubature, DT_(0e3));
        }
        watch_asm_mat.stop();

        // Assemble RHS vector for the velocity
        watch_asm_rhs.start();

        // Compute convection vector for rhs assembly

        // Interpolate the mesh velocity to Q2
        if(time_disc.ale_handling != Control::Time::TermHandling::off)
        {
          Assembly::FEInterpolator<SpaceVeloType, ExtrudedTrafoFESpace>::interpolate(
            mesh_velocity_q2.local(), extruded_mesh_velocity.local(),
            the_domain_level.space_velo, the_domain_level.space_pres);
        }

        vec_conv.format();
        if(time_disc.conv_handling == Control::Time::TermHandling::expl)
        {
          if((time_step > Index(2)))
          {
            // linear extrapolation of solution in time
            vec_conv.scale(vec_sol_v.at(1), DataType(2));
            vec_conv.axpy(vec_sol_v.at(2), vec_conv, -DataType(1));
          }
          else
          {
            // constant extrapolation of solution in time
            vec_conv.copy(vec_sol_v.at(1));
          }
        }

        // Use ALE velocity if necessary
        if(time_disc.ale_handling == Control::Time::TermHandling::expl)
        {
          vec_conv.axpy(mesh_velocity_q2, vec_conv, -DataType(1));
        }

        // Assemble rhs for NVS
        vec_rhs_v.format();

        // Add the terms depending on the last solution to the rhs
        burgers_rhs_1.assemble_vector(
          vec_rhs_v.local(), vec_conv.local(), vec_sol_v.at(Index(1)).local(),
          the_domain_level.space_velo, cubature);

        vec_conv.format();

        // Add the terms depending on the previous solution to the rhs
        if(time_disc.get_num_steps() == Index(2))
        {
          // Does the old velocity have to be filtered with the new BCs? - I think not!
          //filter_v.filter_sol(vec_sol_v.at(step+Index(1)));
          burgers_rhs_2.assemble_vector(
            vec_rhs_v.local(), vec_conv.local(), vec_sol_v.at(Index(2)).local(),
            the_domain_level.space_velo, cubature);
        }

#ifdef ANALYTIC_SOLUTION
        // Add analytic rhs as force functional
        Assembly::Common::ForceFunctional<decltype(velo_rhs)> velo_rhs_func(velo_rhs);
        Assembly::LinearFunctionalAssembler::assemble_vector(
          vec_rhs_v.local(), velo_rhs_func, extruded_dom_ctrl.front()->space_velo, cubature, time_disc.coeff_rhs_f);
#endif

        // Add phi gradient on rhs
        for(Index step(0); step < time_disc.get_num_steps(); ++step)
        {
          Assembly::GradOperatorAssembler::assemble(
            vec_rhs_v.local(), vec_sol_phi.at(step+Index(1)).local(),
            the_domain_level.space_velo, the_domain_level.space_pres,
            cubature, time_disc.coeff_rhs_v_phi.at(step));
        }

        // Synchronise the rhs vector after the local assemblers are finished
        vec_rhs_v.sync_0();

        //// Add phi gradient on rhs - superseeded by the GradOperatorAssembler above
        //for(Index step(0); step < time_disc.get_num_steps(); ++step)
        //{
        //  matrix_b.apply(vec_rhs_v, vec_sol_phi.at(step+Index(1)), vec_rhs_v, time_disc.coeff_rhs_v_phi.at(step));
        //}

        // Add pressure gradient on rhs if it got extrapolated at all. This does not have to be synchronised because
        // it is a multiplication with a Global::Matrix
        if(time_disc.get_p_extrapolation_steps() > Index(0))
        {
          matrix_b.apply(vec_rhs_v, vec_extrapolated_p, vec_rhs_v, time_disc.coeff_rhs_v_p);
        }

        // Apply RHS filter
        filter_v.filter_rhs(vec_rhs_v);

        watch_asm_rhs.stop();

        // Start of matrix assembly
        watch_asm_mat.start();
        vec_conv.format();
        // Compute convection vector for matrix assemly
        if(time_disc.conv_handling == Control::Time::TermHandling::impl)
        {
          if((time_step > Index(2)))
          {
            // linear extrapolation of solution in time
            vec_conv.scale(vec_sol_v.at(1), DataType(2));
            vec_conv.axpy(vec_sol_v.at(2), vec_conv, -DataType(1));
          }
          else
          {
            // constant extrapolation of solution in time
            vec_conv.copy(vec_sol_v.at(1));
          }
        }

        // Use ALE velocity if necessary
        if(time_disc.ale_handling == Control::Time::TermHandling::impl)
        {
          vec_conv.axpy(mesh_velocity_q2, vec_conv, -DataType(1));
        }

        // Phase 2: loop over all levels and assemble the burgers matrices

        // assemble burgers matrices on all levels
        for(std::size_t i(0); i < system_levels.size(); ++i)
        {
          auto& loc_mat_a = system_levels.at(i)->matrix_a.local();
          typename GlobalVeloVector::LocalVectorType vec_cv(vec_conv.local(), loc_mat_a.rows(), IndexType(0));
          loc_mat_a.format();
          burgers_lhs.assemble_matrix(loc_mat_a, vec_cv, extruded_dom_ctrl.at(i)->space_velo, cubature);
        }
        watch_asm_mat.stop();

        // Phase 3: initialise linear solvers
        watch_sol_init.start();
        matrix_stock_velo.hierarchy_init_numeric();
        solver_a->init_numeric();
        watch_sol_init.stop();

        // Solve velocity system
        FEAT::Statistics::expression_target = "solver_a";
        watch_solver_a.start();
        Solver::Status status_a = solver_a->correct(vec_sol_v.at(0), vec_rhs_v);
        watch_solver_a.stop();

        // release the solver
        solver_a->done_numeric();
        matrix_stock_velo.hierarchy_done_numeric();

        // check the solver's output
        if(!Solver::status_success(status_a))
        {
          comm.print("\n\nERROR: velocity solver broke down!\n");
          failure = true;
          ++failed_checks;
          break;
        }

        if(test)
        {
          if(status_a != Solver::Status::success)
          {
            comm.print("FAILED: V-Bur solver status was "+stringify(status_a));
            ++failed_checks;
          }
        }

        watch_asm_rhs.start();
        // Compute the weak divergence of the new velocity
        matrix_d.apply(vec_D_v, vec_sol_v.at(0));

        // Assemble rhs for the auxillary pressure poisson problem
        vec_rhs_phi.format();
        vec_rhs_phi.axpy(vec_D_v, vec_rhs_phi, time_disc.coeff_rhs_phi);
        filter_p.filter_rhs(vec_rhs_phi);
        watch_asm_rhs.stop();

        // initialise pressure poisson solver
        watch_sol_init.start();
        matrix_stock_pres.hierarchy_init_numeric();
        solver_s->init_numeric();
        watch_sol_init.stop();

        // Solve auxillary pressure Poisson
        FEAT::Statistics::expression_target = "solver_s";
        watch_solver_s.start();
        Solver::Status status_s = solver_s->correct(vec_sol_phi.at(0), vec_rhs_phi);
        watch_solver_s.stop();

        // release pressure poisson solver
        solver_s->done_numeric();
        matrix_stock_pres.hierarchy_done_numeric();

        // check solver output
        if(!Solver::status_success(status_s))
        {
          comm.print("\n\nERROR: pressure poisson solver broke down!\n");
          failure = true;
          ++failed_checks;
          break;
        }

        // Assemble rhs for projection step - no filtering here
        watch_asm_rhs.start();
        vec_rhs_p.format();
        // rhs_p = M (p + phi) + 1/reynolds D v
        matrix_m_p.apply(vec_rhs_p, vec_sol_phi.at(0), vec_rhs_p, DataType(1));
        // Add pressure on rhs if it got extrapolated at all
        if(time_disc.get_p_extrapolation_steps() > Index(0))
        {
          matrix_m_p.apply(vec_rhs_p, vec_extrapolated_p, vec_rhs_p, DataType(1));
        }
        // Add rotational part (or not)
        if(time_disc.is_rotational())
        {
          vec_rhs_p.axpy(vec_D_v, vec_rhs_p, time_disc.coeff_rhs_proj_D_v);
        }
        watch_asm_rhs.stop();

        // Update the pressure mass matrix and initialise the solver
        watch_sol_init.start();
        ms_mass_p.refresh();
        solver_m_p->init_numeric();
        watch_sol_init.stop();

        // Solve pressure projection problem: M_p p[k+1] = M_p p[k] + M_p phi[k+1] + 1/Re D u[k+1]
        FEAT::Statistics::expression_target = "solver_m_p";
        watch_solver_m_p.start();
        Solver::Status status_m_p = solver_m_p->correct(vec_sol_p.at(0), vec_rhs_p);
        watch_solver_m_p.stop();

        // release solver
        solver_m_p->done_numeric();

        // Is this necessary?
        filter_p.filter_sol(vec_sol_p.at(0));
        if(!Solver::status_success(status_m_p))
        {
          comm.print("\n\nERROR: pressure projection solver broke down!\n");
          failure = true;
          ++failed_checks;
          break;
        }

#ifdef ANALYTIC_SOLUTION
        {
          comm.print("");
          // compute local errors
          Assembly::VectorErrorInfo<DataType> error_v = Assembly::VectorErrorComputer<1>::compute(
            vec_sol_v.at(0).local(), velo_sol, the_domain_level.space_velo, cubature);

          Assembly::ScalarErrorInfo<DataType> error_p = Assembly::ScalarErrorComputer<0>::compute(
            vec_sol_p.at(0).local(), pres_sol, the_domain_level.space_pres, cubature);

          // synchronise all local errors
          error_v.synchronise(comm);
          error_p.synchronise(comm);

          err_v_L2_max = Math::max(err_v_L2_max, error_v.norm_h0);
          err_v_H1_max = Math::max(err_v_H1_max, error_v.norm_h1);
          err_p_L2_max = Math::max(err_p_L2_max, error_p.norm_h0);

          err_v_l2_L2 += Math::sqr(error_v.norm_h0)*time_disc.delta_t();
          err_v_l2_H1 += Math::sqr(error_v.norm_h1)*time_disc.delta_t();
          err_p_l2_L2 += Math::sqr(error_p.norm_h0)*time_disc.delta_t();

          // print errors
          if (comm.rank() == 0)
          {
            std::cout << "||u - u_ex||_L2: " << stringify_fp_sci(error_v.norm_h0) << ", time-l2: " <<
              stringify_fp_sci(Math::sqrt(err_v_l2_L2)) << std::endl;
            std::cout << " |u - u_ex|_H1 : " << stringify_fp_sci(error_v.norm_h1) << ", time-l2: " <<
              stringify_fp_sci(Math::sqrt(err_v_l2_H1)) << std::endl;
            std::cout << "||p - p_ex||_L2: " << stringify_fp_sci(error_p.norm_h0) << ", time-l2: " <<
              stringify_fp_sci(Math::sqrt(err_p_l2_L2)) << std::endl;
          }
        }
#endif
      } // solve_flow

      // Write output on the finest level only if scheduled or something went wrong
      if(write_vtk && ( (time_step%vtk_freq == 0) || failure))
      {

        watch_vtk.start();
        // 2d export
        //if(solve_mesh_optimisation)
        //{
        //  String vtk_name = String("mesh_post_n"+stringify(comm.size())+"_"+stringify(time_step));
        //  comm.print("Writing "+vtk_name+".vtk");
        //  auto& dom_lvl = dom_ctrl.front();

        //  // Compute mesh quality on the finest
        //  meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise.data(), qi_cellwise.data());

        //  // Create a VTK exporter for our mesh
        //  Geometry::ExportVTK<MeshType> exporter(dom_lvl->get_mesh());

        //  exporter.add_cell_scalar("Worst angle", edge_angle_cellwise.data());
        //  exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise.data());

        //  meshopt_ctrl->add_to_vtk_exporter(exporter, (*dom_lvl).get_level_index());

        //  exporter.add_vertex_vector("mesh velocity", mesh_velocity.local());

        //  exporter.write(vtk_name, comm.rank(), comm.size());
        //}

        // 3d export
        if(solve_flow)
        {
          auto& dom_lvl = extruded_dom_ctrl.front();

          String vtk_name = String("flow_n"+stringify(comm.size())+"_"+stringify(time_step).pad_front(6,'0'));
          comm.print("Writing "+vtk_name+".vtk");

          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<ExtrudedMeshType> exporter(dom_lvl->get_mesh());

          exporter.add_vertex_vector("mesh_velocity_q2", mesh_velocity_q2.local());

          exporter.add_vertex_vector("v", vec_sol_v.at(0).local());
          //for(Index step(0); step < time_disc.get_num_steps(); ++step)
          //{
          //  exporter.add_vertex_vector("v_"+stringify(step+1), vec_sol_v.at(step+1).local());
          //}

          //exporter.add_vertex_vector("vec_conv", vec_conv.local());
          //exporter.add_vertex_scalar("div_v", vec_D_v.local().elements());

          exporter.add_vertex_scalar("p", vec_sol_p.at(0).local().elements());
          //exporter.add_vertex_scalar("phi", vec_sol_phi.at(0).local().elements());

          //// compute and write time-derivatives
          //vec_der_v.axpy(vec_sol_v.at(1), vec_sol_v.at(0), -DataType(1));
          //vec_der_p.axpy(vec_sol_p.at(1), vec_sol_p.at(0), -DataType(1));

          //vec_der_v.scale(vec_der_v, DataType(1) / time_disc.delta_t());
          //vec_der_p.scale(vec_der_p, DataType(1) / time_disc.delta_t());

          //exporter.add_vertex_vector("v_dt", vec_der_v.local());
          //exporter.add_vertex_scalar("p_dt", vec_der_p.local().elements());

          //exporter.add_vertex_vector("v_rhs", vec_rhs_v.local());

#ifdef ANALYTIC_SOLUTION
          Assembly::AnalyticVertexProjector::project(vec_sol_v_analytic.local(), velo_rhs, the_domain_level.trafo);
          exporter.add_vertex_vector("v_rhs_analytic", vec_sol_v_analytic.local());

          Assembly::AnalyticVertexProjector::project(vec_sol_v_analytic.local(), velo_sol, the_domain_level.trafo);
          exporter.add_vertex_vector("v_analytic", vec_sol_v_analytic.local());

          Assembly::AnalyticVertexProjector::project(vec_sol_p_analytic.local(), pres_sol, the_domain_level.trafo);
          exporter.add_vertex_scalar("p_analytic", vec_sol_p_analytic.local().elements());
#endif

          exporter.write(vtk_name, comm.rank(), comm.size());
        }
        watch_vtk.stop();
        comm.print("\n");
      }

      // Write xml
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
          String xml_name("mesh_lvl_"+stringify(lvl_index)+"_"+stringify(time_step)+".xml");

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
        comm.print("");
      }

      // Compute and print quality indicators on the finest level only
      {
        watch_quality.start();
        DT_ lambda_min(Math::huge<DT_>());
        DT_ lambda_max(0);
        DT_ vol(0);
        DT_ vol_min(Math::huge<DT_>());
        DT_ vol_max(0);

        cell_size_defect = meshopt_ctrl->compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max, vol);

        meshopt_ctrl->compute_mesh_quality(edge_angle, qi_min, qi_mean, edge_angle_cellwise.data(), qi_cellwise.data());

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
        watch_quality.stop();
      }

      if(edge_angle < DT_(1))
      {
        comm.print("FAILED: Mesh deteriorated, stopping");
        ++failed_checks;
        failure = true;
      }

      // compress all statistics from the current timestep for further analysis after the solution is finished
      FEAT::Statistics::compress_solver_expressions();

      // write timings
      if(comm.rank() == 0)
      {
        double t_total = watch_total.elapsed();

        double t_asm_mat = watch_asm_mat.elapsed();
        double t_asm_rhs = watch_asm_rhs.elapsed();
        double t_asm_fil = watch_asm_fil.elapsed();
        double t_meshopt = watch_meshopt.elapsed();
        double t_meshopt_preproc = watch_meshopt_preproc.elapsed();
        double t_mesh_quality = watch_quality.elapsed();
        double t_sol_init = watch_sol_init.elapsed();
        double t_solver_a = watch_solver_a.elapsed();
        double t_solver_m_p = watch_solver_m_p.elapsed();
        double t_solver_s = watch_solver_s.elapsed();
        double t_vtk = watch_vtk.elapsed();

        double t_sum =
          t_asm_mat + t_asm_rhs + t_asm_fil +
          t_meshopt + t_meshopt_preproc + t_mesh_quality +
          t_sol_init + t_solver_a + t_solver_s + t_solver_m_p + t_vtk;

        dump_time(comm, "Total solver time", t_total, t_total);
        dump_time(comm, "Matrix assembly time", t_asm_mat, t_total);
        dump_time(comm, "Vector assembly time", t_asm_rhs, t_total);
        dump_time(comm, "Filter assembly time", t_asm_fil, t_total);
        dump_time(comm, "Solver init time", t_sol_init, t_total);
        dump_time(comm, "Solver-A time", t_solver_a, t_total);
        dump_time(comm, "Solver-S time", t_solver_s, t_total);
        dump_time(comm, "Solver-M_p time", t_solver_m_p, t_total);
        dump_time(comm, "Mesh optimisation time", t_meshopt, t_total);
        dump_time(comm, "Mesh opt preproc time", t_meshopt_preproc, t_total);
        dump_time(comm, "Mesh quality computation time", t_mesh_quality, t_total);
        dump_time(comm, "VTK write time", t_vtk, t_total);
        dump_time(comm, "Other time", t_total-t_sum, t_total);
      }

      if(failure)
      {
        break;
      }

      // continue with next time-step

    } // time loop

    // Check for the hard coded settings for test mode
    if(test)
    {
      if( edge_angle < DT_(27.8))
      {
        comm.print("FAILED: Final worst edge angle should be >= "+stringify_fp_fix(27.8)
            + " but is "+stringify_fp_fix(edge_angle));
        ++failed_checks;
      }
      if( qi_min < DT_(3e-2))
      {
        comm.print("FAILED: Final minimal quality indicator should be >= "+stringify_fp_fix(3.07e-2)
            + " but is "+stringify_fp_fix(qi_min));
        ++failed_checks;
      }
      if( cell_size_defect > DT_(6.68e-1))
      {
        comm.print("FAILED: Final cell size defect should be <= "+stringify_fp_fix(6.68e-1)
            + " but is "+stringify_fp_fix(cell_size_defect));
        ++failed_checks;
      }
    }

#ifdef ANALYTIC_SOLUTION
    err_v_l2_L2 = Math::sqrt(err_v_l2_L2);
    err_v_l2_H1 = Math::sqrt(err_v_l2_H1);
    err_p_l2_L2 = Math::sqrt(err_p_l2_L2);

    // print errors
    if (comm.rank() == 0)
    {
      std::cout << "Velocity error: " << stringify_fp_sci(err_v_L2_max) << " (l_infty-L2) " <<
        stringify_fp_sci(err_v_l2_L2) << " (l2-L2)" << std::endl;
      std::cout << "Velocity error: " << stringify_fp_sci(err_v_H1_max) << " (l_infty-H1) " <<
        stringify_fp_sci(err_v_l2_H1) << " (l2-H1)"<< std::endl;
      std::cout << "Pressure error: " << stringify_fp_sci(err_p_L2_max) << " (l_infty-L2) " <<
        stringify_fp_sci(err_p_l2_L2) << " (l2-L2)" << std::endl;
    }
#endif

    solver_a->done_symbolic();
    matrix_stock_velo.hierarchy_done_symbolic();

    solver_s->done_symbolic();
    matrix_stock_pres.hierarchy_done_symbolic();

    solver_m_p->done_symbolic();
    ms_mass_p.hierarchy_done_symbolic();

    //// Write final vtk output
    //if(write_vtk)
    //{
    //  watch_vtk.start();
    //  for(size_t l(0); l < dom_ctrl.size_physical(); ++l)
    //  {
    //    auto& dom_lvl = dom_ctrl.at(l);
    //    int lvl_index(dom_lvl->get_level_index());

    //    String vtk_name = String(file_basename+"_final_lvl_"+stringify(lvl_index));
    //    comm.print("Writing "+vtk_name+".vtk");

    //    // Create a VTK exporter for our mesh
    //    Geometry::ExportVTK<MeshType> exporter(dom_lvl->get_mesh());

    //    exporter.add_cell_scalar("Worst angle", edge_angle_cellwise);
    //    exporter.add_cell_scalar("Shape quality heuristic", qi_cellwise);

    //    meshopt_ctrl->add_to_vtk_exporter(exporter, lvl_index);
    //    exporter.add_vertex_vector("v", vec_sol_v.at(0).local());
    //    exporter.add_vertex_scalar("p", vec_sol_p.at(0).local().elements());

    //    //Assembly::AnalyticVertexProjector::project(vec_sol_v_analytic.local(), velo_sol, the_domain_level.trafo);
    //    //exporter.add_vertex_vector("v_analytic", vec_sol_v_analytic.local());

    //    //Assembly::AnalyticVertexProjector::project(vec_sol_p_analytic.local(), pres_sol, the_domain_level.trafo);
    //    //exporter.add_vertex_scalar("p_analytic", vec_sol_p_analytic.local().elements());

    //    exporter.write(vtk_name, comm.rank(), comm.size());
    //  }
    //  watch_vtk.stop();
    //}

    // Print success or not
    if(failed_checks == 0)
    {
      comm.print("\nFinished successfully!");
    }
    else
    {
      String msg("\nFAILED: "+stringify(failed_checks) + " check");
      if(failed_checks > 1)
      {
        msg+="s";
      }
      comm.print(msg);
    }

    watch_total.stop();

    double t_total = watch_total.elapsed();
    double t_asm_mat = watch_asm_mat.elapsed();
    double t_asm_rhs = watch_asm_rhs.elapsed();
    double t_asm_fil = watch_asm_fil.elapsed();
    double t_meshopt = watch_meshopt.elapsed();
    double t_meshopt_preproc = watch_meshopt_preproc.elapsed();
    double t_mesh_quality = watch_quality.elapsed();
    double t_sol_init = watch_sol_init.elapsed();
    double t_solver_a = watch_solver_a.elapsed();
    double t_solver_m_p = watch_solver_m_p.elapsed();
    double t_solver_s = watch_solver_s.elapsed();
    double t_vtk = watch_vtk.elapsed();

    double t_sum =
      t_asm_mat + t_asm_rhs + t_asm_fil +
      t_meshopt + t_meshopt_preproc + t_mesh_quality +
      t_sol_init + t_solver_a + t_solver_s + t_solver_m_p + t_vtk;

    // write timings
    if(comm.rank() == 0)
    {
      comm.print("");
      comm.print(FEAT::Statistics::get_formatted_solver_tree("solver_a").trim());
      comm.print(FEAT::Statistics::get_formatted_solver_tree("solver_s").trim());
      comm.print("");
      dump_time(comm, "Total solver time", t_total, t_total);
      dump_time(comm, "Matrix assembly time", t_asm_mat, t_total);
      dump_time(comm, "Vector assembly time", t_asm_rhs, t_total);
      dump_time(comm, "Solver init time", t_sol_init, t_total);
      dump_time(comm, "Solver-A time", t_solver_a, t_total);
      dump_time(comm, "Solver-S time", t_solver_s, t_total);
      dump_time(comm, "Solver-M_p time", t_solver_m_p, t_total);
      dump_time(comm, "Mesh optimisation time", t_meshopt, t_total);
      dump_time(comm, "Mesh quality computation time", t_mesh_quality, t_total);
      dump_time(comm, "VTK write time", t_vtk, t_total);
      dump_time(comm, "Other time", t_total-t_sum, t_total);
    }

    TimeStamp stamp_end;
    comm.print("Elapsed time: "+ stringify(stamp_end.elapsed(stamp_start)));

    return failed_checks;

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
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;
  //typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
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
  String mesh_path;
  std::deque<String> mesh_files;
  // String containing the mesh type, read from the header of the mesh file
  String mesh_type("");
  // Is the application running as a test? Read from the command line arguments
  bool test(false);

  // Streams for synchronising information read from files
  std::stringstream synchstream_app_config;
  std::stringstream synchstream_meshopt_config;
  std::stringstream synchstream_solver_config;

  // Create a parser for command line arguments.
  SimpleArgParser args(argc, argv);
  args.support("application_config");
  args.support("help");
  args.support("mesh-path");
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

  if( args.check("test") >= 0 )
  {
    test = true;
  }

  // query mesh directory path
  args.parse("mesh-path", mesh_path);
  if(!mesh_path.empty() && !mesh_path.ends_with('/'))
    mesh_path.push_back('/');

  // create a mesh file reader
  Geometry::MeshFileReader mesh_file_reader;

  // Application settings, has to be created here because it gets filled differently according to test
  PropertyMap application_config;
  PropertyMap meshopt_config;
  PropertyMap solver_config;

  // If we are not in test mode, parse command line arguments, read files, synchronise streams
  if(!test)
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

    // Parse the application config from the (synchronised) stream
    application_config.read(synchstream_app_config, true);

    // Get the application settings section
    auto app_settings_section = application_config.query_section("ApplicationSettings");
    XASSERTM(app_settings_section != nullptr,
    "Application config is missing the mandatory ApplicationSettings section!");

    auto mesh_files_p = app_settings_section->query("mesh_files");
    mesh_files = mesh_files_p.first.split_by_whitespaces();

    // Read configuration for mesh optimisation to stream
    auto meshopt_config_filename_p = app_settings_section->query("meshopt_config_file");

    XASSERTM(meshopt_config_filename_p.second,
    "ApplicationConfig section is missing the mandatory meshopt_config_file entry!");

    comm.print("Reading mesh optimisation config from file "+meshopt_config_filename_p.first);
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
    read_test_application_config(synchstream_app_config);
    application_config.read(synchstream_app_config, true);

    read_test_meshopt_config(synchstream_meshopt_config);
    meshopt_config.read(synchstream_meshopt_config, true);

    read_test_solver_config(synchstream_solver_config);
    solver_config.read(synchstream_solver_config, true);

    read_test_mesh_file_names(mesh_files);
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
    ret = NavierStokesScrewsApp<MemType, DataType, IndexType, H2M2D>::run(
      args, comm, application_config, meshopt_config, solver_config, mesh_file_reader);
  }
  //else if(mesh_type == "conformal:simplex:2:2")
  //{
  //  ret = NavierStokesScrewsApp<MemType, DataType, IndexType, S2M2D>::run(
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
  FEAT::Runtime::initialise(argc, argv);
  int ret = run_app(argc, argv);
  FEAT::Runtime::finalise();
  return ret;
}

static void read_test_application_config(std::stringstream& iss)
{
  iss << "[ApplicationSettings]" << std::endl;
  iss << "mesh_optimiser = HyperelasticityDefault" << std::endl;
  iss << "t_end = 1e-4" << std::endl;
  iss << "solve_flow = 1" << std::endl;
  iss << "solve_mesh_optimisation = 1" << std::endl;
  iss << "reynolds = 1e1" << std::endl;
  iss << "use_deformation = 0" << std::endl;

  iss << "[DomainControlSettings]" << std::endl;
  //iss << "parti-type = fallback parmetis" << std::endl;
  //iss << "parti-rank-elems = 4" << std::endl;
  iss << "lvl_min = 0" << std::endl;
  iss << "lvl_max = 1" << std::endl;
  iss << "z_min = 0.0" << std::endl;
  iss << "z_max = 1.0" << std::endl;
  iss << "slices = 1" << std::endl;

  iss << "[TimeDiscretisation]" << std::endl;
  iss << "delta_t = 1e-4" << std::endl;
  iss << "num_steps = 2" << std::endl;
  iss << "p_extrapolation_steps = 1" << std::endl;
  iss << "use_rotational_form = 1" << std::endl;
  iss << "ALE = impl" << std::endl;
  iss << "convection = impl" << std::endl;
  iss << "viscous = impl" << std::endl;
}

static void read_test_meshopt_config(std::stringstream& iss)
{
  iss << "[HyperElasticityDefault]" << std::endl;
  iss << "type = Hyperelasticity" << std::endl;
  iss << "config_section = HyperelasticityDefaultParameters" << std::endl;
  iss << "slip_boundaries = bnd:i bnd:o" << std::endl;
  iss << "meshopt_lvl = 0" << std::endl;

  iss << "[DuDvPreproc]" << std::endl;
  iss << "type = DuDv" << std::endl;
  iss << "config_section = DuDvDefaultParameters" << std::endl;
  iss << "dirichlet_boundaries = bnd:i bnd:o" << std::endl;
  iss << "meshopt_lvl = -1" << std::endl;

  iss << "[HyperelasticityDefaultParameters]" << std::endl;
  iss << "global_functional = HyperelasticityFunctional" << std::endl;
  iss << "cell_functional = RumpfFunctional" << std::endl;
  iss << "solver_config = NLCG" << std::endl;
  iss << "fac_norm = 1e0" << std::endl;
  iss << "fac_det = 1.0" << std::endl;
  iss << "fac_cof = 0.0" << std::endl;
  iss << "fac_reg = 1e-8" << std::endl;
  iss << "exponent_det = 2" << std::endl;
  iss << "scale_computation = iter_concentration" << std::endl;
  iss << "conc_function = GapWidth" << std::endl;

  iss << "[GapWidth]" << std::endl;
  iss << "type = ChartDistance" << std::endl;
  iss << "operation = add" << std::endl;
  iss << "chart_list = screw:i screw:o" << std::endl;
  iss << "function_type = default" << std::endl;
}

static void read_test_solver_config(std::stringstream& iss)
{
  iss << "[linsolver_a]" << std::endl;
  iss << "type = fgmres" << std::endl;
  iss << "max_iter = 1000" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "precon = mgv_a" << std::endl;
  iss << "precon_variant = left" << std::endl;
  iss << "krylov_dim = 7" << std::endl;
  iss << "plot_mode = all" << std::endl;

  iss << "[mgv_a]" << std::endl;
  iss << "type = mg" << std::endl;
  iss << "hierarchy = h1_a" << std::endl;
  iss << "lvl_min = -1" << std::endl;
  iss << "lvl_max = 0" << std::endl;
  iss << "cycle = f" << std::endl;

  iss << "[h1_a]" << std::endl;
  iss << "type = hierarchy" << std::endl;
  iss << "smoother = smoother_a" << std::endl;
  iss << "coarse = scale" << std::endl;

  iss << "[scale]" << std::endl;
  iss << "type = scale" << std::endl;
  iss << "omega = 0.7" << std::endl;

  iss << "[smoother_a]" << std::endl;
  iss << "type = fgmres" << std::endl;
  iss << "min_iter = 16" << std::endl;
  iss << "max_iter = 16" << std::endl;
  iss << "krylov_dim = 16" << std::endl;
  iss << "precon = jac" << std::endl;

  iss << "[linsolver_s]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "max_iter = 1000" << std::endl;
  iss << "#tol_abs = 1e-4" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "precon = mgv_s" << std::endl;
  iss << "min_stag_iter = 3" << std::endl;
  iss << "plot_mode = summary" << std::endl;

  iss << "[mgv_s]" << std::endl;
  iss << "type = mg" << std::endl;
  iss << "hierarchy = h1_s" << std::endl;
  iss << "lvl_min = 0" << std::endl;
  iss << "lvl_max = -1" << std::endl;
  iss << "cycle = v" << std::endl;

  iss << "[h1_s]" << std::endl;
  iss << "type = hierarchy" << std::endl;
  iss << "smoother = cg" << std::endl;
  iss << "coarse = Coarse-S" << std::endl;

  iss << "[Coarse-S]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "plot_mode = none" << std::endl;
  iss << "max_iter = 1000" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "precon = jac" << std::endl;

  iss << "[solver_m_p]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "max_iter = 100" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "tol_abs = 1e-4" << std::endl;
  iss << "precon = jac" << std::endl;
  iss << "plot_mode = summary" << std::endl;

  iss << "[NLCG]" << std::endl;
  iss << "type = NLCG" << std::endl;
  iss << "precon = none" << std::endl;
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

  iss << "[Meshopt-PCG]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "max_iter = 500" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "plot_mode = summary" << std::endl;
  iss << "precon = Meshopt-MG" << std::endl;

  iss << "[Meshopt-MG]" << std::endl;
  iss << "type = mg" << std::endl;
  iss << "hierarchy = h_cg" << std::endl;
  iss << "lvl_min = -1" << std::endl;
  iss << "lvl_max = 0" << std::endl;
  iss << "cycle = v" << std::endl;

  iss << "[h_cg]" << std::endl;
  iss << "type = hierarchy" << std::endl;
  iss << "smoother = cg" << std::endl;
  iss << "coarse = PCG-Jacobi" << std::endl;

  iss << "[cg]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "min_iter = 6" << std::endl;
  iss << "max_iter = 6" << std::endl;

  iss << "[PCG-Jacobi]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "plot_mode = none" << std::endl;
  iss << "max_iter = 1000" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "precon = jac" << std::endl;

  iss << "[jac]" << std::endl;
  iss << "type = jacobi" << std::endl;
  iss << "omega = 0.5" << std::endl;

  iss << "[PCG-JAC]" << std::endl;
  iss << "type = pcg" << std::endl;
  iss << "max_iter = 1000" << std::endl;
  iss << "tol_rel = 1e-8" << std::endl;
  iss << "precon = jac" << std::endl;

}

static void read_test_mesh_file_names(std::deque<String>& mesh_files)
{
  mesh_files.push_back("screws_2d_chart_bezier_24_28.xml");
  mesh_files.push_back("screws_2d_mesh_quad_opt_ss_360_1.xml");
}

static void display_help(const Dist::Comm& comm)
{
  if(comm.rank() == 0)
  {
    std::cout << "Navier Stokes Screws Application" << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << " --application_config: Path to the application configuration file" << std::endl;
    std::cout << "Optional arguments:" << std::endl;
    std::cout << " --test: Run as a test. Ignores configuration files and uses hard coded settings." << std::endl;
    std::cout << " --vtk <FREQ>: If this is set, vtk files are written every <FREQ> time steps." << std::endl;
    std::cout << " --help: Displays this text" << std::endl;
  }
}
