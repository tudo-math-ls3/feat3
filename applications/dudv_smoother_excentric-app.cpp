#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>
#include <kernel/meshopt/dudv_smoother.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/mesh_streamer.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/tiny_algebra.hpp>

using namespace FEAST;

  template<typename PointType, typename DataType>
  void centre_point_outer(PointType& my_point, DataType time)
  {
    my_point.v[0] = DataType(0.5) - DataType(0.125)*Math::cos(DataType(2)*Math::pi<DataType>()*DataType(time));
    my_point.v[1] = DataType(0.5) - DataType(0.125)*Math::sin(DataType(2)*Math::pi<DataType>()*DataType(time));
  }

  template<typename PointType, typename DataType>
  void centre_point_inner(PointType& my_point, DataType time)
  {
    my_point.v[0] = DataType(0.5) - DataType(0.1875)*Math::cos(DataType(2)*Math::pi<DataType>()*DataType(time));
    my_point.v[1] = DataType(0.5) - DataType(0.1875)*Math::sin(DataType(2)*Math::pi<DataType>()*DataType(time));
  }

/**
 * Application that shows the limitation of linear variational mesh optimisers, like the Du:Dv smoother
 */
/**
 * \brief Wrapper struct as functions do not seem to agree with template template parameters
 */
template
<
  typename Mem_,
  typename DT_,
  typename IT_,
  typename MeshType_
>
struct DuDvSmootherExcentricApp
{
  /// Precision for meshes etc, everything else uses the same data type
  typedef DT_ DataType;
  /// Memory architecture for solvers
  typedef Mem_ MemType;
  /// So we use Index
  typedef IT_ IndexType;
  /// The type of the mesh
  typedef MeshType_ MeshType;
  /// Shape of the mesh cells
  typedef typename MeshType::ShapeType ShapeType;
  /// Shape of mesh facets
  typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension - 1>::ShapeType FacetShapeType;
  /// Type of a surface mesh of facets
  typedef typename Geometry::ConformalMesh
  <FacetShapeType, MeshType::world_dim, MeshType::world_dim, typename MeshType::CoordType> SurfaceMeshType;

  /// The corresponding transformation
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  /// The smoother
  typedef Meshopt::DuDvSmoother<Mem_, DT_, IT_, TrafoType> SmootherType;
  /// Type for points in the mesh
  typedef Tiny::Vector<DataType, MeshType::world_dim> ImgPointType;

  /**
   * \brief Routine that does the actual work
   *
   * \param[in] my_streamer
   * MeshStreamer that contains the data from the mesh file.
   *
   * \param[in] level
   * Number of refines.
   */
  static int run(MeshStreamer& my_streamer, Index lvl_max, DT_ deltat)
  {
    // Filename for writing .vtu output
    String filename("");

    std::deque<String> dirichlet_list;
    dirichlet_list.push_back("inner");
    dirichlet_list.push_back("outer");
    std::deque<String> slip_list;

    // Read mesh from the MeshStreamer and create the MeshAtlas
    std::cout << "Creating mesh atlas..." << std::endl;
    Geometry::MeshAtlas<MeshType_>* atlas = nullptr;
    try
    {
      atlas = new Geometry::MeshAtlas<MeshType_>(my_streamer);
    }
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return 1;
    }

    // Create mesh node
    std::cout << "Creating mesh node..." << std::endl;
    Geometry::RootMeshNode<MeshType_>* rmn = nullptr;
    try
    {
      rmn = new Geometry::RootMeshNode<MeshType_>(my_streamer, atlas);
      rmn ->adapt();
    }
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return 1;
    }

    // refine
    for(Index lvl(1); lvl <= lvl_max; ++lvl)
    {
      std::cout << "Refining up to level " << lvl << "..." << std::endl;

      SmootherType refinement_smoother(rmn, dirichlet_list, slip_list);
      refinement_smoother.init();
      refinement_smoother.prepare();

      // Write initial state to file
      filename = "refinement_pre_"+stringify(lvl);
      Geometry::ExportVTK<MeshType> writer_refinement_pre(*(rmn->get_mesh()));
      writer_refinement_pre.write(filename);

      refinement_smoother.optimise();

      // Compute initial functional gradient
      filename = "refinement_post_"+stringify(lvl);
      Geometry::ExportVTK<MeshType> writer_refinement_post(*(rmn->get_mesh()));
      writer_refinement_post.write(filename);

      auto* old = rmn;
      rmn = old->refine();
      delete old;

    }

    auto* mesh = rmn->get_mesh();

    // This is the centre reference point
    ImgPointType x_0(DataType(0));

    // This is the centre point of the rotation of the inner screw
    ImgPointType x_inner(DataType(0));
    DataType excentricity_inner(DataType(0.2833));
    x_inner.v[0] = -excentricity_inner;
    // The indices for the inner screw
    auto& inner_indices = rmn->find_mesh_part("inner")->template get_target_set<0>();

    // This is the centre point of the rotation of the outer screw
    ImgPointType x_outer(DataType(0));
    // The indices for the outer screw
    auto& outer_indices = rmn->find_mesh_part("outer")->template get_target_set<0>();

    // The smoother in all its template glory
    SmootherType mr_dudv(rmn, dirichlet_list, slip_list);
    mr_dudv.init();
    mr_dudv.prepare();

    // Write initial state to file
    Geometry::ExportVTK<MeshType> writer_initial_pre(*mesh);
    writer_initial_pre.write("pre_initial");

    // Smooth the mesh
    mr_dudv.optimise();

    // Write optimised initial mesh
    Geometry::ExportVTK<MeshType> writer_initial_post(*mesh);
    writer_initial_post.write("post_initial");

    // For saving the old coordinates
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      coords_old(mesh->get_num_entities(0),DataType(0));
    // For writing out the pre optimisation mesh
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      coords_intermediate(mesh->get_num_entities(0),DataType(0));
    // For computing the mesh velocity
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      mesh_velocity(mesh->get_num_entities(0), DataType(0));

    // Initial time
    DataType time(0);
    // Timestep size
    std::cout << "deltat = " << stringify_fp_sci(deltat) << std::endl;

    // Counter for timesteps
    Index n(0);

    // This is the absolute turning angle of the screws
    DataType alpha(0);
    // Need some pi for all the angles
    DataType pi(Math::pi<DataType>());

    while(time < DataType(1))
    {
      std::cout << "timestep " << n << std::endl;
      time+= deltat;

      // Save old vertex coordinates
      coords_old.clone(mr_dudv._coords);

      // Compute new target scales
      //mr_dudv.compute_h_uniform();

      DataType alpha_old = alpha;
      alpha = -DataType(2)*pi*time;

      DataType delta_alpha = alpha - alpha_old;

      // Update boundary of the inner screw
      // This is the 2x2 matrix representing the turning by the angle delta_alpha of the inner screw
      Tiny::Matrix<DataType, 2, 2> rot(DataType(0));

      rot(0,0) = Math::cos(delta_alpha*DataType(7)/DataType(6));
      rot(0,1) = - Math::sin(delta_alpha*DataType(7)/DataType(6));
      rot(1,0) = -rot(0,1);
      rot(1,1) = rot(0,0);

      ImgPointType tmp(DataType(0));
      ImgPointType tmp2(DataType(0));
      for(Index i(0); i < inner_indices.get_num_entities(); ++i)
      {
        // Index of boundary vertex i in the mesh
        Index j(inner_indices[i]);
        // Translate the point to the centre of rotation
        tmp = mr_dudv._coords(j) - x_inner;
        // Rotate
        tmp2.set_vec_mat_mult(tmp, rot);
        // Translate the point by the new centre of rotation
        mr_dudv._coords(j, x_inner + tmp2);
      }


      // Rotate the mesh in the discrete chart. This has to use an evil downcast for now
      auto* inner_chart = reinterpret_cast< Geometry::Atlas::DiscreteChart<MeshType, SurfaceMeshType>*>
        (atlas->find_mesh_chart("inner"));

      auto& vtx_inner = inner_chart->_surface_mesh->get_vertex_set();

      for(Index i(0); i < inner_chart->_surface_mesh->get_num_entities(0); ++i)
      {
        tmp = vtx_inner[i] - x_inner;
        // Rotate
        tmp2.set_vec_mat_mult(tmp, rot);
        // Translate the point by the new centre of rotation
        vtx_inner[i] = x_inner + tmp2;
      }

      //filename = "chart_inner_" + stringify(n);
      //Geometry::ExportVTK<SurfaceMeshType> writer_chart_inner(*(inner_chart->_surface_mesh));
      //writer_chart_inner.write(filename);

      // The outer screw has 7 teeth as opposed to the inner screw with 6, and it rotates at 6/7 of the speed
      rot(0,0) = Math::cos(delta_alpha);
      rot(0,1) = - Math::sin(delta_alpha);
      rot(1,0) = -rot(0,1);
      rot(1,1) = rot(0,0);

      // The outer screw rotates centrically, so x_outer remains the same at all times

      for(Index i(0); i < outer_indices.get_num_entities(); ++i)
      {
        // Index of boundary vertex i in the mesh
        Index j(outer_indices[i]);
        tmp = mr_dudv._coords(j) - x_outer;

        tmp2.set_vec_mat_mult(tmp, rot);

        mr_dudv._coords(j, x_outer + tmp2);
      }

      // Rotate the mesh in the discrete chart. This has to use an evil downcast for now
      auto* outer_chart = reinterpret_cast<Geometry::Atlas::DiscreteChart<MeshType, SurfaceMeshType>*>
        (atlas->find_mesh_chart("outer"));

      auto& vtx_outer = outer_chart->_surface_mesh->get_vertex_set();

      for(Index i(0); i < outer_chart->_surface_mesh->get_num_entities(0); ++i)
      {
        tmp = vtx_outer[i] - x_outer;
        // Rotate
        tmp2.set_vec_mat_mult(tmp, rot);
        vtx_outer[i] = x_outer + tmp2;
      }

      //filename = "chart_outer_" + stringify(n);
      //Geometry::ExportVTK<SurfaceMeshType> writer_chart_outer(*(outer_chart->_surface_mesh));
      //writer_chart_outer.write(filename);
      //
      coords_intermediate.clone(mr_dudv._coords);

      // Write new boundary to mesh
      mr_dudv.set_coords();

      // Write pre-optimisation mesh
      filename = "pre_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_field_vertex_blocked_vector("mesh_velocity", mesh_velocity);
      std::cout << "Writing " << filename << std::endl;
      writer_pre.write(filename);

      mr_dudv._coords.clone(coords_intermediate);
      mr_dudv.set_coords();

      // Optimise the mesh
      mr_dudv.optimise();

      // Compute max. mesh velocity
      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;

      for(Index i(0); i < mesh->get_num_entities(0); ++i)
      {
        mesh_velocity(i, ideltat*(mr_dudv._coords(i) - coords_old(i)));

        DataType my_mesh_velocity(mesh_velocity(i).norm_euclid());

        if(my_mesh_velocity > max_mesh_velocity)
          max_mesh_velocity = my_mesh_velocity;
      }
      std::cout << "max mesh velocity = " << stringify_fp_sci(max_mesh_velocity) << std::endl;

      // Write post-optimisation mesh
      filename = "post_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_field_vertex_blocked_vector("mesh_velocity", mesh_velocity);
      std::cout << "Writing " << filename << std::endl;
      writer_post.write(filename);

      n++;
    }

    // Clean up
    delete rmn;
    if(atlas != nullptr)
      delete atlas;

    return 0;

  }


}; // struct LevelsetApp

/**
 * \cond internal
 *
 * Mesh Streamer Application
 *
 */
int main(int argc, char* argv[])
{
  FEAST::Runtime::initialise(argc, argv);
  // Creata a parser for command line arguments.
  SimpleArgParser args(argc, argv);

  if( args.check("help") > -1 || args.num_args()==1)
  {
    std::cout << "Du:Dv Smoother Application for Excentric Screws usage: " << std::endl;
    std::cout << "Required arguments: --filename [String]: Path to a FEAST mesh file." << std::endl;
    std::cout << "Optional arguments: --level [unsigned int]: Number of refines, defaults to 0." << std::endl;
    exit(1);
  }
  // Specify supported command line switches
  args.support("level");
  args.support("filename");
  args.support("help");
  // Refinement level
  Index lvl_max(0);

  // Input file name, required
  String filename;
  // Get unsupported command line arguments
  std::deque<std::pair<int,String> > unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
  }

  // Check and parse --filename
  if(args.check("filename") != 1 )
  {
    std::cout << "You need to specify a mesh file with --filename.";
    throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --filename");
  }
  else
  {
    args.parse("filename", filename);
    std::cout << "Reading mesh from file " << filename << std::endl;
  }

  // Check and parse --level
  if(args.check("level") != 1)
    std::cout << "No refinement level specified, defaulting to 0." << std::endl;
  else
  {
    args.parse("level", lvl_max);
    std::cout << "Refinement level " << lvl_max << std::endl;
  }

  // Create a MeshStreamer and read the mesh file
  MeshStreamer my_streamer;
  my_streamer.parse_mesh_file(filename);

  // This is the raw mesh data my_streamer read from filename
  auto& mesh_data = my_streamer.get_root_mesh_node()->mesh_data;
  // Marker int for the MeshType
  int mesh_type = mesh_data.mesh_type;
  // Marker int for the ShapeType
  int shape_type = mesh_data.shape_type;

  ASSERT(mesh_type == mesh_data.mt_conformal, "This application only works for conformal meshes!");

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  DataType deltat(DataType(1e-3));

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, DataType> Simplex2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, DataType> Hypercube2Mesh_2d;

  int ret(1);

  // Call the run() method of the appropriate wrapper class
  if(shape_type == mesh_data.st_tria)
  {
    ret = DuDvSmootherExcentricApp<MemType, DataType, IndexType, Simplex2Mesh_2d>::
      run(my_streamer, lvl_max, deltat);
  }

  if(shape_type == mesh_data.st_quad)
  {
    ret = DuDvSmootherExcentricApp<MemType, DataType, IndexType, Hypercube2Mesh_2d>::
      run(my_streamer, lvl_max, deltat);
  }

  ret = ret | FEAST::Runtime::finalise();
  // If no MeshType from the list was in the file, return 1
  return ret;
}
/// \endcond
