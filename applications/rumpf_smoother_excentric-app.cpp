#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/meshopt/rumpf_smoother.hpp>
#include <kernel/meshopt/rumpf_smoother_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1hack.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <kernel/util/runtime.hpp>

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
 * \brief This application demonstrates the usage of some of the RumpfSmoother classes for boundary deformations
 *
 * \note Because an application of the (nonlinear) Rumpf smoother requires operations similar to a matrix assembly,
 * Rumpf smoothers are implemented for Mem::Main only.
 *
 * In this application, a mesh with two excentric screws is read from a mesh. The screws rotate with different
 * angular velocities, so large mesh deformations occur.
 *
 * \author Jordi Paul
 *
 * \tparam DT_
 * The precision of the mesh etc.
 *
 * \tparam MeshType
 * The mesh type, has to be known because we stream the mesh from a file
 *
 * \tparam FunctionalType
 * The Rumpf functional variant to use
 *
 * \tparam RumpfSmootherType_
 * The Rumpf smoother variant to use
 *
 **/
/**
 * \brief Wrapper struct as functions do not seem to agree with template template parameters
 **/
template
<
  typename DT_,
  typename MeshType_,
  template<typename, typename> class FunctionalType_,
  template<typename ... > class RumpfSmootherType_
  >
  struct RumpfSmootherExcentricApp
{
  /// Precision for meshes etc, everything else uses the same data type
  typedef DT_ DataType;
  /// Rumpf Smoothers are implemented for Mem::Main only
  typedef Mem::Main MemType;
  /// So we use Index
  typedef Index IndexType;
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
  /// Our functional type
  typedef FunctionalType_<DataType, ShapeType> FunctionalType;
  /// The Rumpf smoother
  typedef RumpfSmootherType_<TrafoType, FunctionalType> RumpfSmootherType;
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
  static int run(Geometry::MeshFileReader& file_reader, Geometry::MeshFileReader& chart_reader, Index lvl_max, DT_ deltat)
  {
    // Filename for writing .vtu output
    String filename("");
    DataType fval(0);

    // Parameters for the Rumpf functional
    DataType fac_norm = DataType(1e-1),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

    std::deque<String> dirichlet_list;
    std::deque<String> slip_list;
    slip_list.push_back("outer");
    slip_list.push_back("inner");

    // create an empty atlas and a root mesh node
    Geometry::MeshAtlas<MeshType>* atlas = new Geometry::MeshAtlas<MeshType>();
    Geometry::RootMeshNode<MeshType>* rmn = new Geometry::RootMeshNode<MeshType>(nullptr, atlas);

    // try to parse the mesh files
#ifndef DEBUG
    try
#endif
    {
      chart_reader.parse(*rmn, *atlas);
      file_reader.parse(*rmn, *atlas);
    }
#ifndef DEBUG
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return 1;
    }
    catch(...)
    {
      std::cerr << "ERROR: unknown exception" << std::endl;
      return 1;
    }
#endif

    // refine
    for(Index lvl(1); lvl <= lvl_max; ++lvl)
    {
      std::cout << "Refining up to level " << lvl << "..." << std::endl;
      // Arrays for saving the contributions of the different Rumpf functional parts
      DataType* func_norm(new DataType[rmn->get_mesh()->get_num_entities(MeshType::shape_dim)]);
      DataType* func_det(new DataType[rmn->get_mesh()->get_num_entities(MeshType::shape_dim)]);
      DataType* func_rec_det(new DataType[rmn->get_mesh()->get_num_entities(MeshType::shape_dim)]);

      RumpfSmootherType refinement_smoother(rmn, dirichlet_list, slip_list, my_functional);
      refinement_smoother.init();
      //refinement_smoother.compute_lambda_uniform();
      //refinement_smoother.compute_h();
      refinement_smoother.prepare();

      fval = refinement_smoother.compute_functional(func_norm, func_det, func_rec_det);
      std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << std::endl;

      // Compute initial functional gradient
      refinement_smoother.compute_gradient();

      // Write initial state to file
      filename = "refinement_pre_"+stringify(lvl);
      Geometry::ExportVTK<MeshType> writer_refinement_pre(*(rmn->get_mesh()));
      writer_refinement_pre.add_cell_vector("h", refinement_smoother._h);
      writer_refinement_pre.add_cell_vector("fval", func_norm, func_det, func_rec_det);
      writer_refinement_pre.add_vertex_vector("grad", refinement_smoother._grad);
      writer_refinement_pre.write(filename);

      refinement_smoother.optimise();

      fval = refinement_smoother.compute_functional(func_norm, func_det, func_rec_det);
      std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << std::endl;

      // Compute initial functional gradient
      refinement_smoother.compute_gradient();
      filename = "refinement_post_"+stringify(lvl);
      Geometry::ExportVTK<MeshType> writer_refinement_post(*(rmn->get_mesh()));
      writer_refinement_post.add_cell_vector("h", refinement_smoother._h);
      writer_refinement_post.add_cell_vector("fval", func_norm, func_det, func_rec_det);
      writer_refinement_post.add_vertex_vector("grad", refinement_smoother._grad);
      writer_refinement_post.write(filename);

      auto* old = rmn;
      rmn = old->refine();
      delete old;

      delete[] func_norm;
      delete[] func_det;
      delete[] func_rec_det;

    }

    rmn->adapt();

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
    RumpfSmootherType rumpflpumpfl(rmn, dirichlet_list, slip_list, my_functional);
    rumpflpumpfl.init();
    rumpflpumpfl.print();

    rumpflpumpfl.prepare();

    // Arrays for saving the contributions of the different Rumpf functional parts
    DataType* func_norm(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_rec_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);

    // Compute initial functional value
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
    std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << std::endl;

    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    // Write initial state to file
    Geometry::ExportVTK<MeshType> writer_initial_pre(*mesh);
    writer_initial_pre.add_cell_vector("h", rumpflpumpfl._h);
    writer_initial_pre.add_cell_vector("fval", func_norm, func_det, func_rec_det);
    writer_initial_pre.add_vertex_vector("grad", rumpflpumpfl._grad);
    writer_initial_pre.write("pre_initial");

    // Smooth the mesh
    rumpflpumpfl.optimise();

    // Call prepare() again because the mesh changed due to the optimisation and it was not called again after the
    // last iteration
    rumpflpumpfl.prepare();
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
    rumpflpumpfl.compute_gradient();

    std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << std::endl;

    // Write optimised initial mesh
    Geometry::ExportVTK<MeshType> writer_initial_post(*mesh);
    writer_initial_post.add_cell_vector("h", rumpflpumpfl._h);
    writer_initial_post.add_cell_vector("fval", func_norm, func_det, func_rec_det);
    writer_initial_post.add_vertex_vector("grad", rumpflpumpfl._grad);
    writer_initial_post.write("post_initial");

    // For saving the old coordinates
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      coords_old(mesh->get_num_entities(0),DataType(0));
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
      coords_old.clone(rumpflpumpfl._coords);

      // Compute new target scales
      //rumpflpumpfl.compute_h_uniform();

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
        tmp = rumpflpumpfl._coords(j) - x_inner;
        // Rotate
        tmp2.set_vec_mat_mult(tmp, rot);
        // Translate the point by the new centre of rotation
        rumpflpumpfl._coords(j, x_inner + tmp2);
      }

      // Rotate the chart. This has to use an evil downcast for now
      auto* inner_chart = reinterpret_cast< Geometry::Atlas::Polyline<MeshType>*>
        (atlas->find_mesh_chart("inner"));

      auto& vtx_inner = inner_chart->get_world_points();

      for(auto& it : vtx_inner)
      {
        tmp = it - x_inner;
        // Rotate
        tmp2.set_vec_mat_mult(tmp, rot);
        // Translate the point by the new centre of rotation
        it = x_inner + tmp2;
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
        tmp = rumpflpumpfl._coords(j) - x_outer;

        tmp2.set_vec_mat_mult(tmp, rot);

        rumpflpumpfl._coords(j, x_outer + tmp2);
      }

      // Rotate the outer chart. This has to use an evil downcast for now
      auto* outer_chart = reinterpret_cast<Geometry::Atlas::Polyline<MeshType>*>
        (atlas->find_mesh_chart("outer"));

      auto& vtx_outer = outer_chart->get_world_points();

      for(auto& it :vtx_outer)
      {
        tmp = it - x_outer;
        // Rotate
        tmp2.set_vec_mat_mult(tmp, rot);
        it = x_outer + tmp2;
      }

      const auto& idx = mesh->template get_index_set<ShapeType::dimension, 0>();
      const auto& vtx = mesh->get_vertex_set();

      typename LAFEM::DenseVector<Mem::Main, DataType, Index> tmp_lambda(mesh->get_num_entities(ShapeType::dimension));
      DataType sum_lambda(0);

      for(Index cell(0); cell < mesh->get_num_entities(ShapeType::dimension); ++cell)
      {
        // Compute midpoint of current cell
        ImgPointType midpoint(DataType(0));
        for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
        {
          Index i(idx(cell, Index(j)));
          midpoint += vtx[i];
        }
        midpoint *= (DataType(1))/DataType(Shape::FaceTraits<ShapeType,0>::count);

        DataType dist_inner(Math::Limits<DataType>::max());
        for(Index j(0); j < inner_indices.get_num_entities(); ++j)
        {
          Index i(inner_indices[j]);
          DataType my_dist = (midpoint - vtx[i]).norm_euclid();
          if(my_dist < dist_inner)
            dist_inner = my_dist;
        }

        DataType dist_outer(Math::Limits<DataType>::max());
        for(Index j(0); j < outer_indices.get_num_entities(); ++j)
        {
          Index i(outer_indices[j]);
          DataType my_dist = (midpoint - vtx[i]).norm_euclid();
          if(my_dist < dist_outer)
            dist_outer = my_dist;
        }

        tmp_lambda(cell, dist_inner+dist_outer);
        sum_lambda+=tmp_lambda(cell);

      }
      tmp_lambda.scale(tmp_lambda, DataType(1)/sum_lambda);
      rumpflpumpfl._lambda.convert(tmp_lambda);
      rumpflpumpfl.compute_h();

      //filename = "chart_outer_" + stringify(n);
      //Geometry::ExportVTK<SurfaceMeshType> writer_chart_outer(*(outer_chart->_surface_mesh));
      //writer_chart_outer.write(filename);

      // Write new boundary to mesh
      rumpflpumpfl.set_coords();

      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << std::endl;

      // Write pre-optimisation mesh
      filename = "pre_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_cell_vector("h", rumpflpumpfl._h);
      writer_pre.add_cell_vector("fval", func_norm, func_det, func_rec_det);
      writer_pre.add_vertex_vector("grad", rumpflpumpfl._grad);
      writer_pre.add_vertex_vector("mesh_velocity", mesh_velocity);
      std::cout << "Writing " << filename << std::endl;
      writer_pre.write(filename);

      // Optimise the mesh
      rumpflpumpfl.optimise();

      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << std::endl;

      // Compute max. mesh velocity
      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;

      for(Index i(0); i < mesh->get_num_entities(0); ++i)
      {
        mesh_velocity(i, ideltat*(rumpflpumpfl._coords(i) - coords_old(i)));

        DataType my_mesh_velocity(mesh_velocity(i).norm_euclid());

        if(my_mesh_velocity > max_mesh_velocity)
          max_mesh_velocity = my_mesh_velocity;
      }
      std::cout << "max mesh velocity = " << stringify_fp_sci(max_mesh_velocity) << std::endl;

      // Write post-optimisation mesh
      filename = "post_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_cell_vector("h", rumpflpumpfl._h);
      writer_post.add_cell_vector("fval", func_norm, func_det, func_rec_det);
      writer_post.add_vertex_vector("grad", rumpflpumpfl._grad);
      writer_post.add_vertex_vector("mesh_velocity", mesh_velocity);
      std::cout << "Writing " << filename << std::endl;
      writer_post.write(filename);

      n++;
    }

    // Clean up
    delete rmn;
    if(atlas != nullptr)
      delete atlas;

    delete[] func_norm;
    delete[] func_det;
    delete[] func_rec_det;

    return 0;

  }


}; // struct LevelsetApp

template<typename A, typename B>
using MyFunctional= Meshopt::RumpfFunctional_D2<A, B>;

template<typename A, typename B>
using MyFunctionalQ1Hack = Meshopt::RumpfFunctionalQ1Hack<A, B, Meshopt::RumpfFunctional_D2>;

template<typename A, typename B>
using MySmoother = Meshopt::RumpfSmoother<A, B>;

template<typename A, typename B>
using MySmootherQ1Hack = Meshopt::RumpfSmootherQ1Hack<A, B>;


/**
 * \cond internal
 *
 * Mesh Streamer Application
 *
 */
int main(int argc, char* argv[])
{
  typedef double DataType;
  FEAST::Runtime::initialise(argc, argv);
  // Creata a parser for command line arguments.
  SimpleArgParser args(argc, argv);

  if( args.check("help") > -1 || args.num_args()==1)
  {
    std::cout << "Rumpf Smoother Application for Excentric Screws usage: " << std::endl;
    std::cout << "Required arguments: --filename [String]: Path to a FEAST mesh file." << std::endl;
    std::cout << "Optional arguments: --level [unsigned int]: Number of refines, defaults to 0." << std::endl;
    std::cout << "                    --q1hack: Use Q1Hack functionals for hypercube meshes." << std::endl;
    exit(1);
  }
  // Specify supported command line switches
  args.support("level");
  args.support("meshfile");
  args.support("chartfile");
  args.support("help");
  args.support("q1hack");
  // Refinement level
  Index lvl_max(0);
  // Switch for the Q1Hack
  bool use_q1hack(false);

  // Input mesh file name, required
  String meshfile;
  // Input chart file name, required
  String chartfile;
  // Get unsupported command line arguments
  std::deque<std::pair<int,String> > unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
  }

  // Check and parse --meshfile
  if(args.check("meshfile") != 1 )
  {
    std::cout << "You need to specify a mesh file with --meshfile.";
    throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --filename");
  }
  else
  {
    args.parse("meshfile", meshfile);
    std::cout << "Reading mesh from file " << meshfile << std::endl;
  }

  // Check and parse --chartfile
  if(args.check("chartfile") != 1 )
  {
    std::cout << "You need to specify a chart file with --chartfile.";
    throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --filename");
  }
  else
  {
    args.parse("chartfile", chartfile);
    std::cout << "Reading chart from file " << chartfile << std::endl;
  }

  // Check and parse --level
  if(args.check("level") != 1)
    std::cout << "No refinement level specified, defaulting to 0." << std::endl;
  else
  {
    args.parse("level", lvl_max);
    std::cout << "Refinement level " << lvl_max << std::endl;
  }

  // Check and parse --q1hack
  if(args.check("q1hack") != 0)
  {
    std::cout << "Not using the Q1Hack for hypercube meshes." << std::endl;
  }
  else
  {
    use_q1hack = true;
    std::cout << "Using the Q1Hack for hypercube meshes." << std::endl;
  }

  // Create mesh input file stream
  std::ifstream ifs_mesh(meshfile, std::ios_base::in);
  if(!ifs_mesh.is_open() || !ifs_mesh.good())
  {
    std::cerr << "ERROR: Failed to open '" << meshfile << "'" << std::endl;
    return 1;
  }
  // Create a reader and read the root markup
  Geometry::MeshFileReader mesh_reader(ifs_mesh);
  mesh_reader.read_root_markup();

  // Create chart input file stream
  std::ifstream ifs_chart(chartfile, std::ios_base::in);
  if(!ifs_chart.is_open() || !ifs_chart.good())
  {
    std::cerr << "ERROR: Failed to open '" << chartfile << "'" << std::endl;
    return 1;
  }

  // Create a reader for the chart file
  Geometry::MeshFileReader chart_reader(ifs_chart);

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, DataType> Simplex2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, DataType> Hypercube2Mesh_2d;

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  int ret(1);

  DataType deltat(DataType(1e-4));

  // Call the run() method of the appropriate wrapper class
  if(mtype == "conformal:simplex:2:2")
    ret = RumpfSmootherExcentricApp<DataType, Simplex2Mesh_2d, MyFunctional, MySmoother>::
      run(mesh_reader, chart_reader, lvl_max, deltat);

  if(mtype == "conformal:hypercube:2:2")
  {
    if(use_q1hack)
    {
      ret = RumpfSmootherExcentricApp<DataType, Hypercube2Mesh_2d, MyFunctionalQ1Hack, MySmootherQ1Hack>::
        run(mesh_reader, chart_reader, lvl_max, deltat);
    }
    else
    {
      ret = RumpfSmootherExcentricApp<DataType, Hypercube2Mesh_2d, MyFunctional, MySmoother>::
        run(mesh_reader, chart_reader, lvl_max, deltat);
    }
  }

  ret = ret | FEAST::Runtime::finalise();
  // If no MeshType from the list was in the file, return 1
  return ret;
}
/// \endcond
