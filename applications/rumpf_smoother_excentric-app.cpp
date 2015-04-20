#include <iostream>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother_q1hack.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1_d2.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1_d2.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1hack.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>
#include <kernel/util/mesh_streamer.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/assembly/common_functions.hpp>

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

using namespace FEAST::Geometry;

/**
 * \brief Wrapper struct as functions do not seem to agree with template template parameters
 **/
template
<
  typename DataType_,
  typename MemType_,
  typename MeshType_,
  template<typename ... > class RumpfSmootherType_,
  template<typename, typename> class FunctionalType_
  > struct RumpfSmootherExcentricApp
{
  typedef DataType_ DataType;
  typedef MemType_ MemType;
  typedef MeshType_ MeshType;
  typedef typename MeshType::ShapeType ShapeType;

  /**
   * \brief Routine that does the actual work
   *
   * \param[in] my_streamer
   * MeshStreamer that contains the data from the mesh file.
   *
   * \param[in] level
   * Number of refines.
   */
  static int run(MeshStreamer& my_streamer, const Index lvl_max)
  {
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef FunctionalType_<DataType, ShapeType> FunctionalType;

    // create atlas
    std::cout << "Creating mesh atlas..." << std::endl;
    MeshAtlas<MeshType_>* atlas = nullptr;
    try
    {
      atlas = new MeshAtlas<MeshType_>(my_streamer);
    }
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return 1;
    }

    // create mesh node
    std::cout << "Creating mesh node..." << std::endl;
    RootMeshNode<MeshType_>* rmn = nullptr;
    try
    {
      rmn = new RootMeshNode<MeshType_>(my_streamer, atlas);
      rmn ->adapt();
    }
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return 1;
    }

    // get all mesh part names
    //std::deque<String> part_names = node->get_mesh_part_names();

    // refine
    for(Index lvl(1); lvl <= lvl_max; ++lvl)
    {
      std::cout << "Refining up to level " << lvl << "..." << std::endl;
      auto* old = rmn;
      rmn = old->refine();
      delete old;
    }

    MeshType* mesh = rmn->get_mesh();
    TrafoType trafo(*mesh);

    auto& screw_1_indices = rmn->find_mesh_part("inner")->template get_target_set<0>();
    auto& screw_2_indices = rmn->find_mesh_part("outer")->template get_target_set<0>();

    typedef Tiny::Vector<DataType,2,2> ImgPointType;

    DataType excentricity_screw_1(DataType(0.2833));

    ImgPointType x_0(DataType(0));

    ImgPointType x_1(DataType(0));
    x_1.v[0] = -excentricity_screw_1;

    ImgPointType x_2(DataType(0));

    typedef RumpfSmootherType_
    <
      DataType,
      MemType,
      TrafoType,
      FunctionalType
    > RumpfSmootherType;

    DataType deltat(DataType(1e-4));
    DataType pi(Math::pi<DataType>());

    DataType fac_norm = DataType(1e-0),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

    // The smoother in all its template glory
    RumpfSmootherType rumpflpumpfl(trafo, my_functional);

    // Print lotsa information
    std::cout << __func__ << " at refinement level " << lvl_max << std::endl;
    std::cout << "deltat = " << scientify(deltat) << std::endl;

    rumpflpumpfl.init();
    rumpflpumpfl.print();

    DataType* func_norm(new DataType[mesh->get_num_entities(2)]);
    DataType* func_det(new DataType[mesh->get_num_entities(2)]);
    DataType* func_rec_det(new DataType[mesh->get_num_entities(2)]);

    // Evaluates the levelset function and its gradient
    rumpflpumpfl.prepare();
    // Compute initial functional value
    DataType fval(0);

    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
    std::cout << "fval pre optimisation = " << scientify(fval) << " cell size quality indicator: " << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;
    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    // Filename for vtk files
    std::string filename;

    Geometry::ExportVTK<MeshType> writer_pre_initial(*mesh);
    writer_pre_initial.add_scalar_cell("norm", func_norm);
    writer_pre_initial.add_scalar_cell("det", func_det);
    writer_pre_initial.add_scalar_cell("rec_det", func_rec_det);
    writer_pre_initial.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_pre_initial.add_scalar_cell("h_0", rumpflpumpfl._h[0].elements() );
    writer_pre_initial.add_scalar_cell("h_1", rumpflpumpfl._h[1].elements() );
    writer_pre_initial.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
    writer_pre_initial.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh->get_num_entities(0)]);
    writer_pre_initial.write("pre_initial.vtk");

    rumpflpumpfl.optimise();

    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
    std::cout << "fval post optimisation = " << scientify(fval) << " cell size quality indicator: " << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;

    Geometry::ExportVTK<MeshType> writer_post_initial(*mesh);
    writer_post_initial.add_scalar_cell("norm", func_norm);
    writer_post_initial.add_scalar_cell("det", func_det);
    writer_post_initial.add_scalar_cell("rec_det", func_rec_det);
    writer_post_initial.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_post_initial.add_scalar_cell("h_0", rumpflpumpfl._h[0].elements() );
    writer_post_initial.add_scalar_cell("h_1", rumpflpumpfl._h[1].elements() );
    writer_post_initial.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
    writer_post_initial.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh->get_num_entities(0)]);
    writer_post_initial.write("post_initial.vtk");

    DataType time(0);
    Index n(0);

    DataType* mesh_velocity(new DataType[mesh->get_num_entities(0)]);

    // Old mesh coordinates for computing the mesh velocity
    LAFEM::DenseVector<MemType, DataType_> coords_old[MeshType::world_dim];
    for(Index d = 0; d < MeshType::world_dim; ++d)
      coords_old[d]= std::move(LAFEM::DenseVector<MemType, DataType>(mesh->get_num_entities(0)));

    DataType alpha(0);

    while(time < DataType(1))
    {
      std::cout << "timestep " << n << std::endl;
      time+= deltat;
      // Save old vertex coordinates
      for(Index d(0); d < MeshType::world_dim; ++d)
      {
        for(Index i(0); i < mesh->get_num_entities(0); ++i)
          coords_old[d](i, rumpflpumpfl._coords[d](i));
      }

      DataType alpha_old = alpha;
      alpha = -DataType(2)*pi*time;

      DataType delta_alpha = alpha - alpha_old;

      Tiny::Matrix<DataType, 2, 2> rot(DataType(0));

      rot(0,0) = Math::cos(delta_alpha);
      rot(0,1) = - Math::sin(delta_alpha);
      rot(1,0) = -rot(0,1);
      rot(1,1) = rot(0,0);

      ImgPointType x_1_old(x_1);

      x_1.v[0] = x_0.v[0] - excentricity_screw_1*Math::cos(alpha);
      x_1.v[1] = x_0.v[1] - excentricity_screw_1*Math::sin(alpha);

      // Boundary update case
      ImgPointType tmp(DataType(0));
      ImgPointType tmp2(DataType(0));
      for(Index i(0); i < screw_1_indices.get_num_entities(); ++i)
      {
        Index j(screw_1_indices[i]);
        tmp.v[0] = (rumpflpumpfl._coords[0](j) - x_1_old.v[0]);
        tmp.v[1] = (rumpflpumpfl._coords[1](j) - x_1_old.v[1]);

        tmp2.set_vec_mat_mult(tmp, rot);

        rumpflpumpfl._coords[0](j, tmp2.v[0] + x_1.v[0]);
        rumpflpumpfl._coords[1](j, tmp2.v[1] + x_1.v[1]);
      }

      rot(0,0) = Math::cos(delta_alpha*DataType(6)/DataType(7));
      rot(0,1) = - Math::sin(delta_alpha*DataType(6)/DataType(7));
      rot(1,0) = -rot(0,1);
      rot(1,1) = rot(0,0);

      for(Index i(0); i < screw_2_indices.get_num_entities(); ++i)
      {
        Index j(screw_2_indices[i]);
        tmp.v[0] = (rumpflpumpfl._coords[0](j) - x_2.v[0]);
        tmp.v[1] = (rumpflpumpfl._coords[1](j) - x_2.v[1]);

        tmp2.set_vec_mat_mult(tmp, rot);

        rumpflpumpfl._coords[0](j, tmp2.v[0] + x_2.v[0]);
        rumpflpumpfl._coords[1](j, tmp2.v[1] + x_2.v[1]);
      }

      rumpflpumpfl.set_coords();

      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
      std::cout << "fval pre optimisation = " << scientify(fval) << " cell size quality indicator: " << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;
      rumpflpumpfl.compute_gradient();

      filename = "pre_" + stringify(n) + ".vtk";
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
      writer_pre.add_scalar_cell("h_0", rumpflpumpfl._h[0].elements() );
      writer_pre.add_scalar_cell("h_1", rumpflpumpfl._h[1].elements() );
      writer_pre.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
      writer_pre.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh->get_num_entities(0)]);
      writer_pre.add_scalar_cell("norm", func_norm);
      writer_pre.add_scalar_cell("det", func_det);
      writer_pre.add_scalar_cell("rec_det", func_rec_det);
      std::cout << "Writing " << filename << std::endl;
      writer_pre.write(filename);

      rumpflpumpfl.optimise();

      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;
      // Compute grid velocity
      for(Index i(0); i < mesh->get_num_entities(0); ++i)
      {
        mesh_velocity[i] = DataType(0);
        for(Index d(0); d < MeshType::world_dim; ++d)
          mesh_velocity[i] += Math::sqr(ideltat*(coords_old[d](i) - rumpflpumpfl._coords[d](i)));

        mesh_velocity[i] = Math::sqrt(mesh_velocity[i]);
        if(mesh_velocity[i] > max_mesh_velocity)
          max_mesh_velocity = mesh_velocity[i];
      }
      std::cout << "max mesh velocity = " << scientify(max_mesh_velocity) << std::endl;

      fval = rumpflpumpfl.compute_functional(func_norm,func_det,func_rec_det);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval post optimisation = " << scientify(fval) << " cell size quality indicator: " << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;

      filename = "post_" + stringify(n) + ".vtk";
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
      writer_post.add_scalar_cell("h_0", rumpflpumpfl._h[0].elements() );
      writer_post.add_scalar_cell("h_1", rumpflpumpfl._h[1].elements() );
      writer_post.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
      writer_post.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh->get_num_entities(0)]);
      writer_post.add_scalar_cell("norm", func_norm);
      writer_post.add_scalar_cell("det", func_det);
      writer_post.add_scalar_cell("rec_det", func_rec_det);
      writer_post.add_scalar_vertex("mesh_velocity", mesh_velocity);
      writer_post.write(filename);

      n++;
    }
    delete mesh_velocity;

    delete func_norm;
    delete func_det;
    delete func_rec_det;

    return 0;

  }


}; // struct LevelsetApp

template<typename A, typename B>
using MyFunctional= Geometry::RumpfFunctional_D2<A, B>;

template<typename A, typename B>
using MyFunctionalQ1Hack = Geometry::RumpfFunctionalQ1Hack<A, B, Geometry::RumpfFunctional_D2>;

template<typename A, typename B, typename C, typename D>
using MySmoother = Geometry::RumpfSmoother<A, B, C, D>;

template<typename A, typename B, typename C, typename D>
using MySmootherQ1Hack = Geometry::RumpfSmootherQ1Hack<A, B, C, D>;


/**
 * \cond internal
 *
 * Mesh Streamer Application
 *
 */
int main(int argc, char* argv[])
{
  // Creata a parser for command line arguments.
  SimpleArgParser args(argc, argv);

  if( args.check("help") > -1 || args.num_args()==1)
  {
    std::cout << "Rumpf Smoother Application for Excentric Screws usage: " << std::endl;
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
  FEAST::String filename;
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
    throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --filename");
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

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<1>, 1, 1, Real> Simplex1Mesh_1d;
  typedef Geometry::ConformalMesh<Shape::Simplex<1>, 2, 2, Real> Simplex1Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Simplex<1>, 3, 3, Real> Simplex1Mesh_3d;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> Simplex2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, 3, Real> Simplex2Mesh_3d;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, Real> Simplex3Mesh;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> Hypercube2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, 3, Real> Hypercube2Mesh_3d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, Real> Hypercube3Mesh;
  // Call the run() method of the appropriate wrapper class
  if(shape_type == mesh_data.st_tria)
    return RumpfSmootherExcentricApp<double, Mem::Main, Simplex2Mesh_2d, MySmoother, RumpfFunctional_D2>::run(my_streamer, lvl_max);
  if(shape_type == mesh_data.st_quad)
    return RumpfSmootherExcentricApp<double, Mem::Main, Hypercube2Mesh_2d, MySmoother, MyFunctional>::run(my_streamer, lvl_max);

  // If no MeshType from the list was in the file, return 1
  return 1;
}
/// \endcond
