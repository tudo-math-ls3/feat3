#include <kernel/base_header.hpp>
#ifdef FEAST_HAVE_ALGLIB
#include <kernel/archs.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/meshopt/rumpf_smoother_conc.hpp>
#include <kernel/meshopt/rumpf_smoother_lvlset.hpp>
#include <kernel/meshopt/rumpf_smoother_lvlset_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/lvlset.hpp>
#include <kernel/meshopt/rumpf_functionals/conc_2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/conc_2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/conc_2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/conc_2d_q1_d2.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/simple_arg_parser.hpp>

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
 * \brief Wrapper struct as functions do not seem to agree with template template parameters
 **/
template
<
  typename DT_,
  typename MeshType_,
  template<typename ... > class RumpfSmootherType_,
  template<typename, typename> class FunctionalType_,
  template<typename, typename> class LevelsetFunctionalType_
> struct LevelsetApp
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
  /// The corresponding transformation
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  /// The Rumpf functional variant
  typedef FunctionalType_<DataType, ShapeType> FunctionalType;
  /// The lvl_maxset part of the Rumpf functional
  typedef LevelsetFunctionalType_<DataType, ShapeType> LevelsetFunctionalType;
  /// Type for points in the mesh
  typedef Tiny::Vector<DataType,2,2> ImgPointType;

  /**
   * \brief Routine that does the actual work
   *
   * \param[in] lvl_max
   * Number of refines.
   */
  static int run(Geometry::MeshFileReader& file_reader, Index lvl_max, DT_ deltat)
  {

    // create an empty atlas and a root mesh node
    Geometry::MeshAtlas<MeshType>* atlas = new Geometry::MeshAtlas<MeshType>();
    Geometry::RootMeshNode<MeshType>* rmn = new Geometry::RootMeshNode<MeshType>(nullptr, atlas);

    // try to parse the mesh files
#ifndef DEBUG
    try
#endif
    {
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
      std::cout << "Refining up to lvl_max " << lvl << "..." << std::endl;
      auto* old = rmn;
      rmn = old->refine();
      delete old;
    }

    MeshType* mesh = rmn->get_mesh();

    typedef Analytic::Common::DistanceFunctionSD<2, DataType> AnalyticFunctionType1;
    typedef Analytic::Common::MinOfTwoFunctions<AnalyticFunctionType1, AnalyticFunctionType1> AnalyticFunctionType;

    // Reference point for the rotation
    ImgPointType x0(DataType(0));

    centre_point_outer(x0,DataType(0));
    DataType scaling_outer(-1);
    DataType radius_outer(0.35);

    // Analytic function for the distance to the outer circle's boundary
    AnalyticFunctionType1 outer(x0, radius_outer, scaling_outer);

    DataType scaling_inner(-scaling_outer);
    DataType radius_inner(-scaling_inner*0.275);
    centre_point_inner(x0,DataType(0));

    // Analytic function for the distance to the inner circle's boundary
    AnalyticFunctionType1 inner(x0, radius_inner, scaling_inner);

    AnalyticFunctionType analytic_lvlset(inner, outer);

    // Now we can define the whole smoother
    typedef RumpfSmootherType_
    <
      AnalyticFunctionType,
      TrafoType,
      FunctionalType,
      LevelsetFunctionalType
    > RumpfSmootherType;

    // Parameters for the functional
    DataType fac_norm = DataType(1e0),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
    // Create the functional
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);
    // Parameters for the Rumpf smoother
    bool align_to_lvlset(false);
    bool r_adaptivity(true);
    DataType r_adapt_reg = DataType(1e-2), r_adapt_pow = DataType(0.5);
    // Create lvl_maxset part of the functional
    LevelsetFunctionalType my_lvl_maxset_functional;

    // Set slip boundary conditions at these parts of the boundary
    std::deque<String> slip_list;
    slip_list.push_back("left");
    slip_list.push_back("right");

    // Set Dirichlet boundary conditions at these parts of the boundary
    std::deque<String> dirichlet_list;
    dirichlet_list.push_back("bottom");
    dirichlet_list.push_back("top");

    // The smoother in all its template glory
    RumpfSmootherType rumpflpumpfl(rmn, dirichlet_list, slip_list, my_functional, my_lvl_maxset_functional,
    align_to_lvlset, r_adaptivity, analytic_lvlset);

    // Set the r-adaptivity parameters
    rumpflpumpfl.set_r_adapt_params(r_adapt_reg,r_adapt_pow);

    rumpflpumpfl.init();
    rumpflpumpfl.print();

    // Arrays for saving the contributions of the different Rumpf functional parts
    DataType* func_norm(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_rec_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_lvlset(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);

    // Evaluates the lvl_maxset function and its gradient
    rumpflpumpfl.prepare();
    // Compute initial functional value
    DataType fval(0);
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
    std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << " cell size quality indicator: " << stringify_fp_sci(rumpflpumpfl.cell_size_quality()) << std::endl;
    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    // Filename for vtk files
    std::string filename;

    // Write initial state to file
    Geometry::ExportVTK<MeshType> writer_pre_initial(*mesh);
    writer_pre_initial.add_cell_scalar("lambda", rumpflpumpfl._lambda.elements() );
    writer_pre_initial.add_cell_vector("h", rumpflpumpfl._h );
    writer_pre_initial.add_vertex_scalar("lvl_maxset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_pre_initial.add_vertex_vector("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec);
    writer_pre_initial.add_cell_vector("fval", func_norm, func_det, func_rec_det);
    writer_pre_initial.add_cell_scalar("lvl_maxset_constraint", func_lvlset );
    writer_pre_initial.add_vertex_vector("grad", rumpflpumpfl._grad);
    writer_pre_initial.write("pre_initial");

    // Optimise the mesh
    rumpflpumpfl.optimise();

    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
    std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << " cell size quality indicator: " << stringify_fp_sci(rumpflpumpfl.cell_size_quality()) << std::endl;

    // Call prepare() again because the mesh changed due to the optimisation and it was not called again after the
    // last iteration
    rumpflpumpfl.prepare();
    rumpflpumpfl.compute_gradient();

    // Write optimised initial mesh
    Geometry::ExportVTK<MeshType> writer_post_initial(*mesh);
    writer_post_initial.add_cell_scalar("lambda", rumpflpumpfl._lambda.elements() );
    writer_post_initial.add_cell_vector("h", rumpflpumpfl._h );
    writer_post_initial.add_vertex_scalar("lvl_maxset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_post_initial.add_vertex_vector("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec);
    writer_post_initial.add_cell_vector("fval", func_norm, func_det, func_rec_det);
    writer_post_initial.add_cell_scalar("lvl_maxset_constraint", func_lvlset );
    writer_post_initial.add_vertex_vector("grad", rumpflpumpfl._grad);
    writer_post_initial.write("post_initial");

    DataType time(0);
    Index n(0);

    // Old mesh coordinates for computing the mesh velocity
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      coords_old(mesh->get_num_entities(0),DataType(0));
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      mesh_velocity(mesh->get_num_entities(0), DataType(0));

    while(time < DataType(1))
    {
      std::cout << "timestep " << n << std::endl;
      time+= deltat;

      // Save old vertex coordinates
      coords_old.clone(rumpflpumpfl._coords);

      // Update lvl_maxset function
      centre_point_outer(x0,time);
      outer.set_point(x0);

      // Update reference point
      centre_point_inner(x0,time);
      inner.set_point(x0);

      // Compute functional value and gradient
      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << " cell size quality indicator: " <<
        stringify_fp_sci(rumpflpumpfl.cell_size_quality()) << std::endl;

      // Write pre-optimisation mesh
      filename = "pre_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_cell_scalar("lambda", rumpflpumpfl._lambda.elements() );
      writer_pre.add_cell_vector("h", rumpflpumpfl._h);
      writer_pre.add_vertex_scalar("lvl_maxset", rumpflpumpfl._lvlset_vtx_vec.elements());
      writer_pre.add_vertex_vector("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec);
      writer_pre.add_cell_vector("fval", func_norm, func_det, func_rec_det);
      writer_pre.add_cell_scalar("lvl_maxset_constraint", func_lvlset );
      writer_pre.add_vertex_vector("grad", rumpflpumpfl._grad);
      writer_pre.write(filename);

      // Optimise the mesh
      rumpflpumpfl.optimise();

      // Compute grid velocity
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

      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << " cell size quality indicator: "
      << stringify_fp_sci(rumpflpumpfl.cell_size_quality()) << std::endl;

      // Write post-optimisation mesh
      filename = "post_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_cell_scalar("lambda", rumpflpumpfl._lambda.elements() );
      writer_post.add_cell_vector("h", rumpflpumpfl._h);
      writer_post.add_vertex_scalar("lvl_maxset", rumpflpumpfl._lvlset_vtx_vec.elements());
      writer_post.add_vertex_vector("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec);
      writer_post.add_cell_vector("fval", func_norm, func_det, func_rec_det);
      writer_post.add_cell_scalar("lvl_maxset_constraint", func_lvlset );
      writer_post.add_vertex_vector("grad", rumpflpumpfl._grad);
      writer_post.add_vertex_vector("mesh_velocity", mesh_velocity);
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
    delete[] func_lvlset;

    return 0;

  }

}; // struct LevelsetApp

template<typename A, typename B, typename C, typename D>
using MySmoother = Meshopt::RumpfSmootherLevelsetConcAnalytic<A, B, C, D>;

template<typename A, typename B, typename C, typename D, typename E, typename F>
using MySmootherQ1Hack = Meshopt::RumpfSmootherLevelsetAnalyticQ1Hack<A, B, C, D>;

template<typename A, typename B>
using MyFunctional= Meshopt::RumpfFunctionalConc<A, B>;

template<typename A, typename B>
using MyFunctionalQ1Hack = Meshopt::RumpfFunctionalQ1Hack<A, B, Meshopt::RumpfFunctional>;

int main(int argc, char* argv[])
{
  typedef double DataType;

  FEAST::Runtime::initialise(argc, argv);
  // Creata a parser for command line arguments.
  SimpleArgParser args(argc, argv);

  if( args.check("help") > -1 || args.num_args()==1)
  {
    std::cout << "Rumpf Smoother Application for Levelset usage: " << std::endl;
    std::cout << "Required arguments: --filename [String]: Path to a FEAST mesh file." << std::endl;
    std::cout << "Optional arguments: --lvl_max [unsigned int]: Number of refines, defaults to 0." << std::endl;
    exit(1);
  }
  // Specify supported command line switches
  args.support("lvl_max");
  args.support("filename");
  args.support("help");
  // Refinement lvl_max
  Index lvl_max(0);
  // Input file name, required
  FEAST::String meshfile;
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

  // Check and parse --lvl_max
  if(args.check("level") != 1)
    std::cout << "No refinement level specified, defaulting to 0." << std::endl;
  else
  {
    args.parse("level", lvl_max);
    std::cout << "Refinement level" << lvl_max << std::endl;
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

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  DataType deltat(DataType(1e-2));

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> Simplex2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> Hypercube2Mesh_2d;

  int ret(1);

  // Call the run() method of the appropriate wrapper class
  if(mtype == "conformal:simplex:2:2")
    ret = LevelsetApp<DataType, Simplex2Mesh_2d, MySmoother, MyFunctional, Meshopt::RumpfFunctionalLevelset>::
      run(mesh_reader, lvl_max, deltat);

  if(mtype == "conformal:hypercube:2:2")
    ret = LevelsetApp<DataType, Hypercube2Mesh_2d, MySmoother, MyFunctional, Meshopt::RumpfFunctionalLevelset>::
      run(mesh_reader, lvl_max, deltat);

  ret = ret | FEAST::Runtime::finalise();
  // If no MeshType from the list was in the file, ret is 1
  return ret;
}
/// \endcond
#endif // FEAST_HAVE_ALGLIB
