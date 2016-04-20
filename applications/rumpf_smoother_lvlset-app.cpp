#include <kernel/base_header.hpp>
#ifdef FEAST_HAVE_ALGLIB
#include <kernel/analytic/common.hpp>
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
#include <kernel/util/simple_arg_parser.hpp>

#include <kernel/util/runtime.hpp>

using namespace FEAST;

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
  /// The levelset part of the Rumpf functional
  typedef LevelsetFunctionalType_<DataType, ShapeType> LevelsetFunctionalType;
  /// Type for points in the mesh
  typedef Tiny::Vector<DataType,2,2> ImgPointType;

  /**
   * \brief Routine that does the actual work
   *
   * \param[in] my_streamer
   * MeshStreamer that contains the data from the mesh file.
   *
   * \param[in] level
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
      std::cout << "Refining up to level " << lvl << "..." << std::endl;

      auto* old = rmn;
      rmn = old->refine();
      delete old;
    }

    MeshType* mesh = rmn->get_mesh();

    typedef Analytic::Common::DistanceFunctionSD<2, DataType> AnalyticFunctionType;
    ImgPointType x0(DataType(0));
    x0.v[0] = DataType(0.25) *(DataType(2) + Math::cos(DataType(0)));
    x0.v[1] = DataType(0.25) *(DataType(2) + Math::sin(DataType(0)));
    DataType displacement(0.15);
    DataType scaling(-1);
    AnalyticFunctionType analytic_lvlset(x0, displacement, scaling);

    //const int plane = 0;
    //typedef Assembly::Common::PlaneDistanceFunctionSD<plane, ImgPointType> AnalyticFunctionType;
    //typedef Assembly::Common::PlaneDistanceFunctionSD_grad<plane, 0, ImgPointType> AnalyticFunctionGrad0Type;
    //typedef Assembly::Common::PlaneDistanceFunctionSD_grad<plane, 1, ImgPointType> AnalyticFunctionGrad1Type;
    //ImgPointType x0(DataType(0));
    //x0.v[0] = DataType(0.5);
    //x0.v[1] = DataType(0.5);
    //DataType scaling(1);
    //AnalyticFunctionType analytic_lvlset(x0, scaling);
    //AnalyticFunctionGrad0Type analytic_lvlset_grad0(x0, scaling);
    //AnalyticFunctionGrad1Type analytic_lvlset_grad1(x0, scaling);

    // The full smoother type
    typedef RumpfSmootherType_
    <
      AnalyticFunctionType,
      TrafoType,
      FunctionalType,
      LevelsetFunctionalType
    > RumpfSmootherType;

    // Parameters for the funcitional
    DataType fac_norm = DataType(1e-2),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
    // Create the functional
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);
    // Parameters for the Rumpf smoother
    bool align_to_lvlset(true);
    bool r_adaptivity(true);
    DataType r_adapt_reg = DataType(1e-2), r_adapt_pow = DataType(0.5);
    // Create levelset part of the functional
    LevelsetFunctionalType my_levelset_functional;

    std::deque<String> slip_list;
    slip_list.push_back("bottom");
    slip_list.push_back("top");

    std::deque<String> dirichlet_list;
    dirichlet_list.push_back("left");
    dirichlet_list.push_back("right");

    // The smoother in all its template glory
    RumpfSmootherType rumpflpumpfl(rmn, slip_list, dirichlet_list, my_functional, my_levelset_functional,
    align_to_lvlset, r_adaptivity, analytic_lvlset);
    rumpflpumpfl.init();

    //rumpflpumpfl.lvlset_constraint_tol = Math::pow(Math::eps<DataType>(),DataType(0.5) );

    //rumpflpumpfl._update_h = true;

    // Set the r-adaptivity parameters
    rumpflpumpfl.set_r_adapt_params(r_adapt_reg,r_adapt_pow);

    // Print lotsa information
    rumpflpumpfl.print();

    // Arrays for saving the contributions of the different Rumpf functional parts
    DataType* func_norm(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_rec_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_lvlset(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);

    // Evaluates the levelset function and its gradient
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
    writer_pre_initial.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_pre_initial.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
    writer_pre_initial.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_pre_initial.add_field_vertex_blocked_vector("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec);
    writer_pre_initial.add_field_cell("fval", func_norm, func_det, func_rec_det);
    writer_pre_initial.add_scalar_cell("levelset_constraint", func_lvlset );
    writer_pre_initial.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
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
    writer_post_initial.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_post_initial.add_field_cell_blocked_vector("h", rumpflpumpfl._h );
    writer_post_initial.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_post_initial.add_field_vertex_blocked_vector("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec);
    writer_post_initial.add_field_cell("fval", func_norm, func_det, func_rec_det);
    writer_post_initial.add_scalar_cell("levelset_constraint", func_lvlset );
    writer_post_initial.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_post_initial.write("post_initial");

    DataType time(0);
    Index n(0);

    // Old mesh coordinates for computing the mesh velocity
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim> coords_old(mesh->get_num_entities(0),DataType(0));
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim> mesh_velocity(mesh->get_num_entities(0), DataType(0));

    std::cout << "deltat = " << stringify_fp_sci(deltat) << std::endl;
    while(time < DataType(2)*Math::pi<DataType>())
    {
      std::cout << "timestep " << n << std::endl;
      time += deltat;

      // Save old vertex coordinates
      coords_old.clone(rumpflpumpfl._coords);

      // update levelset function
      x0.v[0] = DataType(0.25) *(DataType(2) + Math::cos(time));
      x0.v[1] = DataType(0.25) *(DataType(2) + Math::sin(DataType(3)*time));
      //x0.v[0] = DataType(0.25) *(DataType(2) + Math::cos(time*Math::pi<DataType>()));
      //x0.v[1] = DataType(0.25) *(DataType(2) + Math::sin(DataType(3)*time*Math::pi<DataType>()));

      // Update the reference point of the distance function
      rumpflpumpfl._analytic_lvlset.set_point(x0);

      // Compute functional value and gradient
      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << " cell size quality indicator: " <<
        stringify_fp_sci(rumpflpumpfl.cell_size_quality()) << std::endl;

      // Write pre-optimisation mesh
      filename = "pre_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
      writer_pre.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
      writer_pre.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
      writer_pre.add_field_vertex_blocked_vector("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec);
      writer_pre.add_field_cell("fval", func_norm, func_det, func_rec_det);
      writer_pre.add_scalar_cell("levelset_constraint", func_lvlset );
      writer_pre.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
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

      // Write post-optimisation mesh
      fval = rumpflpumpfl.compute_functional(func_norm,func_det,func_rec_det,func_lvlset);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << " cell size quality indicator: "
      << stringify_fp_sci(rumpflpumpfl.cell_size_quality()) << std::endl;

      filename = "post_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_scalar_cell("norm", func_norm);
      writer_post.add_scalar_cell("det", func_det);
      writer_post.add_scalar_cell("rec_det", func_rec_det);
      writer_post.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
      writer_post.add_field_cell_blocked_vector("h", rumpflpumpfl._h );
      writer_post.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
      writer_post.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
      writer_post.add_field_vertex_blocked_vector("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec);
      writer_post.add_scalar_cell("levelset_constraint", func_lvlset );
      writer_post.add_field_vertex_blocked_vector("mesh_velocity", mesh_velocity);
      writer_post.write(filename);

      n++;
    } // time loop

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

template<typename A, typename B, typename C, typename D>
using MySmootherQ1Hack = Meshopt::RumpfSmootherLevelsetAnalyticQ1Hack<A, B, C, D>;

template<typename A, typename B>
using MyFunctional= Meshopt::RumpfFunctionalConc_D2<A, B>;

template<typename A, typename B>
using MyFunctionalQ1Hack = Meshopt::RumpfFunctionalQ1Hack<A, B, Meshopt::RumpfFunctional>;

/**
 * \cond internal
 *
 * RumpfSmootherLevelset demo application
 *
 */
int main(int argc, char* argv[])
{
  typedef double DataType;

  // Creata a parser for command line arguments.
  SimpleArgParser args(argc, argv);

  if( args.check("help") > -1 || args.num_args()==1)
  {
    std::cout << "Rumpf Smoother Application for Levelset usage: " << std::endl;
    std::cout << "Required arguments: --filename [String]: Path to a FEAST mesh file." << std::endl;
    std::cout << "Optional arguments: --level [unsigned int]: Number of refines, defaults to 0." << std::endl;
    exit(1);
  }
  // Specify supported command line switches
  args.support("level");
  args.support("meshfile");
  args.support("help");

  // Refinement level
  Index lvl_max(0);
  // Input mesh file name, required
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

  // Check and parse --level
  if(args.check("level") != 1)
    std::cout << "No refinement level specified, defaulting to 0." << std::endl;
  else
  {
    args.parse("level", lvl_max);
    std::cout << "Refinement level " << lvl_max << std::endl;
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

  DataType deltat(DataType(1e-2));

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> Simplex2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> Hypercube2Mesh_2d;

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  int ret(1);

  // Call the run() method of the appropriate wrapper class
  if(mtype == "conformal:simplex:2:2")
    return LevelsetApp<DataType, Simplex2Mesh_2d, MySmoother, MyFunctional, Meshopt::RumpfFunctionalLevelset>::
      run(mesh_reader, lvl_max, deltat);
  if(mtype == "conformal:hypercube:2:2")
    return LevelsetApp<DataType, Hypercube2Mesh_2d, MySmoother, MyFunctional, Meshopt::RumpfFunctionalLevelset>::
      run(mesh_reader, lvl_max, deltat);

  ret = ret | FEAST::Runtime::finalise();
  // If no MeshType from the list was in the file, return 1
  return ret;
}
/// \endcond
#endif
