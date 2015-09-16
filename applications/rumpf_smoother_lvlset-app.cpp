#include <kernel/base_header.hpp>
#ifdef FEAST_HAVE_ALGLIB
#include <kernel/assembly/common_functions.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
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

using namespace FEAST;

namespace FEAST
{
  namespace Assembly
  {
    namespace Common
    {
      /**
       * @brief Class Template for the gradient of DistanceFunctionSD
       *
       * \tparam ImgPointType_
       * Type for the point in the DistanceFunction
       *
       * \tparam component
       * Component of the gradient
       *
       **/
      template<typename ImgPointType_, int component>
      class DistanceFunctionSD_grad :
        public AnalyticFunction
      {
        public:
          /// Datatype from the point
          typedef typename ImgPointType_::DataType DataType;
          /** \copydoc AnalyticFunction::FunctionCapabilities */
          enum FunctionCapabilities
          {
            can_value = 1,
            can_grad = 0,
            can_hess = 0
          };

          /** \copydoc AnalyticFunction::ConfigTraits */
          template<typename Config_>
          struct ConfigTraits
          {
            /**
             * \brief Trafo configuration tag class
             *
             * \see Trafo::ConfigBase
             */
            struct TrafoConfig :
              public Trafo::ConfigBase
            {
              enum
              {
                need_img_point = 1
              };
            };
          };

          /** \copydoc AnalyticFunction::Evaluator */
          template<typename EvalTraits_>
          class Evaluator :
            public AnalyticFunction::Evaluator<EvalTraits_>
        {
          public:
            /// trafo evaluator data
            typedef typename EvalTraits_::TrafoEvaluator TrafoEvaluator;
            /// trafo data type
            typedef typename EvalTraits_::TrafoData TrafoData;
            /// coefficient data type
            typedef typename EvalTraits_::DataType DataType;
            /// value type
            typedef typename EvalTraits_::ValueType ValueType;
            /// gradient type
            typedef typename EvalTraits_::GradientType GradientType;
            /// hessian type
            typedef typename EvalTraits_::HessianType HessianType;
            /// type for the points the analytic function is evaluated at
            typedef typename EvalTraits_::ImagePointType ImgPointType;

          private:
            /// Function to evaluate
            const DistanceFunctionSD_grad& _function;

          public:
            /// Constructor
            explicit Evaluator(const DistanceFunctionSD_grad& function) :
              _function(function)
              {
              }

            ValueType value(const TrafoData& tau) const
            {
              ImgPointType tmp (tau.img_point - _function._point);
              DataType norm = tmp.norm_euclid();
              if(norm <= Math::eps<DataType>())
              {
                return DataType(0);
              }
              else
              {
                return _function._b*( tau.img_point[component] - _function._point(component))/norm;
              }
            }
        }; // class DistanceFunctionSD::Evaluator<...>

        private:
          /// The point to which the distance to is calculated
          ImgPointType_& _point;
          /// Displacement of the function
          DataType _a;
          /// Scaling factor
          DataType _b;

        public:
          /// Constructor
          explicit DistanceFunctionSD_grad(ImgPointType_& x0_, const DataType a_, const DataType b_) :
            _point(x0_),
            _a(a_),
            _b(b_)
            {
            }

      }; // class DistanceFunctionSD_grad

      template<int component, int derivative_component, typename ImgPointType_>
      class PlaneDistanceFunctionSD_grad :
        public AnalyticFunction
      {
        public:
          /// Datatype from the point
          typedef typename ImgPointType_::DataType DataType;
          /** \copydoc AnalyticFunction::FunctionCapabilities */
          enum FunctionCapabilities
          {
            can_value = 1,
            can_grad = 0,
            can_hess = 0
          };

          /** \copydoc AnalyticFunction::ConfigTraits */
          template<typename Config_>
          struct ConfigTraits
          {
            /**
             * \brief Trafo configuration tag class
             *
             * \see Trafo::ConfigBase
             */
            struct TrafoConfig :
              public Trafo::ConfigBase
            {
              enum
              {
                need_img_point = 1
              };
            };
          };

          /** \copydoc AnalyticFunction::Evaluator */
          template<typename EvalTraits_>
          class Evaluator :
            public AnalyticFunction::Evaluator<EvalTraits_>
        {
          public:
            /// trafo evaluator data
            typedef typename EvalTraits_::TrafoEvaluator TrafoEvaluator;
            /// trafo data type
            typedef typename EvalTraits_::TrafoData TrafoData;
            /// coefficient data type
            typedef typename EvalTraits_::DataType DataType;
            /// value type
            typedef typename EvalTraits_::ValueType ValueType;
            /// gradient type
            typedef typename EvalTraits_::GradientType GradientType;
            /// hessian type
            typedef typename EvalTraits_::HessianType HessianType;
            /// type for the points the analytic function is evaluated at
            typedef typename EvalTraits_::ImagePointType ImgPointType;

          private:
            /// Function to evaluate
            const PlaneDistanceFunctionSD_grad& _function;

          public:
            /// Constructor
            explicit Evaluator(const PlaneDistanceFunctionSD_grad& function) :
              _function(function)
              {
              }

            ValueType value(const TrafoData& DOXY(tau)) const
            {
              return _function._b*ValueType(component == derivative_component);
            }
        }; // class PlaneDistanceFunctionSD::Evaluator<...>

        private:
          /// The point to which the distance to is calculated
          ImgPointType_& _point;
          /// Scaling factor
          DataType _b;

        public:
          /// Constructor
          explicit PlaneDistanceFunctionSD_grad(ImgPointType_& x0_, const DataType b_) :
            _point(x0_),
            _b(b_)
            {
            }

      }; // class DistanceFunctionSD_grad

    } // namespace Common
  } // namespace Assembly
} // namespace FEAST

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
  static int run(MeshStreamer& my_streamer, Index level, DT_ deltat)
  {
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
    for(Index lvl(1); lvl <= level; ++lvl)
    {
      std::cout << "Refining up to level " << lvl << "..." << std::endl;
      auto* old = rmn;
      rmn = old->refine();
      delete old;
    }

    MeshType* mesh = rmn->get_mesh();

    typedef Assembly::Common::DistanceFunctionSD<ImgPointType> AnalyticFunctionType;
    typedef Assembly::Common::DistanceFunctionSD_grad<ImgPointType, 0> AnalyticFunctionGrad0Type;
    typedef Assembly::Common::DistanceFunctionSD_grad<ImgPointType, 1> AnalyticFunctionGrad1Type;
    ImgPointType x0(DataType(0));
    x0.v[0] = DataType(0.25) *(DataType(2) + Math::cos(DataType(0)));
    x0.v[1] = DataType(0.25) *(DataType(2) + Math::sin(DataType(0)));
    DataType displacement(0.15);
    DataType scaling(-1);
    AnalyticFunctionType analytic_lvlset(x0, displacement, scaling);
    AnalyticFunctionGrad0Type analytic_lvlset_grad0(x0, displacement, scaling);
    AnalyticFunctionGrad1Type analytic_lvlset_grad1(x0, displacement, scaling);

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
      AnalyticFunctionGrad0Type,
      AnalyticFunctionGrad1Type,
      TrafoType,
      FunctionalType,
      LevelsetFunctionalType
    > RumpfSmootherType;

    // Parameters for the funcitional
    DataType fac_norm = DataType(1e0),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
    // Create the functional
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);
    // Parameters for the Rumpf smoother
    bool align_to_lvlset(false);
    bool r_adaptivity(true);
    DataType r_adapt_reg = DataType(1e-2), r_adapt_pow = DataType(0.5);
    // Create levelset part of the functional
    LevelsetFunctionalType my_levelset_functional;

    std::deque<String> slip_list;
    slip_list.push_back("left");
    slip_list.push_back("right");

    std::deque<String> dirichlet_list;
    dirichlet_list.push_back("bottom");
    dirichlet_list.push_back("top");

    // The smoother in all its template glory
    RumpfSmootherType rumpflpumpfl(rmn, slip_list, dirichlet_list, my_functional, my_levelset_functional,
    align_to_lvlset, r_adaptivity, analytic_lvlset, analytic_lvlset_grad0, analytic_lvlset_grad1);
    rumpflpumpfl.init();

    //rumpflpumpfl._update_h = true;

    // Set the r-adaptivity parameters
    rumpflpumpfl.set_r_adapt_params(r_adapt_reg,r_adapt_pow);

    // Print lotsa information
    std::cout << __func__ << " at refinement level " << level << std::endl;
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
    std::cout << "fval pre optimisation = " << scientify(fval) << " cell size quality indicator: " << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;
    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    // Filename for vtk files
    std::string filename;

    // Write initial state to file
    Geometry::ExportVTK<MeshType> writer_pre_initial(*mesh);
    writer_pre_initial.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_pre_initial.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
    writer_pre_initial.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_pre_initial.add_field_vertex("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements(), rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
    writer_pre_initial.add_field_cell("fval", func_norm, func_det, func_rec_det);
    writer_pre_initial.add_scalar_cell("levelset_constraint", func_lvlset );
    writer_pre_initial.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_pre_initial.write("pre_initial");

    // Optimise the mesh
    rumpflpumpfl.optimise();

    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
    std::cout << "fval post optimisation = " << scientify(fval) << " cell size quality indicator: " << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;

    // Call prepare() again because the mesh changed due to the optimisation and it was not called again after the
    // last iteration
    rumpflpumpfl.prepare();
    rumpflpumpfl.compute_gradient();

    // Write optimised initial mesh
    Geometry::ExportVTK<MeshType> writer_post_initial(*mesh);
    writer_post_initial.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_post_initial.add_field_cell_blocked_vector("h", rumpflpumpfl._h );
    writer_post_initial.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_post_initial.add_field_vertex("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements(), rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
    writer_post_initial.add_field_cell("fval", func_norm, func_det, func_rec_det);
    writer_post_initial.add_scalar_cell("levelset_constraint", func_lvlset );
    writer_post_initial.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_post_initial.write("post_initial");

    DataType time(0);
    Index n(0);

    // Old mesh coordinates for computing the mesh velocity
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim> coords_old(mesh->get_num_entities(0),DataType(0));
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim> mesh_velocity(mesh->get_num_entities(0), DataType(0));

    std::cout << "deltat = " << scientify(deltat) << std::endl;
    while(time < DataType(deltat))
    {
      std::cout << "timestep " << n << std::endl;
      time+= deltat;

      // Save old vertex coordinates
      coords_old.clone(rumpflpumpfl._coords);

      // update levelset function
      x0.v[0] = DataType(0.25) *(DataType(2) + Math::cos(time));
      x0.v[1] = DataType(0.25) *(DataType(2) + Math::sin(DataType(3)*time));

      // Update the reference point of the distance function
      analytic_lvlset.set_point(x0);

      // Compute functional value and gradient
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval pre optimisation = " << scientify(fval) << " cell size quality indicator: " <<
        scientify(rumpflpumpfl.cell_size_quality()) << std::endl;

      // Write pre-optimisation mesh
      filename = "pre_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
      writer_pre.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
      writer_pre.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
      writer_pre.add_field_vertex("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements(),
      rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
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
      std::cout << "max mesh velocity = " << scientify(max_mesh_velocity) << std::endl;

      // Write post-optimisation mesh
      fval = rumpflpumpfl.compute_functional(func_norm,func_det,func_rec_det,func_lvlset);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval post optimisation = " << scientify(fval) << " cell size quality indicator: "
      << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;

      filename = "post_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_scalar_cell("norm", func_norm);
      writer_post.add_scalar_cell("det", func_det);
      writer_post.add_scalar_cell("rec_det", func_rec_det);
      writer_post.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
      writer_post.add_field_cell_blocked_vector("h", rumpflpumpfl._h );
      writer_post.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
      writer_post.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
      writer_post.add_field_vertex("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements(),
      rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
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

template<typename A, typename B, typename C, typename D, typename E, typename F>
using MySmoother = Meshopt::RumpfSmootherLevelsetConcAnalytic<A, B, C, D, E, F>;

template<typename A, typename B, typename C, typename D, typename E, typename F>
using MySmootherQ1Hack = Meshopt::RumpfSmootherLevelsetAnalyticQ1Hack<A, B, C, D, E, F>;

template<typename A, typename B>
using MyFunctional= Meshopt::RumpfFunctionalConc_D2<A, B>;

template<typename A, typename B>
using MyFunctionalQ1Hack = Meshopt::RumpfFunctionalQ1Hack<A, B, Meshopt::RumpfFunctionalConc>;

/**
 * \cond internal
 *
 * RumpfSmootherLevelset demo application
 *
 */
int main(int argc, char* argv[])
{
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

  typedef double DataType;

  DataType deltat(DataType(1e-2));

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> Simplex2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> Hypercube2Mesh_2d;

  // Call the run() method of the appropriate wrapper class
  if(shape_type == mesh_data.st_tria)
    return LevelsetApp<DataType, Simplex2Mesh_2d, MySmoother, MyFunctional, Meshopt::RumpfFunctionalLevelset>::
      run(my_streamer, lvl_max, deltat);
  if(shape_type == mesh_data.st_quad)
    return LevelsetApp<DataType, Hypercube2Mesh_2d, MySmoother, MyFunctional, Meshopt::RumpfFunctionalLevelset>::
      run(my_streamer, lvl_max, deltat);

  // If no MeshType from the list was in the file, return 1
  return 1;
}
/// \endcond
#endif
