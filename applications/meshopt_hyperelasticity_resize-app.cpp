#include <kernel/base_header.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/meshopt/hyperelasticity_functional.hpp>
#include <kernel/meshopt/rumpf_functionals/q1.hpp>
#include <kernel/meshopt/rumpf_functionals/p1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_unrolled.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_unrolled.hpp>
#include <kernel/meshopt/rumpf_functionals/3d_p1_unrolled.hpp>
//#include <kernel/meshopt/rumpf_functionals/3d_q1_unrolled.hpp>
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAT;

/// \cond internal
// Helper class for resizing tests
template<typename ShapeType_>
struct helperclass;

/// \endcond

/**
 * \brief This application demonstrates the usage of some of the HyperelasticityFunctional classes to resize reference cells
 *
 * This is excellent for error checking in smoothers and functionals, as a correct implementation should resize
 * a Rumpf reference cell to the desired volume and achieve a global minimum of the functional value.
 *
 * \note Because an application of the (nonlinear) Rumpf smoother requires operations similar to a matrix assembly,
 * Rumpf smoothers are implemented for Mem::Main only.
 *
 * \author Jordi Paul
 *
 * \tparam DT_
 * The precision of the mesh etc.
 *
 * \tparam ShapeType
 * The shape of the mesh's cells
 *
 * \tparam FunctionalType
 * The Rumpf functional variant to use
 *
 * \tparam HyperelasticityFunctionalType_
 * The Rumpf smoother variant to use
 *
 */
template
<
  typename DT_,
  typename ShapeType_,
  template<typename, typename> class FunctionalType_,
  template<typename ... > class HyperelasticityFunctionalType_
  > struct ResizeApp
{
  /**
   * \brief Runs mesh optimiser stuff
   *
   */
  static void run()
  {
    /// Precision for meshes etc, everything else uses the same data type
    typedef DT_ DataType;
    // Rumpf Smoothers are implemented for Mem::Main only
    typedef Mem::Main MemType;
    // So we use Index
    typedef Index IndexType;
    /// Shape of the mesh cells
    typedef ShapeType_ ShapeType;
    /// The complete mesh type
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension,ShapeType::dimension, DataType> MeshType;
    /// The corresponding transformation
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    /// The FE space for the transformation
    typedef typename FEAT::Meshopt::Intern::TrafoFE<TrafoType>::Space TrafoSpace;
    /// Our functional type
    typedef FunctionalType_<DataType, TrafoType> FunctionalType;
    /// The Rumpf smoother
    typedef HyperelasticityFunctionalType_<MemType, DataType, IndexType, TrafoType, FunctionalType> HyperelasticityFunctionalType;
    /// Filter for Dirichlet boundary conditions
    typedef LAFEM::UnitFilterBlocked<MemType, DataType, IndexType, MeshType::world_dim> DirichletFilterType;
    /// Filter for slip boundary conditions
    typedef LAFEM::SlipFilter<MemType, DataType, IndexType, MeshType::world_dim> SlipFilterType;
    /// Combined filter
    typedef LAFEM::FilterChain
    <
      LAFEM::FilterSequence<SlipFilterType>,
      LAFEM::FilterSequence<DirichletFilterType>
    > FilterType;


    // Create a single reference cell of the shape type
    Geometry::ReferenceCellFactory<ShapeType, DataType> mesh_factory;
    // Create the mesh
    MeshType* mesh(new MeshType(mesh_factory));
    // The filters will be empty, but we still need the filter and assembler objects
    std::map<String, std::shared_ptr<Assembly::UnitFilterAssembler<MeshType>>> dirichlet_asm;
    std::map<String, std::shared_ptr<Assembly::SlipFilterAssembler<MeshType>>> slip_asm;
    FilterType my_filter;

    // Create the root mesh node
    Geometry::RootMeshNode<MeshType>* rmn(new Geometry::RootMeshNode<MeshType>(mesh, nullptr));

    // Parameters for the Rumpf functional
    DataType fac_norm(2);
    DataType fac_det(2.5);
    DataType fac_cof(3*(MeshType::world_dim == 3));
    DataType fac_reg(DataType(1e-8));

    // Create the functional with these parameters
    auto my_functional = std::make_shared<FunctionalType>(fac_norm, fac_det, fac_cof, fac_reg, 1);

    // Set optimal scale
    DataType target_scaling(DataType(2.5));
    helperclass<ShapeType>::set_coords(*(rmn->get_mesh()), target_scaling);

    // For checking if the method generates a mesh only consisting of optimal cells, the mesh can be refined here.
    // Note that this will not work for Simplex<3> meshes, as a refinement of the optimal cell creates cells of
    // different shapes
    int lvl_max(0);
    for(int lvl(0); lvl < lvl_max; ++lvl)
    {
      Geometry::RootMeshNode<MeshType>* coarse_node = rmn;
      rmn = coarse_node->refine(Geometry::AdaptMode::none);
      delete coarse_node;
    }

    // Trafo and trafo FE space
    TrafoType trafo(*(rmn->get_mesh()));
    TrafoSpace trafo_space(trafo);

    // Create the mesh quality functional
    HyperelasticityFunctionalType rumpflpumpfl(rmn, trafo_space, dirichlet_asm, slip_asm, my_functional, Meshopt::ScaleComputation::current_uniform);
    // Print information
    rumpflpumpfl.print();

    // init() sets the coordinates in the mesh and computes h
    rumpflpumpfl.init();

    // Now set the initial state
    auto& coords_buffer = rumpflpumpfl.get_coords();
    coords_buffer.scale(coords_buffer, DataType(6)/DataType(2.5));
    rumpflpumpfl.buffer_to_mesh();

    // Arrays for saving the contributions of the different Rumpf functional parts
    DataType* fval_norm(new DataType[rmn->get_mesh()->get_num_entities(MeshType::shape_dim)]);
    DataType* fval_det(new DataType[rmn->get_mesh()->get_num_entities(MeshType::shape_dim)]);
    DataType* fval_rec_det(new DataType[rmn->get_mesh()->get_num_entities(MeshType::shape_dim)]);

    // Dummy vector for rhs
    auto rhs = rumpflpumpfl.create_vector_r();
    // Vector to save the gradient in for visualisation
    auto grad = rumpflpumpfl.create_vector_r();
    // Solver vector containing the initial guess
    auto new_coords = rumpflpumpfl.get_coords().clone(LAFEM::CloneMode::Deep);

    // Compute initial functional value
    DataType fval(0);
    rumpflpumpfl.eval_fval_cellwise(fval, fval_norm, fval_det, fval_rec_det);
    std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << std::endl;

    std::string filename;
    // Write initial state to file
    filename = "pre_" + helperclass<ShapeType>::print_typename();
    Geometry::ExportVTK<MeshType> writer_initial_pre(*(rmn->get_mesh()));
    rumpflpumpfl.add_to_vtk_exporter(writer_initial_pre);
    writer_initial_pre.write(filename);

    // Create a solver
    auto linesearch = Solver::new_strong_wolfe_linesearch(rumpflpumpfl, my_filter);

    auto solver = Solver::new_nlcg(
      rumpflpumpfl, // operator
      my_filter, // filter
      linesearch, // linesearch
      Solver::NLCGDirectionUpdate::DYHSHybrid, // search direction update
      false /*, // do not keep iterates
      nullptr*/); // no preconditioner

    solver->init();
    solver->set_plot(true);
    solver->set_tol_rel(Math::pow(Math::eps<DataType>(),DataType(0.9)));
    solver->set_max_iter(30);
    solver->correct(new_coords, rhs);
    solver->done();
    std::cout << "Solver used: " << FEAT::Statistics::get_formatted_solver_tree().trim() <<std::endl;

    rumpflpumpfl.eval_fval_cellwise(fval, fval_norm, fval_det, fval_rec_det);
    std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << std::endl;

    // Write optimised initial mesh
    filename = "post_" + helperclass<ShapeType>::print_typename();
    Geometry::ExportVTK<MeshType> writer_initial_post(*(rmn->get_mesh()));
    rumpflpumpfl.add_to_vtk_exporter(writer_initial_post);
    writer_initial_post.write(filename);

    // Clean up
    delete rmn;
    delete[] fval_norm;
    delete[] fval_det;
    delete[] fval_rec_det;

  }

};

// Template aliases to easier switch between variants

template<typename A, typename B>
using MyLocalFunctional = Meshopt::RumpfFunctionalUnrolled<A, B>;

// For using the Q1 split functional, the functional is a bit more complicated
template<typename A, typename B>
using MyLocalFunctionalQ1Split = Meshopt::RumpfFunctionalQ1Split<A, B, Meshopt::RumpfFunctional>;

// Vanilla Rumpf smoother
template<typename A, typename B, typename C, typename D, typename E>
using MyQualityFunctional = Meshopt::HyperelasticityFunctional<A, B, C, D, E>;

int main(int argc, char** argv)
{
  FEAT::Runtime::initialise(argc, argv);

  ResizeApp<double, Shape::Hypercube<2>, MyLocalFunctional, MyQualityFunctional>::run();

  return FEAT::Runtime::finalise();
}

/// \brief Specialisation for hypercubes
template<int shape_dim>
struct helperclass< FEAT::Shape::Hypercube<shape_dim> >
{
  typedef FEAT::Shape::Hypercube<shape_dim> ShapeType;
  /// \brief Sets coordinates so we deal the the reference element
  template<typename DT_>
  static void set_coords(Geometry::ConformalMesh<ShapeType, shape_dim, shape_dim, DT_>& mesh_, const DT_& scaling)
  {
    auto& coords = mesh_.get_vertex_set();
    Tiny::Vector<DT_, shape_dim > tmp;

    Tiny::Matrix<DT_,shape_dim, shape_dim> trafo(DT_(0));

    trafo(0,0) = DT_(2);
    trafo(0,1) = DT_(1);

    trafo(1,1) = DT_(1.5);
    trafo(1,0) = DT_(0.5);

    for(Index i(0); i < Index(1 << shape_dim); ++i)
    {
      for(int d(0); d < shape_dim; ++d)
      {
        tmp(d) = (DT_(((i >> d) & 1) << 1) - DT_(1)) * scaling;
      }

      //coords[i] = trafo*tmp;
      coords[i] = tmp;
    }


  }

  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Hypercube_" + stringify(shape_dim);
  }
};

/// \brief Specialisation for 2d simplices
template<>
struct helperclass< FEAT::Shape::Simplex<2> >
{
  /// \brief Sets coordinates so we deal the the Rumpf reference element
  template<typename DT_>
  static void set_coords(Geometry::ConformalMesh<FEAT::Shape::Simplex<2>, 2, 2, DT_>& mesh_, const DT_& scaling)
  {
    auto& coords = mesh_.get_vertex_set();
    Tiny::Vector<DT_, 2> tmp(0);
    coords[0] = tmp;

    tmp(0) = scaling;
    coords[1] = tmp;

    tmp(0) = DT_(0.5) * scaling;
    //tmp(0) = DT_(0.75) * scaling;
    tmp(1) = DT_(0.5)*Math::sqrt(DT_(3))*scaling;
    coords[2] = tmp;
  }

  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Simplex_2";
  }
};

/// \brief Specialisation for 3d simplices
template<>
struct helperclass< FEAT::Shape::Simplex<3> >
{
  /// \brief Sets coordinates so we deal the the Rumpf reference element
  template<typename DT_>
  static void set_coords(Geometry::ConformalMesh<FEAT::Shape::Simplex<3>,3,3, DT_>& mesh_, const DT_& scaling)
  {

    auto& coords = mesh_.get_vertex_set();
    Tiny::Vector<DT_, 3> tmp(0);
    coords[0] = tmp;

    tmp(0) = scaling;
    tmp(1) = DT_(0.1);
    tmp(2) = DT_(0.2);
    coords[1] = tmp;

    tmp(0) = DT_(0.5) * scaling;
    tmp(1) = DT_(0.5)*Math::sqrt(DT_(3))*scaling;
    tmp(2) = DT_(0);
    coords[2] = tmp;

    tmp(0) = DT_(0.5) * scaling;
    tmp(1) = Math::sqrt(DT_(3))/DT_(6)*scaling;
    tmp(2) = Math::sqrt(DT_(6))/DT_(3)*scaling;
    coords[3] = tmp;
  }

  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Simplex_3";
  }
};
