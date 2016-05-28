#include <kernel/base_header.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/meshopt/hyperelasticity_functional.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1split.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAST;

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
 **/
template
<
  typename DT_,
  typename ShapeType_,
  template<typename, typename> class FunctionalType_,
  template<typename ... > class HyperelasticityFunctionalType_
  > struct ResizeApp
{
  /**
   * @brief Runs mesh smoother stuff
   *
   **/
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
    typedef typename FEAST::Meshopt::Intern::TrafoFE<TrafoType>::Space TrafoSpace;
    /// Our functional type
    typedef FunctionalType_<DataType, ShapeType> FunctionalType;
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
    // Trafo and trafo FE space
    TrafoType trafo(*mesh);
    TrafoSpace trafo_space(trafo);
    // The filters will be empty, but we still need the filter and assembler objects
    std::map<String, std::shared_ptr<Assembly::UnitFilterAssembler<MeshType>>> dirichlet_asm;
    std::map<String, std::shared_ptr<Assembly::SlipFilterAssembler<MeshType>>> slip_asm;
    FilterType my_filter;

    // Create the root mesh node
    Geometry::RootMeshNode<MeshType>* rmn(new Geometry::RootMeshNode<MeshType>(mesh, nullptr));

    // Parameters for the Rumpf functional
    DataType fac_norm = DataType(1e0),fac_det = DataType(1),fac_cof = DataType(0), fac_reg(DataType(0e0));
    // Create the functional with these parameters
    auto my_functional = std::make_shared<FunctionalType>(fac_norm, fac_det, fac_cof, fac_reg);

    // Create the mesh quality functional
    HyperelasticityFunctionalType rumpflpumpfl(rmn, trafo_space, dirichlet_asm, slip_asm, my_functional);
    // Print information
    rumpflpumpfl.print();

    // Vector for setting coordinates
    auto new_coords = rumpflpumpfl.create_vector_r();

    // First we scale the reference element and call compute_h(). This sets the optimal scales. Then we rescale the
    // cell again WITHOUT calling compute_h(): this sets our initial state from which we start the optimisation.
    // After optimisation, the original optimal scale should have been recovered.

    // Set optimal scale
    DataType target_scaling(DataType(2.5));
    helperclass<ShapeType>::set_coords(new_coords, target_scaling);
    // init() sets the coordinates in the mesh and computes h
    rumpflpumpfl.init();
    // Now set
    rumpflpumpfl.prepare(new_coords, my_filter);
    rumpflpumpfl.compute_h();

    // Transform the cell to the initial state
    DataType scaling(DataType(5.5));
    helperclass<ShapeType>::set_coords(new_coords, scaling);
    rumpflpumpfl.prepare(new_coords, my_filter);

    // Arrays for saving the contributions of the different Rumpf functional parts
    DataType* func_norm(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_rec_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);

    // Dummy vector for rhs
    auto rhs = rumpflpumpfl.create_vector_r();
    // Vector to save the gradient in for visualisation
    auto grad = rumpflpumpfl.create_vector_r();

    // Compute initial functional value
    DataType fval(0);
    fval = rumpflpumpfl.compute_func(func_norm, func_det, func_rec_det);
    std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << std::endl;

    // Compute initial functional gradient
    rumpflpumpfl.compute_grad(grad);

    std::string filename;
    // Write initial state to file
    filename = "pre_" + helperclass<ShapeType>::print_typename();
    Geometry::ExportVTK<MeshType> writer_initial_pre(*mesh);
    writer_initial_pre.add_cell_vector("h", rumpflpumpfl._h);
    writer_initial_pre.add_cell_vector("fval", func_norm, func_det, func_rec_det);
    writer_initial_pre.add_vertex_vector("grad", grad);
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
    solver->correct(new_coords, rhs);
    solver->done();

    fval = rumpflpumpfl.compute_func(func_norm, func_det, func_rec_det);
    std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << std::endl;

    rumpflpumpfl.compute_grad(grad);

    // Write optimised initial mesh
    filename = "post_" + helperclass<ShapeType>::print_typename();
    Geometry::ExportVTK<MeshType> writer_initial_post(*mesh);
    writer_initial_post.add_cell_vector("h", rumpflpumpfl._h);
    writer_initial_post.add_cell_vector("fval", func_norm, func_det, func_rec_det);
    writer_initial_post.add_vertex_vector("grad", grad);
    writer_initial_post.write(filename);

    // Clean up
    delete rmn;
    delete[] func_norm;
    delete[] func_det;
    delete[] func_rec_det;

  }

};

// Template aliases to easier switch between variants

template<typename A, typename B>
using MyLocalFunctional = Meshopt::RumpfFunctional<A, B>;

// For using the Q1 split functional, the functional is a bit more complicated
template<typename A, typename B>
using MyLocalFunctionalQ1Split = Meshopt::RumpfFunctionalQ1Split<A, B, Meshopt::RumpfFunctional_D2>;

// Vanilla Rumpf smoother
template<typename A, typename B, typename C, typename D, typename E>
using MyQualityFunctional = Meshopt::HyperelasticityFunctional<A, B, C, D, E>;

int main(int argc, char** argv)
{
  FEAST::Runtime::initialise(argc, argv);

  ResizeApp<double, Shape::Hypercube<2>, MyLocalFunctional, MyQualityFunctional>::run();

  return FEAST::Runtime::finalise();
}

/// \brief Specialisation for hypercubes
template<int shape_dim_>
struct helperclass< FEAST::Shape::Hypercube<shape_dim_> >
{
  /// \brief Sets coordinates so we deal the the reference element
  template<typename VectorType_, typename DataType_>
  static void set_coords(VectorType_& coords_, const DataType_& scaling)
  {
    for(Index i(0); i < Index(1 << shape_dim_); ++i)
    {
      Tiny::Vector<DataType_, VectorType_::BlockSize, VectorType_::BlockSize> tmp;
      for(int d(0); d < shape_dim_; ++d)
        tmp(d) = (DataType_(((i >> d) & 1) << 1) - DataType_(1)) * scaling ;

      coords_(i, tmp);
    }
  }

  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Hypercube_" + stringify(shape_dim_);
  }
};

/// \brief Specialisation for 2d simplices
template<>
struct helperclass< FEAST::Shape::Simplex<2> >
{
  /// \brief Sets coordinates so we deal the the Rumpf reference element
  template<typename VectorType_, typename DataType_>
  static void set_coords(VectorType_& coords_, const DataType_& scaling)
  {
    Tiny::Vector<DataType_, 2, 2> tmp(0);
    coords_(0, tmp);

    tmp(0) = scaling;
    coords_(1, tmp);

    tmp(0) = DataType_(0.5) * scaling;
    tmp(1) = DataType_(0.5)*Math::sqrt(DataType_(3))*scaling;
    coords_(2, tmp);
  }

  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Simplex_2";
  }
};
