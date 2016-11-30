#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/meshopt/hyperelasticity_functional.hpp>

#include <kernel/meshopt/rumpf_functionals/p1.hpp>
#include <kernel/meshopt/rumpf_functionals/q1.hpp>

#include <kernel/meshopt/rumpf_functionals/2d_p1_unrolled.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_unrolled.hpp>

#include <kernel/meshopt/rumpf_functionals/3d_p1_unrolled.hpp>
// Not implemented yet due to too bloody much stuff for Maple to compute
//#include <kernel/meshopt/rumpf_functionals/3d_q1_unrolled.hpp>

#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlcg.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/// \cond internal

/// \brief Helper class for resizing tests
template<typename ShapeType_>
struct helperclass;

/// \cond internal

/**
 * \brief Test for Hyperelasticity-based mesh optimisation
 *
 * The input mesh consists of a single Rumpf reference cell of some target scaling. This is then rescaled and the
 * mesh optimiser is supposed to scale it back to the original scaling.
 *
 * If the resulting cell is optimal in the defined sense, the Frobenius norm term should be zero and the determinant
 * should be 1 (mind the scaling from fac_det etc. in the functional!).
 *
 * \author Jordi Paul
 **/
template
<
  typename DT_,
  typename ShapeType_,
  template<typename, typename> class FunctionalType_,
  template<typename ... > class MeshQualityFunctional_
  >
  class HyperelasticityFunctionalTest
  : public TestSystem::FullTaggedTest<Mem::Main, DT_, Index>
{
  public:
    typedef Mem::Main MemType;
    typedef Index IndexType;
    typedef DT_ DataType;

    typedef ShapeType_ ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    /// The FE space for the transformation
    typedef typename FEAT::Meshopt::Intern::TrafoFE<TrafoType>::Space TrafoSpace;
    /// Our functional type
    typedef FunctionalType_<DataType, TrafoType> FunctionalType;
    /// The Rumpf smoother
    typedef MeshQualityFunctional_<MemType, DataType, IndexType, TrafoType, FunctionalType> MeshQualityFunctional;
    /// Filter for Dirichlet boundary conditions
    typedef LAFEM::UnitFilterBlocked<MemType, DataType, IndexType, MeshType::world_dim> DirichletFilterType;
    /// Filter for slip boundary conditions
    typedef LAFEM::SlipFilter<MemType, DataType, IndexType, MeshType::world_dim> SlipFilterType;
    /// Combined filter
    typedef LAFEM::FilterChain<LAFEM::FilterSequence<SlipFilterType>, LAFEM::FilterSequence<DirichletFilterType>> FilterType;

  private:
    const int _exponent_det;

  public:
    explicit HyperelasticityFunctionalTest(int exponent_det) : TestSystem::FullTaggedTest<MemType, DataType, IndexType>
    ("hyperelasticity_functional_test-"+FunctionalType::name()), _exponent_det(exponent_det)
      {
      }

    virtual ~HyperelasticityFunctionalTest()
    {
    }

    virtual void run() const override
    {
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
      DataType fac_norm(1e-1);
      DataType fac_det(2.5);
      DataType fac_cof(MeshType::world_dim == 3);
      DataType fac_reg(DataType(1e-8));
      // Create the functional with these parameters
      auto my_functional = std::make_shared<FunctionalType>(fac_norm, fac_det, fac_cof, fac_reg, _exponent_det);

      // Set optimal scale
      DataType target_scaling(DataType(2.5));
      helperclass<ShapeType>::set_coords(*(rmn->get_mesh()), target_scaling);

      // Trafo and trafo FE space
      TrafoType trafo(*(rmn->get_mesh()));
      TrafoSpace trafo_space(trafo);

      // Create the mesh quality functional
      MeshQualityFunctional rumpflpumpfl(rmn, trafo_space, dirichlet_asm, slip_asm, my_functional, Meshopt::ScaleComputation::current_uniform);

      // init() sets the coordinates in the mesh and computes h
      rumpflpumpfl.init();

      // Now set the initial state
      auto& coords_buffer = rumpflpumpfl.get_coords();
      coords_buffer.scale(coords_buffer, DataType(2));
      rumpflpumpfl.buffer_to_mesh();

      // Arrays for saving the contributions of the different Rumpf functional parts
      DataType func_norm;
      DataType func_cof;
      DataType func_det;

      // Dummy vector for rhs
      auto rhs = rumpflpumpfl.create_vector_r();
      // Vector to save the gradient in for visualisation
      auto grad = rumpflpumpfl.create_vector_r();
      // Solver vector containing the initial guess
      auto new_coords = rumpflpumpfl.get_coords().clone(LAFEM::CloneMode::Deep);

      // Compute initial functional value
      DataType fval_pre(0);
      rumpflpumpfl.eval_fval_cellwise(fval_pre, &func_norm, &func_cof, &func_det);

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
      solver->set_tol_rel(Math::pow(Math::eps<DataType>(), DataType(0.9)));
      solver->correct(new_coords, rhs);
      solver->done();
      // Compute functional value post optimisation
      DataType fval_post(0);
      rumpflpumpfl.eval_fval_cellwise(fval_post, &func_norm, &func_cof, &func_det);

      const DataType eps = Math::pow(Math::eps<DataType>(),DataType(0.5));

      // Check functional value contributions
      TEST_CHECK(fval_pre > fval_post);
      // Both of these should evaluate to 0
      TEST_CHECK_EQUAL_WITHIN_EPS(func_norm, DataType(0), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(func_cof, DataType(0), eps);
      // Because of how the stability term depend on 1/det and det=1 for the identity mapping, this should be 2*fac_det
      TEST_CHECK_EQUAL_WITHIN_EPS(func_det, fac_det*DataType(2), eps);

      // Now do the negative test: Change the functional in a nonsensical manner. Calling the optimiser should NOT
      // give the correctly scaled element
      my_functional->_fac_rec_det = DataType(0.6676);

      // Compute initial functional value
      rumpflpumpfl.eval_fval_cellwise(fval_pre, &func_norm, &func_cof, &func_det);

      // Optimise again
      rumpflpumpfl.init();
      rumpflpumpfl.prepare(new_coords, my_filter);

      solver->init();
      solver->set_plot(true);
      solver->set_tol_rel(Math::eps<DataType>());
      solver->correct(new_coords, rhs);
      solver->done();

      // Compute new functional value
      rumpflpumpfl.eval_fval_cellwise(fval_post, &func_norm, &func_cof, &func_det);

      // With the new functional, the functional value should still have decreased
      TEST_CHECK(fval_pre > fval_post);
      // These differences should all be greater than eps
      TEST_CHECK(Math::abs(func_norm - DataType(0)) > eps);
      TEST_CHECK(Math::abs(func_det - my_functional->_fac_det*DataType(1)) > eps);

      delete rmn;

    }
};

// Use template alias to make everything more readable
template<typename A, typename B, typename C, typename D, typename E>
using MyQualityFunctional = Meshopt::HyperelasticityFunctional<A, B, C, D, E>;

HyperelasticityFunctionalTest<double, Shape::Hypercube<2>, Meshopt::RumpfFunctional, MyQualityFunctional> test_hc_1(1);
HyperelasticityFunctionalTest<double, Shape::Hypercube<2>, Meshopt::RumpfFunctionalUnrolled, MyQualityFunctional> test_hc_1_u(1);

HyperelasticityFunctionalTest<double, Shape::Hypercube<2>, Meshopt::RumpfFunctional, MyQualityFunctional> test_hc_2(2);
HyperelasticityFunctionalTest<double, Shape::Hypercube<2>, Meshopt::RumpfFunctionalUnrolled, MyQualityFunctional> test_hc_2_u(2);

HyperelasticityFunctionalTest<double, Shape::Hypercube<3>, Meshopt::RumpfFunctional, MyQualityFunctional> test_hc3_1(1);
HyperelasticityFunctionalTest<double, Shape::Hypercube<3>, Meshopt::RumpfFunctional, MyQualityFunctional> test_hc3_2(2);

HyperelasticityFunctionalTest<double, Shape::Simplex<2>, Meshopt::RumpfFunctional, MyQualityFunctional> test_s_1(1);
HyperelasticityFunctionalTest<double, Shape::Simplex<2>, Meshopt::RumpfFunctionalUnrolled, MyQualityFunctional> test_s_1_u(1);

HyperelasticityFunctionalTest<double, Shape::Simplex<2>, Meshopt::RumpfFunctional, MyQualityFunctional> test_s_2(2);
HyperelasticityFunctionalTest<double, Shape::Simplex<2>, Meshopt::RumpfFunctionalUnrolled, MyQualityFunctional> test_s_s_2_u(2);

HyperelasticityFunctionalTest<double, Shape::Simplex<3>, Meshopt::RumpfFunctional, MyQualityFunctional> test_s3_1(1);
HyperelasticityFunctionalTest<double, Shape::Simplex<3>, Meshopt::RumpfFunctionalUnrolled, MyQualityFunctional> test_s3_1_u(1);
HyperelasticityFunctionalTest<double, Shape::Simplex<3>, Meshopt::RumpfFunctional, MyQualityFunctional> test_s3_2(2);
HyperelasticityFunctionalTest<double, Shape::Simplex<3>, Meshopt::RumpfFunctionalUnrolled, MyQualityFunctional> test_s3_2_u(2);

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

    for(Index i(0); i < Index(1 << shape_dim); ++i)
    {
      for(int d(0); d < shape_dim; ++d)
        tmp(d) = (DT_(((i >> d) & 1) << 1) - DT_(1)) * scaling ;

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
    tmp(1) = DT_(0);
    tmp(2) = DT_(0);
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
