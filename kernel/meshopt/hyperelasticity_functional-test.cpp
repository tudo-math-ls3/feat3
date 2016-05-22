#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/meshopt/hyperelasticity_functional.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1split.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>

#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlcg.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/// \cond internal

/// \brief Helper class for resizing tests
template<typename ShapeType_>
struct helperclass;

/// \cond internal

/**
 * \brief Test for Rumpf smoothers and functionals
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
  class HyperelasticityFunctionalTest_2d
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
    typedef typename FEAST::Meshopt::Intern::TrafoFE<TrafoType>::Space TrafoSpace;
    /// Our functional type
    typedef FunctionalType_<DataType, ShapeType> FunctionalType;
    /// The Rumpf smoother
    typedef MeshQualityFunctional_<MemType, DataType, IndexType, TrafoType, FunctionalType> MeshQualityFunctional;
    /// Filter for Dirichlet boundary conditions
    typedef LAFEM::UnitFilterBlocked<MemType, DataType, IndexType, MeshType::world_dim> DirichletFilterType;
    /// Filter for slip boundary conditions
    typedef LAFEM::SlipFilter<MemType, DataType, IndexType, MeshType::world_dim> SlipFilterType;
    /// Combined filter
    typedef LAFEM::FilterChain<LAFEM::FilterSequence<SlipFilterType>, LAFEM::FilterSequence<DirichletFilterType>> FilterType;


    HyperelasticityFunctionalTest_2d() :
      TestSystem::FullTaggedTest<MemType, DataType, IndexType>("hyperelasticity_functional_test")
      {
      }

    virtual void run() const override
    {

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
      DataType fac_norm = DataType(1e0),fac_det = DataType(1), fac_cof = DataType(0), fac_reg(DataType(0e0));
      // Create the functional with these parameters
      auto my_functional = std::make_shared<FunctionalType>(fac_norm, fac_det, fac_cof, fac_reg);

      // Create the smoother
      MeshQualityFunctional rumpflpumpfl(rmn, trafo_space, dirichlet_asm, slip_asm, my_functional);

      // Vector for setting coordinates
      auto new_coords = rumpflpumpfl.create_vector_r();

      // First we scale the reference element and call compute_h(). This sets the optimal scales. Then we rescale the
      // cell again WITHOUT calling compute_h(): this sets our initial state from which we start the optimisation.
      // After optimisation, the original optimal scale should have been recovered.

      // Set optimal scale
      DataType target_scaling(DataType(5.5));
      helperclass<ShapeType>::buffer_to_mesh(new_coords, target_scaling);
      // init() sets the coordinates in the mesh and computes h
      rumpflpumpfl.init();
      // Now set
      rumpflpumpfl.prepare(new_coords, my_filter);
      rumpflpumpfl.compute_h();

      // Transform the cell to the initial state
      DataType scaling(DataType(1.25));
      helperclass<ShapeType>::buffer_to_mesh(new_coords, scaling);
      rumpflpumpfl.prepare(new_coords, my_filter);

      // For saving parts of the functional value
      DataType func_norm;
      DataType func_det;
      DataType func_rec_det;

      DataType fval_pre = rumpflpumpfl.compute_func();

      // Dummy vector for rhs
      auto rhs = rumpflpumpfl.create_vector_r();
      // Create a solver
      auto linesearch = Solver::new_strong_wolfe_linesearch(rumpflpumpfl, my_filter);

      auto solver = Solver::new_nlcg(
        rumpflpumpfl, // operator
        my_filter, // filter
        linesearch, // linesearch
        Solver::NLCGDirectionUpdate::DYHSHybrid, // search direction update
        false); // do not keep iterates
        //nullptr); // no preconditioner

      solver->init();
      // Set very low relative tolerance, it's just one cell
      solver->set_tol_rel(Math::eps<DataType>());
      solver->set_plot(false);
      solver->correct(new_coords, rhs);

      // Compute functional value post optimisation
      DataType fval_post = rumpflpumpfl.compute_func(&func_norm, &func_det, &func_rec_det);

      const DataType eps = Math::pow(Math::eps<DataType>(),DataType(0.5));

      // Only check func_norm and func_det. Because of the different factors fac_rec_det depending on the
      // functionals, func_rec_det is not the same in every case. If func_det==1, we have the correct volume anyway.
      TEST_CHECK(fval_pre > fval_post);
      TEST_CHECK_EQUAL_WITHIN_EPS(func_norm, DataType(0), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(func_det, my_functional->_fac_det*DataType(1), eps);

      // Now do the negative test: Change the functional in a nonsensical manner. Calling the optimiser should NOT
      // give the correctly scaled element
      my_functional->_fac_rec_det = DataType(0.6676);

      // Compute initial functional value
      fval_pre = rumpflpumpfl.compute_func();

      // Optimise again
      rumpflpumpfl.prepare(new_coords, my_filter);
      solver->correct(new_coords, rhs);

      // Compute new functional value
      fval_post = rumpflpumpfl.compute_func(&func_norm, &func_det, &func_rec_det);

      // With the new functional, the functional value should still have decreased
      TEST_CHECK(fval_pre > fval_post);
      // These differences should all be greater than eps
      TEST_CHECK(Math::abs(func_norm - DataType(0)) > eps);
      TEST_CHECK(Math::abs(func_det - my_functional->_fac_det*DataType(1)) > eps);

      solver->done();
      delete rmn;

    }
};
// Vanilla Rumpf smoother
template<typename A, typename B, typename C, typename D, typename E>
using MyQualityFunctional = Meshopt::HyperelasticityFunctional<A, B, C, D, E>;

HyperelasticityFunctionalTest_2d<float, Shape::Hypercube<2>, Meshopt::RumpfFunctional, MyQualityFunctional> test_hc_1;
HyperelasticityFunctionalTest_2d<double, Shape::Hypercube<2>, Meshopt::RumpfFunctional_D2, MyQualityFunctional> test_hc_2;
HyperelasticityFunctionalTest_2d<double, Shape::Simplex<2>, Meshopt::RumpfFunctional, MyQualityFunctional> test_s_1;
HyperelasticityFunctionalTest_2d<float, Shape::Simplex<2>, Meshopt::RumpfFunctional_D2, MyQualityFunctional> test_s_2;

template<typename A, typename B>
using MyFunctionalQ1Split = Meshopt::RumpfFunctionalQ1Split<A, B, Meshopt::RumpfFunctional>;

template<typename A, typename B>
using MyFunctionalQ1Split_D2 = Meshopt::RumpfFunctionalQ1Split<A, B, Meshopt::RumpfFunctional_D2>;

HyperelasticityFunctionalTest_2d<float, Shape::Hypercube<2>, MyFunctionalQ1Split, MyQualityFunctional> test_q1hack_f_1;
HyperelasticityFunctionalTest_2d<double, Shape::Hypercube<2>, MyFunctionalQ1Split_D2, MyQualityFunctional> test_q1hack_d_2;
HyperelasticityFunctionalTest_2d<double, Shape::Simplex<2>, MyFunctionalQ1Split, MyQualityFunctional> test_q1hack_ds_1;

/// \brief Specialisation for hypercubes
template<int shape_dim_>
struct helperclass< FEAST::Shape::Hypercube<shape_dim_> >
{
  /// \brief Sets coordinates so we deal the the reference element
  template<typename VectorType_, typename DataType_>
  static void buffer_to_mesh(VectorType_& coords_, const DataType_& scaling)
  {
    for(Index i(0); i < Index(1 << shape_dim_); ++i)
    {
      Tiny::Vector<DataType_, VectorType_::BlockSize, VectorType_::BlockSize> tmp;
      for(int d(0); d < shape_dim_; ++d)
        tmp(d) = (DataType_(((i >> d) & 1) << 1) - DataType_(1)) * scaling ;

      coords_(i, tmp);
    }
  }

};

/// \brief Specialisation for 2d simplices
template<>
struct helperclass< FEAST::Shape::Simplex<2> >
{
  /// \brief Sets coordinates so we deal the the Rumpf reference element
  template<typename VectorType_, typename DataType_>
  static void buffer_to_mesh(VectorType_& coords_, const DataType_& scaling)
  {
    Tiny::Vector<DataType_, 2, 2> tmp(0);
    coords_(0, tmp);

    tmp(0) = scaling;
    coords_(1, tmp);

    tmp(0) = DataType_(0.5) * scaling;
    tmp(1) = DataType_(0.5)*Math::sqrt(DataType_(3))*scaling;
    coords_(2, tmp);
  }
};
