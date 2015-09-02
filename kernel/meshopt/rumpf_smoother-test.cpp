#include <kernel/base_header.hpp>
#ifdef FEAST_HAVE_ALGLIB
#include <test_system/test_system.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/meshopt/rumpf_smoother.hpp>
#include <kernel/meshopt/rumpf_smoother_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>

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
 * When the input mesh consists of a single Rumpf reference element, the mesh optimiser is supposed not to change it,
 * as:
 * 1) The functional value is minimal for the reference element and
 * 2) the functional gradient evaluates to zero meaning that
 * 3) the optimiser stops after one evaluation of the gradient.
 *
 * This means the functional value pre and post optimisation should be the same.
 *
 * \author Jordi Paul
 **/
template
<
  typename DT_,
  typename ShapeType_,
  template<typename, typename> class FunctionalType_,
  template<typename ... > class RumpfSmootherType_
>
class RumpfSmootherTest_2d
: public TestSystem::FullTaggedTest<Mem::Main, DT_, Index>
{
  public:
    typedef Mem::Main MemType;
    typedef Index IndexType;
    typedef DT_ DataType;

    typedef ShapeType_ ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    typedef FunctionalType_<DataType, ShapeType> FunctionalType;
    typedef RumpfSmootherType_<TrafoType, FunctionalType> RumpfSmootherType;

    RumpfSmootherTest_2d() :
      TestSystem::FullTaggedTest<MemType, DataType, IndexType>("rumpf_smoother_test")
      {
      }

    virtual void run() const
    {
      // Mesh and trafo
      Geometry::ReferenceCellFactory<ShapeType, DataType> mesh_factory;
      MeshType mesh(mesh_factory);
      TrafoType trafo(mesh);

      // In 2d, the cofactor matrix is not used
      DataType fac_norm = DataType(1e0),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
      FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

      // The smoother in all its template glory
      RumpfSmootherType rumpflpumpfl(trafo, my_functional);

      // It is possible to scale the reference element, but we do not want that here
      DataType scaling(1);
      // This transforms the unit element to the Rumpf reference element
      helperclass<ShapeType>::set_coords(rumpflpumpfl._coords, scaling);

      // Since we changed the internal _coords, they have to be copied back to the mesh
      rumpflpumpfl.set_coords();
      rumpflpumpfl.init();

      // Compute initial functional value
      DataType fval_pre = rumpflpumpfl.compute_functional();
      // Optimise the mesh
      rumpflpumpfl.optimise();
      // Compute new functional value
      DataType fval_post = rumpflpumpfl.compute_functional();

      const DataType eps = Math::pow(Math::eps<DataType>(),DataType(0.8));
      TEST_CHECK_EQUAL_WITHIN_EPS(fval_pre, fval_post, eps);

    }
};

template<typename A, typename B>
using MySmoother = Meshopt::RumpfSmoother<A, B>;

template<typename A, typename B>
using MySmootherQ1Hack = Meshopt::RumpfSmootherQ1Hack<A, B>;

RumpfSmootherTest_2d<float, Shape::Hypercube<2>, Meshopt::RumpfFunctional, MySmoother> test_hc_f_1;
RumpfSmootherTest_2d<float, Shape::Hypercube<2>, Meshopt::RumpfFunctional_D2, MySmoother> test_hc_f_2;
RumpfSmootherTest_2d<float, Shape::Simplex<2>, Meshopt::RumpfFunctional, MySmoother> test_s_f_1;
RumpfSmootherTest_2d<float, Shape::Simplex<2>, Meshopt::RumpfFunctional_D2, MySmoother> test_s_f_2;

RumpfSmootherTest_2d<double, Shape::Hypercube<2>, Meshopt::RumpfFunctional, MySmoother> test_hc_d_1;
RumpfSmootherTest_2d<double, Shape::Hypercube<2>, Meshopt::RumpfFunctional_D2, MySmoother> test_hc_d_2;
RumpfSmootherTest_2d<double, Shape::Simplex<2>, Meshopt::RumpfFunctional, MySmoother> test_s_d_1;
RumpfSmootherTest_2d<double, Shape::Simplex<2>, Meshopt::RumpfFunctional_D2, MySmoother> test_s_d_2;

template<typename A, typename B>
using MyFunctionalQ1Hack = Meshopt::RumpfFunctionalQ1Hack<A, B, Meshopt::RumpfFunctional>;

template<typename A, typename B>
using MyFunctionalQ1Hack_D2 = Meshopt::RumpfFunctionalQ1Hack<A, B, Meshopt::RumpfFunctional_D2>;

RumpfSmootherTest_2d<float, Shape::Hypercube<2>, MyFunctionalQ1Hack, MySmootherQ1Hack> test_q1hack_f_1;
RumpfSmootherTest_2d<double, Shape::Hypercube<2>, MyFunctionalQ1Hack, MySmootherQ1Hack> test_q1hack_d_1;

/// \compilerhack
// The MSVC CTP_Nov2013 contains a bug that causes the compiler to skip template aliases when generating its
// COMDATs (compiler generated decorated function names) so it cannot distinguish between these function names and
// the ones above. The bug has been fixed in the regular version of the compiler (which we cannot use due to missing
// C++11 features), but not in the CTP_Nov2013 version.
//
// See http://blogs.msdn.com/b/vcblog/archive/2014/08/04/bugs-fixed-in-visual-studio-2013-update-3.aspx
//
// Compiler hack to be removed once the new (Visual Studio 201{>13}) compiler is deployed.
//#ifndef FEAST_COMPILER_MICROSOFT
RumpfSmootherTest_2d<double, Shape::Hypercube<2>, MyFunctionalQ1Hack_D2, MySmootherQ1Hack> test_q1hack_d_2;
RumpfSmootherTest_2d<float, Shape::Hypercube<2>, MyFunctionalQ1Hack_D2, MySmootherQ1Hack> test_q1hack_f_2;
//#endif

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
};

#endif
