#include <test_system/test_system.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/rew_projector.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief RewProjector test class template
 *
 * \test Tests the Assembly::RewProjector class
 *
 * \tparam DataType_
 * The data type for the test. Shall be either double or float.
 *
 * \author Peter Zajac
 */
template<typename DataType_>
class RewProjectorTest :
  public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef LAFEM::DenseVector<Mem::Main, DataType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
  typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceRT;

public:
  RewProjectorTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("RewProjectorTest")
  {
  }

  virtual void run() const
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

    // create coarse mesh
    Geometry::RefinedUnitCubeFactory<QuadMesh> unit_factory(3);
    QuadMesh mesh(unit_factory);

    // create trafo
    QuadTrafo trafo(mesh);

    // project Q0
    DataType_ q0_err = project<QuadSpaceQ0>(trafo, "gauss-legendre:2");
    q0_err = q0_err/DataType_(7.9731492672E-2) - DataType_(1);
    TEST_CHECK_EQUAL_WITHIN_EPS(q0_err, DataType_(0), eps);
    //TEST_CHECK_EQUAL_WITHIN_EPS(q0_err, DataType_(7.97315E-2), eps);

    // project Q1
    DataType_ q1_err = project<QuadSpaceQ1>(trafo, "gauss-legendre:3");
    q1_err = q1_err/DataType_(4.1430308042E-3) - DataType_(1);
    TEST_CHECK_EQUAL_WITHIN_EPS(q1_err, DataType_(0), eps);
    //TEST_CHECK_EQUAL_WITHIN_EPS(q1_err, DataType_(4.14303E-3), eps);

    // project RT
    DataType_ rt_err = project<QuadSpaceRT>(trafo, "gauss-legendre:3");
    rt_err = rt_err/DataType_(7.5903151744E-3) - DataType_(1);
    TEST_CHECK_EQUAL_WITHIN_EPS(rt_err, DataType_(0), eps);
    //TEST_CHECK_EQUAL_WITHIN_EPS(rt_err, DataType_(7.59032E-3), eps);
  }

  template<typename Space_>
  DataType_ project(QuadTrafo& trafo, String cubature_name) const
  {
    // create space
    Space_ space(trafo);

    // define functor
    Assembly::Common::SineBubbleFunction function;

    // define a cubature factory
    Cubature::DynamicFactory cubature_factory(cubature_name);

    // project functor into FE space
    VectorType vector;
    Assembly::RewProjector::project(vector, function, space, cubature_factory);

    // compute L2-Error
    return //std::make_pair(
      Assembly::ScalarErrorComputerL2::compute(vector, function, space, cubature_factory)/*,
      Assembly::ScalarErrorComputerH1::compute(vector, function, space, cubature_factory))*/;
  }

};

RewProjectorTest<float> rew_projector_test_float;
RewProjectorTest<double> rew_projector_test_double;
