#include <test_system/test_system.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename DataType_>
class SineBubble;


template<typename DataType_>
class InterpolatorTest :
  public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef LAFEM::DenseVector<Mem::Main, DataType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  InterpolatorTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("InterpolatorTest")
  {
  }

  virtual void run() const
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    // create coarse mesh
    Geometry::RefinedUnitCubeFactory<QuadMesh> unit_factory(3);
    QuadMesh mesh(unit_factory);

    // run tests
    test_unit_2d_q1(mesh);
    test_unit_2d_q0(mesh);
  }

  void test_unit_2d_q1(const QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ1 space(trafo);

    // define functor
    typedef Analytic::StaticWrapperFunctor<SineBubble, true, true> FunctorType;
    FunctorType functor;

    // interpolate functor into FE space
    VectorType vector;
    Assembly::Interpolator::project(vector, functor, space);

    // get mesh vertex count
    Index num_verts = mesh.get_num_entities(0);

    // validate vector size
    TEST_CHECK_EQUAL(vector.size(), num_verts);

    // get the vertex set of the mesh
    const typename QuadMesh::VertexSetType vertex_set(mesh.get_vertex_set());

    // loop over all vertices
    for(Index i(0); i < num_verts; ++i)
    {
      // compute sine-bubble value in vertex position
      DataType_ s = SineBubble<DataType_>::eval(DataType_(vertex_set[i][0]), DataType_(vertex_set[i][1]));

      // validate vector data
      TEST_CHECK_EQUAL_WITHIN_EPS(vector(i), s, eps);
    }
  }

  void test_unit_2d_q0(const QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps =Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ0 space(trafo);

    // define functor
    typedef Analytic::StaticWrapperFunctor<SineBubble, true, true> FunctorType;
    FunctorType functor;

    // interpolate functor into FE space
    VectorType vector;
    Assembly::Interpolator::project(vector, functor, space);

    // get mesh element count
    Index num_quads = mesh.get_num_entities(2);

    // validate vector size
    TEST_CHECK_EQUAL(vector.size(), num_quads);

    // get the vertex set of the mesh
    const typename QuadMesh::VertexSetType vertex_set(mesh.get_vertex_set());

    // get the index set of the mesh
    const Geometry::IndexSet<4> index_set(mesh.get_index_set<2,0>());

    // loop over all quads
    for(Index i(0); i < num_quads; ++i)
    {
      // compute quad center
      DataType_ x(0), y(0);
      for(int j(0); j < 4; ++j)
      {
        x += DataType_(vertex_set[index_set[i][j]][0]);
        y += DataType_(vertex_set[index_set[i][j]][1]);
      }

      // compute sine-bubble value in quad center
      DataType_ s = SineBubble<DataType_>::eval(DataType_(0.25)*x, DataType_(0.25)*y);

      // validate vector data
      TEST_CHECK_EQUAL_WITHIN_EPS(vector(i), s, eps);
    }
  }
};

InterpolatorTest<float> interpolator_test_float;
InterpolatorTest<double> interpolator_test_double;


template<typename DataType_>
class SineBubble :
  public Analytic::StaticFunction<DataType_>
{
public:
  /// returns the constant pi
  static DataType_ pi()
  {
    return Math::Limits<DataType_>::pi();
  }

  /// 1D: function value
  static DataType_ eval(DataType_ x)
  {
    return Math::sin(pi() * x);
  }

  /// 2D: function value
  static DataType_ eval(DataType_ x, DataType_ y)
  {
    return Math::sin(pi() * x) * Math::sin(pi() * y);
  }

  /// 3D: function value
  static DataType_ eval(DataType_ x, DataType_ y, DataType_ z)
  {
    return Math::sin(pi() * x) * Math::sin(pi() * y) * Math::sin(pi() * z);
  }

  /// 1D: X-derivative
  static DataType_ der_x(DataType_ x)
  {
    return pi() * Math::cos(pi() * x);
  }

  /// 2D: X-derivative
  static DataType_ der_x(DataType_ x, DataType_ y)
  {
    return pi() * Math::cos(pi() * x) * Math::sin(pi() * y);
  }

  /// 2D: Y-derivative
  static DataType_ der_y(DataType_ x, DataType_ y)
  {
    return pi() * Math::sin(pi() * x) * Math::cos(pi() * y);
  }

  /// 3D: X-derivative
  static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z)
  {
    return pi() * Math::cos(pi() * x) * Math::sin(pi() * y) * Math::sin(pi() * z);
  }

  /// 3D: Y-derivative
  static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z)
  {
    return pi() * Math::sin(pi() * x) * Math::cos(pi() * y) * Math::sin(pi() * z);
  }

  /// 3D: Z-derivative
  static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z)
  {
    return pi() * Math::sin(pi() * x) * Math::sin(pi() * y) * Math::cos(pi() * z);
  }
}; // class SineBubble<...>
