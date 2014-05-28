#include <test_system/test_system.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename DataType_>
class LinearFunctionalTest :
  public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef LAFEM::DenseVector<Mem::Main, DataType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  LinearFunctionalTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("LinearFunctionalTest")
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

  void test_unit_2d_q1(QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ1 space(trafo);

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("gauss-legendre:2");

    // assemble the linear functional into a vector
    VectorType vector(space.get_num_dofs(), DataType_(0));
    Assembly::Common::ConstantFunction function(DataType_(1));
    Assembly::Common::ForceFunctional<Assembly::Common::ConstantFunction> functional(function);
    Assembly::LinearFunctionalAssembler::assemble_vector(vector, functional, space, cubature_factory);

    // get mesh element count
    Index num_verts = mesh.get_num_entities(0);
    Index num_quads = mesh.get_num_entities(2);

    // create another vector
    VectorType vector2(num_verts, DataType_(0));

    // compute squared mesh width
    const DataType_ weight = DataType_(1) / DataType_(4*num_quads);

    // get the vertex set of the mesh
    const typename QuadMesh::VertexSetType vertex_set(mesh.get_vertex_set());

    // get the index set of the mesh
    const Geometry::IndexSet<4> index_set(mesh.get_index_set<2,0>());

    // loop over all quads
    for(Index i(0); i < num_quads; ++i)
    {
      for(int j(0); j < 4; ++j)
        vector2(index_set[i][j], vector2(index_set[i][j]) + weight);
    }

    // loop over all vertices
    for(Index i(0); i < num_verts; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vector(i), vector2(i), eps);
    }
  }

  void test_unit_2d_q0(QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ0 space(trafo);

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("gauss-legendre:5");

    // assemble the linear functional into a vector
    VectorType vector(space.get_num_dofs(), DataType_(0));
    Assembly::Common::SineBubbleFunction function;
    Assembly::Common::ForceFunctional<Assembly::Common::SineBubbleFunction> functional(function);
    Assembly::LinearFunctionalAssembler::assemble_vector(vector, functional, space, cubature_factory);

    // get mesh element count
    Index num_quads = mesh.get_num_entities(2);

    // validate vector size
    TEST_CHECK_EQUAL(vector.size(), num_quads);

    // get the vertex set of the mesh
    const typename QuadMesh::VertexSetType vertex_set(mesh.get_vertex_set());

    // get the index set of the mesh
    const Geometry::IndexSet<4> index_set(mesh.get_index_set<2,0>());

    // get the constant pi
    const DataType_ pi = Math::pi<DataType_>();

    // loop over all quads
    for(Index i(0); i < num_quads; ++i)
    {
      // get quad dimensions
      DataType_ x0 = DataType_(vertex_set[index_set[i][0]][0]);
      DataType_ x1 = DataType_(vertex_set[index_set[i][3]][0]);
      DataType_ y0 = DataType_(vertex_set[index_set[i][0]][1]);
      DataType_ y1 = DataType_(vertex_set[index_set[i][3]][1]);

      // compute analytic integral:
      // int_[x0,x1]x[y0,y1] sin(Pi*x)*sin(pi*y) dxy = (cos(pi*x0) - cos(pi*x1))*(cos(pi*y0) - cos(pi*y1)) / pi^2
      DataType_ s = ((Math::cos(pi*x0) - Math::cos(pi*x1)) * (Math::cos(pi*y0)-Math::cos(pi*y1))) / (pi*pi);

      // validate vector data
      TEST_CHECK_EQUAL_WITHIN_EPS(vector(i), s, eps);
    }
  }
};

LinearFunctionalTest<float> linear_functional_test_float;
LinearFunctionalTest<double> linear_functional_test_double;
