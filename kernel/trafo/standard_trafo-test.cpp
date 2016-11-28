#include <test_system/test_system.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/geometry/shape_convert_factory.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/standard/inverse_mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Standard Trafo test
 *
 * \test Tests the Standard Trafo class template space
 *
 * \tparam DataType_
 * The data type for the test. Shall be either double or float.
 *
 * \author Stefan Wahlers
 */
template<typename DataType_>
class StandardTrafoTest
  : public TestSystem::TaggedTest<Archs::None, DataType_>
{
  static constexpr TrafoTags unit_config = static_cast<TrafoTags>(~0);

public:
  StandardTrafoTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("StandardTrafoTest")
    {
    }

  virtual ~StandardTrafoTest()
  {
  }

  virtual void run() const override
  {
    // test assembly on unit quad
    test_unit_quad();
    test_quad();
    test_unit_tria();
    test_tria();
    test_unit_cube();
    test_cube();
    test_unit_tetra();
    test_tetra();
  }

  void test_unit_quad() const
  {
    typedef Shape::Quadrilateral ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a quad mesh
    Geometry::UnitCubeFactory<MeshType> mesh_factory;
    MeshType mesh(mesh_factory);

    // create a quad-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);

    // create a trafo evaluation data
    typename TrafoEvaluator::template ConfigTraits<unit_config>::EvalDataType trafo_data;

    // create a domain point
    typename TrafoEvaluator::DomainPointType dom_point;

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(0);
    dom_point[1] = DataType_(0);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqr(trafo_data.img_point[0] - DataType_(1)/DataType_(2))
      + Math::sqr(trafo_data.img_point[1] - DataType_(1)/DataType_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_(1)/DataType_(2), eps);

    // check hessian tensor
    for(int i(0); i < 2; ++i)
    {
      for(int j(0); j < 2; ++j)
      {
        for(int k(0); k < 2; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0), eps);
        }
      }
    }

    trafo_eval.finish();
  } // test_unit_quad()

  void test_quad() const
  {
    typedef Shape::Quadrilateral ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a quad mesh
    Geometry::UnitCubeFactory<MeshType> mesh_factory;
    MeshType mesh(mesh_factory);

    // create a quad-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);

    // create a trafo evaluation data
    typename TrafoEvaluator::template ConfigTraits<unit_config>::EvalDataType trafo_data;

    // create a domain point
    typename TrafoEvaluator::DomainPointType dom_point;

    mesh.get_vertex_set()[0][0] = -DataType_(1)/DataType_(4);
    mesh.get_vertex_set()[0][1] =  DataType_(1)/DataType_(4);

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(0);
    dom_point[1] = DataType_(0);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqr(trafo_data.img_point[0] - DataType_(7)/DataType_(16))
      + Math::sqr(trafo_data.img_point[1] - DataType_(9)/DataType_(16));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0],  DataType_(9)/DataType_(16), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1],  DataType_(1)/DataType_(16), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], -DataType_(1)/DataType_(16), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1],  DataType_(7)/DataType_(16), eps);

    // check hessian tensor
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][1], -DataType_(1)/DataType_(16), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][0], -DataType_(1)/DataType_(16), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][1], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][1], DataType_(1)/DataType_(16), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][0], DataType_(1)/DataType_(16), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][1], DataType_(0), eps);

    trafo_eval.finish();
  }// test_quad()

  void test_unit_cube() const
  {
    typedef Shape::Hexahedron ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a cube mesh
    Geometry::UnitCubeFactory<MeshType> mesh_factory;
    MeshType mesh(mesh_factory);

    // create a cube-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);

    // create a trafo evaluation data
    typename TrafoEvaluator::template ConfigTraits<unit_config>::EvalDataType trafo_data;

    // create a domain point
    typename TrafoEvaluator::DomainPointType dom_point;

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(0);
    dom_point[1] = DataType_(0);
    dom_point[2] = DataType_(0);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqr(trafo_data.img_point[0] - DataType_(1)/DataType_(2))
      + Math::sqr(trafo_data.img_point[1] - DataType_(1)/DataType_(2))
      + Math::sqr(trafo_data.img_point[2] - DataType_(1)/DataType_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][2], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][2], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][1], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][2], DataType_(1)/DataType_(2), eps);

    // check hessian tensor
    for(int i(0); i < 3; ++i)
    {
      for(int j(0); j < 3; ++j)
      {
        for(int k(0); k < 3; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0), eps);
        }
      }
    }

    trafo_eval.finish();
  }// test_unit_cube()

  void test_cube() const
  {
    typedef Shape::Hexahedron ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a cube mesh
    Geometry::UnitCubeFactory<MeshType> mesh_factory;
    MeshType mesh(mesh_factory);

    // create a cube-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);

    // create a trafo evaluation data
    typename TrafoEvaluator::template ConfigTraits<unit_config>::EvalDataType trafo_data;

    // create a domain point
    typename TrafoEvaluator::DomainPointType dom_point;

    mesh.get_vertex_set()[7][0] = DataType_(9)/DataType_(10);
    mesh.get_vertex_set()[7][1] = DataType_(8)/DataType_(10);
    mesh.get_vertex_set()[7][2] = DataType_(3)/DataType_(4);

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(0);
    dom_point[1] = DataType_(0);
    dom_point[2] = DataType_(0);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqr(trafo_data.img_point[0] - DataType_(39)/DataType_(80))
      + Math::sqr(trafo_data.img_point[1] - DataType_(19)/DataType_(40))
      + Math::sqr(trafo_data.img_point[2] - DataType_(15)/DataType_(32));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_(39)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], -DataType_(1)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][2], -DataType_(1)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], -DataType_(1)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_(19)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][2], -DataType_(1)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][0], -DataType_(1)/DataType_(32), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][1], -DataType_(1)/DataType_(32), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][2], DataType_(15)/DataType_(32), eps);

    // check hessian tensor
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][1], -DataType_(1)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][2], -DataType_(1)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][0], -DataType_(1)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][1], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][2], -DataType_(1)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][2][0], -DataType_(1)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][2][1], -DataType_(1)/DataType_(80), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][2][2], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][1], -DataType_(1)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][2], -DataType_(1)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][0], -DataType_(1)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][1], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][2], -DataType_(1)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][2][0], -DataType_(1)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][2][1], -DataType_(1)/DataType_(40), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][2][2], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][0][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][0][1], -DataType_(1)/DataType_(32), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][0][2], -DataType_(1)/DataType_(32), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][1][0], -DataType_(1)/DataType_(32), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][1][1], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][1][2], -DataType_(1)/DataType_(32), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][2][0], -DataType_(1)/DataType_(32), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][2][1], -DataType_(1)/DataType_(32), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][2][2], DataType_(0), eps);

    trafo_eval.finish();
  } // test_cube()

  void test_unit_tria() const
  {
    typedef Shape::Triangle ShapeType;
    typedef Shape::Quadrilateral QuadType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::ConformalMesh<QuadType> QuadMeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a quad mesh
    Geometry::UnitCubeFactory<QuadMeshType> quad_mesh_factory;
    QuadMeshType quad_mesh(quad_mesh_factory);

    // convert to tria mesh
    Geometry::ShapeConvertFactory<MeshType> tria_mesh_factory(quad_mesh);
    MeshType mesh(tria_mesh_factory);

    // create a tria-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);

    // create a trafo evaluation data
    typename TrafoEvaluator::template ConfigTraits<unit_config>::EvalDataType trafo_data;

    // create a domain point
    typename TrafoEvaluator::DomainPointType dom_point;

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(1)/DataType_(2);
    dom_point[1] = DataType_(1)/DataType_(2);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqr(trafo_data.img_point[0] - DataType_(1)/DataType_(2))
      + Math::sqr(trafo_data.img_point[1] - DataType_(0));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], -DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1],  DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], -DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], -DataType_(1)/DataType_(2), eps);

    // check hessian tensor
    for(int i(0); i < 2; ++i)
    {
      for(int j(0); j < 2; ++j)
      {
        for(int k(0); k < 2; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0), eps);
        }
      }
    }

    trafo_eval.finish();
  } // test_unit_tria()

  void test_tria() const
  {
    typedef Shape::Triangle ShapeType;
    typedef Shape::Quadrilateral QuadType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::ConformalMesh<QuadType> QuadMeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a quad mesh
    Geometry::UnitCubeFactory<QuadMeshType> quad_mesh_factory;
    QuadMeshType quad_mesh(quad_mesh_factory);

    // convert to tria mesh
    Geometry::ShapeConvertFactory<MeshType> tria_mesh_factory(quad_mesh);
    MeshType mesh(tria_mesh_factory);

    // create a tria-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);

    // create a trafo evaluation data
    typename TrafoEvaluator::template ConfigTraits<unit_config>::EvalDataType trafo_data;

    // create a domain point
    typename TrafoEvaluator::DomainPointType dom_point;

    mesh.get_vertex_set()[1][0] =  DataType_(1)/DataType_(10);
    mesh.get_vertex_set()[1][1] = -DataType_(1)/DataType_(10);

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(1)/DataType_(4);
    dom_point[1] = DataType_(1)/DataType_(4);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqr(trafo_data.img_point[0] - DataType_(11)/DataType_(40))
      + Math::sqr(trafo_data.img_point[1] - DataType_(9)/DataType_(40));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], -DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], -DataType_(4)/DataType_(10), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], -DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], -DataType_(6)/DataType_(10), eps);

    // check hessian tensor
    for(int i(0); i < 2; ++i)
    {
      for(int j(0); j < 2; ++j)
      {
        for(int k(0); k < 2; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0), eps);
        }
      }
    }

    trafo_eval.finish();
  } // test_tria()

  void test_unit_tetra() const
  {
    typedef Shape::Tetrahedron ShapeType;
    typedef Shape::Hexahedron QuadType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::ConformalMesh<QuadType> QuadMeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a cube mesh
    Geometry::UnitCubeFactory<QuadMeshType> quad_mesh_factory;
    QuadMeshType quad_mesh(quad_mesh_factory);

    // convert to tetra mesh
    Geometry::ShapeConvertFactory<MeshType> tria_mesh_factory(quad_mesh);
    MeshType mesh(tria_mesh_factory);

    // create a quad-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);

    // create a trafo evaluation data
    typename TrafoEvaluator::template ConfigTraits<unit_config>::EvalDataType trafo_data;

    // create a domain point
    typename TrafoEvaluator::DomainPointType dom_point;

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(1)/DataType_(4);
    dom_point[1] = DataType_(1)/DataType_(4);
    dom_point[2] = DataType_(1)/DataType_(4);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqr(trafo_data.img_point[0] - DataType_(1)/DataType_(2))
      + Math::sqr(trafo_data.img_point[1] - DataType_(1)/DataType_(4))
      + Math::sqr(trafo_data.img_point[2] - DataType_(1)/DataType_(8));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_(1), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][2], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][2], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][1], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][2], DataType_(1)/DataType_(2), eps);

    // check hessian tensor
    for(int i(0); i < 3; ++i)
    {
      for(int j(0); j < 3; ++j)
      {
        for(int k(0); k < 3; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0), eps);
        }
      }
    }

    trafo_eval.finish();
  }// test_unit_tetra()

  void test_tetra() const
  {
    typedef Shape::Tetrahedron ShapeType;
    typedef Shape::Hexahedron QuadType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::ConformalMesh<QuadType> QuadMeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a cube mesh
    Geometry::UnitCubeFactory<QuadMeshType> quad_mesh_factory;
    QuadMeshType quad_mesh(quad_mesh_factory);

    // convert to tetra mesh
    Geometry::ShapeConvertFactory<MeshType> tria_mesh_factory(quad_mesh);
    MeshType mesh(tria_mesh_factory);

    // create a quad-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);

    // create a trafo evaluation data
    typename TrafoEvaluator::template ConfigTraits<unit_config>::EvalDataType trafo_data;

    // create a domain point
    typename TrafoEvaluator::DomainPointType dom_point;

    mesh.get_vertex_set()[8][0] = DataType_(4)/DataType_(10);
    mesh.get_vertex_set()[8][1] = DataType_(3)/DataType_(10);
    mesh.get_vertex_set()[8][2] = DataType_(1)/DataType_(4);

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(1)/DataType_(4);
    dom_point[1] = DataType_(1)/DataType_(4);
    dom_point[2] = DataType_(1)/DataType_(4);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqr(trafo_data.img_point[0] - DataType_(19)/DataType_(40))
      + Math::sqr(trafo_data.img_point[1] - DataType_(1)/DataType_(5))
      + Math::sqr(trafo_data.img_point[2] - DataType_(3)/DataType_(16));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_(1), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_(4)/DataType_(10), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][2], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_(3)/DataType_(10), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][2], DataType_(1)/DataType_(2), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][0], DataType_(0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][1], DataType_(1)/DataType_(4), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][2], DataType_(1)/DataType_(2), eps);

    // check hessian tensor
    for(int i(0); i < 3; ++i)
    {
      for(int j(0); j < 3; ++j)
      {
        for(int k(0); k < 3; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0), eps);
        }
      }
    }

    trafo_eval.finish();
  } // test_tetra()
};

StandardTrafoTest<double> standard_trafo_test_double;
StandardTrafoTest<float> standard_trafo_test_float;

/**
 * \brief Test for the inverse mapping for simplices in 1d, 2d, 3d.
 *
 * This also tests the variants with shape_dim < world_dim
 *
 * \author Jordi Paul
 */
template<typename DT_>
class InverseMappingTest
: public TestSystem::TaggedTest<Mem::Main, DT_>
{
  public:
    InverseMappingTest() :
      TestSystem::TaggedTest<Mem::Main, DT_>("InverseMappingTest")
      {
      }

    virtual ~InverseMappingTest()
    {
    }

    virtual void run() const override
    {
      run_1d();
      run_2d();
      run_3d();
      run_2d_edge();
      run_3d_edge();
      run_3d_tria();
    }

    /// Runs tests for a Simplex<1> in 1d
    void run_1d() const
    {
      // Number of points we check
      static constexpr int num_points = 4;

      static constexpr int dim = 1;
      // Tolerance
      DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.9)));

      // 2d test
      Tiny::Vector<DT_, dim> coeffs(DT_(0));
      Tiny::Matrix<DT_, dim+1, dim > coords;

      coords(0,0) = DT_(-1);
      coords(1,0) = DT_(0.142);

      Tiny::Matrix<DT_, num_points, dim> points(DT_(0));
      // Point inside
      points[0] = DT_(0.53)*coords[0] + DT_(0.47)*coords[1];
      // One of the vertices
      points[1] = coords[0];
      // Point outside, coefficients are not barycentric coordinates because sum != 1
      points[2] = DT_(-1.5)*coords[0] + DT_(0.25)*coords[1];
      // Nearly one of the vertices
      points[3] = (DT_(1)+DT_(10)*tol)*coords[0];

      for(int np(0); np < num_points; ++np)
      {
        Trafo::Standard::inverse_mapping(coeffs, points[np], coords);

        Tiny::Vector<DT_, dim> test(coords[0]);
        for(int i(0); i < dim; ++i)
          test+= coeffs(i)*(coords[i+1]-coords[0]);

        for(int i(0); i < dim; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(test(i), points[np](i), tol);
      }
    }

    /// Runs tests for a Simplex<2> in 2d
    void run_2d() const
    {
      // Number of points we check
      static constexpr int num_points = 5;

      static constexpr int dim = 2;
      // Tolerance
      DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.9)));

      // 2d test
      Tiny::Vector<DT_, dim> coeffs(DT_(0));
      Tiny::Matrix<DT_, dim+1, dim > coords;

      coords(0,0) = DT_(-1);
      coords(0,1) = DT_(-1);

      coords(1,0) = DT_(0);
      coords(1,1) = DT_(1);

      coords(2,0) = DT_(1);
      coords(2,1) = DT_(0.2);

      Tiny::Matrix<DT_, num_points, dim> points(DT_(0));
      // Point inside
      points[0] = DT_(0.3)*coords[0] + DT_(0.47)*coords[1] + DT_(0.23)*coords[2];
      // One of the vertices
      points[1] = coords[0];
      // Point on one edge
      points[2] = DT_(0.5)*coords[0] + DT_(0.5)*coords[2];
      // Point outside, coefficients are not barycentric coordinates because sum != 1
      points[3] = DT_(-2)*coords[0] + DT_(1.25)*coords[1] + DT_(0.75)*coords[2];
      // Nearly one of the vertices
      points[4] = (DT_(1)+DT_(10)*tol)*coords[0];

      for(int np(0); np < num_points; ++np)
      {
        Trafo::Standard::inverse_mapping(coeffs, points[np], coords);

        Tiny::Vector<DT_, dim> test(coords[0]);
        for(int i(0); i < dim; ++i)
          test+= coeffs(i)*(coords[i+1]-coords[0]);

        for(int i(0); i < dim; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(test(i), points[np](i), tol);
      }
    }

    // Runs tests for a Simplex<3> in 3d
    void run_3d() const
    {
      static constexpr int num_points = 5;
      static constexpr int dim = 3;
      // Tolerance
      DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.75)));

      // 2d test
      Tiny::Vector<DT_, dim> coeffs(DT_(0));
      Tiny::Matrix<DT_, dim+1, dim > coords;

      coords(0,0) = -DT_(1);
      coords(0,1) = -DT_(1);
      coords(0,2) = -DT_(6);

      coords(1,0) = DT_(1);
      coords(1,1) = DT_(0);
      coords(1,2) = DT_(0);

      coords(2,0) = DT_(0);
      coords(2,1) = DT_(1);
      coords(2,2) = DT_(0.2);

      coords(3,0) = -DT_(0.1);
      coords(3,1) = DT_(0);
      coords(3,2) = DT_(2);

      Tiny::Matrix<DT_, num_points, dim> points(DT_(0));
      // Point inside
      points[0] = DT_(0.15)*coords[0] + DT_(0.17)*coords[1] + DT_(0.53)*coords[2] + DT_(0.15)*coords[3];
      // One of the vertices
      points[1] = coords[2];
      // Point on one edge
      points[2] = DT_(0.5)*coords[0] + DT_(0.5)*coords[3];
      // Point outside, coefficients are not barycentric coordinates because sum != 1
      points[3] = DT_(-2)*coords[0] + DT_(1.25)*coords[1] + DT_(0.75)*coords[2] + DT_(0.5)*coords[3];
      // Nearly one of the vertices
      points[4] = (DT_(1)+DT_(10)*tol)*coords[2];

      for(int np(0); np < num_points; ++np)
      {
        Trafo::Standard::inverse_mapping(coeffs, points[np], coords);

        Tiny::Vector<DT_, dim> test(coords[0]);
        for(int i(0); i < dim; ++i)
          test+= coeffs(i)*(coords[i+1] - coords[0]);

        for(int i(0); i < dim; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(test(i), points[np](i), tol);
      }
    }

    // Runs tests for a Simplex<1> in 2d
    void run_2d_edge() const
    {
      // Number of points we check
      static constexpr int num_points = 5;
      static constexpr int world_dim = 2;
      static constexpr int shape_dim = 1;
      // Tolerance
      DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.8)));

      Tiny::Vector<DT_, world_dim> coeffs(DT_(0));
      // The last column of coords will contain the orthogonal to the edge immersed in 2d
      Tiny::Matrix<DT_, world_dim+1, world_dim > coords;

      coords(0,0) = -DT_(1);
      coords(0,1) = DT_(0.3);

      coords(1,0) = DT_(2);
      coords(1,1) = -DT_(0.2);

      Tiny::Matrix<DT_, world_dim, world_dim-1> tmp_coords(DT_(0));
      for(int i(0); i < world_dim; ++i)
      {
        for(int j(0); j < world_dim-1; ++j)
          tmp_coords(i,j) = coords(j+1,i) - coords(0,i);
      }

      coords[world_dim] = Tiny::orthogonal(tmp_coords);
      coords[world_dim].normalise();

      Tiny::Matrix<DT_, num_points, world_dim> points(DT_(0));
      // Point inside
      points[0] = DT_(0.53)*coords[0] + DT_(0.47)*coords[1];
      // One of the vertices
      points[1] = coords[0];
      // Point outside, coefficients are not barycentric coordinates because sum != 1
      points[2] = coords[0] + DT_(1.25)*(coords[1] - coords[0]);
      // Nearly one of the vertices
      points[3] = (DT_(1)+ DT_(10)*tol)*coords[0] + ( - DT_(10)*tol )*coords[1] ;
      // Point with an orthogonal component, coordinates do not fulfill sum == 1
      points[4] = -DT_(0.2)*coords[0] + DT_(0.11)*coords[1] + DT_(0.69)*coords[2];

      for(int np(0); np < num_points; ++np)
      {
        Trafo::Standard::inverse_mapping(coeffs, points[np], coords.template size_cast<shape_dim+1, world_dim>());

        Tiny::Vector<DT_, world_dim> test(coords[0]);
        for(int i(0); i < world_dim-1; ++i)
          test += coeffs(i)*(coords[i+1] - coords[0]);

        test += coeffs(world_dim-1)*(coords[world_dim]);

        for(int i(0); i < world_dim; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(test(i), points[np](i), tol);

      }
    }

    // Runs tests for a Simplex<1> in 3d
    void run_3d_edge() const
    {
      static constexpr int num_points = 5;
      static constexpr int world_dim = 3;
      static constexpr int shape_dim = 1;
      // Tolerance
      DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.5)));

      Tiny::Vector<DT_, shape_dim+1> coeffs(DT_(0));
      // The first two colums in coords are the points defining the edge. The third is something to define a plane,
      // and the last will be the orthonormal to that plane, so we can check the distance of the point to the edge.
      Tiny::Matrix<DT_, world_dim+1, world_dim > coords;

      coords(0,0) = -DT_(1);
      coords(0,1) = -DT_(0.3);
      coords(0,2) = -DT_(0.3);

      coords(1,0) = DT_(2);
      coords(1,1) = -DT_(0.2);
      coords(1,2) = DT_(0.3);

      // This point is not part of the edge and is used solely for testing purposes
      coords(2,0) = DT_(0);
      coords(2,1) = DT_(0);
      coords(2,2) = DT_(1);

      Tiny::Matrix<DT_, world_dim, world_dim-1> tmp_coords(DT_(0));

      for(int i(0); i < world_dim; ++i)
      {
        for(int j(0); j < world_dim-1; ++j)
          tmp_coords(i,j) = coords(j+1,i) - coords(0,i);
      }

      coords[world_dim] = Tiny::orthogonal(tmp_coords);
      coords[world_dim].normalise();

      Tiny::Matrix<DT_, num_points, world_dim> points(DT_(0));
      // Point inside
      points[0] = DT_(0.11)*coords[0] + DT_(0.89)*coords[1];
      // One of the vertices
      points[1] = coords[1];
      // Point outside, coefficients are not barycentric coordinates because sum != 1
      points[2] = coords[0] + DT_(1.25)*(coords[1] - coords[0]);
      // Nearly one of the vertices
      points[3] = (DT_(1)+DT_(10)*tol)*coords[0] + (-DT_(10)*tol)*coords[1];
      // Something with an orthogonal component
      points[4] = DT_(1)*coords[0] - DT_(0.1)*tol*coords[1] + DT_(0.18)*coords[3];

      for(int np(0); np < num_points-1; ++np)
      {
        // We size-cast away the last column since it is just for result checking, so the correct version gets called
        Trafo::Standard::inverse_mapping(coeffs, points[np], coords.template size_cast<shape_dim+1, world_dim>());

        Tiny::Vector<DT_, world_dim> test(coords[0]);
        test += coeffs(0)*(coords[1] - coords[0]);

        for(int i(0); i < world_dim; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(test(i), points[np](i), tol);
      }

      // The last point is different because it is not in the edge and we cannot reconstruct it, but we can check if
      // the distance is correct
      Trafo::Standard::inverse_mapping(coeffs, points[4], coords.template size_cast<shape_dim+1, world_dim>());

      TEST_CHECK_EQUAL_WITHIN_EPS(coeffs(shape_dim), DT_(0.18), tol);

    }

    // Runs tests for a Simplex<2> in 3d
    void run_3d_tria() const
    {
      static constexpr int num_points = 5;
      static constexpr int world_dim = 3;
      static constexpr int shape_dim = 2;
      // Tolerance
      DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.8)));

      Tiny::Vector<DT_, world_dim> coeffs(DT_(0));
      // The last column of coords will contain the orthogonal to the triangle immersed in 3d
      Tiny::Matrix<DT_, world_dim+1, world_dim > coords;

      coords(0,0) = -DT_(1);
      coords(0,1) = -DT_(0.3);
      coords(0,2) = -DT_(0.3);

      coords(1,0) = DT_(2);
      coords(1,1) = -DT_(0.2);
      coords(1,2) = DT_(0.3);

      coords(2,0) = DT_(2);
      coords(2,1) = DT_(6);
      coords(2,2) = DT_(0.3);

      Tiny::Matrix<DT_, world_dim, world_dim-1> tmp_coords(DT_(0));

      for(int i(0); i < world_dim; ++i)
      {
        for(int j(0); j < world_dim-1; ++j)
          tmp_coords(i,j) = coords(j+1,i) - coords(0,i);
      }

      coords[world_dim] = Tiny::orthogonal(tmp_coords);
      coords[world_dim].normalise();

      Tiny::Matrix<DT_, num_points, world_dim> points(DT_(0));
      // Point inside
      points[0] = DT_(0.53)*coords[0] + DT_(0.17)*coords[1] + DT_(0.3)*coords[2];
      // One of the vertices
      points[1] = coords[2];
      // Point outside, coefficients are not barycentric coordinates because sum != 1
      points[2] = DT_(-1)*coords[0] + DT_(1.25)*coords[1] + DT_(0.75)*coords[2];
      // Nearly one of the vertices
      points[3] = (DT_(1)+DT_(10)*tol)*coords[1] + (-DT_(10)*tol)*coords[2];
      // With component in orthogonal direction
      points[4] = DT_(0.1)*coords[0] - DT_(2.5)*coords[1] + DT_(0.2)*coords[2] + DT_(0.45)*coords[3];

      for(int np(0); np < num_points; ++np)
      {
        // We size-cast away the last column since it is just for result checking, so the correct version gets called
        Trafo::Standard::inverse_mapping(coeffs, points[np], coords.template size_cast<shape_dim+1, world_dim>());

        Tiny::Vector<DT_, world_dim> test(coords[0]);
        for(int i(0); i < world_dim-1; ++i)
          test+= coeffs(i)*(coords[i+1] - coords[0]);

        // The last column is different
        test += coeffs(world_dim-1)*coords[world_dim];

        for(int i(0); i < world_dim; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(test(i), points[np](i), tol);

      }
    }
};

InverseMappingTest<float> btf;
InverseMappingTest<double> btd;
#ifdef FEAT_HAVE_QUADMATH
InverseMappingTest<__float128> btq;
#endif // FEAT_HAVE_QUADMATH

/*
 * \brief Test for cell volume computations
 *
 * \test Verfies the routine's output against pre computed values
 *
 * \author Jordi Paul
 * */
template<typename DataType_>
class StandardTrafoVolumeTest
: public TestSystem::TaggedTest<Archs::None, DataType_>
{
  public:
    StandardTrafoVolumeTest() :
      TestSystem::TaggedTest<Archs::None, DataType_>("standard_trafo_volume_test")
      {
      }

    virtual ~StandardTrafoVolumeTest()
    {
    }

    virtual void run() const override
    {
      test_1d_quad();
      test_1d_simplex();
      test_2d_quad();
      test_2d_simplex();
      test_3d_quad();
      test_3d_simplex();
    }
    void test_1d_quad() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Hypercube<1> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      static const CoordType vtx[2*1] = { CoordType(1), Math::sqrt(CoordType(2))};
      Geometry::TestAux::copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      DataType_ vol = trafo.compute_vol<ShapeType,DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(DataType_(2)) - DataType_(1), eps);
    }

    void test_1d_simplex() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Simplex<1> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      static const CoordType vtx[2*1] = { CoordType(1), Math::sqrt(CoordType(2))};
      Geometry::TestAux::copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      DataType_ vol = trafo.compute_vol<ShapeType,DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(DataType_(2)) - DataType_(1), eps);
    }

    void test_2d_simplex() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Simplex<2> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      // Sufficiently wild simplex
      static const CoordType vtx[3*2] =
      { -CoordType(1), -CoordType(1),
        CoordType(1), CoordType(1)/CoordType(2),
        -CoordType(1)/CoordType(2), CoordType(3)/CoordType(4)};
      Geometry::TestAux::copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // Everything checked against has been computed by hand
      DataType_ vol = trafo.compute_vol<ShapeType, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(11)/DataType_(8), eps);
      // Check volume of sub simplices
      vol = trafo.compute_vol<Shape::Simplex<1>, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(DataType_(37)/DataType_(16)), eps);
      vol = trafo.compute_vol<Shape::Simplex<1>, DataType_>(Index(1));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(DataType_(53)/DataType_(16)), eps);
      vol = trafo.compute_vol<Shape::Simplex<1>, DataType_>(Index(2));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(5)/DataType_(2), eps);
    }

    void test_2d_quad() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Hypercube<2> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      // Sufficiently wild hypercube
      static const CoordType vtx[4*2] =
      { -Math::sqrt(CoordType(2)), -CoordType(3),
        CoordType(2), -CoordType(1),
        -CoordType(1), CoordType(1),
        CoordType(5)/CoordType(2), Math::sqrt(CoordType(3))};
      Geometry::TestAux::copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // Everything checked against has been computed by hand
      DataType_ vol = trafo.compute_vol<ShapeType, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(17025)/DataType_(1546), eps);

      static const DataType_ l[4] = {
        DataType_(10735)/DataType_(2713), DataType_(10788)/DataType_(3017),
        DataType_(22749)/DataType_(5657), DataType_(26405)/DataType_(9507)};

      // Check volume of sub elements
      for(Index i(0); i < 4; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<1>,DataType_>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, l[i], eps);
      }
    }

    void test_3d_quad() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Hypercube<3> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      /* In 3d. the volume for nontrivial hypercubes (especially those with nonplanar faces) is difficult to compute
       * by hand, so this is just a parallel piped
       */
      static const CoordType vtx[8*3] =
      { -CoordType(2), -CoordType(2), -CoordType(2),
        CoordType(2), -CoordType(2), -CoordType(2),
        -CoordType(2), CoordType(2), -CoordType(2),
        CoordType(2), CoordType(2), -CoordType(2),
        -CoordType(1), -CoordType(2), CoordType(2),
        CoordType(3), -CoordType(2), CoordType(2),
        -CoordType(1), CoordType(2), CoordType(2),
        CoordType(3), CoordType(2), CoordType(2) };
      Geometry::TestAux::copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // Everything checked against has been computed by hand
      DataType_ vol = trafo.compute_vol<ShapeType, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(64), eps);

      // Check the 2d volume of the faces
      static const DataType_ f[6] =
      { DataType_(16), DataType_(16),
        DataType_(16), DataType_(16),
        DataType_(4)*Math::sqrt(DataType_(17)), DataType_(4)*Math::sqrt(DataType_(17))};

      for(Index i(0); i < 6; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<2>, DataType_>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, f[i], eps);
      }

      // Check the 1d volume of the edges
      static const DataType_ l[12] = {
        DataType_(4), DataType_(4), DataType_(4), DataType_(4),
        DataType_(4), DataType_(4), DataType_(4), DataType_(4),
        Math::sqrt(DataType_(17)), Math::sqrt(DataType_(17)),
        Math::sqrt(DataType_(17)), Math::sqrt(DataType_(17))};

      for(Index i(0); i < 12; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<1>, DataType_>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, l[i], eps);
      }
    }

    void test_3d_simplex() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));
      typedef Shape::Simplex<3> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      /* Our simplex is the standard simplex scaled by a factor of 2 in x_3 direction, and vertex 3 has been moved
         which does not alter the volume of 2*1/d! = 1/3
         */
      static const CoordType vtx[4*3] =
      { CoordType(0), CoordType(0), CoordType(0),
        CoordType(1), CoordType(0), CoordType(0),
        CoordType(0), CoordType(1), CoordType(0),
        CoordType(0), CoordType(0.5), CoordType(2), };
      Geometry::TestAux::copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      DataType_ vol = trafo.compute_vol<ShapeType, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(1)/DataType_(3), eps);

      /* Edge lengths computed by hand */
      static const DataType_ l[6] =
      {DataType_(1), DataType_(1),
        Math::sqrt(DataType_(17)/DataType_(4)), Math::sqrt(DataType_(2)),
        Math::sqrt(DataType_(21)/DataType_(4)), Math::sqrt(DataType_(17)/DataType_(4))};
      for(Index i=0; i<6; i++)
      {
        vol = trafo.compute_vol<Shape::Simplex<1>, DataType_>(Index(i));
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, l[i], eps);
      }

      /* With the edgelengths, check the volume of the sub simplices via Heron's formula */
      // s will be half the circumference of a sub simplex
      DataType_ s(0);

      // Face 0 consists of edges 3, 4, 5
      s = DataType_(0.5)*(l[3] + l[4] + l[5]);
      vol = trafo.compute_vol<Shape::Simplex<2>, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(s * (s - l[3]) * (s - l[4]) * (s - l[5]) ), eps);
      // Face 1 consists of edges 1, 2, 5
      s = DataType_(0.5)*(l[1] + l[2] + l[5]);
      vol = trafo.compute_vol<Shape::Simplex<2>, DataType_>(Index(1));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(s * (s - l[1]) * (s - l[2]) * (s - l[5]) ), eps);
      // Face 2 consists of edges 0, 2, 4
      s = DataType_(0.5)*(l[0] + l[2] + l[4]);
      vol = trafo.compute_vol<Shape::Simplex<2>, DataType_>(Index(2));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(s * (s - l[0]) * (s - l[2]) * (s - l[4]) ), eps);
      // Face 3 consists of edges 0, 1, 3
      s = DataType_(0.5)*(l[0] + l[1] + l[3]);
      vol = trafo.compute_vol<Shape::Simplex<2>, DataType_>(Index(3));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(s * (s - l[0]) * (s - l[1]) * (s - l[3]) ), eps);
    }
};

StandardTrafoVolumeTest<float> standard_trafo_volume_test_float;
StandardTrafoVolumeTest<double> standard_trafo_volume_test_double;
