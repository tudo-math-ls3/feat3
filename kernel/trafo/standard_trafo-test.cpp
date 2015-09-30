#include <test_system/test_system.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/standard/inverse_mapping.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/shape_convert_factory.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

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

  struct UnitTrafoConfig : public Trafo::ConfigBase
  {
    enum
    {
      need_dom_point = 1,
      need_img_point = 1,
      need_jac_mat = 1,
      need_jac_inv = 1,
      need_jac_det = 1,
      need_hess_ten = 1,
      need_hess_inv = 1
    };
  };


public:
  StandardTrafoTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("StandardTrafoTest")
  {
  }

  virtual void run() const
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
    typename TrafoEvaluator::template ConfigTraits<UnitTrafoConfig>::EvalDataType trafo_data;

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
    typename TrafoEvaluator::template ConfigTraits<UnitTrafoConfig>::EvalDataType trafo_data;

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
    typename TrafoEvaluator::template ConfigTraits<UnitTrafoConfig>::EvalDataType trafo_data;

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
    typename TrafoEvaluator::template ConfigTraits<UnitTrafoConfig>::EvalDataType trafo_data;

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
    typename TrafoEvaluator::template ConfigTraits<UnitTrafoConfig>::EvalDataType trafo_data;

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
    typename TrafoEvaluator::template ConfigTraits<UnitTrafoConfig>::EvalDataType trafo_data;

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
    typename TrafoEvaluator::template ConfigTraits<UnitTrafoConfig>::EvalDataType trafo_data;

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
    typename TrafoEvaluator::template ConfigTraits<UnitTrafoConfig>::EvalDataType trafo_data;

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

    virtual void run() const
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

        //std::cout << "coeffs = ";
        //for(int i(0); i < world_dim; ++i)
        //  std::cout << " " << scientify(coeffs(i));
        //std::cout << std::endl;

        //std::cout << "org = ";
        //for(int i(0); i < world_dim; ++i)
        //  std::cout << " " << scientify(points[np](i));
        //std::cout << std::endl;

        //std::cout << "new = ";
        //for(int i(0); i < world_dim; ++i)
        //  std::cout << " " << scientify(test(i));
        //std::cout << std::endl;

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
      DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.8)));

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

        //std::cout << "coeffs = ";
        //for(int i(0); i < shape_dim+1; ++i)
        //  std::cout << " " << scientify(coeffs(i));
        //std::cout << std::endl;

        //std::cout << "org = ";
        //for(int i(0); i < world_dim; ++i)
        //  std::cout << " " << scientify(points[np](i));
        //std::cout << std::endl;

        //std::cout << "new = ";
        //for(int i(0); i < world_dim; ++i)
        //  std::cout << " " << scientify(test(i));
        //std::cout << std::endl;

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

        //std::cout << "coeffs = ";
        //for(int i(0); i < world_dim; ++i)
        //  std::cout << " " << scientify(coeffs(i));
        //std::cout << std::endl;

        //std::cout << "org = ";
        //for(int i(0); i < world_dim; ++i)
        //  std::cout << " " << scientify(points[np](i));
        //std::cout << std::endl;

        //std::cout << "new = ";
        //for(int i(0); i < world_dim; ++i)
        //  std::cout << " " << scientify(test(i));
        //std::cout << std::endl;

        for(int i(0); i < world_dim; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(test(i), points[np](i), tol);

      }
    }


};

InverseMappingTest<float> btf;
InverseMappingTest<double> btd;
#ifdef FEAST_HAVE_QUADMATH
InverseMappingTest<__float128> btq;
#endif // FEAST_HAVE_QUADMATH
