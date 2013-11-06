#include <test_system/test_system.hpp>
#include <kernel/trafo/standard/mapping.hpp>
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
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

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

    dom_point[0] = DataType_(0.0);
    dom_point[1] = DataType_(0.0);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqrt(Math::sqr(trafo_data.img_point[0] - DataType_(0.5))
                           + Math::sqr(trafo_data.img_point[1] - DataType_(0.5)));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_(0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_(0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_(0.5), eps);

    // check hessian tensor
    for(Index i(0); i < 2; ++i)
    {
      for(Index j(0); j < 2; ++j)
      {
        for(Index k(0); k < 2; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0.0), eps);
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
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

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

    mesh.get_vertex_set()[0][0] = DataType_(-0.25);
    mesh.get_vertex_set()[0][1] = DataType_( 0.25);

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(0.0);
    dom_point[1] = DataType_(0.0);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqrt(Math::sqr(trafo_data.img_point[0] - DataType_(0.4375))
                           + Math::sqr(trafo_data.img_point[1] - DataType_(0.5625)));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_( 0.5625), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_( 0.0625), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(-0.0625), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_( 0.4375), eps);

    // check hessian tensor
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][0], DataType_(0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][1], DataType_(-0.0625), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][0], DataType_(-0.0625), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][1], DataType_(0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][0], DataType_(0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][1], DataType_(0.0625), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][0], DataType_(0.0625), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][1], DataType_(0.0), eps);

    trafo_eval.finish();
  }// test_quad()

  void test_unit_cube() const
  {
    typedef Shape::Hexahedron ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute eps
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

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

    dom_point[0] = DataType_(0.0);
    dom_point[1] = DataType_(0.0);
    dom_point[2] = DataType_(0.0);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqrt(Math::sqr(trafo_data.img_point[0] - DataType_(0.5))
                           + Math::sqr(trafo_data.img_point[1] - DataType_(0.5))
                           + Math::sqr(trafo_data.img_point[2] - DataType_(0.5)));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][2], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][2], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][1], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][2], DataType_( 0.5), eps);

    // check hessian tensor
    for(Index i(0); i < 3; ++i)
    {
      for(Index j(0); j < 3; ++j)
      {
        for(Index k(0); k < 3; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0.0), eps);
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
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

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

    mesh.get_vertex_set()[7][0] = DataType_(0.9);
    mesh.get_vertex_set()[7][1] = DataType_(0.8);
    mesh.get_vertex_set()[7][2] = DataType_(0.75);

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(0.0);
    dom_point[1] = DataType_(0.0);
    dom_point[2] = DataType_(0.0);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqrt(Math::sqr(trafo_data.img_point[0] - DataType_(0.4875))
                           + Math::sqr(trafo_data.img_point[1] - DataType_(0.475))
                           + Math::sqr(trafo_data.img_point[2] - DataType_(0.46875)));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_( 0.4875), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_(-0.0125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][2], DataType_(-0.0125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(-0.025), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_( 0.475), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][2], DataType_(-0.025), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][0], DataType_(-0.03125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][1], DataType_(-0.03125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][2], DataType_( 0.46875), eps);

    // check hessian tensor
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][1], DataType_(-0.0125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][0][2], DataType_(-0.0125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][0], DataType_(-0.0125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][1], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][1][2], DataType_(-0.0125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][2][0], DataType_(-0.0125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][2][1], DataType_(-0.0125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[0][2][2], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][1], DataType_(-0.025), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][0][2], DataType_(-0.025), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][0], DataType_(-0.025), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][1], DataType_(0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][1][2], DataType_(-0.025), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][2][0], DataType_(-0.025), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][2][1], DataType_(-0.025), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[1][2][2], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][0][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][0][1], DataType_(-0.03125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][0][2], DataType_(-0.03125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][1][0], DataType_(-0.03125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][1][1], DataType_(0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][1][2], DataType_(-0.03125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][2][0], DataType_(-0.03125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][2][1], DataType_(-0.03125), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[2][2][2], DataType_(0.0), eps);

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
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

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

    dom_point[0] = DataType_(0.5);
    dom_point[1] = DataType_(0.5);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqrt(Math::sqr(trafo_data.img_point[0] - DataType_(0.5))
                           + Math::sqr(trafo_data.img_point[1] - DataType_(0.0)));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_(-0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(-0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_(-0.5), eps);

    // check hessian tensor
    for(Index i(0); i < 2; ++i)
    {
      for(Index j(0); j < 2; ++j)
      {
        for(Index k(0); k < 2; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0.0), eps);
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
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

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

    mesh.get_vertex_set()[1][0] = DataType_( 0.1);
    mesh.get_vertex_set()[1][1] = DataType_(-0.1);

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(0.25);
    dom_point[1] = DataType_(0.25);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqrt(Math::sqr(trafo_data.img_point[0] - DataType_(0.275))
                           + Math::sqr(trafo_data.img_point[1] - DataType_(0.225)));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_(-0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_(-0.4), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_(-0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_(-0.6), eps);

    // check hessian tensor
    for(Index i(0); i < 2; ++i)
    {
      for(Index j(0); j < 2; ++j)
      {
        for(Index k(0); k < 2; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0.0), eps);
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
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

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

    dom_point[0] = DataType_(0.25);
    dom_point[1] = DataType_(0.25);
    dom_point[2] = DataType_(0.25);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqrt(Math::sqr(trafo_data.img_point[0] - DataType_(0.5))
                           + Math::sqr(trafo_data.img_point[1] - DataType_(0.25))
                           + Math::sqr(trafo_data.img_point[2] - DataType_(0.125)));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_( 1.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][2], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][2], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][1], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][2], DataType_( 0.5), eps);

    // check hessian tensor
    for(Index i(0); i < 3; ++i)
    {
      for(Index j(0); j < 3; ++j)
      {
        for(Index k(0); k < 3; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0.0), eps);
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
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

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

    mesh.get_vertex_set()[8][0] = DataType_(0.4);
    mesh.get_vertex_set()[8][1] = DataType_(0.3);
    mesh.get_vertex_set()[8][2] = DataType_(0.25);

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    dom_point[0] = DataType_(0.25);
    dom_point[1] = DataType_(0.25);
    dom_point[2] = DataType_(0.25);

    // evaluate trafo
    trafo_eval(trafo_data, dom_point);

    // check image point coordinates
    DataType_ x = Math::sqrt(Math::sqr(trafo_data.img_point[0] - DataType_(0.475))
                           + Math::sqr(trafo_data.img_point[1] - DataType_(0.2))
                           + Math::sqr(trafo_data.img_point[2] - DataType_(0.1875)));
    TEST_CHECK_EQUAL_WITHIN_EPS(x, DataType_(0), eps);

    // check jacobian matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][0], DataType_( 1.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][1], DataType_( 0.4), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[0][2], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][1], DataType_( 0.3), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[1][2], DataType_( 0.5), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][0], DataType_( 0.0), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][1], DataType_( 0.25), eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.jac_mat[2][2], DataType_( 0.5), eps);

    // check hessian tensor
    for(Index i(0); i < 3; ++i)
    {
      for(Index j(0); j < 3; ++j)
      {
        for(Index k(0); k < 3; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(trafo_data.hess_ten[i][j][k], DataType_(0.0), eps);
        }
      }
    }

    trafo_eval.finish();
  } // test_tetra()
};

StandardTrafoTest<double> standard_trafo_test_double;
StandardTrafoTest<float> standard_trafo_test_float;
