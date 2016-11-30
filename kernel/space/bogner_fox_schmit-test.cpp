#include <test_system/test_system.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/bogner_fox_schmit/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/adjacency/graph.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename DataType_>
class BognerFoxSchmitTest
  : public TestSystem::TaggedTest<Archs::None, DataType_>
{
  static constexpr TrafoTags unit_trafo_config = TrafoTags::jac_det | TrafoTags::jac_inv;

  static constexpr SpaceTags unit_space_config = SpaceTags::value | SpaceTags::grad;

public:
  BognerFoxSchmitTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("BognerFoxSchmit Test")
  {
  }

  virtual ~BognerFoxSchmitTest()
  {
  }

  virtual void run() const override
  {
    // test assembly on 1D unit interval
    asm_unit_interval();
  }

  void asm_unit_interval() const
  {
    typedef Shape::Hypercube<1> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::BognerFoxSchmit::Element<TrafoType> SpaceType;
    typedef Cubature::Rule<ShapeType, DataType_, DataType_, Tiny::Vector<DataType_, 1> > CubatureRule;

    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a single line mesh over [-1,1]
    Geometry::UnitCubeFactory<MeshType> mesh_factory;
    MeshType mesh(mesh_factory);
    mesh.get_vertex_set()[0][0] = -1.0;

    // create a quad-trafo
    TrafoType trafo(mesh);

    // create a BFS space
    SpaceType space(trafo);

    // create a trafo evaluator
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);
    typename TrafoEvaluator::template ConfigTraits<unit_trafo_config>::EvalDataType trafo_data;

    // create a space evaluator
    typedef typename SpaceType::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;
    SpaceEvaluator space_eval(space);
    typename SpaceEvaluator::template ConfigTraits<unit_space_config>::EvalDataType space_data;

    // create a 4x4 Gauss-Legendre cubature formula
    CubatureRule cubature_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:4"));

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    // prepare space evaluator
   space_eval.prepare(trafo_eval);

    // check the number of local DOFs
    int num_loc_dofs = space_eval.get_num_local_dofs();
    TEST_CHECK_EQUAL(num_loc_dofs, Index(4));

    // create local matrix assembly data
    Tiny::Matrix<DataType_, 4, 4> L, M, Lr, Mr;
    L = DataType_(0);
    M = DataType_(0);

    // loop over all 16 quadrature points and integrate
    for(int k(0); k < cubature_rule.get_num_points(); ++k)
    {
      // compute trafo data
      trafo_eval(trafo_data, cubature_rule.get_point(k));

      // compute space data
      space_eval(space_data, trafo_data);

      // test function loop
      for(int i(0); i < num_loc_dofs; ++i)
      {
        // trial function loop
        for(int j(0); j < num_loc_dofs; ++j)
        {
          // mass matrix entry
          M(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) * (
            space_data.phi[i].value * space_data.phi[j].value);

          // laplace matrix entry
          L(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) * (
            space_data.phi[i].grad[0] * space_data.phi[j].grad[0]);
          // continue with next trial function
        }
        // continue with next test function
      }
      // continue with next cubature point
    }

    // finish evaluators
    space_eval.finish();
    trafo_eval.finish();

    // set up analytic matrices
    Mr(0,0) =  DataType_(26) / DataType_( 35);
    Mr(0,1) =  DataType_(22) / DataType_(105);
    Mr(0,2) =  DataType_( 9) / DataType_( 35);
    Mr(0,3) = -DataType_(13) / DataType_(105);
    Mr(1,0) =  DataType_(22) / DataType_(105);
    Mr(1,1) =  DataType_( 8) / DataType_(105);
    Mr(1,2) =  DataType_(13) / DataType_(105);
    Mr(1,3) = -DataType_( 2) / DataType_( 35);
    Mr(2,0) =  DataType_( 9) / DataType_( 35);
    Mr(2,1) =  DataType_(13) / DataType_(105);
    Mr(2,2) =  DataType_(26) / DataType_( 35);
    Mr(2,3) = -DataType_(22) / DataType_(105);
    Mr(3,0) = -DataType_(13) / DataType_(105);
    Mr(3,1) = -DataType_( 2) / DataType_( 35);
    Mr(3,2) = -DataType_(22) / DataType_(105);
    Mr(3,3) =  DataType_( 8) / DataType_(105);
    Lr(0,0) =  DataType_(3) / DataType_( 5);
    Lr(0,1) =  DataType_(1) / DataType_(10);
    Lr(0,2) = -DataType_(3) / DataType_( 5);
    Lr(0,3) =  DataType_(1) / DataType_(10);
    Lr(1,0) =  DataType_(1) / DataType_(10);
    Lr(1,1) =  DataType_(4) / DataType_(15);
    Lr(1,2) = -DataType_(1) / DataType_(10);
    Lr(1,3) = -DataType_(1) / DataType_(15);
    Lr(2,0) = -DataType_(3) / DataType_( 5);
    Lr(2,1) = -DataType_(1) / DataType_(10);
    Lr(2,2) =  DataType_(3) / DataType_( 5);
    Lr(2,3) = -DataType_(1) / DataType_(10);
    Lr(3,0) =  DataType_(1) / DataType_(10);
    Lr(3,1) = -DataType_(1) / DataType_(15);
    Lr(3,2) = -DataType_(1) / DataType_(10);
    Lr(3,3) =  DataType_(4) / DataType_(15);

    // test function loop
    for(int i(0); i < num_loc_dofs; ++i)
    {
      // trial function loop
      for(int j(0); j < num_loc_dofs; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), Mr(i,j), eps);
        TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), Lr(i,j), eps);
        // continue with next trial function
      }
      // continue with next test function
    }
  }
};

BognerFoxSchmitTest<double> bogner_fox_schmit_test_double;
BognerFoxSchmitTest<float> bogner_fox_schmit_test_float;
