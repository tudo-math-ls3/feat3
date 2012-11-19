#include <test_system/test_system.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>

#include <limits>
#include <cmath>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief Discontinuous Element test
 *
 * \test Tests the Discontinuous Finite-Element space
 *
 * \tparam DataType_
 * The data type for the test. Shall be either double or float.
 *
 * \author Peter Zajac
 */
template<typename DataType_>
class DiscontinuousTest
  : public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;

  typedef Cubature::Rule<ShapeType, DataType_, DataType_, Tiny::Vector<DataType_, 2> > CubatureRule;
  typedef typename Cubature::DynamicFactorySelect<CubatureRule>::Type CubatureFactory;

  struct UnitSpaceConfig : public Space::ConfigBase
  {
    enum
    {
      need_value = 1
    };
  };

  struct UnitTrafoConfig : public QuadSpaceQ0::TrafoConfig<UnitSpaceConfig>
  {
    enum
    {
      need_jac_det = 1
    };
  };


public:
  DiscontinuousTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("Discontinuous Test")
  {
  }

  virtual void run() const
  {
    // test assembly on unit quad
    asm_unit_quad_q0();
  }

  void asm_unit_quad_q0() const
  {
    // compute eps
    const DataType_ eps = std::pow(std::numeric_limits<DataType_>::epsilon(), DataType_(0.8));

    // create a quad mesh
    Geometry::UnitCubeFactory<QuadMesh> mesh_factory;
    QuadMesh mesh(mesh_factory);

    // create a quad-trafo
    QuadTrafo trafo(mesh);

    // create a Q0 space
    QuadSpaceQ0 space(trafo);

    // create a trafo evaluator
    typedef typename QuadTrafo::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);
    Trafo::EvalData<TrafoEvaluator, UnitTrafoConfig> trafo_data;

    // create a space evaluator
    typedef typename QuadSpaceQ0::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;
    SpaceEvaluator space_eval(space);
    Space::EvalData<SpaceEvaluator, UnitSpaceConfig> space_data;

    // create a 2x2 Gauss-Legendre cubature formula
    CubatureRule cubature_rule(CubatureFactory::create("gauss-legendre:2"));

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    // prepare space evaluator
    space_eval.prepare(trafo_eval);

    // check the number of local DOFs
    Index num_loc_dofs = space_eval.get_num_local_dofs();
    TEST_CHECK_EQUAL(num_loc_dofs, 1u);

    // create local matrix assembly data
    DataType_ M(0);

    // loop over all 4 quadrature points and integrate
    for(Index k(0); k < cubature_rule.get_num_points(); ++k)
    {
      // compute trafo data
      trafo_data(trafo_eval, cubature_rule.get_point(k));

      // compute space data
      space_data(space_eval, trafo_data);

      // mass matrix entry
      M += trafo_data.jac_det * cubature_rule.get_weight(k) * space_data.values[0] * space_data.values[0];
    }

    // finish evaluators
    space_eval.finish();
    trafo_eval.finish();

    TEST_CHECK_EQUAL_WITHIN_EPS(M, DataType_(1), eps);
  }
};

DiscontinuousTest<double> discontinuous_test_double;
DiscontinuousTest<float> discontinuous_test_float;
