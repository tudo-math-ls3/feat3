#include <test_system/test_system.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief Ranancher-Turek Element test
 *
 * \test Tests the Ranancher-Turek Finite-Element space
 *
 * \tparam DataType_
 * The data type for the test. Shall be either double or float.
 *
 * \author Peter Zajac
 */
template<typename DataType_>
class RannacherTurekTest
  : public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceStdNonPar;

  typedef Cubature::Rule<ShapeType, DataType_, DataType_, Tiny::Vector<DataType_, 2> > CubatureRule;

  struct UnitTrafoConfig : public Trafo::ConfigBase
  {
    enum
    {
      need_dom_point = 1,
      need_img_point = 1,
      need_jac_mat = 1,
      need_jac_det = 1
    };
  };

  struct UnitSpaceConfig : public Space::ConfigBase
  {
    enum
    {
      need_value = 1,
      need_grad = 1
    };
  };

public:
  RannacherTurekTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("Rannacher-Turek Test")
  {
  }

  virtual void run() const
  {
    // test assembly on unit quad
    asm_unit_quad_std_non_par();
  }

  void asm_unit_quad_std_non_par() const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

    // create a quad mesh
    Geometry::UnitCubeFactory<QuadMesh> mesh_factory;
    QuadMesh mesh(mesh_factory);

    // create a quad-trafo
    QuadTrafo trafo(mesh);

    // create a Q1~ space, standard non-parametric
    QuadSpaceStdNonPar space(trafo);

    // create a trafo evaluator
    typedef typename QuadTrafo::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);
    Trafo::EvalData<typename TrafoEvaluator::EvalTraits, UnitTrafoConfig> trafo_data;

    // create a space evaluator
    typedef typename QuadSpaceStdNonPar::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;
    SpaceEvaluator space_eval(space);
    Space::EvalData<typename SpaceEvaluator::SpaceEvalTraits, UnitSpaceConfig> space_data;

    // create a 3x3 Gauss-Legendre cubature formula
    //CubatureRule cubature_rule;
    //CubatureFactory::create(cubature_rule, "gauss-legendre:3");
    CubatureRule cubature_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:3"));

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    // prepare space evaluator
    space_eval.prepare(trafo_eval);

    // check the number of local DOFs
    Index num_loc_dofs = space_eval.get_num_local_dofs();
    TEST_CHECK_EQUAL(num_loc_dofs, 4u);

    // create local matrix assembly data
    Tiny::Matrix<DataType_, 4, 4> L, M;
    L = DataType_(0);
    M = DataType_(0);

    // loop over all 4 quadrature points and integrate
    for(Index k(0); k < cubature_rule.get_num_points(); ++k)
    {
      // compute trafo data
      trafo_eval(trafo_data, cubature_rule.get_point(k));

      // compute space data
      space_eval(space_data, trafo_data);

      // test function loop
      for(Index i(0); i < num_loc_dofs; ++i)
      {
        // trial function loop
        for(Index j(0); j < num_loc_dofs; ++j)
        {
          // mass matrix entry
          M(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) * (
            space_data.phi[i].value * space_data.phi[j].value);

          // laplace matrix entry
          L(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) * (
            space_data.phi[i].grad[0] * space_data.phi[j].grad[0] +
            space_data.phi[i].grad[1] * space_data.phi[j].grad[1]);
          // continue with next trial function
        }
        // continue with next test function
      }
      // continue with next cubature point
    }

    // finish evaluators
    space_eval.finish();
    trafo_eval.finish();

    // test function loop
    for(Index i(0); i < num_loc_dofs; ++i)
    {
      // trial function loop
      for(Index j(0); j < num_loc_dofs; ++j)
      {
        // check entries
        if(i == j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(41) / DataType_(240), eps); // = 1.708p3e-1
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(5) / DataType_(2), eps); // = 2.5
        }
        else if((int(i >> 1)) == int(j >> 1)) // i-j-pairs: (0-1) and (2-3)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(5) / DataType_(1200), eps); // = 4.1p6e-3
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(1) / DataType_(2), eps); // = 0.5
        }
        else
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(3) / DataType_(80), eps); // = 3.75e-2
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(-3) / DataType_(2), eps); // = -1.5
        }
        // continue with next trial function
      }
      // continue with next test function
    }
  }
};

RannacherTurekTest<double> rannacher_turek_test_double;
RannacherTurekTest<float> rannacher_turek_test_float;
