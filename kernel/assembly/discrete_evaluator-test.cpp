// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/analytic/function.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/util/random.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename MT_, typename DT_, typename IT_>
class DiscreteEvaluatorTest :
  public TestSystem::FullTaggedTest<MT_, DT_, IT_>
{
public:
  static constexpr int dim = 2;
  typedef MT_ MemType;
  typedef DT_ DataType;
  typedef IT_ IndexType;

  typedef LAFEM::DenseVector<MemType, DataType, IndexType> ScalarVectorType;
  typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> BlockedVectorType;


  DiscreteEvaluatorTest() :
    TestSystem::FullTaggedTest<MT_, DT_, IT_>("DiscreteEvaluatorTest")
  {
  }

  virtual void run() const override
  {
    // Note:
    // The tolerances need to be large, as we are checking analytic
    // function values against discrete function values!
    test_fe_eval(3, 1E-3, 1E-2);
  }

  void test_fe_eval(Index level, DataType val_tol, DataType grad_tol) const
  {
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, dim, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Trafo::InverseMapping<TrafoType, DataType> InvMappingType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

    // create mesh, trafo and space
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);
    TrafoType trafo(mesh);
    SpaceType space(trafo);

    // create an inverse mapping
    InvMappingType inv_mapping(trafo);

    // create and interpolate a scalar function
    typedef Analytic::Common::SineBubbleFunction<dim> ScalarFuncType;
    typedef Analytic::Gradient<ScalarFuncType> VectorFuncType;

    ScalarFuncType scalar_func;
    VectorFuncType vector_func(scalar_func);

    typedef Analytic::EvalTraits<DataType, ScalarFuncType> ScalarEvalTraits;
    typedef Analytic::EvalTraits<DataType, VectorFuncType> VectorEvalTraits;

    typedef typename ScalarEvalTraits::PointType PointType;
    typedef typename ScalarEvalTraits::ValueType ScalarValueType;
    typedef typename VectorEvalTraits::ValueType VectorValueType;

    typename ScalarFuncType::template Evaluator<ScalarEvalTraits> scalar_func_eval(scalar_func);
    typename VectorFuncType::template Evaluator<VectorEvalTraits> vector_func_eval(vector_func);

    ScalarVectorType scalar_vec;
    BlockedVectorType blocked_vec;

    Assembly::Interpolator::project(scalar_vec, scalar_func, space);
    Assembly::Interpolator::project(blocked_vec, vector_func, space);

    // create an RNG
    Random rng;

    // test a few points
    for(int k(0); k < 10; ++k)
    {
      // create a random point on the unit-square
      PointType point;
      for(int i(0); i < dim; ++i)
        point[i] = rng(DataType(0.0001), DataType(0.9999));

      // unmap point
      auto inv_point = inv_mapping.unmap_point(point);

      // ensure that the point was unmapped
      TEST_CHECK(!inv_point.empty());

      // evaluate both functions
      auto scalar_eval_data = Assembly::DiscreteEvaluator::eval_fe_function(inv_point, scalar_vec, space);
      auto vector_eval_data = Assembly::DiscreteEvaluator::eval_fe_function(inv_point, blocked_vec, space);

      // ensure that the data structures are not empty
      TEST_CHECK(!scalar_eval_data.empty());
      TEST_CHECK(!vector_eval_data.empty());

      // compute mean values
      ScalarValueType scalar_value_fe = scalar_eval_data.mean_value();
      VectorValueType vector_value_fe = vector_eval_data.mean_value();

      // evaluate the analytic functions
      ScalarValueType scalar_value_ana = scalar_func_eval.value(point);
      VectorValueType vector_value_ana = vector_func_eval.value(point);

      // compute differences
      ScalarValueType scalar_diff = scalar_value_fe - scalar_value_ana;
      VectorValueType vector_diff = vector_value_fe - vector_value_ana;

      // compute errors
      DataType scalar_err = Math::abs(scalar_diff);
      DataType vector_err = vector_diff.norm_euclid();

      //std::cout << point << " > " << scalar_err << " | " << vector_err << std::endl;

      // check errors against tolerances
      TEST_CHECK(scalar_err < val_tol);
      TEST_CHECK(vector_err < grad_tol);
    }
  }
}; // class DiscreteEvaluatorTest<...>

DiscreteEvaluatorTest<Mem::Main, double, Index> discrete_evaluator_test_main_double_index;
