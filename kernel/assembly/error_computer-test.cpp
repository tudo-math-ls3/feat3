// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
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

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename DT_, typename IT_>
class ErrorComputerTest :
  public UnitTest
{
public:
  static constexpr int dim = 2;
  typedef DT_ DataType;
  typedef IT_ IndexType;

  typedef LAFEM::DenseVector<DataType, IndexType> ScalarVectorType;
  typedef LAFEM::DenseVectorBlocked<DataType, IndexType, dim> BlockedVectorType;


  ErrorComputerTest(PreferredBackend backend) :
    UnitTest("ErrorComputerTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    test_fe_error(IT_(3), Math::pow(Math::eps<DT_>(), DT_(0.7)), Math::pow(Math::eps<DT_>(), DT_(0.7)), Math::pow(Math::eps<DT_>(), DT_(0.7)));
  }

  void test_fe_error(Index level, DataType val_tol, DataType grad_tol, DataType hess_tol) const
  {
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, dim, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

    String cubature("gauss-legendre:5");

    // create mesh, trafo and space
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);
    TrafoType trafo(mesh);
    SpaceType space(trafo);

    // create and interpolate a scalar function
    typedef Analytic::Common::SineBubbleFunction<dim> ScalarFuncType;
    typedef Analytic::Gradient<ScalarFuncType> VectorFuncType;

    ScalarFuncType scalar_func;
    VectorFuncType vector_func(scalar_func);

    ScalarVectorType scalar_vec;
    BlockedVectorType blocked_vec;

    Assembly::Interpolator::project(scalar_vec, scalar_func, space);
    Assembly::Interpolator::project(blocked_vec, vector_func, space);

    Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
    domain_assembler.compile_all_elements();
    {
      // Now create our intgeral jobs
      Assembly::DomainAssemblyErrorFunctionIntegralJob<ScalarFuncType, ScalarVectorType, SpaceType, 0> base_err_job(scalar_func, scalar_vec, space, cubature);
      Assembly::CellErrorFunctionIntegralJob<ScalarFuncType, ScalarVectorType, SpaceType, 0> cell_err_job(scalar_func, scalar_vec, space, cubature);

      domain_assembler.assemble(base_err_job);
      domain_assembler.assemble(cell_err_job);

      auto result_base = base_err_job.result();
      auto results_cell = cell_err_job.result();

      TEST_CHECK_EQUAL_WITHIN_EPS(result_base.value, results_cell.integral_info.value, val_tol);

      DataType summed_h0{DataType(0)};
      for(std::size_t i = 0; i < results_cell.vec.size(); ++i)
      {
        auto loc_vec = results_cell.vec(i);
        summed_h0 += loc_vec;
      }

      TEST_CHECK_EQUAL_WITHIN_EPS(result_base.norm_h0_sqr, summed_h0, val_tol);
    }
    {
      // Now create our intgeral jobs
      Assembly::DomainAssemblyErrorFunctionIntegralJob<ScalarFuncType, ScalarVectorType, SpaceType, 2> base_err_job(scalar_func, scalar_vec, space, cubature);
      Assembly::CellErrorFunctionIntegralJob<ScalarFuncType, ScalarVectorType, SpaceType, 2> cell_err_job(scalar_func, scalar_vec, space, cubature);

      domain_assembler.assemble(base_err_job);
      domain_assembler.assemble(cell_err_job);

      auto result_base = base_err_job.result();
      auto results_cell = cell_err_job.result();

      TEST_CHECK_EQUAL_WITHIN_EPS(result_base.value, results_cell.integral_info.value, val_tol);

      DataType summed_h0{DataType(0)}, summed_h1{DataType(0)}, summed_h2{DataType(0)};
      for(std::size_t i = 0; i < results_cell.vec.size(); ++i)
      {
        auto loc_vec = results_cell.vec(i);
        summed_h0 += loc_vec[0];
        summed_h1 += loc_vec[1];
        summed_h2 += loc_vec[2];
      }

      TEST_CHECK_EQUAL_WITHIN_EPS(result_base.norm_h0_sqr, summed_h0, val_tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(result_base.norm_h1_sqr, summed_h1, grad_tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(result_base.norm_h2_sqr, summed_h2, hess_tol);
    }
    {
      // Now create our intgeral jobs
      Assembly::DomainAssemblyErrorFunctionIntegralJob<VectorFuncType, BlockedVectorType, SpaceType, 0> base_err_job(vector_func, blocked_vec, space, cubature);
      Assembly::CellErrorFunctionIntegralJob<VectorFuncType, BlockedVectorType, SpaceType, 0> cell_err_job(vector_func, blocked_vec, space, cubature);

      domain_assembler.assemble(base_err_job);
      domain_assembler.assemble(cell_err_job);

      auto result_base = base_err_job.result();
      auto results_cell = cell_err_job.result();

      TEST_CHECK_EQUAL_WITHIN_EPS(result_base.value[0], results_cell.integral_info.value[0], val_tol);

      DataType summed_h0{DataType(0)};
      for(std::size_t i = 0; i < results_cell.vec.size(); ++i)
      {
        auto loc_vec = results_cell.vec(i);
        summed_h0 += loc_vec;
      }

      TEST_CHECK_EQUAL_WITHIN_EPS(result_base.norm_h0_sqr, summed_h0, val_tol);
    }
    // {
    //   // Now create our intgeral jobs
    //   Assembly::ErrorFunctionIntegralJob<VectorFuncType, BlockedVectorType, SpaceType, 2> base_err_job(vector_func, blocked_vec, space, cubature);
    //   Assembly::CellErrorFunctionIntegralJob<VectorFuncType, BlockedVectorType, SpaceType, 2> cell_err_job(vector_func, blocked_vec, space, cubature);

    //   domain_assembler.assemble(base_err_job);
    //   // domain_assembler.assemble(cell_err_job);

    //   auto result_base = base_err_job.result();
    //   auto results_cell = cell_err_job.result();

    //   TEST_CHECK_EQUAL_WITHIN_EPS(result_base.value[0], results_cell.integral_info.value[0], val_tol);

    //   DataType summed_h0{DataType(0)}, summed_h1{DataType(0)}, summed_h2{DataType(0)};
    //   for(std::size_t i = 0; i < results_cell.vec.size(); ++i)
    //   {
    //     auto loc_vec = results_cell.vec(i);
    //     summed_h0 += loc_vec[0];
    //     summed_h1 += loc_vec[1];
    //     summed_h2 += loc_vec[2];
    //   }

    //   TEST_CHECK_EQUAL_WITHIN_EPS(result_base.norm_h0_sqr, summed_h0, val_tol);
    //   TEST_CHECK_EQUAL_WITHIN_EPS(result_base.norm_h1_sqr, summed_h1, grad_tol);
    //   TEST_CHECK_EQUAL_WITHIN_EPS(result_base.norm_h2_sqr, summed_h2, hess_tol);
    // }

  }
}; // class DiscreteEvaluatorTest<...>

ErrorComputerTest <double, std::uint32_t> discrete_evaluator_test_main_double_uint32(PreferredBackend::generic);
ErrorComputerTest <float, std::uint32_t> discrete_evaluator_test_main_float_uint32(PreferredBackend::generic);
ErrorComputerTest <double, std::uint64_t> discrete_evaluator_test_main_double_uint64(PreferredBackend::generic);
ErrorComputerTest <float, std::uint64_t> discrete_evaluator_test_main_float_uint64(PreferredBackend::generic);
