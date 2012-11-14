#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/solver_functor.hpp>
#include <kernel/scarc/solver_pattern.hpp>
#include <kernel/util/cpp11_smart_pointer.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;
using namespace FEAST::ScaRC;


template<typename Tag_, typename DataType_>
class SolverPatternTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    SolverPatternTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("SolverPatternTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      std::shared_ptr<MatrixData> A(new MatrixData("A_omega_i"));
      std::shared_ptr<MatrixData> P(new MatrixData("P_omega_i"));
      std::shared_ptr<VectorData> x(new VectorData("x_omega_i"));
      std::shared_ptr<VectorData> b(new VectorData("b_omega_i"));
      std::shared_ptr<VectorData> c(new VectorData("c_omega_i"));

      ///suppose our solver is a defect correction scheme
      std::shared_ptr<ProxyDefect<VectorData, MatrixData, VectorData> > defect(new ProxyDefect<VectorData, MatrixData, VectorData>(b, A, x));

      ///suppose our solver is preconditioned locally by Approximate Inverse multiplies
      std::shared_ptr<ProxyMatrixVectorProduct<MatrixData, VectorData> > precon(SolverOperatorGeneration<ApproximateInverseMultiply>::execute(P, c));

      CompoundFunctor<> solver_layers;

      solver_layers.add_functor(std::shared_ptr<FunctorBase>());
      solver_layers.add_functor(std::shared_ptr<FunctorBase>());

      std::shared_ptr<FunctorBase> scarc(SolverPatternGeneration<Richardson>::execute(x, solver_layers.get_functors().at(0)));
      ((ProxyPreconApply*)solver_layers.get_functors().at(0).get())->get() = SolverPatternGeneration<Richardson>::execute(A, c, defect, precon, solver_layers.get_functors().at(1));

      TEST_CHECK_EQUAL(scarc.get()->type_name(), "[PRECON([PRECON([RICHARDSON([MatrixData A_omega_i], [VectorData c_omega_i], [DEFECT([VectorData b_omega_i],[MatrixData A_omega_i],[VectorData x_omega_i])], [MatrixData P_omega_i] * [VectorData c_omega_i])])])] + [VectorData x_omega_i]");
      TEST_CHECK_EQUAL(solver_layers.get_functors().at(0).get()->type_name(), "[PRECON([PRECON([RICHARDSON([MatrixData A_omega_i], [VectorData c_omega_i], [DEFECT([VectorData b_omega_i],[MatrixData A_omega_i],[VectorData x_omega_i])], [MatrixData P_omega_i] * [VectorData c_omega_i])])])]");
      TEST_CHECK_EQUAL(solver_layers.get_functors().at(1).get()->type_name(), "[MatrixData P_omega_i] * [VectorData c_omega_i]");
    }
};
SolverPatternTest<Mem::Main, double> sf_cpu_double("StorageType: std::vector, DataType: double");
