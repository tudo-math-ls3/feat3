#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/solver_functor.hpp>
#include <kernel/scarc/solver_pattern.hpp>
#include <kernel/util/cpp11_smart_pointer.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;


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
      std::shared_ptr<ScaRC::MatrixData> A3(new ScaRC::MatrixData);
      std::shared_ptr<ScaRC::VectorData> x3(new ScaRC::VectorData);
      std::shared_ptr<ScaRC::VectorData> b3(new ScaRC::VectorData);
      std::shared_ptr<FunctorBase> p3_master;
      std::shared_ptr<FunctorBase> p3_local;

      std::shared_ptr<FunctorBase> scarc(ScaRC::SolverPatternGeneration<ScaRC::Richardson>::execute(x3, p3_master));
      ((ScaRC::ProxyPreconApply*)p3_master.get())->get() = ScaRC::SolverPatternGeneration<ScaRC::Richardson>::execute(A3, x3, b3, p3_local);

      TEST_CHECK_EQUAL(scarc.get()->type_name(), "__precon__(__precon__(__defect__(ProxyVector,ProxyMatrix,ProxyVector)) + ProxyVector) + ProxyVector");
    }
};
SolverPatternTest<Archs::CPU, double> sf_cpu_double("StorageType: std::vector, DataType: double");
