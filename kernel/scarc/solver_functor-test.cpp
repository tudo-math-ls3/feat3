#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/solver_functor.hpp>
#include <kernel/util/cpp11_smart_pointer.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_, typename DataType_>
class SolverFunctorTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    SolverFunctorTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("SolverFunctorTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      ScaRC::MatrixData A;
      ScaRC::VectorData b, x;

      ScaRC::ProxyDefect<ScaRC::VectorData, ScaRC::MatrixData, ScaRC::VectorData> r(b, A, x);
      ScaRC::ProxyPreconApply<ScaRC::ProxyDefect<ScaRC::VectorData, ScaRC::MatrixData, ScaRC::VectorData> > p(r);
      ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply<ScaRC::ProxyDefect<ScaRC::VectorData, ScaRC::MatrixData, ScaRC::VectorData> > > richardson(x, p);

      TEST_CHECK_EQUAL(richardson.type_name(), "ProxyVector + __precon__(__defect__(ProxyVector,ProxyMatrix,ProxyVector))");
    }
};
SolverFunctorTest<Archs::CPU, double> sf_cpu_double("StorageType: std::vector, DataType: double");
