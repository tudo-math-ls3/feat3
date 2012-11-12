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
      //-----------------------------------
      //create layer n-1
      ScaRC::MatrixData A;
      ScaRC::VectorData b, x;

      std::shared_ptr<FunctorBase> pr(new ScaRC::ProxyDefect<ScaRC::VectorData, ScaRC::MatrixData, ScaRC::VectorData>(b, A, x));
      ScaRC::ProxyPreconApply p(pr);

      ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply> richardson(x, p);

      TEST_CHECK_EQUAL(richardson.type_name(), "ProxyVector + __precon__(__defect__(ProxyVector,ProxyMatrix,ProxyVector))");

      //create layer 0
      ScaRC::MatrixData A1;
      ScaRC::VectorData b1, x1;

      std::shared_ptr<FunctorBase> pr1(new ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply>(richardson));
      ScaRC::ProxyPreconApply p1(pr1);

      ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply> richardson1(x1, p1);

      TEST_CHECK_EQUAL(richardson1.type_name(), "ProxyVector + __precon__(ProxyVector + __precon__(__defect__(ProxyVector,ProxyMatrix,ProxyVector)))");

      //----------------------------------
      //create layer 0
      ScaRC::MatrixData A2;
      ScaRC::VectorData b2, x2;

      ScaRC::ProxyPreconApply p2;

      ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply> richardson2(x2, p2);

      TEST_CHECK_EQUAL(richardson2.type_name(), "ProxyVector + __precon__(__UNINITIALIZED_PRECONDITIONER__())");

      //create layer n-1
      std::shared_ptr<FunctorBase> pr2(new ScaRC::ProxyDefect<ScaRC::VectorData, ScaRC::MatrixData, ScaRC::VectorData>(b2, A2, x2));
      ScaRC::ProxyPreconApply p3(pr2);

      p2.get() = std::shared_ptr<FunctorBase>(new ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply>(x2, p3));
      TEST_CHECK_EQUAL(richardson2.type_name(), "ProxyVector + __precon__(ProxyVector + __precon__(__defect__(ProxyVector,ProxyMatrix,ProxyVector)))");
    }
};
SolverFunctorTest<Archs::CPU, double> sf_cpu_double("StorageType: std::vector, DataType: double");
