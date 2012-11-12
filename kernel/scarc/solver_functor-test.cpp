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
      //-----------------------------------
      //create layer n-1
      std::shared_ptr<ScaRC::MatrixData> A(new ScaRC::MatrixData);
      std::shared_ptr<ScaRC::VectorData> b(new ScaRC::VectorData);
      std::shared_ptr<ScaRC::VectorData> x(new ScaRC::VectorData);

      std::shared_ptr<FunctorBase> pr(new ScaRC::ProxyDefect<ScaRC::VectorData, ScaRC::MatrixData, ScaRC::VectorData>(b, A, x));
      std::shared_ptr<ScaRC::ProxyPreconApply> p(new ScaRC::ProxyPreconApply(pr));

      std::shared_ptr<ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply> > richardson(new ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply>(x, p));

      TEST_CHECK_EQUAL(richardson.get()->type_name(), "ProxyVector + (__defect__(ProxyVector,ProxyMatrix,ProxyVector))");

      //create layer 0
      std::shared_ptr<ScaRC::MatrixData> A1(new ScaRC::MatrixData);
      std::shared_ptr<ScaRC::VectorData> b1(new ScaRC::VectorData);
      std::shared_ptr<ScaRC::VectorData> x1(new ScaRC::VectorData);

      ///only in case, where the master layer is created last, we need to reinterpret the shared ptr because we MUST use its copy-ctor
      std::shared_ptr<FunctorBase> pr1(*(reinterpret_cast<std::shared_ptr<FunctorBase>* >(&richardson)));

      std::shared_ptr<ScaRC::ProxyPreconApply> p1(new ScaRC::ProxyPreconApply(pr1));

      std::shared_ptr<ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply> > richardson1(new ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply>(x1, p1));

      TEST_CHECK_EQUAL(richardson1.get()->type_name(), "ProxyVector + (ProxyVector + (__defect__(ProxyVector,ProxyMatrix,ProxyVector)))");

      //----------------------------------
      //create layer 0
      ///this is the common case
      std::shared_ptr<ScaRC::MatrixData> A2(new ScaRC::MatrixData);
      std::shared_ptr<ScaRC::VectorData> b2(new ScaRC::VectorData);
      std::shared_ptr<ScaRC::VectorData> x2(new ScaRC::VectorData);

      std::shared_ptr<ScaRC::ProxyPreconApply> p2(new ScaRC::ProxyPreconApply);

      std::shared_ptr<ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply> > richardson2(new ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply>(x2, p2));

      TEST_CHECK_EQUAL(richardson2.get()->type_name(), "ProxyVector + (__UNINITIALIZED_PRECONDITIONER_APPLICATION__())");

      //create layer n-1
      std::shared_ptr<FunctorBase> pr2(new ScaRC::ProxyDefect<ScaRC::VectorData, ScaRC::MatrixData, ScaRC::VectorData>(b2, A2, x2));
      std::shared_ptr<ScaRC::ProxyPreconApply> p3(new ScaRC::ProxyPreconApply(pr2));

      p2.get()->get() = std::shared_ptr<FunctorBase>(new ScaRC::ProxyVectorSum<ScaRC::VectorData, ScaRC::ProxyPreconApply>(x2, p3));
      TEST_CHECK_EQUAL(richardson2.get()->type_name(), "ProxyVector + (ProxyVector + (__defect__(ProxyVector,ProxyMatrix,ProxyVector)))");
    }
};
SolverFunctorTest<Archs::CPU, double> sf_cpu_double("StorageType: std::vector, DataType: double");
