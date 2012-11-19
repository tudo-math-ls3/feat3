#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/solver_expression.hpp>
#include <kernel/archs.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::ScaRC;


template<typename Tag_, typename DataType_>
class SolverExpressionTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    SolverExpressionTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("SolverExpressionTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      VGVector a("a");
      VGVector b("b");
      VGVector c(a + b);
      TEST_CHECK_EQUAL(c.get_id(), "a + b");

      VGScalar s1("s1"), s2("s2");
      VGScalar s3(s1 + s2);
      TEST_CHECK_EQUAL(s3.get_id(), "s1 + s2");

      VGBool b1("b1"), b2("b2");
      VGBool b3(b1 < b2);
      TEST_CHECK_EQUAL(b3.get_id(), "b1 < b2");

      VGMatrix A("A");
      VGVector d(defect_expr(b, A, a));
      TEST_CHECK_EQUAL(d.get_id(), "DEFECT(b, A, a)");

      d = defect_expr(b, A, c);
      TEST_CHECK_EQUAL(d.get_id(), "DEFECT(b, A, a + b)");

      s3 = normvec_expr(d);
      TEST_CHECK_EQUAL(s3.get_id(), "NORMVEC(DEFECT(b, A, a + b))");

      Vector a_i(a);
      TEST_CHECK_EQUAL(a_i.get_id(), "[a]_chunk");

      a = vec_synch_expr(a_i);
      TEST_CHECK_EQUAL(a.get_id(), "SYNCHVEC([a]_chunk)");

      a = vec_iterate_expr(a, b3);
      TEST_CHECK_EQUAL(a.get_id(), "ITERATE(SYNCHVEC([a]_chunk) UNTIL b1 < b2)");
    }
};
SolverExpressionTest<Mem::Main, double> sf_cpu_double("StorageType: std::vector, DataType: double");
