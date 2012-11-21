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
      VGBool b3(s1 < s2);
      TEST_CHECK_EQUAL(b3.get_id(), "s1 < s2");

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
      TEST_CHECK_EQUAL(a.get_id(), "ITERATE(SYNCHVEC([a]_chunk) UNTIL s1 < s2)");

      a = vec_preconapply_expr(a);
      TEST_CHECK_EQUAL(a.get_id(), "PRECONAPPLY(ITERATE(SYNCHVEC([a]_chunk) UNTIL s1 < s2))");

      a = A * d;
      TEST_CHECK_EQUAL(a.get_id(), "A * DEFECT(b, A, a + b)");

      //--------------------------------------------------------------------------------------

      VGVector gx("gx"), gb("gb");
      VGMatrix gA("gA"), gP("gP");

      Matrix lP(gP), lA(gA);
      Vector lx(gx), lb(gb);

      VGBool converged_outer("converged_outer");
      VGBool converged_inner("converged_inner");

      //global solver
      gx = vec_iterate_expr(vec_preconapply_expr(vec_synch_expr(defect_expr(lb, lA, lx))) + lx, converged_outer);
      TEST_CHECK_EQUAL(gx.get_id(), "ITERATE(PRECONAPPLY(SYNCHVEC(DEFECT([gb]_chunk, [gA]_chunk, [gx]_chunk))) + [gx]_chunk UNTIL converged_outer)");

      //local jacobi
      lx = vec_iterate_expr(lx + lP * defect_expr(lb, lA, lx), converged_inner);
      TEST_CHECK_EQUAL(lx.get_id(), "ITERATE([gx]_chunk + [gP]_chunk * DEFECT([gb]_chunk, [gA]_chunk, [gx]_chunk) UNTIL converged_inner)");

      //gx stays the old value, solvers not connected
      TEST_CHECK_EQUAL(gx.get_id(), "ITERATE(PRECONAPPLY(SYNCHVEC(DEFECT([gb]_chunk, [gA]_chunk, [gx]_chunk))) + [gx]_chunk UNTIL converged_outer)");

      //lx is substituted, solvers connected, do not do this
      gx = vec_iterate_expr(vec_preconapply_expr(vec_synch_expr(defect_expr(lb, lA, lx))) + lx, converged_outer);
      TEST_CHECK_EQUAL(gx.get_id(), "ITERATE(PRECONAPPLY(SYNCHVEC(DEFECT([gb]_chunk, [gA]_chunk, ITERATE([gx]_chunk + [gP]_chunk * DEFECT([gb]_chunk, [gA]_chunk, [gx]_chunk) UNTIL converged_inner)))) + ITERATE([gx]_chunk + [gP]_chunk * DEFECT([gb]_chunk, [gA]_chunk, [gx]_chunk) UNTIL converged_inner) UNTIL converged_outer)");

      //------------------------------------------------------------

      //test inner-outer schemes

      //start over
      VGVector gx1("gx"), gb1("gb"), ld1("ld");
      VGMatrix gA1("gA"), gP1("gP");

      Matrix lP1(gP1), lA1(gA1);
      Vector lx1(gx1), lb1(gb1);

      VGScalar initial_defect_norm_outer("gamma_out");
      VGScalar tolerance("tol");
      converged_outer = (initial_defect_norm_outer / normvec_expr(defect_expr(lb1, lA1, lx1))) < tolerance;

      //generate outermost solver
      ld1 = defect_expr(lb1, lA1, lx1);
      VGVector scarc_precon_0(vec_preconapply_expr(ld1)); //entry-point for inner solver
      lx1 = lx1 + scarc_precon_0;
      lx1 = vec_synch_expr(lx1);
      lx1 = vec_iterate_expr(lx1, converged_outer);
      gx1 = lx1;
      TEST_CHECK_EQUAL(gx1.get_id(), "ITERATE(SYNCHVEC([gx]_chunk + PRECONAPPLY(DEFECT([gb]_chunk, [gA]_chunk, [gx]_chunk))) UNTIL gamma_out / NORMVEC(DEFECT([gb]_chunk, [gA]_chunk, [gx]_chunk)) < tol)");

      //different objects, same handles
      VGVector gx2("gx"), gb2("gb"), ld2("ld");
      VGMatrix gA2("gA"), gP2("gP");

      Matrix lP2(gP2), lA2(gA2);
      Vector lx2(gx2), lb2(gb2);

      //generate innermost solver
      ld2 = defect_expr(lb2, lA2, lx2);
      VGVector scarc_precon_1(lP2 * ld2); //real local preconditioner
      lx2 = lx2 + scarc_precon_1;
      lx2 = vec_iterate_expr(lx2, converged_inner);
      gx2 = lx2;
      TEST_CHECK_EQUAL(gx2.get_id(), "ITERATE([gx]_chunk + [gP]_chunk * DEFECT([gb]_chunk, [gA]_chunk, [gx]_chunk) UNTIL converged_inner)");

    }
};
SolverExpressionTest<Mem::Main, double> sf_cpu_double("StorageType: std::vector, DataType: double");
