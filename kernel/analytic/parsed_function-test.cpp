// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/parsed_function.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Analytic;

#ifdef FEAT_HAVE_FPARSER

template<typename DT_>
class ParsedFunctionTest :
  public UnitTest
{
private:
  const DT_ tol;

public:
  explicit ParsedFunctionTest(DT_ tol_) :
    UnitTest("ParsedFunctionTest", Type::Traits<DT_>::name(), "none", PreferredBackend::generic),
    tol(tol_)
  {
  }

  virtual void run() const override
  {
    test_parsed_scalar_function();
    test_parsed_vector_function();
    test_inline_variables_1d();
    test_inline_variables_2d();
    combined_test_inline_variables_add_variable_add_const();

  }

  void test_parsed_scalar_function() const
  {
    ParsedScalarFunction<2> psf;
    psf.add_variable("t", DT_(2.0));
    psf.add_constant("c", DT_(3.0));
    psf.parse("c*(t*x + y)");

    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(psf, DT_(1.0), DT_(-1.0)), DT_(3.0), tol);

    psf.set_variable("t", DT_(4.0));

    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(psf, DT_(1.0), DT_(-1.0)), DT_(9.0), tol);
  }

  void test_parsed_vector_function()const
  {
    ParsedVectorFunction<3, 2> pvf;
    pvf.add_variable("t", DT_(2.0));
    pvf.add_constant("c", DT_(3.0));
    pvf.parse("[c*(t*x+y) ' t*x+c*z]");

    Tiny::Vector<DT_, 2> v1 = Analytic::eval_value_x(pvf, DT_(1.0), DT_(-1.0), DT_(0.2));
    TEST_CHECK_EQUAL_WITHIN_EPS(v1[0], DT_(3.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v1[1], DT_(2.6), tol);

    pvf.set_variable("t", DT_(4.0));

    Tiny::Vector<DT_, 2> v2 = Analytic::eval_value_x(pvf, DT_(1.0), DT_(-1.0), DT_(0.2));
    TEST_CHECK_EQUAL_WITHIN_EPS(v2[0], DT_(9.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v2[1], DT_(4.6), tol);
  }

  void test_inline_variables_1d()const
  {
    ParsedScalarFunction<2> pf;
    pf.parse("length :=sqrt(x*x + y*y); 2*length*sin(length)");
    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(pf, DT_(3.0), DT_(4.0)), DT_(10.0) * Math::sin(DT_(5.0)), tol);
  }

  void test_inline_variables_2d()const
  {
    //in this test we use the trafo from  Cartesian coordinates to polar coordinates and also the reversed Trafo
    ParsedVectorFunction<2, 2>pvf_2d;
    pvf_2d.parse("[r:= sqrt(x*x+y*y); phi:=atan(y/x);r*cos(phi)'r:= sqrt(x*x+y*y); phi:=atan(y/x);r*sin(phi)]");
    Tiny::Vector<DT_, 2> v3 = Analytic::eval_value_x(pvf_2d, DT_(5.0), DT_(6.5));

    TEST_CHECK_EQUAL_WITHIN_EPS(v3[0], DT_(5.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v3[1], DT_(6.5), tol);
  }

  void combined_test_inline_variables_add_variable_add_const() const
  {
    //in this test we calculate the alteration of the inner energy under the assumption the gas is ideal
    ParsedScalarFunction<2> psf;
    psf.add_constant("R", DT_(8.3145));
    psf.add_variable("n", DT_(3.)/DT_(2.));
    psf.parse("T1:=x+273.15;T2:=y+273.15;n*R*(T2-T1)");

    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(psf, DT_(250.), DT_(270.)), DT_(20.)*DT_(1.5)* DT_(8.3145), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(psf, DT_(270.), DT_(250.)), DT_(-20.) * DT_(1.5) * DT_(8.3145), tol);
  }

};

ParsedFunctionTest <double> parsed_function_test_double(1E-13);
ParsedFunctionTest <float> parsed_function_test_float(1E-5f);
#ifdef FEAT_HAVE_QUADMATH
ParsedFunctionTest <__float128> parsed_function_test_float128(1E-13);
#endif

#endif // FEAT_HAVE_FPARSER
