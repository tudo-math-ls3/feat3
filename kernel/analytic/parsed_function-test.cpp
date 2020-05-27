// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/parsed_function.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Analytic;

#ifdef FEAT_HAVE_FPARSER

template<typename DT_>
class ParsedFunctionTest :
  public FullTaggedTest<Mem::Main, DT_, Index>
{
public:
  ParsedFunctionTest() :
    FullTaggedTest<Mem::Main, DT_, Index>("ParsedFunctionTest")
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));

    ParsedScalarFunction<2> psf;
    psf.add_variable("t", 2.0);
    psf.add_constant("c", 3.0);
    psf.parse("c*(t*x + y)");

    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(psf, 1.0, -1.0), 3.0, tol);

    psf.set_variable("t", 4.0);

    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(psf, 1.0, -1.0), 9.0, tol);


    ParsedVectorFunction<3,2> pvf;
    pvf.add_variable("t", 2.0);
    pvf.add_constant("c", 3.0);
    pvf.parse("[c*(t*x+y) ' t*x+c*z]");

    Tiny::Vector<DT_,2> v1 = Analytic::eval_value_x(pvf, 1.0, -1.0, 0.2);
    TEST_CHECK_EQUAL_WITHIN_EPS(v1[0], 3.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v1[1], 2.6, tol);

    pvf.set_variable("t", 4.0);

    Tiny::Vector<DT_,2> v2 = Analytic::eval_value_x(pvf, 1.0, -1.0, 0.2);
    TEST_CHECK_EQUAL_WITHIN_EPS(v2[0], 9.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v2[1], 4.6, tol);
  }
};

ParsedFunctionTest<double> parsed_function_test_double;

#endif // FEAT_HAVE_FPARSER
