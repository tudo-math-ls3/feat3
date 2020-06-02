// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/common.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Analytic;

template<typename DT_>
class CommonFunctionTest :
  public FullTaggedTest<Mem::Main, DT_, Index>
{
public:
  CommonFunctionTest() :
    FullTaggedTest<Mem::Main, DT_, Index>("CommonFunctionTest")
  {
  }

  void test_par_profile_scalar() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create scalar parabolic profile function
    Analytic::Common::ParProfileScalar<DT_> pprof;
    TEST_CHECK(pprof.parse("(1 2, 3 4, 5)"));

    // evaluate endpoints and midpoint
    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(pprof, DT_(1), DT_(2)), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(pprof, DT_(3), DT_(4)), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(pprof, DT_(2), DT_(3)), DT_(5), tol);
  }

  void test_par_profile_vector() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create vector parabolic profile function
    Analytic::Common::ParProfileVector<DT_> pprof;
    TEST_CHECK(pprof.parse("(1 2, 3 4, 5)"));

    // evaluate endpoints and midpoint
    auto v_0 = Analytic::eval_value_x(pprof, DT_(1), DT_(2));
    auto v_1 = Analytic::eval_value_x(pprof, DT_(3), DT_(4));
    auto v_c = Analytic::eval_value_x(pprof, DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(v_0[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v_0[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v_1[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v_1[1], DT_(0), tol);
    const DT_ x_c = DT_(5) * Math::sqrt(DT_(0.5));
    TEST_CHECK_EQUAL_WITHIN_EPS(v_c[0], +x_c, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v_c[1], -x_c, tol);
  }

  virtual void run() const override
  {
    test_par_profile_scalar();
    test_par_profile_vector();
  }
};

CommonFunctionTest<double> common_function_test_double;
