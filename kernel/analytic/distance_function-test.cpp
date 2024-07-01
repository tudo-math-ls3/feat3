// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/distance_function.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Analytic;

template<typename DT_>
class DistanceFunctionTest :
  public UnitTest
{
public:
  DistanceFunctionTest() :
    UnitTest("DistanceFunctionTest", Type::Traits<DT_>::name(), "none", PreferredBackend::generic)
  {
  }

  void test_distance_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 2d DistanceFunction with origin in (0,2).
    typename Analytic::Distance::DistanceFunction<2, DT_>::PointType orig2d;
    orig2d[0] = DT_(0);
    orig2d[1] = DT_(2);
    Analytic::Distance::DistanceFunction<2, DT_> dfunc2d(orig2d);

    // evaluate in (4, -1)
    auto val2d = Analytic::eval_value_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(val2d, DT_(5), tol);

    auto grad2d = Analytic::eval_gradient_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad2d[0], DT_( 4)/DT_(5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad2d[1], DT_(-3)/DT_(5), tol);

    auto hess2d = Analytic::eval_hessian_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[0][0], DT_( 9)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[0][1], DT_(12)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[1][0], DT_(12)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[1][1], DT_(16)/DT_(125), tol);
  }

  void test_distance_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 3d DistanceFunction with origin in (2,2,1).
    typename Analytic::Distance::DistanceFunction<3, DT_>::PointType orig3d;
    orig3d[0] = DT_(2);
    orig3d[1] = DT_(2);
    orig3d[2] = DT_(1);
    Analytic::Distance::DistanceFunction<3, DT_> dfunc3d(orig3d);

    // evaluate in (4, 0, 2)
    auto val3d = Analytic::eval_value_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val3d, DT_(3), tol);

    auto grad3d = Analytic::eval_gradient_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[0], DT_( 2)/DT_(3), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[1], DT_(-2)/DT_(3), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[2], DT_( 1)/DT_(3), tol);

    auto hess3d = Analytic::eval_hessian_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][0], DT_( 5)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][1], DT_( 4)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][2], DT_(-2)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][0], DT_( 4)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][1], DT_( 5)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][2], DT_( 2)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][0], DT_(-2)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][1], DT_( 2)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][2], DT_( 8)/DT_(27), tol);
  }

  void test_distance_function_sd_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 2d DistanceFunction with origin in (0,2).
    typename Analytic::Distance::DistanceFunctionSD<2, DT_>::PointType orig2d;
    orig2d[0] = DT_(0);
    orig2d[1] = DT_(2);
    Analytic::Distance::DistanceFunctionSD<2, DT_> dfunc2d(orig2d, DT_(0.25), 10);

    // evaluate in (4, -1)
    auto val2d = Analytic::eval_value_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(val2d, DT_(50.25), tol);

    auto grad2d = Analytic::eval_gradient_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad2d[0], DT_( 8), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad2d[1], DT_(-6), tol);

    auto hess2d = Analytic::eval_hessian_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[0][0], DT_( 90)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[0][1], DT_(120)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[1][0], DT_(120)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[1][1], DT_(160)/DT_(125), tol);
  }

  void test_distance_function_sd_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 3d DistanceFunction with origin in (2,2,1).
    typename Analytic::Distance::DistanceFunctionSD<3, DT_>::PointType orig3d;
    orig3d[0] = DT_(2);
    orig3d[1] = DT_(2);
    orig3d[2] = DT_(1);
    Analytic::Distance::DistanceFunctionSD<3, DT_> dfunc3d(orig3d, DT_(4), DT_(3));

    // evaluate in (4, 0, 2)
    auto val3d = Analytic::eval_value_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val3d, DT_(13), tol);

    auto grad3d = Analytic::eval_gradient_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[0], DT_( 2), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[1], DT_(-2), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[2], DT_( 1), tol);

    auto hess3d = Analytic::eval_hessian_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][0], DT_( 5)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][1], DT_( 4)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][2], DT_(-2)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][0], DT_( 4)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][1], DT_( 5)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][2], DT_( 2)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][0], DT_(-2)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][1], DT_( 2)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][2], DT_( 8)/DT_(9), tol);
  }

  void test_plane_distance_function_sd() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 3d DistanceFunction with origin in (2,2,1).
    typename Analytic::Distance::PlaneDistanceFunctionSD<1, 3, DT_>::PointType orig3d;
    orig3d[0] = DT_(2);
    orig3d[1] = DT_(2);
    orig3d[2] = DT_(1);
    Analytic::Distance::PlaneDistanceFunctionSD<1, 3, DT_> dfunc3d(orig3d, DT_(3));

    auto val3d = Analytic::eval_value_x(dfunc3d, DT_(4), DT_(7), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val3d, DT_(15), tol);

    auto grad3d = Analytic::eval_gradient_x(dfunc3d, DT_(4), DT_(-3), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[1], DT_(-3), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[2], DT_(0), tol);
  }

#ifdef FEAT_HAVE_CGAL
  void test_cgal_signed_dist_3d() const
  {
    // CGAL uses double internally, so we cannot expect higher precision even with __float128
    const DT_ tol = Math::pow(Math::max(Math::eps<DT_>(), DT_(1E-16)), DT_(0.7));
    std::stringstream mts;
    mts<<"OFF"<<std::endl;
    mts<<"8 12 18"<<std::endl;
    mts<<"0.0 0.0 0.0"<<std::endl;
    mts<<"1.0 0.0 0.0"<<std::endl;
    mts<<"0.0 1.0 0.0"<<std::endl;
    mts<<"1.0 1.0 0.0"<<std::endl;
    mts<<"0.0 0.0 1.0"<<std::endl;
    mts<<"1.0 0.0 1.0"<<std::endl;
    mts<<"0.0 1.0 1.0"<<std::endl;
    mts<<"1.0 1.0 1.0"<<std::endl;
    mts<<"3 0 1 2"<<std::endl;
    mts<<"3 1 3 2"<<std::endl;
    mts<<"3 0 5 1"<<std::endl;
    mts<<"3 0 4 5"<<std::endl;
    mts<<"3 2 4 0"<<std::endl;
    mts<<"3 4 2 6"<<std::endl;
    mts<<"3 3 1 5"<<std::endl;
    mts<<"3 3 5 7"<<std::endl;
    mts<<"3 2 3 7"<<std::endl;
    mts<<"3 2 7 6"<<std::endl;
    mts<<"3 5 4 7"<<std::endl;
    mts<<"3 7 4 6"<<std::endl;

    Analytic::Distance::CGALDistanceFunction<DT_> func(mts, Geometry::CGALFileMode::fm_off);

    DT_ val_1 = Analytic::eval_value_x(func, DT_(0.0), DT_(0.0), DT_(0.0));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(0), tol);
    DT_ val_2 = Analytic::eval_value_x(func, DT_(0.5), DT_(0.5), DT_(0.5));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(0.5), tol);
    DT_ val_3 = Analytic::eval_value_x(func, DT_(0.1), DT_(0.2), DT_(0.2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_3, DT_(0.1), tol);
    DT_ val_4 = Analytic::eval_value_x(func, DT_(0.3), DT_(0.2), DT_(-0.132));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_4, DT_(-0.132), tol);

  }
#endif // FEAT_HAVE_CGAL

  virtual void run() const override
  {
    test_distance_function_2d();
    test_distance_function_3d();
    test_distance_function_sd_2d();
    test_distance_function_sd_3d();
    test_plane_distance_function_sd();
#ifdef FEAT_HAVE_CGAL
    test_cgal_signed_dist_3d();
#endif // FEAT_HAVE_CGAL
  }
};

DistanceFunctionTest <double> common_function_test_double;
DistanceFunctionTest <float> common_function_test_float;
#ifdef FEAT_HAVE_QUADMATH
DistanceFunctionTest <__float128> common_function_test_float128;
#endif
#ifdef FEAT_HAVE_HALFMATH
DistanceFunctionTest <Half> common_function_test_half;
#endif
#ifdef FEAT_HAVE_CUDA
DistanceFunctionTest <float> cuda_common_function_test_float;
DistanceFunctionTest <double> cuda_common_function_test_double;
#endif
