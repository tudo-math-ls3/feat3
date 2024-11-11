// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/random.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the basic Math functions.
 *
 * \test Tests various math function templates in the Math namespace.
 *
 * \author Peter Zajac
 */
class BasicMathTest :
  public TestSystem::UnitTest
{
public:
  BasicMathTest() :
    TestSystem::UnitTest("BasicMathTest")
  {
  }

  virtual ~BasicMathTest()
  {
  }

  virtual void run() const override
  {
    test_minimaxi();
    test_factorial();
    test_binomial();
    test_ilog10();
  }

  void test_minimaxi() const
  {
    int imin(0), imax(0);

    imin = Math::min(3, 7);
    TEST_CHECK_EQUAL(imin, 3);

    imax = Math::max(3, 7);
    TEST_CHECK_EQUAL(imax, 7);

    Math::mini(imin, 5);
    TEST_CHECK_EQUAL(imin, 3);

    Math::mini(imin, 2);
    TEST_CHECK_EQUAL(imin, 2);

    Math::maxi(imax, 5);
    TEST_CHECK_EQUAL(imax, 7);

    Math::maxi(imax, 8);
    TEST_CHECK_EQUAL(imax, 8);

    Math::minimax(5, imin, imax);
    TEST_CHECK_EQUAL(imin, 2);
    TEST_CHECK_EQUAL(imax, 8);

    Math::minimax(1, imin, imax);
    TEST_CHECK_EQUAL(imin, 1);
    TEST_CHECK_EQUAL(imax, 8);

    Math::minimax(9, imin, imax);
    TEST_CHECK_EQUAL(imin, 1);
    TEST_CHECK_EQUAL(imax, 9);
  }

  void test_factorial() const
  {
    // 0! = 1 (by definition)
    TEST_CHECK_EQUAL(Math::factorial(0), 1);
    // 1! = 1
    TEST_CHECK_EQUAL(Math::factorial(1,0), 1);
    // 5! = 120
    TEST_CHECK_EQUAL(Math::factorial(5), 120);
    // 5*6*7 = 210
    TEST_CHECK_EQUAL(Math::factorial(7,5), 210);
  }

  void test_binomial() const
  {
    TEST_CHECK_EQUAL(Math::binomial(0,0), 1);
    TEST_CHECK_EQUAL(Math::binomial(3,0), 1);
    TEST_CHECK_EQUAL(Math::binomial(3,3), 1);
    TEST_CHECK_EQUAL(Math::binomial(5,2), 10);
    TEST_CHECK_EQUAL(Math::binomial(49,6), 13983816);
  }

  void test_ilog10() const
  {
    TEST_CHECK_EQUAL(Math::ilog10( 0 ), 0 );
    TEST_CHECK_EQUAL(Math::ilog10(10u), 2u);
    TEST_CHECK_EQUAL(Math::ilog10(-9l), 1l);
  }
} basic_math_test;

/**
 * \brief Test class for the Math functions.
 *
 * \test Tests the floating-point math function templates in the Math namespace.
 *
 * \author Peter Zajac
 */
template<typename DT_>
class MathTest :
  public TestSystem::UnitTest
{
public:
  const DT_ tol, tol2;

  /**
   * \brief Constructor
   *
   * \param[in] tol_
   * Tolerance for all functions except sqrt, typically around eps.
   *
   * \param[in] tol2_
   * Tolerance for sqrt, typically around sqrt(eps).
   */
  explicit MathTest(DT_ tol_, DT_ tol2_) :
    TestSystem::UnitTest("MathTest", Type::Traits<DT_>::name()),
    tol(tol_), tol2(tol2_)
  {
  }

  virtual ~MathTest()
  {
  }

  virtual void run() const override
  {
    test_sqrt();
    test_sin();
    test_cos();
    test_sinh();
    test_cosh();
    test_exp();
    test_log();
    test_log10();
    test_pow();
    test_atan();
    test_atan2();
    test_pi();
    test_asin();
    test_acos();
    test_nan();
    test_opening_angle();
  }

  void test_sqrt() const
  {
    // exact root tests
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(0)), DT_(0), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(0.25)), DT_(0.5), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(1)), DT_(1), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(64)), DT_(8), tol2);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(2)), Math::sqrt(DT_(2)), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(0.5)), Math::sqrt(DT_(0.5)), tol2);
  }

  void test_sin() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sin<DT_>(DT_(0)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sin<DT_>(DT_(0.7)), Math::sin(DT_(0.7)), tol);
  }

  void test_cos() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cos<DT_>(DT_(0)), DT_(1), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cos<DT_>(DT_(0.7)), Math::cos(DT_(0.7)), tol);
  }

  void test_sinh() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sinh<DT_>(DT_(0)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sinh<DT_>(DT_(0.7)), Math::sinh(DT_(0.7)), tol);
  }

  void test_cosh() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cosh<DT_>(DT_(0)), DT_(1), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cosh<DT_>(DT_(0.7)), Math::cosh(DT_(0.7)), tol);
  }

  void test_exp() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::exp<DT_>(DT_(0)), DT_(1), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::exp<DT_>( DT_(2)), Math::exp( DT_(2)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::exp<DT_>(-DT_(2)), Math::exp(-DT_(2)), tol);
  }

  void test_log() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log<DT_>(DT_(1)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log<DT_>(DT_(0.5)), Math::log(DT_(0.5)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log<DT_>(DT_(2)), Math::log(DT_(2)), tol);
  }

  void test_log10() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<DT_>(DT_(1)), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<DT_>(DT_(10)), DT_(1), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<DT_>(DT_(0.5)), Math::log10(DT_(0.5)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<DT_>(DT_(2)), Math::log10(DT_(2)), tol);
  }

  void test_pow() const
  {
    // test exact (positive exponents)
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(1), DT_(1)), DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(2), DT_(3)), DT_(8), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(4), DT_(0.5)), DT_(2), tol);

    // test exact (negative exponents)
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(1), DT_(1)), DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(2), DT_(-3)), DT_(0.125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(4), DT_(-0.5)), DT_(0.5), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(0.7), DT_(3.1)), Math::pow(DT_(0.7), DT_(3.1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(3.1), DT_(0.7)), Math::pow(DT_(3.1), DT_(0.7)), tol);
  }

  void test_atan() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<DT_>(DT_(0)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<DT_>( DT_(1.5)), Math::atan( DT_(1.5)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<DT_>(-DT_(1.5)), Math::atan(-DT_(1.5)), tol);
  }

  void test_atan2() const
  {
    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan2<DT_>(DT_(0.7), DT_(3.1)), Math::atan2(DT_(0.7), DT_(3.1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan2<DT_>(DT_(3.1), DT_(0.7)), Math::atan2(DT_(3.1), DT_(0.7)), tol);
  }

  void test_pi() const
  {
    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pi<DT_>(), DT_(2)*Math::acos(DT_(0)), tol);
  }

  void test_asin() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::asin<DT_>(DT_(0)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::asin<DT_>(DT_(0.7)), Math::asin(DT_(0.7)), tol);
  }

  void test_acos() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::acos<DT_>(DT_(1)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::acos<DT_>(DT_(0.7)), Math::acos(DT_(0.7)), tol);
  }

  void test_nan() const
  {
    // generate NaN
    const DT_ q = Math::nan<DT_>();

    // should be NaN
    TEST_CHECK(Math::isnan(q));

    // should not be finite
    TEST_CHECK(!Math::isfinite(q));

    // should not be normal
    TEST_CHECK(!Math::isnormal(q));
  }

  void test_opening_angle() const
  {
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::calc_opening_angle_intern(DT_(FEAT_F128C(-0.8910065241883678623597095714136263127705)), DT_(FEAT_F128C(-0.453990499739546791560408366357871198983))), DT_(FEAT_F128C(1.35))*Math::pi<DT_>(), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::calc_opening_angle_intern(DT_(FEAT_F128C(0.)), DT_(FEAT_F128C(-1.))), Math::pi<DT_>(), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::calc_opening_angle(DT_(FEAT_F128C(1.8)), DT_(FEAT_F128C(0.6)), DT_(FEAT_F128C(-0.00720165995109338994286772454927749004757)), DT_(FEAT_F128C(0.9486559629781224421354834251962896467232))), DT_(FEAT_F128C(0.4))*Math::pi<DT_>(), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::calc_opening_angle(DT_(FEAT_F128C(0.5)), DT_(FEAT_F128C(0.5)), DT_(FEAT_F128C(-0.0442463484164948589867174925087492580525)), DT_(FEAT_F128C(-0.2793604493334841106541998743643783654916))), DT_(FEAT_F128C(1.2))*Math::pi<DT_>(), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::calc_opening_angle(DT_(FEAT_F128C(-13.5)), DT_(FEAT_F128C(0.14)), DT_(FEAT_F128C(-10.83943948874084398729734229831858710701)), DT_(FEAT_F128C(8.048363285160879883151851466033077044309))), DT_(FEAT_F128C(1.8))*Math::pi<DT_>(), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::calc_opening_angle(DT_(FEAT_F128C(0.0)), DT_(FEAT_F128C(2.3)), DT_(FEAT_F128C(-0.6297833374543935798337006890115125282984)), DT_(FEAT_F128C(-0.6706511372276986012949918935308383824628))), DT_(FEAT_F128C(0.76))*Math::pi<DT_>(), tol);
  }
};

MathTest<float> math_test_float(5.8706E-07F, 1.7003E-03F);
MathTest<double> math_test_double(8.1620E-15, 5.4774E-07);
#ifdef FEAT_HAVE_QUADMATH
MathTest<__float128> math_test_float128(4.5308E-31Q, 3.2648E-14Q);
#endif // FEAT_HAVE_QUADMATH


/**
 * \brief Test class for the invert_matrix function.
 *
 * \test Tests the matrix inversion on a set of random matrices.
 *
 * \author Peter Zajac
 */
template<typename DT_, typename IT_>
class MatrixInvertTest :
  public TestSystem::UnitTest
{
public:
  MatrixInvertTest() :
    TestSystem::UnitTest("MatrixInvertTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name())
  {
  }

  virtual ~MatrixInvertTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create an RNG
    Random rng;

    // choose minimum and maximum system size
    static constexpr IT_ n_min = 1;
    static constexpr IT_ n_max = 9;
    static constexpr IT_ k_max = 5;
    static constexpr IT_ nnze = n_max*(n_max+k_max);

    // create four matrices
    DT_ a[nnze], b[nnze];
    // create a pivot array
    IT_ p[n_max];

    // loop over all n
    for(IT_ n(n_min); n <= n_max; ++n)
    {
      // for each dimension, we loop several times, eventually increasing the stride
      for(IT_ k(0); k <= k_max; ++k)
      {
        // compute stride
        const IT_ s = n + (k / 2);

        // generate a random matrix with values in range [-1,+1]
        for(IT_ i(0); i < s*n; ++i)
        {
          a[i] = b[i] = rng(-DT_(1), DT_(1));
        }

        // Note:
        // There is no guarantee that the generated matrix is actually regular.
        // However, for our RNG implementation and the ranges chosen here, it seems
        // that all generated matrices are fit for inversion.

        // invert our matrix
        const DT_ ret = Math::invert_matrix(n, s, b, p);

        // compute xl := ||I - A^{-1}*A||, xr := ||I - A*A^{-1}||
        DT_ xl(0), xr(0);
        for(IT_ i(0); i < n; ++i)
        {
          for(IT_ j(0); j < n; ++j)
          {
            DT_ yr(i == j ? DT_(1) : DT_(0));
            DT_ yl(yr);
            for(IT_ l(0); l < n; ++l)
            {
              yl -= b[i*s+l] * a[l*s+j];
              yr -= a[i*s+l] * b[l*s+j];
            }
            xl = Math::max(xl, Math::abs(yl));
            xr = Math::max(xr, Math::abs(yr));
          }
        }

        // print some numbers to the console
        std::cout << n << "/" << k << ": ";
        std::cout << stringify_fp_sci(Math::abs(ret)) << " : ";
        std::cout << stringify_fp_sci(xl) << " / " << stringify_fp_sci(xr) << "\n";

        // make sure the pivot quotient is not bogus
        TEST_CHECK(Math::isfinite(ret));

        // only check result ifpivot quotient is normal
        if(!Math::isnormal(ret))
        {
          std::cerr << "WARNING: pivot quotient is not normal!" << "\n";
          continue;
        }

        TEST_CHECK_EQUAL_WITHIN_EPS(xl, DT_(0), tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(xr, DT_(0), tol);
      }
    }
  }
};

MatrixInvertTest<double, std::uint32_t> matrix_invert_test_double_uint32;
#ifdef FEAT_HAVE_QUADMATH
MatrixInvertTest<__float128, int> matrix_invert_test_float128_int;
#endif // FEAT_HAVE_QUADMATH
