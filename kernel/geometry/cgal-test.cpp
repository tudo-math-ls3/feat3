// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "kernel/util/math.hpp"
#include <iomanip>
#include <test_system/test_system.hpp>
#include <kernel/util/half.hpp>
#include <kernel/util/random.hpp>

#ifdef FEAT_HAVE_CGAL
#include <kernel/geometry/cgal.hpp>

#include <algorithm>
#ifdef FEAT_HAVE_OMP
#include "omp.h"
#endif

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

template<typename DT_>
class CGALTest
  : public UnitTest
{
public:
  CGALTest() :
    UnitTest("CGALTest", Type::Traits<DT_>::name())
  {
  }

  virtual ~CGALTest()
  {
  }

#if defined(__clang__)
  __attribute__((no_sanitize("undefined")))
#endif
  virtual void run() const override
  {
    // choose tolerance, but not tighter than double precision
    const DT_ tol = DT_(2) * Math::max(Math::eps<DT_>(), DT_(Math::eps<double>()));
    {
      std::stringstream mts;
      mts<<"OFF"<<"\n";
      mts<<"4 4 6"<<"\n";
      mts<<"0.0 0.0 2.0"<<"\n";
      mts<<"1.632993 -0.942809 -0.666667"<<"\n";
      mts<<"0.000000 1.885618 -0.666667"<<"\n";
      mts<<"-1.632993 -0.942809 -0.666667"<<"\n";
      mts<<"3 1 0 3"<<"\n";
      mts<<"3 2 0 1"<<"\n";
      mts<<"3 3 0 2"<<"\n";
      mts<<"3 3 2 1"<<"\n";

      CGALWrapper<DT_> cw2(mts, CGALFileMode::fm_off);
      CGALWrapper<DT_> cw_alt(std::move(cw2));

      TEST_CHECK_EQUAL(cw_alt.point_inside(DT_(0.1), DT_(0.1), DT_(0.1)), true);
      TEST_CHECK_EQUAL(cw_alt.point_inside(DT_(50), DT_(50), DT_(50)), false);
      cw2 = std::move(cw_alt);
      TEST_CHECK_EQUAL(cw2.point_inside(DT_(0.1), DT_(0.1), DT_(0.1)), true);
      TEST_CHECK_EQUAL(cw2.point_inside(DT_(50), DT_(50), DT_(50)), false);

    }

    // threadparallel test
    {
      const std::size_t num_points = 582;
      constexpr DT_ range = 3.;
      // create random points
      std::vector<Tiny::Vector<DT_, 3>> points(num_points);
      // create rng object

      Random rng;
      std::for_each(points.begin(), points.end(), [&rng](Tiny::Vector<DT_,3>& vec){vec = {DT_(rng.next())/DT_(3.14), -DT_(rng.next())/DT_(7.19), DT_(rng.next())/DT_(1.554)};});
      DT_ max_val = (*std::max_element(points.begin(), points.end(), [](Tiny::Vector<DT_, 3>& a, Tiny::Vector<DT_, 3>& b){return a.norm_euclid_sqr() < b.norm_euclid_sqr();})).norm_euclid();
      std::transform(points.begin(), points.end(), points.begin(), [&max_val](const Tiny::Vector<DT_, 3>& a){return a * (range/max_val);});

      // setup cgal
      std::stringstream mts;
      mts<<"OFF"<<"\n";
      mts<<"4 4 6"<<"\n";
      mts<<"0.0 0.0 2.0"<<"\n";
      mts<<"1.632993 -0.942809 -0.666667"<<"\n";
      mts<<"0.000000 1.885618 -0.666667"<<"\n";
      mts<<"-1.632993 -0.942809 -0.666667"<<"\n";
      mts<<"3 1 0 3"<<"\n";
      mts<<"3 2 0 1"<<"\n";
      mts<<"3 3 0 2"<<"\n";
      mts<<"3 3 2 1"<<"\n";

      CGALWrapper<DT_> cw(mts, CGALFileMode::fm_off);

      //evalute points single threaded
      std::vector<int> inside_tests(num_points);
      std::transform(points.begin(), points.end(), inside_tests.begin(), [&cw](const Tiny::Vector<DT_, 3>& a){return int(cw.point_inside(a[0], a[1], a[2]));});

      // now do a omp based loop
      std::vector<int> inside_test2(num_points);
      #pragma omp parallel for
      for(std::size_t i = 0; i < num_points; ++i)
      {
        // #ifdef FEAT_HAVE_OMP
        // if(i == 0)
        //   std::cout << "Num threads are " << omp_get_num_threads() << "\n";
        // #endif
        inside_test2[i] = int(cw.point_inside(points[i][0], points[i][1], points[i][2]));
      }

      for(std::size_t i = 0; i < num_points; ++i)
      {
        TEST_CHECK_EQUAL(inside_tests[i], inside_test2[i]);
      }
    }

    {
      std::stringstream mts;
      mts<<"OFF"<<"\n";
      mts<<"4 4 6"<<"\n";
      mts<<"0.0 0.0 2.0"<<"\n";
      mts<<"1.632993 -0.942809 -0.666667"<<"\n";
      mts<<"0.000000 1.885618 -0.666667"<<"\n";
      mts<<"-1.632993 -0.942809 -0.666667"<<"\n";
      mts<<"3 1 0 3"<<"\n";
      mts<<"3 2 0 1"<<"\n";
      mts<<"3 3 0 2"<<"\n";
      mts<<"3 3 2 1"<<"\n";

      CGALWrapper<DT_> cw2(mts, CGALFileMode::fm_off);
      TEST_CHECK_EQUAL(cw2.point_inside(DT_(50), DT_(50), DT_(50)), false);
      TEST_CHECK_EQUAL(cw2.point_inside(DT_(0.1), DT_(0.1), DT_(0.1)), true);
    }
    {
      std::stringstream mts;
      mts<<"OFF"<<"\n";
      mts<<"4 2 5"<<"\n";
      mts<<"0.0 0.0 0.0"<<"\n";
      mts<<"0.0 0.0 1.0"<<"\n";
      mts<<"0.0 1.0 0.0"<<"\n";
      mts<<"0.0 1.0 1.0"<<"\n";
      mts<<"3 0 1 2"<<"\n";
      mts<<"3 1 3 2"<<"\n";

      FEAT::Tiny::Vector<DT_, 3> spoint{{DT_(0.11), DT_(0.65), DT_(0.234)}};

      CGALWrapper<DT_> cw(mts, CGALFileMode::fm_off);
      auto point = cw.closest_point(spoint);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[0], DT_(0), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[1], spoint[1], tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[2], spoint[2], tol);
      //rotate three times by 2*pi/3 regarding x-axis
      typename CGALWrapper<DT_>::TransformMatrix mat;
      typename CGALWrapper<DT_>::PointType trans{DT_(0),DT_(0),DT_(0)};
      mat.set_rotation_3d(DT_(2.) * Math::pi<DT_>() / DT_(3.), DT_(0), DT_(0));
      cw.transform(mat,trans);
      cw.transform(mat,trans);
      cw.transform(mat,trans);
      auto point2 = cw.closest_point(spoint);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[0], point2[0], tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[1], point2[1], tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[2], point2[2], tol);

      trans[0] = DT_(0.11);
      trans[1] = DT_(0.1);
      trans[2] = DT_(0.334);
      mat.set_identity();
      cw.transform(mat, trans);
      auto point3 = cw.closest_point(spoint);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[0], point3[0] - trans[0], tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[1], point3[1], tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[2], point3[2] - trans[2] + spoint[2], tol);
    }

    {
      std::stringstream mts;
      mts<<"OFF"<<"\n";
      mts<<"4 2 5"<<"\n";
      mts<<"0.0 0.0 0.0"<<"\n";
      mts<<"0.0 0.0 1.0"<<"\n";
      mts<<"0.0 1.0 0.0"<<"\n";
      mts<<"0.0 1.0 1.0"<<"\n";
      mts<<"3 0 1 2"<<"\n";
      mts<<"3 1 3 2"<<"\n";

      FEAT::Tiny::Vector<DT_, 3> spoint{{DT_(0.11), DT_(0.65), DT_(0.234)}};
      typename CGALWrapper<DT_>::TransformMatrix mat;
      typename CGALWrapper<DT_>::PointType trans{DT_(0),DT_(0),DT_(1)};
      mat.set_rotation_3d(DT_(0), Math::pi<DT_>(), DT_(0));

      CGALWrapper<DT_> cw(mts, CGALFileMode::fm_off);
      cw.transform(mat, trans);
      auto point = cw.closest_point(spoint);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[0], DT_(0.), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[1], spoint[1], tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[2], spoint[2], tol);
    }
  }
};

CGALTest<double> cgal_test_double;
CGALTest<float> cgal_test_float;
#ifdef FEAT_HAVE_QUADMATH
CGALTest<__float128> cgal_test_float128;
#endif
#ifdef FEAT_HAVE_HALFMATH
CGALTest<Half> cgal_test_half;
#endif

#endif //FEAT_HAVE_CGAL
