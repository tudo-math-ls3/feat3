// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "kernel/util/math.hpp"
#include <iomanip>
#include <test_system/test_system.hpp>

#ifdef FEAT_HAVE_CGAL
#include <kernel/geometry/cgal.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;


class CGALTest
  : public UnitTest
{
public:
  CGALTest() :
    UnitTest("CGALTest")
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
    {
      std::stringstream mts;
      mts<<"OFF"<<std::endl;
      mts<<"4 4 6"<<std::endl;
      mts<<"0.0 0.0 2.0"<<std::endl;
      mts<<"1.632993 -0.942809 -0.666667"<<std::endl;
      mts<<"0.000000 1.885618 -0.666667"<<std::endl;
      mts<<"-1.632993 -0.942809 -0.666667"<<std::endl;
      mts<<"3 1 0 3"<<std::endl;
      mts<<"3 2 0 1"<<std::endl;
      mts<<"3 3 0 2"<<std::endl;
      mts<<"3 3 2 1"<<std::endl;

      CGALWrapper cw2(mts, CGALFileMode::fm_off);
      TEST_CHECK_EQUAL(cw2.point_inside(50, 50, 50), false);
      TEST_CHECK_EQUAL(cw2.point_inside(0.1, 0.1, 0.1), true);
    }
    {
      std::stringstream mts;
      mts<<"OFF"<<std::endl;
      mts<<"4 2 5"<<std::endl;
      mts<<"0.0 0.0 0.0"<<std::endl;
      mts<<"0.0 0.0 1.0"<<std::endl;
      mts<<"0.0 1.0 0.0"<<std::endl;
      mts<<"0.0 1.0 1.0"<<std::endl;
      mts<<"3 0 1 2"<<std::endl;
      mts<<"3 1 3 2"<<std::endl;

      FEAT::Tiny::Vector<double, 3> spoint{{0.11, 0.65, 0.234}};

      CGALWrapper cw(mts, CGALFileMode::fm_off);
      auto point = cw.closest_point(spoint);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[0], 0., Math::eps<double>());
      TEST_CHECK_EQUAL_WITHIN_EPS(point[1], spoint[1], Math::eps<double>());
      TEST_CHECK_EQUAL_WITHIN_EPS(point[2], spoint[2], Math::eps<double>());
      //rotate three times by 2*pi/3 regarding x-axis
      typename CGALWrapper::TransformMatrix mat;
      typename CGALWrapper::PointType trans{0,0,0};
      mat.set_rotation_3d(double(2.) * Math::pi<double>() / double(3.), double(0), double(0));
      cw.transform(mat,trans);
      cw.transform(mat,trans);
      cw.transform(mat,trans);
      auto point2 = cw.closest_point(spoint);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[0], point2[0], double(2) * Math::eps<double>());
      TEST_CHECK_EQUAL_WITHIN_EPS(point[1], point2[1], double(2) * Math::eps<double>());
      TEST_CHECK_EQUAL_WITHIN_EPS(point[2], point2[2], double(2) * Math::eps<double>());
      trans[0] = double(0.11);
      trans[1] = double(0.1);
      trans[2] = double(0.334);
      mat.set_identity();
      cw.transform(mat, trans);
      auto point3 = cw.closest_point(spoint);

      TEST_CHECK_EQUAL_WITHIN_EPS(point[0], point3[0] - trans[0], double(2) * Math::eps<double>());
      TEST_CHECK_EQUAL_WITHIN_EPS(point[1], point3[1], double(2) * Math::eps<double>());
      TEST_CHECK_EQUAL_WITHIN_EPS(point[2], point3[2] - trans[2] + spoint[2], double(2) * Math::eps<double>());
    }

    {
      std::stringstream mts;
      mts<<"OFF"<<std::endl;
      mts<<"4 2 5"<<std::endl;
      mts<<"0.0 0.0 0.0"<<std::endl;
      mts<<"0.0 0.0 1.0"<<std::endl;
      mts<<"0.0 1.0 0.0"<<std::endl;
      mts<<"0.0 1.0 1.0"<<std::endl;
      mts<<"3 0 1 2"<<std::endl;
      mts<<"3 1 3 2"<<std::endl;

      FEAT::Tiny::Vector<double, 3> spoint{{0.11, 0.65, 0.234}};
      typename CGALWrapper::TransformMatrix mat;
      typename CGALWrapper::PointType trans{0,0,1};
      mat.set_rotation_3d(double(0), Math::pi<double>(), double(0));

      CGALWrapper cw(mts, CGALFileMode::fm_off);
      cw.transform(mat, trans);
      auto point = cw.closest_point(spoint);
      TEST_CHECK_EQUAL_WITHIN_EPS(point[0], 0., Math::eps<double>());
      TEST_CHECK_EQUAL_WITHIN_EPS(point[1], spoint[1], Math::eps<double>());
      TEST_CHECK_EQUAL_WITHIN_EPS(point[2], spoint[2], Math::eps<double>());


    }


  }

} cgal_test;
#endif //FEAT_HAVE_CGAL
