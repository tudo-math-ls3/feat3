// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/half.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/random.hpp>
#include <test_system/test_system.hpp>

#ifdef FEAT_HAVE_CGAL
#include <kernel/geometry/cgal.hpp>

#include <sstream>
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

  /**
   * \brief Returns a std::stringstream containing a OFF file of a tetrahedron
   */
  static std::stringstream tetrahedron_mesh()
  {
    std::stringstream mts;
    mts<<"OFF\n";
    mts<<"4 4 6\n";
    mts<<"0.0 0.0 2.0\n";
    mts<<"1.632993 -0.942809 -0.666667\n";
    mts<<"0.000000 1.885618 -0.666667\n";
    mts<<"-1.632993 -0.942809 -0.666667\n";
    mts<<"3 1 0 3\n";
    mts<<"3 2 0 1\n";
    mts<<"3 3 0 2\n";
    mts<<"3 3 2 1\n";

    return mts;
  }

  /**
   * \brief Returns a std::stringstream containing a OFF file of a unit cube
   */
  static std::stringstream cube_mesh()
  {
    std::stringstream mts;
    mts << "OFF\n";
    mts << "8 12 0\n";
    mts << "-0.5 -0.5 -0.5\n";
    mts << "0.5 -0.5 -0.5\n";
    mts << "-0.5 0.5 -0.5\n";
    mts << "0.5 0.5 -0.5\n";
    mts << "-0.5 -0.5 0.5\n";
    mts << "0.5 -0.5 0.5\n";
    mts << "-0.5 0.5 0.5\n";
    mts << "0.5 0.5 0.5\n";
    mts << "3 1 0 3\n";
    mts << "3 0 2 3\n";
    mts << "3 4 5 6\n";
    mts << "3 5 7 6\n";
    mts << "3 0 1 4\n";
    mts << "3 1 5 4\n";
    mts << "3 3 2 7\n";
    mts << "3 2 6 7\n";
    mts << "3 2 0 6\n";
    mts << "3 0 4 6\n";
    mts << "3 1 3 5\n";
    mts << "3 3 7 5\n";

    return mts;
  }

  /**
   * \brief Returns a std::stringstream containing a OFF file of a single quad
   */
  static std::stringstream quad_mesh()
  {
    std::stringstream mts;
    mts<<"OFF\n";
    mts<<"4 2 5\n";
    mts<<"0.0 0.0 0.0\n";
    mts<<"0.0 0.0 1.0\n";
    mts<<"0.0 1.0 0.0\n";
    mts<<"0.0 1.0 1.0\n";
    mts<<"3 0 1 2\n";
    mts<<"3 1 3 2\n";

    return mts;
  }

#if defined(__clang__)
  __attribute__((no_sanitize("undefined")))
#endif
  virtual void run() const override
  {
    // choose tolerance, but not tighter than double precision
    const DT_ tol = DT_(2) * Math::max(Math::eps<DT_>(), DT_(Math::eps<double>()));

    test_move_constructor();
    test_thread_safety();
    test_point_inside();
    test_closest_point(tol);
    test_transform(tol);
    test_squared_distance();
    test_feature_detection(tol);
    test_mesh_constructor(tol);
  }

  void test_move_constructor() const
  {
    std::stringstream mts = tetrahedron_mesh();
    CGALWrapper<DT_> cw2(mts, CGALFileMode::fm_off);
    CGALWrapper<DT_> cw_alt(std::move(cw2));

    TEST_CHECK_EQUAL(cw_alt.point_inside(DT_(0.1), DT_(0.1), DT_(0.1)), true);
    TEST_CHECK_EQUAL(cw_alt.point_inside(DT_(50), DT_(50), DT_(50)), false);
    cw2 = std::move(cw_alt);
    TEST_CHECK_EQUAL(cw2.point_inside(DT_(0.1), DT_(0.1), DT_(0.1)), true);
    TEST_CHECK_EQUAL(cw2.point_inside(DT_(50), DT_(50), DT_(50)), false);
  }

  void test_thread_safety() const
  {
    const std::size_t num_points = 582;
    constexpr DT_ range = 3.;
    // create random points
    std::vector<Tiny::Vector<DT_, 3>> points(num_points);
    // create rng object

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    std::for_each(points.begin(), points.end(), [&rng](Tiny::Vector<DT_,3>& vec){vec = {DT_(rng.next())/DT_(3.14), -DT_(rng.next())/DT_(7.19), DT_(rng.next())/DT_(1.554)};});
    DT_ max_val = (*std::max_element(points.begin(), points.end(), [](Tiny::Vector<DT_, 3>& a, Tiny::Vector<DT_, 3>& b){return a.norm_euclid_sqr() < b.norm_euclid_sqr();})).norm_euclid();
    std::transform(points.begin(), points.end(), points.begin(), [&max_val](const Tiny::Vector<DT_, 3>& a){return a * (range/max_val);});

    // setup cgal
    std::stringstream mts = tetrahedron_mesh();
    CGALWrapper<DT_> cw(mts, CGALFileMode::fm_off);

    //evalute points single threaded
    std::vector<int> inside_tests(num_points);
    std::transform(points.begin(), points.end(), inside_tests.begin(), [&cw](const Tiny::Vector<DT_, 3>& a){return int(cw.point_inside(a[0], a[1], a[2]));});

    // now do a omp based loop
    std::vector<int> inside_test2(num_points);
    FEAT_PRAGMA_OMP(parallel for)
    for(std::size_t i = 0; i < num_points; ++i)
    {
      inside_test2[i] = int(cw.point_inside(points[i][0], points[i][1], points[i][2]));
    }

    for(std::size_t i = 0; i < num_points; ++i)
    {
      TEST_CHECK_EQUAL(inside_tests[i], inside_test2[i]);
    }
  }

  void test_point_inside() const
  {
    std::stringstream mts = tetrahedron_mesh();
    CGALWrapper<DT_> cw2(mts, CGALFileMode::fm_off);
    TEST_CHECK_EQUAL(cw2.point_inside(DT_(50), DT_(50), DT_(50)), false);
    TEST_CHECK_EQUAL(cw2.point_inside(DT_(0.1), DT_(0.1), DT_(0.1)), true);
  }

  void test_closest_point(DT_ tol) const
  {
    std::stringstream mts = quad_mesh();
    CGALWrapper<DT_> cw(mts, CGALFileMode::fm_off);

    FEAT::Tiny::Vector<DT_, 3> spoint{{DT_(0.11), DT_(0.65), DT_(0.234)}};

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

  void test_transform(DT_ tol) const
  {
    std::stringstream mts = quad_mesh();
    CGALWrapper<DT_> cw(mts, CGALFileMode::fm_off);

    FEAT::Tiny::Vector<DT_, 3> spoint{{DT_(0.11), DT_(0.65), DT_(0.234)}};
    typename CGALWrapper<DT_>::TransformMatrix mat;
    typename CGALWrapper<DT_>::PointType trans{DT_(0),DT_(0),DT_(1)};
    mat.set_rotation_3d(DT_(0), Math::pi<DT_>(), DT_(0));

    cw.transform(mat, trans);
    auto point = cw.closest_point(spoint);
    TEST_CHECK_EQUAL_WITHIN_EPS(point[0], DT_(0.), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(point[1], spoint[1], tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(point[2], spoint[2], tol);
  }

  void test_squared_distance() const
  {
    std::stringstream mts = tetrahedron_mesh();
    CGALWrapper<DT_> cw(mts, CGALFileMode::fm_off);

    // closest point is a vertex
    DT_ d1 = cw.squared_distance(DT_(0.0), DT_(2.0), DT_(-0.666667));
    TEST_CHECK_EQUAL_WITHIN_EPS(d1, DT_(0.013083241924), DT_(1E-3));

    // closest point is on an edge
    DT_ d2 = cw.squared_distance(DT_(0.4), DT_(-1.0), DT_(-1.0));
    TEST_CHECK_EQUAL_WITHIN_EPS(d2, DT_(0.11438169937), DT_(1E-3));

    // closest point is on a face
    DT_ d3 = cw.squared_distance(DT_(0.0), DT_(0.0), DT_(-1.0));
    TEST_CHECK_EQUAL_WITHIN_EPS(d3, DT_(0.111110888889), DT_(1E-3));
  }

  void test_feature_detection(DT_ tol) const
  {
    std::stringstream tet_file = tetrahedron_mesh();
    CGALWrapper<DT_> tet_wrapper(tet_file, CGALFileMode::fm_off);

    std::stringstream cube_file = cube_mesh();
    CGALWrapper<DT_> cube_wrapper(cube_file, CGALFileMode::fm_off);

    std::stringstream quad_file = quad_mesh();
    CGALWrapper<DT_> quad_wrapper(quad_file, CGALFileMode::fm_off);

    // Mark all edges whose adjacent normal vectors differ by more than 45 degrees
    CGALFeatureNetwork tet_network = tet_wrapper.detect_features(45.0);
    CGALFeatureNetwork cube_network = cube_wrapper.detect_features(45.0);
    CGALFeatureNetwork quad_network = quad_wrapper.detect_features(45.0);

    TEST_CHECK_EQUAL(tet_network.size(), 6);
    TEST_CHECK_EQUAL(cube_network.size(), 12);
    TEST_CHECK_EQUAL(quad_network.size(), 1);

    // Tet-Mesh network features should consist of a single edge (two vertices)
    for(const FEAT::Geometry::CGALFeature& f : tet_network)
    {
      TEST_CHECK_EQUAL(f.size(), 2);
    }

    // Cube-Mesh network features should consist of a single edge (two vertices)
    for(const FEAT::Geometry::CGALFeature& f : cube_network)
    {
      TEST_CHECK_EQUAL(f.size(), 2);
    }

    // Quad-mesh feature should consist of one edge loop (4 vertices, start duplicated)
    TEST_CHECK_EQUAL(quad_network[0].size(), 5);
    TEST_CHECK_EQUAL(quad_network[0].front(), quad_network[0].back());

    // Tet-Mesh network should contain all edges of the tetrahedron
    for(Index i(0); i < 4; i++)
    {
      for(Index j(0); j < 4; j++)
      {
        if(i != j)
        {
          bool found = false;
          for(const FEAT::Geometry::CGALFeature& f : tet_network)
          {
            found = found || (f[0] == i && f[1] == j) || (f[0] == j && f[1] == i);
          }
          TEST_CHECK(found);
        }
      }
    }

    // Quad feature should be circular around the quad
    auto it = quad_network[0].begin();
    for(Index i : {1, 0, 2, 3, 1})
    {
      TEST_CHECK_EQUAL(*(it++), i);
    }

    // Check closest points
    using PointType = typename CGALWrapper<DT_>::PointType;

    PointType query_left{{DT_(0), DT_(-1.0), DT_(0.5)}};
    PointType query_top{{DT_(0), DT_(0.5), DT_(5)}};

    // Closest point should be center of left edge
    PointType closest_left = quad_wrapper.closest_point_on_feature(quad_network[0], query_left);
    TEST_CHECK_EQUAL_WITHIN_EPS(closest_left[0], DT_(0.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(closest_left[1], DT_(0.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(closest_left[2], DT_(0.5), tol);

    // Closest point should be center of top edge
    PointType closest_top = quad_wrapper.closest_point_on_feature(quad_network[0], query_top);
    TEST_CHECK_EQUAL_WITHIN_EPS(closest_top[0], DT_(0.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(closest_top[1], DT_(0.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(closest_top[2], DT_(1.0), tol);

    DT_ distance_left = quad_wrapper.squared_distance_to_feature(quad_network[0], query_left);
    DT_ distance_top = quad_wrapper.squared_distance_to_feature(quad_network[0], query_top);

    TEST_CHECK_EQUAL_WITHIN_EPS(distance_left, DT_(1.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(distance_top, DT_(16.0), tol);

    // Squared distance and closest points should agree
    TEST_CHECK_EQUAL_WITHIN_EPS(distance_left, (query_left - closest_left).norm_euclid_sqr(), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(distance_top, (query_top - closest_top).norm_euclid_sqr(), tol);
  }

  void test_mesh_constructor(DT_ tol) const
  {
    using MeshType = ConformalMesh<Shape::Hexahedron, 3, DT_>;
    RefinedUnitCubeFactory<MeshType> factory(1);
    MeshType mesh(factory);
    mesh.reorient_boundary_facets();
    MeshPart<MeshType> boundary = make_boundary_meshpart(mesh);

    CGALWrapper<DT_> wrapper(mesh, boundary);

    TEST_CHECK(wrapper.point_inside(0.5, 0.5, 0.5));
    TEST_CHECK(!wrapper.point_inside(2.0, 2.0, 2.0));

    CGALFeatureNetwork fnet = wrapper.detect_features(45.0);
    TEST_CHECK_EQUAL(fnet.size(), 12);

    TEST_CHECK_EQUAL_WITHIN_EPS(wrapper.squared_distance(5.0, 0.5, 0.5), DT_(16), tol);
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
