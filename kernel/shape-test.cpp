// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/shape.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Shape;

/**
 * \brief Test class for the Shape and FaceTraits class templates.
 *
 * \test Tests the Shape and FaceTraits class templates.
 *
 * \author Peter Zajac
 */
class ShapeTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  ShapeTest() :
    TaggedTest<Archs::None, Archs::None>("shape_test")
  {
  }

  int bla(int dim, int i) const
  {
    int ret = 1 - ((((dim >> 1) ^ 1) ^ i) << 1);
    std::cout << ret << std::endl;
    return ret;
  }

  virtual void run() const override
  {
    // test simplex vertex counts

    // 1-Simplex (line): 2 vertices
    int simplex1_nverts = FaceTraits<Simplex<1>,0>::count;
    TEST_CHECK_EQUAL(simplex1_nverts, 2);

    // 2-Simplex (triangle): 3 vertices
    int simplex2_nverts = FaceTraits<Simplex<2>,0>::count;
    TEST_CHECK_EQUAL(simplex2_nverts, 3);

    // 3-Simplex (tetrahedron): 4 vertices
    int simplex3_nverts = FaceTraits<Simplex<3>,0>::count;
    TEST_CHECK_EQUAL(simplex3_nverts, 4);


    // test hypercube vertex counts

    // 1-Hypercube (line): 2 vertices
    int hcube1_nverts = FaceTraits<Hypercube<1>,0>::count;
    TEST_CHECK_EQUAL(hcube1_nverts, 2);

    // 2-Hypercube (quadrilateral): 4 vertices
    int hcube2_nverts = FaceTraits<Hypercube<2>,0>::count;
    TEST_CHECK_EQUAL(hcube2_nverts, 4);

    // 3-Hypercube (hexahedron): 8 vertices
    int hcube3_nverts = FaceTraits<Hypercube<3>,0>::count;
    TEST_CHECK_EQUAL(hcube3_nverts, 8);


    // test 2-Simplex face counts
    // 2-Simplex: 3x 1-faces
    int simplex2_nface1 = FaceTraits<Simplex<2>,1>::count;
    TEST_CHECK_EQUAL(simplex2_nface1, 3);


    // test 3-Simplex face counts
    // 3-Simplex: 6x 1-faces
    int simplex3_nface1 = FaceTraits<Simplex<3>,1>::count;
    TEST_CHECK_EQUAL(simplex3_nface1, 6);

    // 3-Simplex: 4x 2-faces
    int simplex3_nface2 = FaceTraits<Simplex<3>,2>::count;
    TEST_CHECK_EQUAL(simplex3_nface2, 4);


    // test 2-hypercube face counts
    // 2-Hypercube: 4x 1-faces
    int hcube2_nface1 = FaceTraits<Hypercube<2>,1>::count;
    TEST_CHECK_EQUAL(hcube2_nface1, 4);


    // test 3-hypercube face counts
    // 3-Hypercube: 12x 1-faces
    int hcube3_nface1 = FaceTraits<Hypercube<3>,1>::count;
    TEST_CHECK_EQUAL(hcube3_nface1, 12);

    // 3-Hypercube: 6x 2-faces
    int hcube3_nface2 = FaceTraits<Hypercube<3>,2>::count;
    TEST_CHECK_EQUAL(hcube3_nface2, 6);

    const double eps = 1E-12;

    // test simplex vertex coords

    // 1-Simplex: {[0], [1]}
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<1>>::vertex<double>(0, 0), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<1>>::vertex<double>(1, 0), 1.0, eps);

    // 2-Simplex: {[0,0], [1,0], [0,1]}
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<2>>::vertex<double>(0, 0), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<2>>::vertex<double>(0, 1), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<2>>::vertex<double>(1, 0), 1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<2>>::vertex<double>(1, 1), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<2>>::vertex<double>(2, 0), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<2>>::vertex<double>(2, 1), 1.0, eps);

    // 3-Simplex: {[0,0,0], [1,0,0], [0,1,0], [0,0,1]}
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(0, 0), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(0, 1), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(0, 2), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(1, 0), 1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(1, 1), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(1, 2), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(2, 0), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(2, 1), 1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(2, 2), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(3, 0), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(3, 1), 0.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::vertex<double>(3, 2), 1.0, eps);

    // test hypercube vertex coords

    // 1-Hypercube: {[-1], [+1]}
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<1>>::vertex<double>(0, 0), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<1>>::vertex<double>(1, 0), +1.0, eps);

    // 2-Hypercube
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::vertex<double>(0, 0), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::vertex<double>(0, 1), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::vertex<double>(1, 0), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::vertex<double>(1, 1), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::vertex<double>(2, 0), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::vertex<double>(2, 1), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::vertex<double>(3, 0), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::vertex<double>(3, 1), +1.0, eps);

    // 3-Hypercube
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(0, 0), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(0, 1), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(0, 2), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(1, 0), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(1, 1), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(1, 2), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(2, 0), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(2, 1), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(2, 2), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(3, 0), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(3, 1), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(3, 2), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(4, 0), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(4, 1), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(4, 2), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(5, 0), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(5, 1), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(5, 2), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(6, 0), -1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(6, 1), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(6, 2), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(7, 0), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(7, 1), +1.0, eps);
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::vertex<double>(7, 2), +1.0, eps);

    // test simplex centre coords

    // 1-Simplex: [1/2]
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<1>>::centre<double>(0), 1.0/2.0, eps);

    // 2-Simplex: [1/3, 1/3]
    for(int i(0); i < 2; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<2>>::centre<double>(i), 1.0/3.0, eps);
    }

    // 3-Simplex: [1/4, 1/4, 1/4]
    for(int i(0); i < 3; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::centre<double>(i), 1.0/4.0, eps);
    }

    // test hypercube centre coords

    // 1-Hypercube
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<1>>::centre<double>(0), 0.0, eps);

    // 2-Hypercube
    for(int i(0); i < 2; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::centre<double>(i), 0.0, eps);
    }

    // 3-Hypercube
    for(int i(0); i < 3; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::centre<double>(i), 0.0, eps);
    }

    // test simplex volumes

    // 1-Simplex: 1
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<1>>::volume<double>(), 1.0, eps);

    // 2-Simplex: 1/2
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<2>>::volume<double>(), 1.0/2.0, eps);

    // 3-Simplex: 1/6
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Simplex<3>>::volume<double>(), 1.0/6.0, eps);

    // test hypercube volumes

    // 1-Hypercube: 2
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<1>>::volume<double>(), 2.0, eps);

    // 2-Hypercube: 4
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<2>>::volume<double>(), 4.0, eps);

    // 3-Hypercube: 8
    TEST_CHECK_EQUAL_WITHIN_EPS(ReferenceCell<Hypercube<3>>::volume<double>(), 8.0, eps);


    // Check reference cell facet orientations
    // For Simplex<1>, these are -1, +1
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<1>>::facet_orientation(0), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<1>>::facet_orientation(1), +1);

    // For Simplex<2>, all are +1
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<2>>::facet_orientation(0), +1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<2>>::facet_orientation(1), +1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<2>>::facet_orientation(2), +1);

    // For Simplex<3>, we have +1 -1 +1 -1
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<3>>::facet_orientation(0), +1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<3>>::facet_orientation(1), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<3>>::facet_orientation(2), +1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<3>>::facet_orientation(3), -1);

    // For Hypercube<1>, these are -1, +1
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<1>>::facet_orientation(0), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<1>>::facet_orientation(1), +1);

    // For Hypercube<2>, these are +1, -1, -1, +1
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<2>>::facet_orientation(0), +1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<2>>::facet_orientation(1), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<2>>::facet_orientation(2), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<2>>::facet_orientation(3), +1);

    // For Hypercube<2>, these are +1, -1, -1, +1
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<3>>::facet_orientation(0), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<3>>::facet_orientation(1), +1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<3>>::facet_orientation(2), +1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<3>>::facet_orientation(3), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<3>>::facet_orientation(4), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Hypercube<3>>::facet_orientation(5), +1);
  }
} shape_test;
