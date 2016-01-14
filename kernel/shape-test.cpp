#include <test_system/test_system.hpp>
#include <kernel/shape.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Shape;

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

    // Check reference cell facet orientations
    // For Simplex<1>, these are -1, +1
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<1>>::orientation(0), -1);
    TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<1>>::orientation(1), 1);

    // For all other simplices, all facet orientations are positive
    for(int facet_index(0); facet_index < FaceTraits<Simplex<2>,1>::count; ++facet_index)
      TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<2>>::orientation(facet_index), 1);

    for(int facet_index(0); facet_index < FaceTraits<Simplex<3>,2>::count; ++facet_index)
      TEST_CHECK_EQUAL(ReferenceCell<Shape::Simplex<3>>::orientation(facet_index), 1);
  }
} shape_test;
