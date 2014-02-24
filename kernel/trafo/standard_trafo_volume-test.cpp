#include <test_system/test_system.hpp>
#include <kernel/archs.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/standard/volume.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry::TestAux;

/*
 * \brief Test for cell volume computations
 *
 * \test Verfies the routine's output against pre computed values
 *
 * \author Jordi Paul
 * */
class StandardTrafoVolumeTest
: public TestSystem::TaggedTest<Archs::None, Archs::None>
{
  public:
    StandardTrafoVolumeTest() :
      TestSystem::TaggedTest<Archs::None, Archs::None>("standard_trafo_volume_test")
  {
  }

    virtual void run() const
    {
      test_1d_quad();
      test_1d_simplex();
      test_2d_quad();
      test_2d_simplex();
      test_3d_quad();
      test_3d_simplex();
    }
    void test_1d_quad() const
    {
      typedef Shape::Hypercube<1> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      static const CoordType vtx[2*1] = { 1., sqrt(2.)};
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      CoordType vol = trafo.compute_vol<ShapeType>(Index(0));
      TEST_CHECK_EQUAL(vol, sqrt(2.)-1.);

    }

    void test_1d_simplex() const
    {
      typedef Shape::Simplex<1> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      static const CoordType vtx[2*1] = { 1., sqrt(2.)};
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      CoordType vol = trafo.compute_vol<ShapeType>(Index(0));
      TEST_CHECK_EQUAL(vol, sqrt(2.)-1.);

    }
    void test_2d_simplex() const
    {
      typedef Shape::Simplex<2> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      // Sufficiently wild simplex
      static const CoordType vtx[3*2] = { -1., -1., 1., 0.5, -0.5, 0.75};
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // Everything checked agains has been computed by hand
      CoordType vol = trafo.compute_vol<ShapeType>(Index(0));
      TEST_CHECK_EQUAL(vol, 11./8.);
      // Check volume of sub simplices
      vol = trafo.compute_vol<Shape::Simplex<1>>(Index(0));
      TEST_CHECK_EQUAL(vol, sqrt(37./16.));
      vol = trafo.compute_vol<Shape::Simplex<1>>(Index(1));
      TEST_CHECK_EQUAL(vol, sqrt(53./16.));
      vol = trafo.compute_vol<Shape::Simplex<1>>(Index(2));
      TEST_CHECK_EQUAL(vol, 5./2.);

    }

    void test_2d_quad() const
    {
      typedef Shape::Hypercube<2> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      // Sufficiently wild hypercube
      static const CoordType vtx[4*2] =
      { -Math::sqrt(CoordType(2.)), -CoordType(3.),
        CoordType(2.), -CoordType(1.),
        -CoordType(1.), CoordType(1.),
        CoordType(2.5), Math::sqrt(CoordType(3.))};
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // This is just sqrt(bloody machine precision) because the computations contain a sqrt() in the end
      CoordType eps = CoordType(1.e-7);
      // Everything checked against has been computed by hand
      CoordType vol = trafo.compute_vol<ShapeType>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, CoordType(17025./1546.), eps);
      static const CoordType l[4] = {10735./2713., 10788./3017., 22749./5657., 26405./9507.};

      // Check volume of sub elements
      for(Index i(0); i < 4; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<1>>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, l[i], eps);
      }

    }

    void test_3d_quad() const
    {
      typedef Shape::Hypercube<3> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      /* In 3d. the volume for nontrivial hypercubes (especially those with nonplanar faces) is difficult to compute
       * by hand, so this is just a parallel piped
       */
      static const CoordType vtx[8*3] =
      { -2., -2., -2.,
        2., -2., -2.,
        -2., 2., -2.,
        2., 2., -2.,
        -1., -2., 2.,
        3., -2., 2.,
        -1., 2., 2.,
        3., 2., 2. };
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // This is just sqrt(bloody machine precision) because the computations contain a sqrt() in the end
      CoordType eps = CoordType(1.e-7);
      // Everything checked against has been computed by hand
      CoordType vol = trafo.compute_vol<ShapeType>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, CoordType(64.), eps);

      // Check the 2d volume of the faces
      static const CoordType f[6] = {16., 16., 16., 16., 4.*Math::sqrt(17.), 4.*Math::sqrt(17.)};

      for(Index i(0); i < 6; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<2>>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, f[i], eps);
      }

      // Check the 1d volume of the edges
      static const CoordType l[12] = {
        4., 4., 4., 4.,
        4., 4., 4., 4.,
        Math::sqrt(17.), Math::sqrt(17.), Math::sqrt(17.), Math::sqrt(17.)};

      for(Index i(0); i < 12; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<1>>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, l[i], eps);
      }

    }

    void test_3d_simplex() const
    {
      typedef Shape::Simplex<3> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      /* Our simplex is the standard simplex scaled by a factor of 2 in x_3 direction, and vertex 3 has been moved
         which does not alter the volume of 2*1/d! = 1/3
         */
      static const CoordType vtx[4*3] =
      { 0., 0., 0.,
        1., 0., 0.,
        0., 1., 0.,
        0., 0.5, 2., };
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      CoordType vol = trafo.compute_vol<ShapeType>(Index(0));
      TEST_CHECK_EQUAL(vol, 1./3.);

      /* Edge lengths computed by hand */
      static const CoordType l[6] = {1., 1., sqrt(17./4.), sqrt(2), sqrt(21./4.), sqrt(17./4.)};
      for(Index i=0; i<6; i++)
      {
        vol = trafo.compute_vol<Shape::Simplex<1>>(Index(i));
        TEST_CHECK_EQUAL(vol, l[i]);
      }

      /* With the edgelengths, check the volume of the sub simplices via Heron's formula */
      // s will be half the circumference of a sub simplex
      CoordType s = 0.;

      // Bloody machine precision
      CoordType eps = 1.e-14;
      // Face 0 consists of edges 3, 4, 5
      s = 0.5*(l[3] + l[4] + l[5]);
      vol = trafo.compute_vol<Shape::Simplex<2>>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, sqrt(s * (s - l[3]) * (s - l[4]) * (s - l[5]) ), eps);
      // Face 1 consists of edges 1, 2, 5
      s = 0.5*(l[1] + l[2] + l[5]);
      vol = trafo.compute_vol<Shape::Simplex<2>>(Index(1));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, sqrt(s * (s - l[1]) * (s - l[2]) * (s - l[5]) ), eps);
      // Face 2 consists of edges 0, 2, 4
      s = 0.5*(l[0] + l[2] + l[4]);
      vol = trafo.compute_vol<Shape::Simplex<2>>(Index(2));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, sqrt(s * (s - l[0]) * (s - l[2]) * (s - l[4]) ), eps);
      // Face 3 consists of edges 0, 1, 3
      s = 0.5*(l[0] + l[1] + l[3]);
      vol = trafo.compute_vol<Shape::Simplex<2>>(Index(3));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, sqrt(s * (s - l[0]) * (s - l[1]) * (s - l[3]) ), eps);
    }
} standard_trafo_volume_test;
