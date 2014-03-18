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
template<typename DataType_>
class StandardTrafoVolumeTest
: public TestSystem::TaggedTest<Archs::None, DataType_>
{
  public:
    StandardTrafoVolumeTest() :
      TestSystem::TaggedTest<Archs::None, DataType_>("standard_trafo_volume_test")
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
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Hypercube<1> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      static const CoordType vtx[2*1] = { CoordType(1), Math::sqrt(CoordType(2))};
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      DataType_ vol = trafo.compute_vol<ShapeType,DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(DataType_(2)) - DataType_(1), eps);
    }

    void test_1d_simplex() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Simplex<1> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      static const CoordType vtx[2*1] = { CoordType(1), Math::sqrt(CoordType(2))};
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      DataType_ vol = trafo.compute_vol<ShapeType,DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(DataType_(2)) - DataType_(1), eps);
    }

    void test_2d_simplex() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Simplex<2> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      // Sufficiently wild simplex
      static const CoordType vtx[3*2] =
      { -CoordType(1), -CoordType(1),
        CoordType(1), CoordType(1)/CoordType(2),
        -CoordType(1)/CoordType(2), CoordType(3)/CoordType(4)};
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // Everything checked agains has been computed by hand
      DataType_ vol = trafo.compute_vol<ShapeType, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(11)/DataType_(8), eps);
      // Check volume of sub simplices
      vol = trafo.compute_vol<Shape::Simplex<1>, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(DataType_(37)/DataType_(16)), eps);
      vol = trafo.compute_vol<Shape::Simplex<1>, DataType_>(Index(1));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(DataType_(53)/DataType_(16)), eps);
      vol = trafo.compute_vol<Shape::Simplex<1>, DataType_>(Index(2));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(5)/DataType_(2), eps);
    }

    void test_2d_quad() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

      typedef Shape::Hypercube<2> ShapeType;
      typedef Geometry::ConformalMesh<ShapeType> MeshType;
      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef MeshType::VertexSetType::CoordType CoordType;

      Geometry::ReferenceCellFactory<ShapeType> factory;

      MeshType mesh(factory);
      // Sufficiently wild hypercube
      static const CoordType vtx[4*2] =
      { -Math::sqrt(CoordType(2)), -CoordType(3),
        CoordType(2), -CoordType(1),
        -CoordType(1), CoordType(1),
        CoordType(5)/CoordType(2), Math::sqrt(CoordType(3))};
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // Everything checked against has been computed by hand
      DataType_ vol = trafo.compute_vol<ShapeType, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(17025)/DataType_(1546), eps);

      static const DataType_ l[4] = {
        DataType_(10735)/DataType_(2713), DataType_(10788)/DataType_(3017),
        DataType_(22749)/DataType_(5657), DataType_(26405)/DataType_(9507)};

      // Check volume of sub elements
      for(Index i(0); i < 4; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<1>,DataType_>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, l[i], eps);
      }
    }

    void test_3d_quad() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

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
      { -CoordType(2), -CoordType(2), -CoordType(2),
        CoordType(2), -CoordType(2), -CoordType(2),
        -CoordType(2), CoordType(2), -CoordType(2),
        CoordType(2), CoordType(2), -CoordType(2),
        -CoordType(1), -CoordType(2), CoordType(2),
        CoordType(3), -CoordType(2), CoordType(2),
        -CoordType(1), CoordType(2), CoordType(2),
        CoordType(3), CoordType(2), CoordType(2) };
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      // Everything checked against has been computed by hand
      DataType_ vol = trafo.compute_vol<ShapeType, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(64), eps);

      // Check the 2d volume of the faces
      static const DataType_ f[6] =
      { DataType_(16), DataType_(16),
        DataType_(16), DataType_(16),
        DataType_(4)*Math::sqrt(DataType_(17)), DataType_(4)*Math::sqrt(DataType_(17))};

      for(Index i(0); i < 6; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<2>, DataType_>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, f[i], eps);
      }

      // Check the 1d volume of the edges
      static const DataType_ l[12] = {
        DataType_(4), DataType_(4), DataType_(4), DataType_(4),
        DataType_(4), DataType_(4), DataType_(4), DataType_(4),
        Math::sqrt(DataType_(17)), Math::sqrt(DataType_(17)),
        Math::sqrt(DataType_(17)), Math::sqrt(DataType_(17))};

      for(Index i(0); i < 12; ++i)
      {
        vol = trafo.compute_vol<Shape::Hypercube<1>, DataType_>(i);
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, l[i], eps);
      }
    }

    void test_3d_simplex() const
    {
      const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));
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
      { CoordType(0), CoordType(0), CoordType(0),
        CoordType(1), CoordType(0), CoordType(0),
        CoordType(0), CoordType(1), CoordType(0),
        CoordType(0), CoordType(0.5), CoordType(2), };
      copy_vtx((&mesh)->get_vertex_set(), vtx);

      TrafoType trafo(mesh);

      DataType_ vol = trafo.compute_vol<ShapeType, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(1)/DataType_(3), eps);

      /* Edge lengths computed by hand */
      static const DataType_ l[6] =
      {DataType_(1), DataType_(1),
        Math::sqrt(DataType_(17)/DataType_(4)), Math::sqrt(DataType_(2)),
        Math::sqrt(DataType_(21)/DataType_(4)), Math::sqrt(DataType_(17)/DataType_(4))};
      for(Index i=0; i<6; i++)
      {
        vol = trafo.compute_vol<Shape::Simplex<1>, DataType_>(Index(i));
        TEST_CHECK_EQUAL_WITHIN_EPS(vol, l[i], eps);
      }

      /* With the edgelengths, check the volume of the sub simplices via Heron's formula */
      // s will be half the circumference of a sub simplex
      DataType_ s = 0.;

      // Face 0 consists of edges 3, 4, 5
      s = DataType_(0.5)*(l[3] + l[4] + l[5]);
      vol = trafo.compute_vol<Shape::Simplex<2>, DataType_>(Index(0));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(s * (s - l[3]) * (s - l[4]) * (s - l[5]) ), eps);
      // Face 1 consists of edges 1, 2, 5
      s = DataType_(0.5)*(l[1] + l[2] + l[5]);
      vol = trafo.compute_vol<Shape::Simplex<2>, DataType_>(Index(1));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(s * (s - l[1]) * (s - l[2]) * (s - l[5]) ), eps);
      // Face 2 consists of edges 0, 2, 4
      s = DataType_(0.5)*(l[0] + l[2] + l[4]);
      vol = trafo.compute_vol<Shape::Simplex<2>, DataType_>(Index(2));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(s * (s - l[0]) * (s - l[2]) * (s - l[4]) ), eps);
      // Face 3 consists of edges 0, 1, 3
      s = DataType_(0.5)*(l[0] + l[1] + l[3]);
      vol = trafo.compute_vol<Shape::Simplex<2>, DataType_>(Index(3));
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, Math::sqrt(s * (s - l[0]) * (s - l[1]) * (s - l[3]) ), eps);
    }
};

StandardTrafoVolumeTest<float> standard_trafo_volume_test_float;
StandardTrafoVolumeTest<double> standard_trafo_volume_test_double;
