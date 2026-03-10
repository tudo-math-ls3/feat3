// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/raycast.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

template<typename DT_>
class RaycastTest : public UnitTest
{
  using Vertex2D = Tiny::Vector<DT_, 2>;

  using Vertex3D = Tiny::Vector<DT_, 3>;

  using RIP = RayIntersectionPrimitives<DT_>;

  using HexMesh = ConformalMesh<Shape::Hexahedron, 3, DT_>;
  using QuadMesh = ConformalMesh<Shape::Quadrilateral, 2, DT_>;
  using TriMesh = ConformalMesh<Shape::Triangle, 2, DT_>;
  using TetMesh = ConformalMesh<Shape::Tetrahedron, 3, DT_>;
  using EdgeMesh = ConformalMesh<Shape::Hypercube<1>, 1, DT_>;

public:
  RaycastTest() : UnitTest("RaycastTest", Type::Traits<DT_>::name())
  {
  }

  ~RaycastTest() override = default;

  void run() const override
  {
    test_segment_intersection();
    test_triangle_intersection();
    test_facet_raycast();
    test_vertex_raycast();
  }

  void test_segment_intersection() const
  {
    // Test hit
    {
      Ray<Vertex2D> r{{0.0, 0.0}, {1.0, 0.0}};

      Vertex2D a{1.0, 1.0};
      Vertex2D b{1.0, -1.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_segment_intersection(r, a, b);
      TEST_CHECK(intersection);
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_entry, DT_(1.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_exit, DT_(1.0), DT_(1e-4));
    }

    // Test miss
    {
      Ray<Vertex2D> r{{0.0, 0.0}, {1.0, 0.0}};

      Vertex2D a{0.0, 1.0};
      Vertex2D b{0.0, 2.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_segment_intersection(r, a, b);

      TEST_CHECK(!intersection);
    }

    // Test miss behind
    {
      Ray<Vertex2D> r{{0.0, 0.0}, {1.0, 0.0}};

      Vertex2D a{-1.0, 1.0};
      Vertex2D b{-1.0, -1.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_segment_intersection(r, a, b);

      TEST_CHECK(!intersection);
    }

    // Test coplanar
    {
      Ray<Vertex2D> r{{0.0, 0.0}, {1.0, 0.0}};

      Vertex2D a{1.0, 0.0};
      Vertex2D b{2.0, 0.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_segment_intersection(r, a, b);
      TEST_CHECK(intersection);
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_entry, DT_(1.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_exit, DT_(2.0), DT_(1e-4));
    }

    // Test partial coplanar
    {
      Ray<Vertex2D> r{{0.0, 0.0}, {1.0, 0.0}};

      Vertex2D a{-2.0, 0.0};
      Vertex2D b{2.0, 0.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_segment_intersection(r, a, b);
      TEST_CHECK(intersection);
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_entry, DT_(0.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_exit, DT_(2.0), DT_(1e-4));
    }
  }

  void test_triangle_intersection() const
  {
    // Test hit
    {
      Ray<Vertex3D> r{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};

      Vertex3D a{-0.5, -0.5, 2.0};
      Vertex3D b{0.5, -0.5, 2.0};
      Vertex3D c{-0.5, 0.5, 2.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_triangle_intersection(r, a, b, c);
      TEST_CHECK(intersection);
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_entry, DT_(2.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_exit, DT_(2.0), DT_(1e-4));
    }

    // Test miss
    {
      Ray<Vertex3D> r{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};

      Vertex3D a{1.0, 1.0, 2.0};
      Vertex3D b{2.0, 1.0, 2.0};
      Vertex3D c{1.0, 2.0, 2.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_triangle_intersection(r, a, b, c);

      TEST_CHECK(!intersection);
    }

    // Test coplanar ray
    {
      Ray<Vertex3D> r{{0.0, 1.0, 2.0}, {1.0, 0.0, 0.0}};

      Vertex3D a{1.0, 1.0, 2.0};
      Vertex3D b{2.0, 1.0, 2.0};
      Vertex3D c{1.0, 2.0, 2.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_triangle_intersection(r, a, b, c);
      TEST_CHECK(intersection);
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_entry, DT_(1.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_exit, DT_(2.0), DT_(1e-4));
    }

    // Test anisotropic coplanar ray
    {
      Ray<Vertex3D> r{{0.0, 0.0, 0.0}, {1.0, 2.0, 0.0}};

      Vertex3D a{1.0, 2.0, 0.0};
      Vertex3D b{2.0, 4.0, 0.0};
      Vertex3D c{1.0, 4.0, 0.0};

      std::optional<IntersectionData<DT_>> intersection = RIP::ray_triangle_intersection(r, a, b, c);

      TEST_CHECK(intersection);
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_entry, DT_(1.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(intersection.value().t_exit, DT_(2.0), DT_(1e-4));
    }
  }

  void test_facet_raycast() const
  {
    // Hex mesh
    {
      UnitCubeFactory<HexMesh> factory;
      HexMesh mesh(factory);

      Ray<Vertex3D> r{{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};

      auto i = facet_raycast(mesh, r, 4);
      TEST_CHECK(i);
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_entry, DT_(1.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(1.0), DT_(1e-4));
    }

    // Hex mesh - diagonal I
    {
      UnitCubeFactory<HexMesh> factory;
      HexMesh mesh(factory);

      Vertex3D origin{-1.0, -1.0, 0.0};
      Vertex3D direction{1.0, 1.0, 0.0};
      direction.normalize();

      Ray<Vertex3D> r{origin, direction};

      auto i = facet_raycast(mesh, r, 0);
      TEST_CHECK(i);
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_entry, DT_(Math::sqrt(2.0)), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(2 * Math::sqrt(2.0)), DT_(1e-4));
    }

    // Hex mesh - diagonal II
    {
      UnitCubeFactory<HexMesh> factory;
      HexMesh mesh(factory);

      Vertex3D origin{2.0, -1.0, 0.0};
      Vertex3D direction{-1.0, 1.0, 0.0};
      direction.normalize();

      Ray<Vertex3D> r{origin, direction};

      auto i = facet_raycast(mesh, r, 0);
      TEST_CHECK(i);
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_entry, DT_(Math::sqrt(2.0)), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(2 * Math::sqrt(2.0)), DT_(1e-4));
    }

    // Tetrahedron
    {
      ReferenceCellFactory<Shape::Tetrahedron, DT_> factory;
      TetMesh mesh(factory);

      Vertex3D origin{0.25, 0.25, -1.0};
      Vertex3D direction{0.0, 0.0, 1.0};

      Ray<Vertex3D> r{origin, direction};

      auto i = facet_raycast(mesh, r, 3);
      TEST_CHECK(i);
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_entry, DT_(1.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(1.0), DT_(1e-4));
    }

    // Quad
    {
      UnitCubeFactory<QuadMesh> factory;
      QuadMesh mesh(factory);

      Vertex2D origin{-1.0, 0.5};
      Vertex2D direction{1.0, 0.0};

      Ray<Vertex2D> r{origin, direction};

      auto i = facet_raycast(mesh, r, 2);
      TEST_CHECK(i);
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_entry, DT_(1.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(1.0), DT_(1e-4));

      i = facet_raycast(mesh, r, 3);
      TEST_CHECK(i);
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_entry, DT_(2.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(2.0), DT_(1e-4));
    }

    // Triangle
    {
      ReferenceCellFactory<Shape::Triangle, DT_> factory;
      TriMesh mesh(factory);

      Vertex2D origin{-1.0, 0.5};
      Vertex2D direction{1.0, 0.0};

      Ray<Vertex2D> r{origin, direction};

      auto i = facet_raycast(mesh, r, 1);
      TEST_CHECK(i);
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_entry, DT_(1.0), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(1.0), DT_(1e-4));

      i = facet_raycast(mesh, r, 0);
      TEST_CHECK(i);
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_entry, DT_(1.5), DT_(1e-4));
      TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(1.5), DT_(1e-4));
    }
  }

  void test_vertex_raycast() const
  {
    RefinedUnitCubeFactory<HexMesh> factory(1);
    HexMesh mesh(factory);

    VertexRaycaster<HexMesh> raycaster(mesh);

    // Cast into cell
    auto i = raycaster.cast(0, {1.0, 1.0, 1.0});

    // Cast should produce a hit
    TEST_CHECK(i);
    // Hit should not stop at start vertex
    TEST_CHECK(i.value().t_exit > DT_(0.0));

    // Cast into outside should produce no hit
    i = raycaster.cast(0, {-1.0, -1.0, -1.0});
    TEST_CHECK(!i);

    // Cast along face in z-direction
    i = raycaster.cast(0, {0.0, 0.0, 1.0});
    // Should hit
    TEST_CHECK(i);
    // Should hit face half a unit above
    TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(0.5), DT_(1e-4));

    // Cast along face in y-direction
    i = raycaster.cast(0, {0.0, 1.0, 0.0});
    // Should hit
    TEST_CHECK(i);
    // Should hit face half a unit above
    TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(0.5), DT_(1e-4));

    // Cast along face in x-direction
    i = raycaster.cast(0, {1.0, 0.0, 0.0});
    // Should hit
    TEST_CHECK(i);
    // Should hit face half a unit above
    TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(0.5), DT_(1e-4));

    // Diagonal cast along face
    Vertex3D dir{1.0, 0.0, 1.0};
    i = raycaster.cast(0, dir.normalize());
    // Cast should hit
    TEST_CHECK(i);
    // Hit should hit middle point of face
    TEST_CHECK_EQUAL_WITHIN_EPS(i.value().t_exit, DT_(Math::sqrt(0.5)), DT_(1e-4));
  }
};

static const RaycastTest<double> raycast_test_double;
static const RaycastTest<float> raycast_test_float;
