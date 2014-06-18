#include <kernel/geometry/hit_test_factory.hpp>
#include <test_system/test_system.hpp>
#include <sstream>
#include <kernel/util/string.hpp>                          // for String

// FEAST-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/cell_sub_set.hpp>                // for CellSubSet
#include <kernel/geometry/conformal_factories.hpp>         // for RefinedUnitCubeFactor
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;

/**
 * \brief Test class for the HitTestFactoryTest class.
 *
 * \test Tests the HitTestFactoryTest class.
 *
 * \author Stefan Wahlers
 */

class HitTestFactoryTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  HitTestFactoryTest() :
    TaggedTest<Archs::None, Archs::None>("hit_test_factory-test")
  {
  }

  void test_0() const
  {
    typedef Shape::Quadrilateral ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::CellSubSet<ShapeType> BoundaryType;
    typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactoryType;

    // some parameters for this special test
    Index level(6);
    double radius(0.3);
    Tiny::Vector<double, 2> mid_point(0.5);

    // create a meshfactory
    MeshFactoryType mesh_factory(level);
    MeshType mesh(mesh_factory);

    // create the hittest
    SphereHitTestFunction<double,Index(2)> hit_func(mid_point,radius);
    HitTestFactory<SphereHitTestFunction<double,Index(2)>, MeshType> hit_test(hit_func, mesh);
    BoundaryType sphere(hit_test);

    // helper vectors to visualize the hittest
    std::vector<double> nodes(mesh.get_num_entities(0),0.0);
    std::vector<double> quads(mesh.get_num_entities(2),0.0);
    // loop over all nodes caught by the hittest
    for (Index i(0); i < sphere.get_num_entities(0); ++i)
    {
      nodes[sphere.get_target_set<0>()[i]] = 1;
    }
    // loop over all quadrilaterals caught by the hittest
    for (Index i(0); i < sphere.get_num_entities(2); ++i)
    {
      quads[sphere.get_target_set<2>()[i]] = 1;
    }

    // export the helper vectors (uncomment this to create the vtk file)
    /*Geometry::ExportVTK<MeshType> exporter(mesh);
    exporter.add_scalar_vertex("nodes", nodes.data());
    exporter.add_scalar_cell("quads", quads.data());
    exporter.write("hit_test.vtk");*/
  } // test_0

  virtual void run() const
  {
    // run test #0
    test_0();
  }
} hit_test_factory_test;
