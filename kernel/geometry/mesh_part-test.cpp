#include <test_system/test_system.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/mesh_attribute.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/target_set.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;

bool check_target_set_consistency(const TargetSet& ts_1, const TargetSet& ts_2)
{
  for(Index i(0); i < ts_1.get_num_entities(); ++i)
  {
    Index parent_reference(ts_1[i]);

    bool found(false);
    Index j(0);
    while(!found && j < ts_2.get_num_entities())
    {
      if(ts_2[j] == parent_reference)
        found = true;

      ++j;
    }
    if(!found )
      return false;
  }
  return true;
}

/**
 * \brief Test class for the MeshPart class template
 *
 * \test Tests the parent set and topology deduction routines. Refinement is already covered in MeshNodeTestConfQuad
 *
 * \author Jordi Paul
 */
class MeshPartTest
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  MeshPartTest() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("mesh_part-test")
  {
  }

  typedef double DataType;
  typedef ConformalMesh<Shape::Hypercube<3>, 3, 3, DataType> MeshType;
  typedef MeshPart<MeshType> MeshPartType;
  typedef MeshAttribute<DataType> MeshAttributeType;

  template<Index dim_>
  using TargetSetType = typename MeshPartType::TargetSet<dim_>::Type;

  /**
   * This constructs a mesh and its boundary MeshPart using the BoundaryFactory. Then some TargetSets get deleted
   * and re-computed
   */
  virtual void run() const
  {
    RefineFactory< MeshType, UnitCubeFactory> my_factory(3);
    MeshType my_mesh(my_factory);

    BoundaryFactory<MeshType>my_boundary_factory(my_mesh);
    MeshPartType mesh_part_reference(my_boundary_factory);
    MeshPartType mesh_part_test(mesh_part_reference);

    // Test bottom to top deduction
    // Delete parent information for all shape types except for edges
    TargetSet tmp0;
    mesh_part_test.get_target_set<0>() = std::move(tmp0);
    mesh_part_test.get_target_set<2>() = std::move(tmp0);
    mesh_part_test.get_target_set<3>() = std::move(tmp0);

    // Get lotsa references for target sets
    TargetSet& ts_test0 = mesh_part_test.get_target_set<0>();
    TargetSet& ts_reference0 = mesh_part_reference.get_target_set<0>();
    TargetSet& ts_test1 = mesh_part_test.get_target_set<1>();
    TargetSet& ts_reference1 = mesh_part_reference.get_target_set<1>();
    TargetSet& ts_test2 = mesh_part_test.get_target_set<2>();
    TargetSet& ts_reference2 = mesh_part_reference.get_target_set<2>();
    TargetSet& ts_test3 = mesh_part_test.get_target_set<3>();
    TargetSet& ts_reference3 = mesh_part_reference.get_target_set<3>();

    // Deduct target sets
    mesh_part_test.compute_target_sets_from_bottom<0>(my_mesh);

    // Check for dimension 0: new one should still be empty
    TEST_CHECK_MSG(ts_test0.get_num_entities() == 0, "num_entities for dimension 0 should be 0!");

    // Check dimension 1: no changes
    TEST_CHECK_MSG(ts_test1.get_num_entities() == ts_reference1.get_num_entities(), "num_entities mismatch for dimension 1");
    TEST_CHECK_MSG(check_target_set_consistency(ts_reference1, ts_test1), "ts_test did not contain all entities from ts_reference for dimension 1!");

    // Check dimension 2
    TEST_CHECK_MSG(ts_test2.get_num_entities() == ts_reference2.get_num_entities(), "num_entities mismatch for dimension 2");
    TEST_CHECK_MSG(check_target_set_consistency(ts_reference2, ts_test2), "ts_test did not contain all entities from ts_reference for dimension 2!");

    // Check dimension 3
    TEST_CHECK_MSG(ts_test3.get_num_entities() == ts_reference3.get_num_entities(), "num_entities mismatch for dimension 3");
    TEST_CHECK_MSG(check_target_set_consistency(ts_reference3, ts_test3), "ts_test did not contain all entities from ts_reference for dimension 3!");

    // Test top to bottom deduction
    // Delete parent information for all shape types except for faces
    mesh_part_test.get_target_set<0>() = std::move(tmp0);
    mesh_part_test.get_target_set<1>() = std::move(tmp0);
    mesh_part_test.get_target_set<3>() = std::move(tmp0);

    mesh_part_test.compute_target_sets_from_top<2>(my_mesh);

    // Check dimension 3: new one should still be empty
    TEST_CHECK_MSG(ts_test3.get_num_entities() == 0, "num_entities for dimension 3 should be 3!");

    // Check dimension 2: No changes
    TEST_CHECK_MSG(ts_test2.get_num_entities() == ts_reference2.get_num_entities(), "num_entities mismatch for dimension 2");
    TEST_CHECK_MSG(check_target_set_consistency(ts_reference2, ts_test2), "ts_test did not contain all entities from ts_reference for dimension 2!");

    // Check dimension 1
    TEST_CHECK_MSG(ts_test1.get_num_entities() == ts_reference1.get_num_entities(), "num_entities mismatch for dimension 1");
    TEST_CHECK_MSG(check_target_set_consistency(ts_reference1, ts_test1), "ts_test did not contain all entities from ts_reference for dimension 1!");

    // Check dimension 0
    TEST_CHECK_MSG(ts_test0.get_num_entities() == ts_reference0.get_num_entities(), "num_entities mismatch for dimension 0");
    TEST_CHECK_MSG(check_target_set_consistency(ts_reference0, ts_test0), "ts_test did not contain all entities from ts_reference for dimension 0!");

    // Check MeshAttribute functionality
    // Create attributes of dimension 0
    MeshAttributeType att_0_1(mesh_part_test.get_num_entities(0),1,1,"SomeName");
    MeshAttributeType att_0_2(mesh_part_test.get_num_entities(0),1,1,"SomeOtherName");
    MeshAttributeType att_0_3(mesh_part_test.get_num_entities(0),1,1,"SomeName");
    // Change some data
    att_0_1[1][0] = DataType(12);
    att_0_3[1][0] = DataType(-5);

    // Add attributes of dimension 0
    TEST_CHECK(mesh_part_test.get_num_attributes() == 0);
    TEST_CHECK(mesh_part_test.add_attribute(att_0_1,0));
    TEST_CHECK(mesh_part_test.add_attribute(att_0_2,0));

    TEST_CHECK(mesh_part_test.get_num_attributes() == 2);

    // As an attribute with the same name is already present, this must return false
    TEST_CHECK(!mesh_part_test.add_attribute(att_0_3,0));
    // There should still be only 2 attributes in the mesh
    TEST_CHECK_EQUAL(mesh_part_test.get_num_attributes(0),2);
    TEST_CHECK_EQUAL((*mesh_part_test.find_attribute("SomeName",0))[1][0], DataType(12));
    // Now we specify replace=true and it should return true
    TEST_CHECK(mesh_part_test.add_attribute(att_0_3,0,true));
    TEST_CHECK_EQUAL((*mesh_part_test.find_attribute("SomeName",0))[1][0], DataType(-5));
    // There should still be only 2 attributes in the mesh
    TEST_CHECK_EQUAL(mesh_part_test.get_num_attributes(0),2);

    // Add attribute of dimension 0 as attribute of dimension 1. This should throw an InternalError.
    TEST_CHECK_THROWS(mesh_part_test.add_attribute(att_0_2,1), InternalError);
  }

} mesh_part_test;
