#include <test_system/test_system.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/mesh_attribute.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/target_set.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

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
  virtual void run() const override
  {
    RefineFactory< MeshType, UnitCubeFactory> my_factory(3);
    MeshType my_mesh(my_factory);

    BoundaryFactory<MeshType>my_boundary_factory(my_mesh);
    MeshPartType mesh_part_reference(my_boundary_factory);
    MeshPartType mesh_part_to_test(my_boundary_factory);

    // Test bottom to top deduction
    // Delete parent information for all shape types except for edges
    TargetSet tmp0;
    mesh_part_to_test.get_target_set<0>() = std::move(tmp0);
    TargetSet tmp2;
    mesh_part_to_test.get_target_set<2>() = std::move(tmp2);
    TargetSet tmp3;
    mesh_part_to_test.get_target_set<3>() = std::move(tmp3);

    // Get lotsa references for target sets
    TargetSet& ts_test0 = mesh_part_to_test.get_target_set<0>();
    TargetSet& ts_reference0 = mesh_part_reference.get_target_set<0>();
    TargetSet& ts_test1 = mesh_part_to_test.get_target_set<1>();
    TargetSet& ts_reference1 = mesh_part_reference.get_target_set<1>();
    TargetSet& ts_test2 = mesh_part_to_test.get_target_set<2>();
    TargetSet& ts_reference2 = mesh_part_reference.get_target_set<2>();
    TargetSet& ts_test3 = mesh_part_to_test.get_target_set<3>();
    TargetSet& ts_reference3 = mesh_part_reference.get_target_set<3>();

    // Deduct target sets
    mesh_part_to_test.deduct_target_sets_from_bottom<0>(my_mesh.get_index_set_holder());

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
    TargetSet tmp00;
    mesh_part_to_test.get_target_set<0>() = std::move(tmp00);
    TargetSet tmp01;
    mesh_part_to_test.get_target_set<1>() = std::move(tmp01);
    TargetSet tmp03;
    mesh_part_to_test.get_target_set<3>() = std::move(tmp03);

    mesh_part_to_test.deduct_target_sets_from_top<2>(my_mesh.get_index_set_holder());

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
    MeshAttributeType* att_0_1 = new MeshAttributeType(mesh_part_to_test.get_num_entities(0),1,1);
    MeshAttributeType* att_0_2 = new MeshAttributeType(mesh_part_to_test.get_num_entities(0),1,1);
    MeshAttributeType* att_0_3 = new MeshAttributeType(mesh_part_to_test.get_num_entities(0),1,1);
    // Change some data
    (*att_0_1)[1][0] = DataType(12);
    (*att_0_3)[1][0] = DataType(-5);

    // Add attributes of dimension 0
    TEST_CHECK(mesh_part_to_test.get_num_attributes() == 0);
    TEST_CHECK(mesh_part_to_test.add_attribute(att_0_1,"SomeName"));
    TEST_CHECK(mesh_part_to_test.add_attribute(att_0_2,"SomeOtherName"));

    TEST_CHECK(mesh_part_to_test.get_num_attributes() == 2);

    // As an attribute with the same name is already present, this must return false
    TEST_CHECK(!mesh_part_to_test.add_attribute(att_0_3,"SomeName"));
    // There should still be only 2 attributes in the mesh
    TEST_CHECK_EQUAL(mesh_part_to_test.get_num_attributes(),2);
    TEST_CHECK_EQUAL((*mesh_part_to_test.find_attribute("SomeName"))[1][0], DataType(12));
    // There should still be only 2 attributes in the mesh
    TEST_CHECK_EQUAL(mesh_part_to_test.get_num_attributes(),2);

    // Since we never successfully inserted this, we have to clean it up
    delete att_0_3;

  }

} mesh_part_test;
