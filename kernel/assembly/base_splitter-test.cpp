#include <test_system/test_system.hpp>
#include <kernel/assembly/base_splitter.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/util/random.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename DT_, typename IT_>
class BaseSplitterTest :
  public TestSystem::FullTaggedTest<Mem::Main, DT_, IT_>
{
public:
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorType;

public:
  BaseSplitterTest() :
    TestSystem::FullTaggedTest<Mem::Main, DT_, IT_>("BaseSplitterTest")
  {
  }

  virtual void run() const override
  {
    // create root mesh node
    RootMeshNodeType* root_mesh_node = nullptr;
    {
      Geometry::RefinedUnitCubeFactory<MeshType> factory(2);
      root_mesh_node = new RootMeshNodeType(new MeshType(factory));
    }

    // create base cell splitting
    root_mesh_node->create_base_splitting();

    // refine mesh node
    for(int i(0); i < 2; ++i)
    {
      auto* old_node = root_mesh_node;
      root_mesh_node = old_node->refine();
      delete old_node;
    }

    {
      // create trafo
      TrafoType trafo(*root_mesh_node->get_mesh());

      // test various spaces
      test_space<Space::Lagrange1::Element<TrafoType>>(trafo, *root_mesh_node);
      test_space<Space::Lagrange2::Element<TrafoType>>(trafo, *root_mesh_node);
      test_space<Space::CroRavRanTur::Element<TrafoType>>(trafo, *root_mesh_node);
    }

    // delete root mesh node
    delete root_mesh_node;
  }

  template<typename Space_>
  void test_space(TrafoType& trafo, const RootMeshNodeType& root_mesh_node) const
  {
    // create space
    Space_ space(trafo);

    // create two vectors
    VectorType vec_ref(space.get_num_dofs());
    VectorType vec_tst(space.get_num_dofs());

    // fill vectors with random values
    Random rng;
    for(Index i(0); i < space.get_num_dofs(); ++i)
    {
      vec_ref(i, rng(-DT_(1), DT_(1)));
      vec_tst(i, rng( DT_(2), DT_(5)));
    }

    // create base-splitter
    Assembly::BaseSplitter<Space_, DT_, IT_> splitter(space, root_mesh_node);

    // split our vector
    std::vector<VectorType> split_vec;
    splitter.split(split_vec, vec_ref);

    // and join it
    splitter.join(vec_tst, split_vec);

    // subtract
    vec_tst.axpy(vec_ref, vec_tst, -DT_(1));

    // compute norm
    DT_ def = vec_tst.norm2();

    // and check
    TEST_CHECK_EQUAL_WITHIN_EPS(def, DT_(0), 1E-8);
  }
};

BaseSplitterTest<double, Index> base_splitter_test_double_index;
