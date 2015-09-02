#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/meshopt/dudv_smoother.hpp>
#include <kernel/meshopt/laplace_smoother.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for linear variational mesh smoothers, like the LaplaceSmoother
 */
template
<
  typename Mem_,
  typename DT_,
  typename IT_,
  typename ShapeType_,
  template<typename, typename, typename, typename> class Smoother_
>
class LinearVariationalSmootherTest
: public TestSystem::FullTaggedTest<Mem_, DT_, IT_>
{
  public:
    typedef ShapeType_ ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DT_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Smoother_<Mem_, DT_, IT_, TrafoType> SmootherType;
    /// Refinementlevel for the RefineFactory
    Index level;

    LinearVariationalSmootherTest(int level_) :
      TestSystem::FullTaggedTest<Mem_, DT_, IT_>("linear_variational_smoother_test"),
      level(Index(level_))
      {
      }

    /**
     * \brief Checks if the smoother reproduces the original mesh
     *
     * A refined unit cube mesh is generated and scaled by 0.75. Solving the appropriate PDE problem on this mesh
     * with the original boundary coordinates as Dirichlet boundary conditions should reproduce the original mesh
     * for the LaplaceSmoother and for the DuDvSmoother
     */
    virtual void run() const
    {
      Geometry::RefineFactory<MeshType, Geometry::UnitCubeFactory> mesh_factory(level);
      MeshType mesh(mesh_factory);
      TrafoType trafo(mesh);

      typedef typename MeshType::VertexSetType::VertexType VertexType;
      VertexType* coords_org(new VertexType[mesh.get_num_entities(0)]);

      auto& vtx = mesh.get_vertex_set();

      // Back up original coordinates and scale the mesh's vertices.
      for(Index i(0); i < mesh.get_num_entities(0); ++i)
      {
        coords_org[i] = vtx[i];
        vtx[i] = DT_(0.75)*vtx[i];
      }

      // Get the boundary so we can set the smoother's _coords field to the original values
      Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
      Geometry::MeshPart<MeshType> boundary(boundary_factory);

      // Create the smoother
      SmootherType my_smoother(trafo);

      // Set _coords to the original values on the boundary
      for(Index i(0); i < boundary.get_num_entities(0); ++i)
      {
        Index j(boundary.template get_target_set<0>()[i]);
        my_smoother._coords(j, coords_org[j]);
      }

      // Smooth the mesh
      my_smoother.optimise();

      // Check results
      const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

      for(Index i(0); i < mesh.get_num_entities(0); ++i)
      {
        for(int d(0); d < MeshType::world_dim; ++d)
          TEST_CHECK_EQUAL_WITHIN_EPS(my_smoother._coords(i)(d), coords_org[i][d], tol);
      }

      // Clean up
      delete[] coords_org;

    }


};

LinearVariationalSmootherTest<Mem::Main, double, unsigned int, Shape::Simplex<2>, Meshopt::LaplaceSmoother> laplace_d_s2(4);
#ifdef FEAST_BACKENDS_CUDA
LinearVariationalSmootherTest<Mem::CUDA, float, Index, Shape::Simplex<3>, Meshopt::LaplaceSmoother> laplace_f_s3(3);
LinearVariationalSmootherTest<Mem::CUDA, float, unsigned int, Shape::Hypercube<2>, Meshopt::LaplaceSmoother> laplace_f_hc2(5);
#endif
LinearVariationalSmootherTest<Mem::Main, double, unsigned int, Shape::Hypercube<3>, Meshopt::LaplaceSmoother> laplace_d_hc3(3);

LinearVariationalSmootherTest<Mem::Main, double, unsigned int, Shape::Simplex<2>, Meshopt::DuDvSmoother> dudv_f_s2(4);
LinearVariationalSmootherTest<Mem::Main, double, unsigned int, Shape::Simplex<3>, Meshopt::DuDvSmoother> dudv_f_s3(2);

LinearVariationalSmootherTest<Mem::Main, float, Index, Shape::Hypercube<2>, Meshopt::DuDvSmoother> dudv_f_hc2(5);
LinearVariationalSmootherTest<Mem::Main, float, Index, Shape::Hypercube<3>, Meshopt::DuDvSmoother> dudv_f_hc3(3);
