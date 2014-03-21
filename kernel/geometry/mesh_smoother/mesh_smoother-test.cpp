#include <test_system/test_system.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1.hpp>

#include <iostream>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename DataType_, typename ShapeType_>
class MeshSmootherTest_2d
: public TestSystem::TaggedTest<Archs::None, DataType_>
{
  public:
    MeshSmootherTest_2d() :
      TestSystem::TaggedTest<Archs::None, DataType_>("mesh_smoother_test")
  {
  }

    virtual void run() const
    {
      typedef DataType_ DataType;
      typedef Mem::Main MemoryType;
      typedef ShapeType_ ShapeType;

      typedef typename Geometry::ConformalMesh<ShapeType, 2, 2, DataType> MeshType;
      //typedef typename Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
      typedef typename Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef typename Geometry::RumpfFunctional<MemoryType, DataType, TrafoType> FunctionalType;

      typedef typename Geometry::CellSubSet<ShapeType> BoundaryType;
      typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

      // Mesh and trafo
      Index level(3);
      Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
      MeshType mesh(mesh_factory);
      TrafoType trafo(mesh);
      // Boundary stuff
      BoundaryFactoryType boundary_factory(mesh);
      BoundaryType boundary(boundary_factory);

      Geometry::TargetSet boundary_set = boundary.template get_target_set<0>();
      int bdry_id[mesh.get_num_entities(0)];
      for(Index i(0); i < mesh.get_num_entities(0); ++i)
        bdry_id[i] = 0;

      for(Index i(0); i < boundary.get_num_entities(0); ++i)
        bdry_id[boundary_set[i]] = -1;

      DataType pi(Math::pi<DataType>());
      DataType deltat(DataType(0.025));

      DataType fac_norm = DataType(1),fac_det = DataType(1), fac_cof = DataType(0), fac_reg(DataType(1e-8));
      FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

      Geometry::RumpfSmoother<FunctionalType, TrafoType, DataType, MemoryType> rumpflpumpfl(trafo, my_functional);
      rumpflpumpfl.init();

      for(Index i(0); i < boundary.get_num_entities(0); ++i)
      {
        Index j = boundary_set[i];
        DataType tmp0 = rumpflpumpfl._coords[0](j);
        DataType tmp1 = rumpflpumpfl._coords[1](j);
        rumpflpumpfl._coords[0](j, tmp0 - ( Math::sin(DataType(2)*pi*tmp1) )/DataType(1 << (level+2)));
        rumpflpumpfl._coords[1](j, tmp1 + ( Math::sin(DataType(2)*pi*tmp0) )/DataType(1 << (level+2)));
      }
      rumpflpumpfl.set_coords();

      rumpflpumpfl.optimise();

      typename MeshType::VertexSetType& vtx = mesh.get_vertex_set();
      for(Index i(0); i < mesh.get_num_entities(0); ++i)
      {
        if(bdry_id[i]==0)
        {
          DataType tmp0 = vtx[i][0];
          DataType tmp1 = vtx[i][1];
          //vtx[i][0]+=deltat*Math::sin(pi*(DataType(0.5)-vtx[i][1]));
          //vtx[i][1]-=deltat*Math::cos(pi*(DataType(0.5)-tmp));
          vtx[i][0]+=deltat*( DataType(1) - Math::cos(DataType(2)*pi*tmp0))*Math::sin(DataType(2)*pi*tmp1);
          vtx[i][1]-=deltat*( DataType(1) - Math::cos(DataType(2)*pi*tmp0))*Math::sin(DataType(2)*pi*tmp1);
        }
      }

    }
};
MeshSmootherTest_2d<float,Shape::Simplex<2>> mesh_smoother_test_float_simplex;
MeshSmootherTest_2d<double,Shape::Simplex<2>> mesh_smoother_test_double_simplex;
MeshSmootherTest_2d<float,Shape::Hypercube<2>> mesh_smoother_test_float_hypercube;
MeshSmootherTest_2d<double,Shape::Hypercube<2>> mesh_smoother_test_double_hypercube;
