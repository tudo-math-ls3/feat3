#include <kernel/runtime.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/atlas/circle.hpp>
#include <kernel/geometry/atlas/sphere.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/geometry/export_vtk.hpp>

namespace DbgSlipAsm
{
  using namespace FEAT;

  void main_2d()
  {
    typedef double DataType;
    static constexpr int dim = 2;
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::MeshPart<MeshType> MeshPartType;
    typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
    typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

    typedef LAFEM::DenseVectorBlocked<DataType, Index, dim> VectorType;
    typedef LAFEM::SlipFilter<DataType, Index, dim> SlipFilterType;

    Geometry::MeshAtlas<MeshType> atlas;
    //atlas.add_mesh_chart("circle", std::make_unique<Geometry::Atlas::Circle<MeshType>>(0.0, 0.0, 1.0));
    atlas.add_mesh_chart("circle", std::unique_ptr<Geometry::Atlas::Circle<MeshType>>(new Geometry::Atlas::Circle<MeshType>(0.0, 0.0, 1.0)));
    const auto* chart = atlas.find_mesh_chart("circle");

    std::unique_ptr<MeshNodeType> mesh_node;
    {
      Geometry::UnitCubeFactory<MeshType> factory;
      mesh_node = MeshNodeType::make_unique(factory.make_unique(), &atlas);
      MeshType* mesh = mesh_node->get_mesh();
      auto& vtx = mesh->get_vertex_set();
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        (vtx[i][0] -= 0.5) *= 2.0;
        (vtx[i][1] -= 0.5) *= 2.0;
      }
      Geometry::BoundaryFactory<MeshType> bnd_factory(*mesh);
      mesh_node->add_mesh_part("bnd", bnd_factory.make_unique(), "circle", chart);
      mesh_node->adapt();
    }

    DataType prev_mean(0.0), prev_max(0.0);
    for(int lvl(0); lvl < 5; ++lvl)
    {
      auto ref_mesh_node = mesh_node->refine_unique();

      const MeshPartType* part = mesh_node->find_mesh_part("bnd");

      TrafoType trafo(*mesh_node->get_mesh());
      trafo.add_meshpart_chart(*mesh_node->find_mesh_part("bnd"), *chart);

      SpaceType space(trafo);

      VectorType vector;
      auto function = Analytic::create_lambda_function_vector_2d([](auto, auto y) {return -y;}, [](auto x, auto) {return x;});
      Assembly::Interpolator::project(vector, function, space);

      VectorType vec_fil = vector.clone(LAFEM::CloneMode::Layout);
      vec_fil.format();

      SlipFilterType filter;
      Assembly::SlipFilterAssembler<TrafoType> slip_asm(trafo);
      slip_asm.add_mesh_part(*part);
      slip_asm.assemble(filter, space);

      const auto* filter_idx = filter.get_indices();
      const auto* filter_val = filter.get_values();
      const auto* vector_val = vector.elements();
      auto* vecfil_val = vec_fil.elements();

      std::vector<DataType> vdot(vector.size(), 0.0);
      DataType rsum(0.0), rmax(0.0);
      for(Index i(0); i < filter.used_elements(); ++i)
      {
        DataType t = Math::abs(Tiny::dot(vector_val[filter_idx[i]], filter_val[i]));
        rsum += t;
        rmax = Math::max(rmax, t);
        vecfil_val[filter_idx[i]] = filter_val[i];
        vdot[filter_idx[i]] = t;
      }
      DataType rmean = rsum / DataType(filter.used_elements());
      std::cout << lvl << ": " << stringify(filter.used_elements()).pad_front(5) << " of " << stringify(filter.size()).pad_front(7) << ": ";
      //std::cout << stringify_fp_fix(rmax, 5, 7) << " / " << stringify_fp_fix(rmean, 5, 7);
      std::cout << stringify_fp_sci(rmax, 5, 10) << " / " << stringify_fp_sci(rmean, 5, 10);
      if(prev_max > 0.0)
        std::cout << " / " << stringify_fp_fix(prev_max / rmax, 3, 7) << " / " << stringify_fp_fix(prev_mean / rmean, 3, 7);
      std::cout << std::endl;
      prev_max = rmax;
      prev_mean = rmean;

      {
        Geometry::ExportVTK<MeshType> writer(*ref_mesh_node->get_mesh());
        writer.add_vertex_vector("vector", vector);
        writer.add_vertex_vector("filter", vec_fil);
        writer.add_vertex_scalar("v_dot_f", vdot.data());
        writer.write("dbg-slip-asm-2d." + stringify(lvl));
      }

      mesh_node = std::move(ref_mesh_node);
    }
  }

  void main_3d()
  {
    typedef double DataType;
    static constexpr int dim = 3;
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::MeshPart<MeshType> MeshPartType;
    typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    //typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    typedef LAFEM::DenseVectorBlocked<DataType, Index, dim> VectorType;
    typedef LAFEM::SlipFilter<DataType, Index, dim> SlipFilterType;

    Geometry::MeshAtlas<MeshType> atlas;
    //atlas.add_mesh_chart("circle", std::make_unique<Geometry::Atlas::Circle<MeshType>>(0.0, 0.0, 1.0));
    //atlas.add_mesh_chart("sphere", std::make_unique<Geometry::Atlas::Sphere<MeshType>>(0.0, 0.0, 0.0, 1.0));
    atlas.add_mesh_chart("sphere", std::unique_ptr<Geometry::Atlas::Sphere<MeshType>>(new Geometry::Atlas::Sphere<MeshType>(0.0, 0.0, 0.0, 1.0)));
    const auto* chart = atlas.find_mesh_chart("sphere");

    std::unique_ptr<MeshNodeType> mesh_node;
    {
      Geometry::UnitCubeFactory<MeshType> factory;
      mesh_node = MeshNodeType::make_unique(factory.make_unique(), &atlas);
      MeshType* mesh = mesh_node->get_mesh();
      auto& vtx = mesh->get_vertex_set();
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        (vtx[i][0] -= 0.5) *= 2.0;
        (vtx[i][1] -= 0.5) *= 2.0;
        (vtx[i][2] -= 0.5) *= 2.0;
      }
      Geometry::BoundaryFactory<MeshType> bnd_factory(*mesh);
      mesh_node->add_mesh_part("bnd", bnd_factory.make_unique(), "sphere", chart);
      mesh_node->adapt();
    }

    DataType prev_mean(0.0), prev_max(0.0);
    for(int lvl(0); lvl < 4; ++lvl)
    {
      auto ref_mesh_node = mesh_node->refine_unique();

      const MeshPartType* part = mesh_node->find_mesh_part("bnd");

      TrafoType trafo(*mesh_node->get_mesh());
      //trafo.add_meshpart_chart(*mesh_node->find_mesh_part("bnd"), *chart);

      SpaceType space(trafo);

      VectorType vector;
      auto function = Analytic::create_lambda_function_vector_3d(
        [](auto, auto y, auto z) {return y - z;},
        [](auto x, auto, auto z) {return z - x;},
        [](auto x, auto y, auto) {return x - y;});
      Assembly::Interpolator::project(vector, function, space);

      VectorType vec_fil = vector.clone(LAFEM::CloneMode::Layout);
      vec_fil.format();

      SlipFilterType filter;
      Assembly::SlipFilterAssembler<TrafoType> slip_asm(trafo);
      slip_asm.add_mesh_part(*part);
      slip_asm.assemble(filter, space);

      const auto* filter_idx = filter.get_indices();
      const auto* filter_val = filter.get_values();
      const auto* vector_val = vector.elements();
      auto* vecfil_val = vec_fil.elements();

      std::vector<DataType> vdot(vector.size(), 0.0);

      DataType rsum(0.0), rmax(0.0);
      for(Index i(0); i < filter.used_elements(); ++i)
      {
        DataType t = Math::abs(Tiny::dot(vector_val[filter_idx[i]], filter_val[i]));
        rsum += t;
        rmax = Math::max(rmax, t);
        vecfil_val[filter_idx[i]] = filter_val[i];
        vdot[filter_idx[i]] = t;
      }
      DataType rmean = rsum / DataType(filter.used_elements());
      std::cout << lvl << ": " << stringify(filter.used_elements()).pad_front(5) << " of " << stringify(filter.size()).pad_front(7) << ": ";
      //std::cout << stringify_fp_fix(rmax, 5, 7) << " / " << stringify_fp_fix(rmean, 5, 7);
      std::cout << stringify_fp_sci(rmax, 5, 10) << " / " << stringify_fp_sci(rmean, 5, 10);
      if(prev_max > 0.0)
        std::cout << " / " << stringify_fp_fix(prev_max / rmax, 3, 7) << " / " << stringify_fp_fix(prev_mean / rmean, 3, 7);
      std::cout << std::endl;
      prev_max = rmax;
      prev_mean = rmean;

      {
        Geometry::ExportVTK<MeshType> writer(*mesh_node->get_mesh());
        writer.add_vertex_vector("vector", vector);
        writer.add_vertex_vector("filter", vec_fil);
        writer.add_vertex_scalar("v_dot_f", vdot.data());
        writer.write("dbg-slip-asm-3d." + stringify(lvl));
      }

      mesh_node = std::move(ref_mesh_node);
    }
  }
}

int main(int argc, char** argv)
{
  FEAT::Runtime::initialize(argc, argv);
  DbgSlipAsm::main_2d();
  DbgSlipAsm::main_3d();
  return FEAT::Runtime::finalize();
}
