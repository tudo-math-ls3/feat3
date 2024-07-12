// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>

using namespace FEAT;

template<typename Mesh_>
void add_mesh_part(Geometry::ExportVTK<Mesh_>& vtk, const String& name, std::vector<double>& vtx_data, const Geometry::MeshPart<Mesh_>& meshpart)
{
  // get the vertex target set
  const Geometry::TargetSet& trg = meshpart.template get_target_set<0>();

  // mark all indexes vertices
  for(Index i(0); i < trg.get_num_entities(); ++i)
    vtx_data[trg[i]] = 1.0;

  // add variable
  vtk.add_vertex_scalar(name, vtx_data.data());

  // unmark vertices
  for(Index i(0); i < trg.get_num_entities(); ++i)
    vtx_data[trg[i]] = 0.0;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  Dist::Comm comm = Dist::Comm::world();
  if(comm.size() != 3)
  {
    comm.print("ERROR: This tool must be run with exactly 3 MPI processes!");
    return 1;
  }
  const int my_rank = comm.rank();

  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

  std::unique_ptr<MeshNodeType> mesh_node;
  {
    Geometry::RefinedUnitCubeFactory<MeshType> factory(2);
    mesh_node.reset(new MeshNodeType(factory.make_unique(), nullptr));
  }

  MeshType& mesh = *mesh_node->get_mesh();

  // we want to create the following domain:
  //
  // +---+
  // | 2 |
  // +---+---+
  // | 0 | 1 |
  // +---+---+

  // adjust vertex coordinates of our rank meshes
  {
    auto& vtx = mesh.get_vertex_set();
    Index nv = vtx.get_num_vertices();
    for(Index i(0); i < nv; ++i)
    {
      vtx[i][0] += double(((my_rank >> 0) & 1) - 1);
      vtx[i][1] += double(((my_rank >> 1) & 1) - 1);
    }
  }

  // generate halos
  switch(my_rank)
  {
  case 0:
    {
      auto hit_lambda_1 = [](const auto& p) {return p[0]+0.001 > 0.0;};
      auto hit_lambda_2 = [](const auto& p) {return p[1]+0.001 > 0.0;};
      Geometry::HitTestFactory<decltype(hit_lambda_1), MeshType> hit_factory_1(hit_lambda_1, mesh);
      Geometry::HitTestFactory<decltype(hit_lambda_2), MeshType> hit_factory_2(hit_lambda_2, mesh);
      mesh_node->add_halo(1, hit_factory_1.make_unique());
      mesh_node->add_halo(2, hit_factory_2.make_unique());
    }
    break;

  case 1:
    {
      auto hit_lambda_0 = [](const auto& p) {return p[0]-0.001 < 0.0;};
      auto hit_lambda_2 = [](const auto& p) {return (p[1]+0.001 > 0.0) && (p[0]-0.001 < 0.0);};
      Geometry::HitTestFactory<decltype(hit_lambda_0), MeshType> hit_factory_0(hit_lambda_0, mesh);
      Geometry::HitTestFactory<decltype(hit_lambda_2), MeshType> hit_factory_2(hit_lambda_2, mesh);
      mesh_node->add_halo(0, hit_factory_0.make_unique());
      mesh_node->add_halo(2, hit_factory_2.make_unique());
    }
    break;

  case 2:
    {
      auto hit_lambda_0 = [](const auto& p) {return p[1]-0.001 < 0.0;};
      auto hit_lambda_1 = [](const auto& p) {return (p[0]+0.001 > 0.0) && (p[1]-0.001 < 0.0);};
      Geometry::HitTestFactory<decltype(hit_lambda_0), MeshType> hit_factory_0(hit_lambda_0, mesh);
      Geometry::HitTestFactory<decltype(hit_lambda_1), MeshType> hit_factory_1(hit_lambda_1, mesh);
      mesh_node->add_halo(0, hit_factory_0.make_unique());
      mesh_node->add_halo(1, hit_factory_1.make_unique());
    }
    break;
  }

  // create normal boundary
  Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
  MeshPartType boundary_part(boundary_factory);

  // create masked boundary factory
  Geometry::MaskedBoundaryFactory<MeshType> masked_factory(mesh);
  for(int i(0); i < 3; ++i)
  {
    if(i != my_rank)
    {
      masked_factory.add_mask_meshpart(*mesh_node->get_halo(i));
    }
  }
  masked_factory.compile();
  MeshPartType masked_part(masked_factory);

  // create masked boundary factory
  Geometry::GlobalMaskedBoundaryFactory<MeshType> global_factory(mesh);
  for(int i(0); i < 3; ++i)
  {
    if(i != my_rank)
    {
      global_factory.add_halo(i, *mesh_node->get_halo(i));
    }
  }
  global_factory.compile(comm);
  MeshPartType global_part(global_factory);

  // write out our VTK
  Geometry::ExportVTK<MeshType> vtk(mesh);
  std::vector<double> vtx_data(mesh.get_num_entities(0), 0.0);
  add_mesh_part(vtk, "bnd", vtx_data, boundary_part);
  add_mesh_part(vtk, "msk", vtx_data, masked_part);
  add_mesh_part(vtk, "gbl", vtx_data, global_part);

  // add halos
  for(int i(0); i < 3; ++i)
  {
    if(i != my_rank)
      add_mesh_part(vtk, "halo:" + stringify(i), vtx_data, *mesh_node->get_halo(i));
    else
      vtk.add_vertex_scalar("halo:" + stringify(i), vtx_data.data());

  }
  vtk.write("dbg-global-masked-boundary", comm);

  return 0;
}
