#ifndef SERIAL
#include <mpi.h>
#endif

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/mesh_control.hpp>
#include<kernel/foundation/dense_data_wrapper.hpp>
#include<kernel/geometry/macro_factory.hpp>
#include<kernel/geometry/patch_factory.hpp>
#include<kernel/lafem/dense_vector.hpp>
#include<kernel/lafem/vector_mirror.hpp>
#include<kernel/foundation/halo_control.hpp>
#include<kernel/foundation/halo.hpp>
#include<kernel/geometry/cell_sub_set.hpp>
#include<kernel/archs.hpp>
#include<deque>
#include<algorithm>
#include<cmath>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/dof_adjacency.hpp>
#include <kernel/space/dof_mirror.hpp>
#include <kernel/assembly/standard_operators.hpp>
#include <kernel/assembly/standard_functionals.hpp>
#include <kernel/assembly/dirichlet_bc.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;
using namespace FEAST::Geometry;

void test_hypercube_2d(int rank, int num_patches)
{
  /*(0,1) (1,1)
   *  *----*
   *  |    |
   *  |    |
   *  *----*
   *(0,0) (1,0)
   */

  //create attributes for vertex coords
  std::vector<std::shared_ptr<Foundation::AttributeBase<> > > attrs;
  attrs.push_back(std::shared_ptr<Foundation::AttributeBase<> >(new Foundation::Attribute<double>)); //vertex x-coords
  attrs.push_back(std::shared_ptr<Foundation::AttributeBase<> >(new Foundation::Attribute<double>)); //vertex y-coords

  ((Foundation::Attribute<double>*)(attrs.at(0).get()))->get_data().push_back(double(0));
  ((Foundation::Attribute<double>*)(attrs.at(1).get()))->get_data().push_back(double(0));

  ((Foundation::Attribute<double>*)(attrs.at(0).get()))->get_data().push_back(double(1));
  ((Foundation::Attribute<double>*)(attrs.at(1).get()))->get_data().push_back(double(0));

  ((Foundation::Attribute<double>*)(attrs.at(0).get()))->get_data().push_back(double(0));
  ((Foundation::Attribute<double>*)(attrs.at(1).get()))->get_data().push_back(double(1));

  ((Foundation::Attribute<double>*)(attrs.at(0).get()))->get_data().push_back(double(1));
  ((Foundation::Attribute<double>*)(attrs.at(1).get()))->get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<> > m(0, &attrs);
  Foundation::MeshAttributeRegistration::execute(m, Foundation::pl_vertex);
  Foundation::MeshAttributeRegistration::execute(m, Foundation::pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);

  ///get a conformal mesh as basemesh
  typedef ConformalMesh<Shape::Hypercube<dim_2D> > BaseMeshType_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(m, size_set);

  BaseMeshType_ basemesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(m, basemesh);
  MeshControl<dim_2D>::fill_vertex_sets(m, basemesh, *((Attribute<double>*)(attrs.at(0).get())), *((Attribute<double>*)(attrs.at(1).get())));

  ///refine basemesh to match process count
  BaseMeshType_* macro_basemesh = new BaseMeshType_(size_set);
  MeshControl<dim_2D>::fill_adjacencies(m, *macro_basemesh);
  MeshControl<dim_2D>::fill_vertex_sets(m, *macro_basemesh, *((Attribute<double>*)(attrs.at(0).get())), *((Attribute<double>*)(attrs.at(1).get())));

  for(Index i(0) ; i < log(num_patches) / log(2) ; ++i)
  {
    Geometry::StandardRefinery<BaseMeshType_> mesh_refinery(*macro_basemesh);
    macro_basemesh = new BaseMeshType_(mesh_refinery);
  }

  delete macro_basemesh;

}

int main(int argc, char* argv[])
{
  int me(0);
  int num_patches(0);

#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &num_patches);
#endif

  test_hypercube_2d(me, num_patches);

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
