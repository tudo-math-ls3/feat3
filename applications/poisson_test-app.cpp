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

template<typename T_>
class RhsFunc
{
public:
  static T_ eval(T_ /*x*/, T_ /*y*/)
  {
    return T_(1);
  }
};

void test_hypercube_2d(int rank, int num_patches, Index desired_refinement_level)
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
  typedef ConformalMesh<Shape::Hypercube<dim_2D> > BaseMeshType;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(m, size_set);

  BaseMeshType basemesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(m, basemesh);
  MeshControl<dim_2D>::fill_vertex_sets(m, basemesh, *((Attribute<double>*)(attrs.at(0).get())), *((Attribute<double>*)(attrs.at(1).get())));

  ///refine basemesh to match process count
  BaseMeshType* macro_basemesh = new BaseMeshType(size_set);
  MeshControl<dim_2D>::fill_adjacencies(m, *macro_basemesh);
  MeshControl<dim_2D>::fill_vertex_sets(m, *macro_basemesh, *((Attribute<double>*)(attrs.at(0).get())), *((Attribute<double>*)(attrs.at(1).get())));

  for(int i(0) ; i < log(num_patches) / log(4) ; ++i)
  {
    BaseMeshType* coarse_mesh(macro_basemesh);
    {
      Geometry::StandardRefinery<BaseMeshType> mesh_refinery(*coarse_mesh);
      macro_basemesh = new BaseMeshType(mesh_refinery);
    }
    delete coarse_mesh;
  }

  ///select macro by rank, create macro subset
  Mesh<rnt_2D, Topology<> > macro_basemesh_found(4711, &attrs);
  MeshControl<dim_2D>::fill_adjacencies(*macro_basemesh, macro_basemesh_found);
  MeshControl<dim_2D>::fill_vertex_sets(*macro_basemesh, macro_basemesh_found, *((Attribute<double>*)(attrs.at(0).get())), *((Attribute<double>*)(attrs.at(1).get())));

  Halo<1, pl_face, Mesh<rnt_2D, Topology<> > > macro_subset(macro_basemesh_found);
  macro_subset.add_element_pair(rank, rank);

  Index* polytopes_in_macrosubset(new Index[3]);
  HaloControl<dim_2D>::fill_sizes(macro_subset, polytopes_in_macrosubset);
  CellSubSet<Shape::Hypercube<2> > macro_subset_geo(polytopes_in_macrosubset);
  HaloControl<dim_2D>::fill_target_set(macro_subset, macro_subset_geo);
  delete polytopes_in_macrosubset;

  ///get a mesh from this
  PatchFactory<BaseMeshType> pf(*macro_basemesh, macro_subset_geo);
  BaseMeshType macro_mesh_geo(pf);

  ///create communication halos
  ///depending on rank: compute adjacent macros to potentially communicate with
  typedef Topology<>::storage_type_ TopologyStorageType;
  TopologyStorageType potential_comm_partners_for_face_rank(macro_basemesh_found.get_adjacent_polytopes(pl_face, pl_face, 0));
  std::vector<std::shared_ptr<HaloBase<Mesh<rnt_2D, Topology<> > > > > macro_comm_halos;
  for(Index i(0) ; i < potential_comm_partners_for_face_rank.size() ; ++i)
  {
    TopologyStorageType comm_intersect_rank_i(macro_basemesh_found.get_comm_intersection(pl_face, pl_edge, rank, i));
    if(comm_intersect_rank_i.size() == 0)
    {
      comm_intersect_rank_i = macro_basemesh_found.get_comm_intersection(pl_face, pl_vertex, rank, i);
      for(Index j(0) ; j < comm_intersect_rank_i.size() ; ++j)
      {
        Halo<0, pl_vertex, Mesh<rnt_2D, Topology<> > > halo(macro_basemesh_found, i);
        halo.add_element_pair(comm_intersect_rank_i.at(j), comm_intersect_rank_i.at(j));

        macro_comm_halos.push_back(std::shared_ptr<HaloBase<Mesh<rnt_2D, Topology<> > > >(new Halo<0, pl_vertex, Mesh<rnt_2D, Topology<> > >(halo)));
      }
    }
    else
    {
      for(Index j(0) ; j < comm_intersect_rank_i.size() ; ++j)
      {
        Halo<0, pl_edge, Mesh<rnt_2D, Topology<> > > halo(macro_basemesh_found, i);
        halo.add_element_pair(comm_intersect_rank_i.at(j), comm_intersect_rank_i.at(j));
        macro_comm_halos.push_back(std::shared_ptr<HaloBase<Mesh<rnt_2D, Topology<> > > >(new Halo<0, pl_edge, Mesh<rnt_2D, Topology<> > >(halo)));
      }
    }
  }
  std::sort(macro_comm_halos.begin(), macro_comm_halos.end(), compare_other<Mesh<rnt_2D, Topology<> >, std::vector, Index>);

  ///get macro boundary components
  Mesh<rnt_2D>::topology_type_::storage_type_ macro_edges(macro_basemesh_found.get_adjacent_polytopes(pl_face, pl_edge, rank));
  std::vector<std::shared_ptr<HaloBase<Mesh<rnt_2D> > > > macro_boundaries_found;
  for(Index i(0) ; i < macro_edges.size() ; ++i)
  {
    Halo<0, pl_edge, Mesh<rnt_2D> > result(macro_basemesh_found);
    result.add_element_pair(macro_edges.at(i), macro_edges.at(i));
    macro_boundaries_found.push_back(std::shared_ptr<HaloBase<Mesh<rnt_2D> > >(new Halo<0, pl_edge, Mesh<rnt_2D> >(result)));
  }

  ///refine everything to desired level of detail
  BaseMeshType* macro_basemesh_fine = new BaseMeshType(*macro_basemesh);
  BaseMeshType* macro_mesh_geo_fine = new BaseMeshType(macro_mesh_geo);
  CellSubSet<Shape::Hypercube<2> >* macro_subset_geo_fine = new CellSubSet<Shape::Hypercube<2> >(macro_subset_geo);

  for(Index i(0) ; i < desired_refinement_level - (log(num_patches) / log(4)) ; ++i)
  {
    BaseMeshType* coarse_macro_basemesh_fine(macro_basemesh_fine);
    BaseMeshType* coarse_macro_mesh_geo_fine(macro_mesh_geo_fine);
    CellSubSet<Shape::Hypercube<2> >* coarse_macro_subset_geo_fine(macro_subset_geo_fine);
    {
      Geometry::StandardRefinery<BaseMeshType> refinery_0(*coarse_macro_basemesh_fine);
      Geometry::StandardRefinery<BaseMeshType> refinery_1(*coarse_macro_mesh_geo_fine);
      Geometry::StandardRefinery<CellSubSet<Shape::Hypercube<2> >, BaseMeshType> refinery_2(*coarse_macro_subset_geo_fine, *macro_basemesh_fine);
      macro_subset_geo_fine = new CellSubSet<Shape::Hypercube<2> >(refinery_2);
      macro_mesh_geo_fine = new BaseMeshType(refinery_1);
      macro_basemesh_fine = new BaseMeshType(refinery_0);
    }
    //delete coarse_macro_basemesh_fine; //TODO why can u not deallocate this
    delete coarse_macro_mesh_geo_fine;
    delete coarse_macro_subset_geo_fine;
  }

  //std::cout << (*macro_mesh_geo_fine).get_index_set<1, 0>().get_num_entities() << std::endl;
  if(rank == 0)
  {
    std::cout << "proc " << rank << " macro_mesh_geo_fine with " << macro_mesh_geo_fine->get_num_entities(0) << " vertices " << std::endl;
    std::cout << "proc " << rank << " macro_mesh_geo_fine with " << macro_mesh_geo_fine->get_num_entities(1) << " edges " << std::endl;
    std::cout << "proc " << rank << " macro_mesh_geo_fine with " << macro_mesh_geo_fine->get_num_entities(2) << " faces " << std::endl;
  }
  std::vector<std::shared_ptr<CellSubSet<Shape::Hypercube<2> > > > macro_comm_halos_fine;
  for(Index i(0) ; i < macro_comm_halos.size() ; ++i)
  {
    Index* polytopes_in_subset = new Index[3];
    if(macro_comm_halos.at(i)->get_level() == pl_vertex)
      HaloControl<dim_1D>::fill_sizes(*((Halo<0, pl_vertex, Mesh<rnt_2D, Topology<> > >*)(macro_comm_halos.at(i).get())), polytopes_in_subset);
    else
      HaloControl<dim_2D>::fill_sizes(*((Halo<0, pl_edge, Mesh<rnt_2D, Topology<> > >*)(macro_comm_halos.at(i).get())), polytopes_in_subset);

    Geometry::CellSubSet<Shape::Hypercube<2> >* cell_sub_set_fine(new Geometry::CellSubSet<Shape::Hypercube<2> >(polytopes_in_subset));

    if(macro_comm_halos.at(i)->get_level() == pl_vertex)
      HaloControl<dim_1D>::fill_target_set(*((Halo<0, pl_vertex, Mesh<rnt_2D, Topology<> > >*)(macro_comm_halos.at(i).get())), *cell_sub_set_fine);
    else
      HaloControl<dim_2D>::fill_target_set(*((Halo<0, pl_edge, Mesh<rnt_2D, Topology<> > >*)(macro_comm_halos.at(i).get())), *cell_sub_set_fine);

    BaseMeshType* macro_basemesh_temp = new BaseMeshType(*macro_basemesh);
    for(Index j(0) ; j < desired_refinement_level - (log(num_patches) / log(4)) ; ++j)
    {
      BaseMeshType* coarse_macro_basemesh_temp(macro_basemesh_temp);
      CellSubSet<Shape::Hypercube<2> >* coarse_cell_sub_set_fine(cell_sub_set_fine);
      //refine
      {
        Geometry::StandardRefinery<BaseMeshType> refinery_0(*coarse_macro_basemesh_temp);
        Geometry::StandardRefinery<Geometry::CellSubSet<Shape::Hypercube<2> >, BaseMeshType> cell_refinery(*coarse_cell_sub_set_fine, *macro_basemesh_temp);
        cell_sub_set_fine = new CellSubSet<Shape::Hypercube<2> >(cell_refinery);
        macro_basemesh_temp = new BaseMeshType(refinery_0);
      }
      delete coarse_macro_basemesh_temp;
      //delete coarse_cell_sub_set_fine;
    }
    delete macro_basemesh_temp;
    //add
    macro_comm_halos_fine.push_back(std::shared_ptr<Geometry::CellSubSet<Shape::Hypercube<2> > >(new Geometry::CellSubSet<Shape::Hypercube<2> >(*cell_sub_set_fine)));
    delete[] polytopes_in_subset;
  }

  std::vector<std::shared_ptr<CellSubSet<Shape::Hypercube<2> > > > macro_boundaries_fine;
  for(Index i(0) ; i < macro_boundaries_found.size() ; ++i)
  {
    Index* polytopes_in_subset = new Index[3];
    if(macro_boundaries_found.at(i)->get_level() == pl_vertex)
      HaloControl<dim_1D>::fill_sizes(*((Halo<0, pl_vertex, Mesh<rnt_2D, Topology<> > >*)(macro_boundaries_found.at(i).get())), polytopes_in_subset);
    else
      HaloControl<dim_2D>::fill_sizes(*((Halo<0, pl_edge, Mesh<rnt_2D, Topology<> > >*)(macro_boundaries_found.at(i).get())), polytopes_in_subset);

    Geometry::CellSubSet<Shape::Hypercube<2> >* cell_sub_set_fine(new Geometry::CellSubSet<Shape::Hypercube<2> >(polytopes_in_subset));

    if(macro_boundaries_found.at(i)->get_level() == pl_vertex)
      HaloControl<dim_1D>::fill_target_set(*((Halo<0, pl_vertex, Mesh<rnt_2D, Topology<> > >*)(macro_boundaries_found.at(i).get())), *cell_sub_set_fine);
    else
      HaloControl<dim_2D>::fill_target_set(*((Halo<0, pl_edge, Mesh<rnt_2D, Topology<> > >*)(macro_boundaries_found.at(i).get())), *cell_sub_set_fine);

    BaseMeshType* macro_basemesh_temp = new BaseMeshType(*macro_basemesh);
    for(Index j(0) ; j < desired_refinement_level - (log(num_patches) / log(4)) ; ++j)
    {
      BaseMeshType* coarse_macro_basemesh_temp(macro_basemesh_temp);
      CellSubSet<Shape::Hypercube<2> >* coarse_cell_sub_set_fine(cell_sub_set_fine);
      //refine
      {
        Geometry::StandardRefinery<BaseMeshType> refinery_0(*coarse_macro_basemesh_temp);
        Geometry::StandardRefinery<Geometry::CellSubSet<Shape::Hypercube<2> >, BaseMeshType> cell_refinery(*coarse_cell_sub_set_fine, *macro_basemesh_temp);
        cell_sub_set_fine = new CellSubSet<Shape::Hypercube<2> >(cell_refinery);
        macro_basemesh_temp = new BaseMeshType(refinery_0);
      }
      delete coarse_macro_basemesh_temp;
      //delete coarse_cell_sub_set_fine;
    }
    delete macro_basemesh_temp;
    //add
    macro_boundaries_fine.push_back(std::shared_ptr<Geometry::CellSubSet<Shape::Hypercube<2> > >(new Geometry::CellSubSet<Shape::Hypercube<2> >(*cell_sub_set_fine)));
    delete[] polytopes_in_subset;
  }

  for(Index i(0) ; i < macro_boundaries_fine.size() ; ++i)
  {
    if(rank == 0)
    {
      std::cout << "proc " << rank << " subset " << i << " with " << macro_boundaries_fine.at(i)->get_num_entities(0) << " vertices." << std::endl;
      std::cout << "proc " << rank << " subset " << i << " with " << macro_boundaries_fine.at(i)->get_num_entities(1) << " edges." << std::endl;
      std::cout << "proc " << rank << " subset " << i << " with " << macro_boundaries_fine.at(i)->get_num_entities(2) << " faces." << std::endl;
      for(Index j(0) ; j < macro_boundaries_fine.at(i)->get_num_entities(0) ; ++j)
        std::cout << "proc " << rank << " subset " << i << " vertex: " << macro_boundaries_fine.at(i)->get_target_set<0>()[j] << std::endl;
      for(Index j(0) ; j < macro_boundaries_fine.at(i)->get_num_entities(1) ; ++j)
        std::cout << "proc " << rank << " subset " << i << " edge: " << macro_boundaries_fine.at(i)->get_target_set<1>()[j] << std::endl;
      for(Index j(0) ; j < macro_boundaries_fine.at(i)->get_num_entities(2) ; ++j)
        std::cout << "proc " << rank << " subset " << i << " face: " << macro_boundaries_fine.at(i)->get_target_set<2>()[j] << std::endl;

    }
  }

  ///assembly
  // create trafo
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > trafo(*macro_mesh_geo_fine);
  // create space
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);

  SparseMatrixCSR<Mem::Main, double> mat_sys(Space::DofAdjacency<>::assemble(space));
  mat_sys.clear();
  Assembly::BilinearScalarLaplaceFunctor::assemble(mat_sys, space, "gauss-legendre:2");

  DenseVector<Mem::Main, double> vec_rhs(space.get_num_dofs(), double(0));
  Assembly::LinearScalarIntegralFunctor<RhsFunc>::assemble(vec_rhs, space, "gauss-legendre:2");

  // assemble homogeneous Dirichlet BCs
  Assembly::DirichletBC<Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > > dirichlet(space);
  for(Index i(0) ; i < macro_boundaries_fine.size() ; ++i)
  {
    std::cout << "Adding cell set on process " << rank << std::endl;
    dirichlet.add_cell_set(*macro_boundaries_fine.at(i).get());
  }
  // allocate solution vector
  DenseVector<Mem::Main, double> vec_sol(space.get_num_dofs(), double(0));

  // assemble filter:
  UnitFilter<Mem::Main, double> filter(dirichlet.assemble<Mem::Main, double>());

  MPI_Barrier(MPI_COMM_WORLD);

  // filter system TODO: check why this segfaults
  filter.filter_mat(mat_sys);
  filter.filter_rhs(vec_rhs);
  filter.filter_sol(vec_sol);


  delete macro_basemesh;
}


int main(int argc, char* argv[])
{
  int me(0);
  int num_patches(0);
  Index desired_refinement_level(4);

#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &num_patches);
#endif

  test_hypercube_2d(me, num_patches, desired_refinement_level);

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
