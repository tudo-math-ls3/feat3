/**
 * \file
 * \brief FEAST milestone 1 test application for the solution of the Poisson problem
 * \author Markus Geveler
 * \date 2012 - 2013
 *
 * This file contains the first application code testing the FEAST kernel w.o. the final application framework.
 */

#include <kernel/base_header.hpp>

#include <kernel/foundation/comm_base.hpp>

#include <test_system/test_system.hpp>

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>
#include <kernel/geometry/patch_factory.hpp>
#include <kernel/geometry/macro_factory.hpp>
#include <kernel/geometry/patch_halo_factory.hpp>
#include <kernel/foundation/attribute.hpp>
#include <kernel/foundation/mesh_control.hpp>
#include <kernel/foundation/dense_data_wrapper.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/matrix_mirror.hpp>
#include <kernel/foundation/halo_control.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/archs.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/dirichlet_assembler.hpp>

#include <kernel/scarc/solver_data.hpp>
#include <kernel/scarc/solver_pattern.hpp>

#include<deque>
#include<algorithm>
#include<cmath>

#ifdef _WIN32
  extern "C" void __cdecl Sleep(unsigned int);
#  define sleep(x) Sleep(1000*(x))
#else
#  include <unistd.h>
#endif

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;
using namespace FEAST::Geometry;
using namespace FEAST::ScaRC;

int pow4(int i)
{
  return (1 << i) * (1 << i); // = 4^i
}

void test_hypercube_2d(Index rank, Index num_patches, Index desired_refinement_level)
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
  Foundation::Mesh<Dim2D, Foundation::Topology<> > m(0);
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
  typedef ConformalMesh<Shape::Hypercube<2> > BaseMeshType;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(m, size_set);

  BaseMeshType basemesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(m, basemesh);
  MeshControl<dim_2D>::fill_vertex_sets(m, basemesh, *((Attribute<double>*)(attrs.at(0).get())), *((Attribute<double>*)(attrs.at(1).get())));

  ///refine basemesh to match process count
  BaseMeshType* macro_basemesh = new BaseMeshType(size_set);
  MeshControl<dim_2D>::fill_adjacencies(m, *macro_basemesh);
  MeshControl<dim_2D>::fill_vertex_sets(m, *macro_basemesh, *((Attribute<double>*)(attrs.at(0).get())), *((Attribute<double>*)(attrs.at(1).get())));

  for(int i(0) ; (Index)pow4((int)i) < num_patches ; ++i)
  {
    BaseMeshType* coarse_mesh(macro_basemesh);
    {
      Geometry::StandardRefinery<BaseMeshType> mesh_refinery(*coarse_mesh);
      macro_basemesh = new BaseMeshType(mesh_refinery);
    }
    delete coarse_mesh;
  }

  ///select macro by rank, create macro subset
  Mesh<Dim2D, Topology<> > macro_basemesh_found(4711);
  MeshControl<dim_2D>::fill_adjacencies(*macro_basemesh, macro_basemesh_found);
  MeshControl<dim_2D>::fill_vertex_sets(*macro_basemesh, macro_basemesh_found, *((Attribute<double>*)(attrs.at(0).get())), *((Attribute<double>*)(attrs.at(1).get())));

  Halo<1, PLFace, Mesh<Dim2D, Topology<> > > macro_subset(macro_basemesh_found);
  macro_subset.push_back(rank);

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
  TopologyStorageType potential_comm_partners_for_face_rank(macro_basemesh_found.get_adjacent_polytopes(pl_face, pl_face, rank));
  for(Index i(0) ; i < potential_comm_partners_for_face_rank.size() ; ++i)
  {
    if(potential_comm_partners_for_face_rank.at(i) == (Index)rank)
    {
      potential_comm_partners_for_face_rank.erase(potential_comm_partners_for_face_rank.begin() + (const TopologyStorageType::difference_type)i);
      break;
    }
  }

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<> > > > > macro_comm_halos;
  for(Index i(0) ; i < potential_comm_partners_for_face_rank.size() ; ++i)
  {
    TopologyStorageType comm_intersect_rank_i(macro_basemesh_found.get_comm_intersection(pl_face, pl_edge, rank, potential_comm_partners_for_face_rank.at(i)));
    if(comm_intersect_rank_i.size() == 0)
    {
      comm_intersect_rank_i = macro_basemesh_found.get_comm_intersection(pl_face, pl_vertex, rank, potential_comm_partners_for_face_rank.at(i));
      for(Index j(0) ; j < comm_intersect_rank_i.size() ; ++j)
      {
        Halo<0, PLVertex, Mesh<Dim2D, Topology<> > > halo(macro_basemesh_found, potential_comm_partners_for_face_rank.at(i));
        halo.push_back(comm_intersect_rank_i.at(j));

        macro_comm_halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<> > > >(new Halo<0, PLVertex, Mesh<Dim2D, Topology<> > >(halo)));
      }
    }
    else
    {
      for(Index j(0) ; j < comm_intersect_rank_i.size() ; ++j)
      {
        Halo<0, PLEdge, Mesh<Dim2D, Topology<> > > halo(macro_basemesh_found, potential_comm_partners_for_face_rank.at(i));
        halo.push_back(comm_intersect_rank_i.at(j));
        macro_comm_halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<> > > >(new Halo<0, PLEdge, Mesh<Dim2D, Topology<> > >(halo)));
      }
    }
  }
  std::sort(macro_comm_halos.begin(), macro_comm_halos.end(), compare_other<Mesh<Dim2D, Topology<> >, std::vector>);

  ///get adjacencies of halos
  std::vector<std::vector<Index> > halo_adjacencies;
  for(Index i(0) ; i < macro_comm_halos.size() ; ++i)
  {
    Index elem(macro_comm_halos.at(i)->get_element(0));
    halo_adjacencies.push_back(std::vector<Index>());

    if(macro_comm_halos.at(i)->get_level() == pl_edge)
    {
      TopologyStorageType adjacent_edges_for_halo_elem(macro_basemesh_found.get_adjacent_polytopes(pl_edge, pl_edge, elem));
      //remove self
      for(Index j(0) ; j < adjacent_edges_for_halo_elem.size() ; ++j)
      {
        if(adjacent_edges_for_halo_elem.at(j) == elem)
        {
          adjacent_edges_for_halo_elem.erase(adjacent_edges_for_halo_elem.begin() + (const TopologyStorageType::difference_type)j);
          break;
        }
      }

      for(Index j(0) ; j < adjacent_edges_for_halo_elem.size() ; ++j)
      {
        for(Index k(0) ; k < macro_comm_halos.size() ; ++k)
        {
          if(macro_comm_halos.at(k)->get_element(0) == adjacent_edges_for_halo_elem.at(j))
          {
            if(macro_comm_halos.at(k)->get_level() == pl_edge)
              halo_adjacencies.at(i).push_back(k);
          }
        }
      }

      TopologyStorageType adjacent_vertices_for_halo_elem(macro_basemesh_found.get_adjacent_polytopes(pl_edge, pl_vertex, elem));
      for(Index j(0) ; j < adjacent_vertices_for_halo_elem.size() ; ++j)
      {
        for(Index k(0) ; k < macro_comm_halos.size() ; ++k)
        {
          if(macro_comm_halos.at(k)->get_element(0) == adjacent_vertices_for_halo_elem.at(j))
          {
            if(macro_comm_halos.at(k)->get_level() == pl_vertex)
              halo_adjacencies.at(i).push_back(k);
          }
        }
      }
    }
    else if(macro_comm_halos.at(i)->get_level() == pl_vertex)
    {
      TopologyStorageType adjacent_edges_for_halo_elem(macro_basemesh_found.get_adjacent_polytopes(pl_vertex, pl_edge, elem));
      for(Index j(0) ; j < adjacent_edges_for_halo_elem.size() ; ++j)
      {
        for(Index k(0) ; k < macro_comm_halos.size() ; ++k)
        {
          if(macro_comm_halos.at(k)->get_element(0) == adjacent_edges_for_halo_elem.at(j))
          {
            if(macro_comm_halos.at(k)->get_level() == pl_edge)
              halo_adjacencies.at(i).push_back(k);
          }
        }
      }
    }
  }
  for(Index i(0) ; i < halo_adjacencies.size() ; ++i)
    for(Index j(0) ; j < halo_adjacencies.at(i).size() ; ++j)
      std::cout << "proc: " << rank << " halo " << i << " adj to " << halo_adjacencies.at(i).at(j) << std::endl;


  ///get (non-inner) macro boundary components
  Mesh<Dim2D>::topology_type_::storage_type_ macro_edges(macro_basemesh_found.get_adjacent_polytopes(pl_face, pl_edge, rank));
  Mesh<Dim2D>::topology_type_::storage_type_ macro_edges_final;
  std::cout << "proc " << rank << " #macro edges initially " << macro_edges.size() << std::endl;
  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D> > > > macro_boundaries_found;
  for(Index i(0) ; i < macro_edges.size() ; ++i)
  {
    bool in(false);
    for(Index j(0) ; j < macro_comm_halos.size() ; ++j)
    {
      if(macro_comm_halos.at(j)->get_level() == 1)
        if(macro_edges.at(i) == macro_comm_halos.at(j)->get_element(0))
          in = true;
    }
    if(!in)
      macro_edges_final.push_back(macro_edges.at(i));
  }
  std::cout << "proc " << rank << " #macro edges after elimination " << macro_edges_final.size() << std::endl;
  for(Index i(0) ; i < macro_edges_final.size() ; ++i)
  {
    Halo<0, PLEdge, Mesh<Dim2D> > result(macro_basemesh_found);
    result.push_back(macro_edges_final.at(i));
    macro_boundaries_found.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D> > >(new Halo<0, PLEdge, Mesh<Dim2D> >(result)));
  }

  ///refine everything to desired level of detail
  BaseMeshType* macro_basemesh_fine = new BaseMeshType(*macro_basemesh);
  BaseMeshType* macro_mesh_geo_fine = new BaseMeshType(macro_mesh_geo);
  CellSubSet<Shape::Hypercube<2> >* macro_subset_geo_fine = new CellSubSet<Shape::Hypercube<2> >(macro_subset_geo);

  for(Index i(0) ; num_patches < (Index)pow4(int(desired_refinement_level - i)) ; ++i)
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
    delete coarse_macro_basemesh_fine;
    delete coarse_macro_mesh_geo_fine;
    delete coarse_macro_subset_geo_fine;
  }

  /*if(rank == 0)
  {
    std::cout << "proc " << rank << " macro_mesh_geo_fine with " << macro_mesh_geo_fine->get_num_entities(0) << " vertices " << std::endl;
    std::cout << "proc " << rank << " macro_mesh_geo_fine with " << macro_mesh_geo_fine->get_num_entities(1) << " edges " << std::endl;
    std::cout << "proc " << rank << " macro_mesh_geo_fine with " << macro_mesh_geo_fine->get_num_entities(2) << " faces " << std::endl;
  }*/

  //refine halos
  std::vector<std::shared_ptr<CellSubSet<Shape::Hypercube<2> > > > macro_comm_halos_fine;
  for(Index i(0) ; i < macro_comm_halos.size() ; ++i)
  {
    Index* polytopes_in_subset = new Index[3];
    if(macro_comm_halos.at(i)->get_level() == pl_vertex)
      HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D, Topology<> > >*)(macro_comm_halos.at(i).get())), polytopes_in_subset);
    else
      HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D, Topology<> > >*)(macro_comm_halos.at(i).get())), polytopes_in_subset);

    Geometry::CellSubSet<Shape::Hypercube<2> >* cell_sub_set_fine(new Geometry::CellSubSet<Shape::Hypercube<2> >(polytopes_in_subset));

    if(macro_comm_halos.at(i)->get_level() == pl_vertex)
      HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D, Topology<> > >*)(macro_comm_halos.at(i).get())), *cell_sub_set_fine);
    else
      HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D, Topology<> > >*)(macro_comm_halos.at(i).get())), *cell_sub_set_fine);

    BaseMeshType* macro_mesh_temp = new BaseMeshType(macro_mesh_geo);

    ///from global to local coordinates
    PatchHaloFactory<CellSubSet<Shape::Hypercube<2> > > map_factory(macro_subset_geo, *cell_sub_set_fine);
    CellSubSet<Shape::Hypercube<2> >* temp(cell_sub_set_fine);
    cell_sub_set_fine = new CellSubSet<Shape::Hypercube<2> >(map_factory);
    delete temp;

    if(!(macro_comm_halos.at(i)->get_level() == pl_vertex))
      for(Index j(0) ; num_patches < (Index)pow4(int(desired_refinement_level - j)) ; ++j)
      {
        BaseMeshType* coarse_macro_mesh_temp(macro_mesh_temp);
        CellSubSet<Shape::Hypercube<2> >* coarse_cell_sub_set_fine(cell_sub_set_fine);
        //refine
        {
          Geometry::StandardRefinery<BaseMeshType> refinery_0(*coarse_macro_mesh_temp);
          Geometry::StandardRefinery<Geometry::CellSubSet<Shape::Hypercube<2> >, BaseMeshType> cell_refinery(*coarse_cell_sub_set_fine, *macro_mesh_temp);
          cell_sub_set_fine = new CellSubSet<Shape::Hypercube<2> >(cell_refinery);
          macro_mesh_temp = new BaseMeshType(refinery_0);
        }
        delete coarse_macro_mesh_temp;
        delete coarse_cell_sub_set_fine;
      }
    delete macro_mesh_temp;
    //add
    macro_comm_halos_fine.push_back(std::shared_ptr<Geometry::CellSubSet<Shape::Hypercube<2> > >(new Geometry::CellSubSet<Shape::Hypercube<2> >(*cell_sub_set_fine)));
    delete[] polytopes_in_subset;
  }

  std::vector<std::shared_ptr<CellSubSet<Shape::Hypercube<2> > > > macro_boundaries_fine;
  for(Index i(0) ; i < macro_boundaries_found.size() ; ++i)
  {
    Index* polytopes_in_subset = new Index[3];

    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D, Topology<> > >*)(macro_boundaries_found.at(i).get())), polytopes_in_subset);
    Geometry::CellSubSet<Shape::Hypercube<2> >* cell_sub_set_fine(new Geometry::CellSubSet<Shape::Hypercube<2> >(polytopes_in_subset));
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D, Topology<> > >*)(macro_boundaries_found.at(i).get())), *cell_sub_set_fine);

    BaseMeshType* macro_mesh_temp = new BaseMeshType(macro_mesh_geo);

    ///from global to local coordinates
    PatchHaloFactory<CellSubSet<Shape::Hypercube<2> > > map_factory(macro_subset_geo, *cell_sub_set_fine);
    CellSubSet<Shape::Hypercube<2> >* temp(cell_sub_set_fine);
    cell_sub_set_fine = new CellSubSet<Shape::Hypercube<2> >(map_factory);
    delete temp;

    for(Index j(0) ; num_patches < (Index)pow4(int(desired_refinement_level - j)) ; ++j)
    {
      BaseMeshType* coarse_macro_mesh_temp(macro_mesh_temp);
      CellSubSet<Shape::Hypercube<2> >* coarse_cell_sub_set_fine(cell_sub_set_fine);
      //refine
      {
        Geometry::StandardRefinery<BaseMeshType> refinery_0(*coarse_macro_mesh_temp);
        Geometry::StandardRefinery<Geometry::CellSubSet<Shape::Hypercube<2> >, BaseMeshType> cell_refinery(*coarse_cell_sub_set_fine, *macro_mesh_temp);
        cell_sub_set_fine = new CellSubSet<Shape::Hypercube<2> >(cell_refinery);
        macro_mesh_temp = new BaseMeshType(refinery_0);
      }
      delete coarse_macro_mesh_temp;
      delete coarse_cell_sub_set_fine;
    }

    delete macro_mesh_temp;
    //add
    macro_boundaries_fine.push_back(std::shared_ptr<Geometry::CellSubSet<Shape::Hypercube<2> > >(new Geometry::CellSubSet<Shape::Hypercube<2> >(*cell_sub_set_fine)));
    delete[] polytopes_in_subset;
  }

  /*for(Index i(0) ; i < macro_boundaries_fine.size() ; ++i)
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
  }*/

  ///assembly
  // create trafo
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > trafo_base(*macro_basemesh_fine); //TODO do we need it?
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > trafo(*macro_mesh_geo_fine);
  // create space
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space_base(trafo_base); // TODO do we need it?
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);

  //assemble dof adjacencies

  SparseMatrixCSR<Mem::Main, double> mat_sys;
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  DenseVector<Mem::Main, double> vec_rhs(space.get_num_dofs(), double(0));
  Assembly::Common::ConstantFunction rhs_func(1.0);
  Assembly::Common::ForceFunctional<Assembly::Common::ConstantFunction> rhs_functional(rhs_func);
  Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, rhs_functional, space, cubature_factory);

  // assemble homogeneous Dirichlet BCs
  Assembly::DirichletAssembler<Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > > dirichlet(space);
  std::cout << "proc " << rank << " #macro boundaries " << macro_boundaries_fine.size() << std::endl;
  for(Index i(0) ; i < macro_boundaries_fine.size() ; ++i)
  {
    //std::cout << "Adding cell set on process " << rank << std::endl;
    dirichlet.add_cell_set(*macro_boundaries_fine.at(i).get());
  }
  // allocate solution vector
  DenseVector<Mem::Main, double> vec_sol(space.get_num_dofs(), double(0));

  ///assemble mirrors
  std::vector<LAFEM::VectorMirror<Mem::Main, double> > mirrors;
  std::vector<LAFEM::DenseVector<Mem::Main, double> > sendbufs;
  std::vector<LAFEM::DenseVector<Mem::Main, double> > recvbufs;
  std::vector<Index> destranks;
  std::vector<Index> sourceranks;

  std::cout << "proc " << rank << " #comm halos " << macro_comm_halos_fine.size() << std::endl;
  for(Index i(0) ; i < macro_comm_halos_fine.size() ; ++i)
  {
    VectorMirror<Mem::Main, double> target_mirror;
    Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, *(macro_comm_halos_fine.at(i).get()));
    DenseVector<Mem::Main, double> sendbuf(target_mirror.size());
    DenseVector<Mem::Main, double> recvbuf(target_mirror.size());

    mirrors.push_back(std::move(target_mirror));
    sendbufs.push_back(sendbuf);
    recvbufs.push_back(recvbuf);

    destranks.push_back(macro_comm_halos.at(i)->get_other());
    sourceranks.push_back(macro_comm_halos.at(i)->get_other());
  }

  ///build up a type-1 matrix for the local solvers
  std::vector<DenseVector<Mem::Main, double> > val_sendbufs;
  std::vector<DenseVector<Mem::Main, double> > val_recvbufs;
  std::vector<DenseVector<Mem::Main, Index> > colind_sendbufs;
  std::vector<DenseVector<Mem::Main, Index> > colind_recvbufs;
  std::vector<DenseVector<Mem::Main, Index> > rp_sendbufs;
  std::vector<DenseVector<Mem::Main, Index> > rp_recvbufs;

  SparseMatrixCSR<Mem::Main, double> mat_localsys(mat_sys.clone());
  for(Index i(0) ; i < macro_comm_halos.size() ; ++i)
  {
    //gather data, exchange
    MatrixMirror<Mem::Main, double> mat_mirror(mirrors.at(i), mirrors.at(i));
    SparseMatrixCSR<Mem::Main, double> buf_mat;
    Assembly::MirrorAssembler::assemble_buffer_matrix(buf_mat, mat_mirror, mat_localsys);
    mat_mirror.gather_op(buf_mat, mat_localsys);

    double* val(buf_mat.val());
    Index* row_ptr(buf_mat.row_ptr());
    Index* col_ind(buf_mat.col_ind());

    DenseVector<Mem::Main, double> val_sendbuf(buf_mat.used_elements());
    DenseVector<Mem::Main, double> val_recvbuf(buf_mat.used_elements());
    DenseVector<Mem::Main, Index> colind_sendbuf(buf_mat.used_elements());
    DenseVector<Mem::Main, Index> colind_recvbuf(buf_mat.used_elements());
    DenseVector<Mem::Main, Index> rp_sendbuf(buf_mat.rows() + 1);
    DenseVector<Mem::Main, Index> rp_recvbuf(buf_mat.rows() + 1);

    for(Index j(0) ; j < buf_mat.used_elements() ; ++j)
    {
      val_sendbuf(j, val[j]);
      colind_sendbuf(j, col_ind[j]);
    }
    for(Index j(0) ; j < buf_mat.rows() + 1 ; ++j)
    {
      rp_sendbuf(j, row_ptr[j]);
    }

#ifndef SERIAL
    Comm::send_recv(val_sendbuf.elements(),
                              val_sendbuf.size(),
                              macro_comm_halos.at(i)->get_other(),
                              val_recvbuf.elements(),
                              val_recvbuf.size(),
                              macro_comm_halos.at(i)->get_other());
    Comm::send_recv(colind_sendbuf.elements(),
                              colind_sendbuf.size(),
                              macro_comm_halos.at(i)->get_other(),
                              colind_recvbuf.elements(),
                              colind_recvbuf.size(),
                              macro_comm_halos.at(i)->get_other());
    Comm::send_recv(rp_sendbuf.elements(),
                              rp_sendbuf.size(),
                              macro_comm_halos.at(i)->get_other(),
                              rp_recvbuf.elements(),
                              rp_recvbuf.size(),
                              macro_comm_halos.at(i)->get_other());
#endif

    val_sendbufs.push_back(val_sendbuf);
    val_recvbufs.push_back(val_recvbuf);
    colind_sendbufs.push_back(colind_sendbuf);
    colind_recvbufs.push_back(colind_recvbuf);
    rp_sendbufs.push_back(rp_sendbuf);
    rp_recvbufs.push_back(rp_recvbuf);
  }
#ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  for(Index i(0) ; i < macro_comm_halos.size() ; ++i)
  {
    MatrixMirror<Mem::Main, double> mat_mirror(mirrors.at(i), mirrors.at(i));
    SparseMatrixCSR<Mem::Main, double> buf_mat;
    Assembly::MirrorAssembler::assemble_buffer_matrix(buf_mat, mat_mirror, mat_localsys);

    mat_mirror.gather_op(buf_mat, mat_localsys);

    SparseMatrixCSR<Mem::Main, double> other_buf_mat(buf_mat.rows(),
        buf_mat.columns(),
        colind_recvbufs.at(i),
        val_recvbufs.at(i),
        rp_recvbufs.at(i));

    std::cout << "proc " << rank << "A0_omega_i " << buf_mat;
    std::cout << "proc " << rank << "A0_omega_j " << other_buf_mat;
    buf_mat.axpy<Algo::Generic>(buf_mat, other_buf_mat);
    std::cout << "proc " << rank << "A1_omega_i " << buf_mat;
    mat_mirror.scatter_op(mat_localsys, buf_mat);
  }

  ///bring up a local preconditioning matrix TODO use a product wrapper and DV
  SparseMatrixCOO<Mem::Main, double> mat_precon_temp(mat_localsys.rows(), mat_localsys.columns());
  for(Index i(0) ; i < mat_localsys.rows() ; ++i)
    mat_precon_temp(i, i, double(0.3) * (double(1)/mat_localsys(i, i)));

  SparseMatrixCSR<Mem::Main, double> mat_precon(mat_precon_temp);

  ///(type-1 to type-0 conversion) -- NO: mats are type 0 at the beginning
  /*std::cout << "proc " << rank << " #mirrors " << mirrors.size() << std::endl;
  for(Index i(0) ; i < mirrors.size() ; ++i)
  {
    if(macro_comm_halos.at(i)->get_level() != pl_vertex)
    {
      MatrixMirror<Mem::Main, double> mat_mirror(mirrors.at(i), mirrors.at(i));
      SparseMatrixCSR<Mem::Main, double> buf_mat(mat_mirror.create_buffer(mat_sys));
      mat_mirror.gather_op(buf_mat, mat_sys);
      DenseVector<Mem::Main, double> t(buf_mat.used_elements(), buf_mat.val());
      t.template scale<Algo::Generic>(t, 0.5);
      mat_mirror.scatter_op(mat_sys, buf_mat);
    }
  }*/

  // assemble filter:
  UnitFilter<Mem::Main, double> filter(dirichlet.assemble<Mem::Main, double>());

  ///filter system
  filter.filter_mat<Algo::Generic>(mat_sys);
  filter.filter_mat<Algo::Generic>(mat_localsys);
  filter.filter_rhs<Algo::Generic>(vec_rhs);
  filter.filter_sol<Algo::Generic>(vec_sol);
  //filter.filter_mat(mat_precon); //NO! we do this in the solver program when applying the correction filter after preconditioning

  std::cout << "proc " << rank << " A0 " << mat_sys << std::endl;
  std::cout << "proc " << rank << " A1 " << mat_localsys;
  std::cout << "proc " << rank << " P " << mat_precon;

  ///bring up solver data
  SynchronisedPreconditionedFilteredSolverData<double,
    Mem::Main,
    DenseVector,
    VectorMirror,
    SparseMatrixCSR,
    SparseMatrixCSR,
    UnitFilter> data(mat_sys, mat_precon, vec_sol, vec_rhs, filter,
                     std::max(SolverPatternGeneration<ScaRCBlockSmoother, Algo::Generic>::min_num_temp_vectors(), SolverPatternGeneration<RichardsonLayer, Algo::Generic>::min_num_temp_vectors()),
                     std::max(SolverPatternGeneration<ScaRCBlockSmoother, Algo::Generic>::min_num_temp_scalars(), SolverPatternGeneration<RichardsonLayer, Algo::Generic>::min_num_temp_scalars()),
                     std::max(SolverPatternGeneration<ScaRCBlockSmoother, Algo::Generic>::min_num_temp_indices(), SolverPatternGeneration<RichardsonLayer, Algo::Generic>::min_num_temp_indices()));

  data.vector_mirrors() = std::move(mirrors);
  data.vector_mirror_sendbufs() = sendbufs;
  data.vector_mirror_recvbufs() = recvbufs;
  data.dest_ranks() = destranks;
  data.source_ranks() = sourceranks;
  data.localsys() = mat_localsys;

  std::shared_ptr<SolverFunctorBase<DenseVector<Mem::Main, double> > > solver(SolverPatternGeneration<ScaRCBlockSmoother, Algo::Generic>::execute(data, 1000, 1e-8));

  DenseVector<Mem::Main, double> dummy;
  std::shared_ptr<SolverFunctorBase<DenseVector<Mem::Main, double> > > block_solver(SolverPatternGeneration<RichardsonLayer, Algo::Generic>::execute(data, dummy, 20, 1e-1));

  solver->set_preconditioner(block_solver);

#ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  std::cout << "proc " << rank << " rhs " << data.rhs();

  solver->execute();

#ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  std::cout << "proc " << rank << " sol " << data.sol();
  std::cout << "proc " << rank << " iters used " << data.used_iters() << std::endl;

#ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  sleep(1);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (rank==0)
    std::cout<<SolverFunctorBase<DenseVector<Mem::Main, double> >::pretty_printer(solver->type_name());

  ///gather solution vector on process 0
  //local mirror creation
  VectorMirror<Mem::Main, double> source_mirror_local;
  Assembly::MirrorAssembler::assemble_mirror(source_mirror_local, space_base, *macro_subset_geo_fine);
  ///now we must send the mirror to rank 0 -> VM must implement foundation's sendable/bufferable interfaces -> TODO
  ///alternatively, we can create all needed subsets again on rank 0 and only communicate sol
  ///TODO send data.sol to rank 0
  if(rank == 0)
  {
    ///TODO receive sol and mirror, scatter into global solution
  }

  delete macro_basemesh;
}

int main(int argc, char* argv[])
{
  int me(0);
  int num_patches(0);
  Index desired_refinement_level(4);

  std::cout<<"CTEST_FULL_OUTPUT"<<std::endl;

#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &num_patches);
#else
  (void)argc;
  (void)argv;
  std::cout << "Parallel tests unavailable on sole process " << me << " with " << num_patches << " patches and l = " << desired_refinement_level << std::endl;
#endif

#ifndef SERIAL
  test_hypercube_2d((Index)me, (Index)num_patches, desired_refinement_level);
#endif

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
