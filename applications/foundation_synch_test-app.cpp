//#define SERIAL
#ifndef SERIAL
#include <mpi.h>
#endif
#include <iostream>
#include <limits>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/foundation/synch.hpp>

#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/space/dof_adjacency.hpp>
#include <kernel/space/dof_mirror.hpp>

#include <kernel/trafo/standard/mapping.hpp>

using namespace FEAST;
using namespace Foundation;

template<typename DT_>
struct TestResult
{
  TestResult(DT_ l, DT_ r, DT_ eps) :
    left(l),
    right(r),
    epsilon(eps)
  {
    passed = abs(l - r) < eps ? true : false;
  }

  TestResult()
  {
  }

  DT_ left;
  DT_ right;
  DT_ epsilon;
  bool passed;
};

typedef Geometry::ConformalMesh<Shape::Quadrilateral> quad_mesh_type_;
typedef Geometry::CellSubSet<Shape::Quadrilateral> quad_cell_subset_type_;
typedef Trafo::Standard::Mapping<quad_mesh_type_> quad_trafo_type_;
typedef Space::Lagrange1::Element<quad_trafo_type_> quad_space_type_;

void fill_quad_mesh_2d(quad_mesh_type_& mesh, double x = 0.0, double y = 0.0)
{
  // set up vertex coordinates array
  const double vtx0[4*2] =
  {
    0.0+x, 0.0+y,
    1.0+x, 0.0+y,
    0.0+x, 1.0+y,
    1.0+x, 1.0+y
  };

  // set up vertices-at-edge array
  static const Index v_e0[4*2] =
  {
    0, 1,
    2, 3,
    0, 2,
    1, 3
  };

  // set up vertices-at-quad array
  static const Index v_q0[1*4] =
  {
    0, 1, 2, 3
  };

  // set up edges-at-quad array
  static const Index e_q0[1*4] =
  {
    0,  1,  2,  3
  };

  Geometry::TestAux::copy_vtx(mesh.get_vertex_set(), vtx0);
  Geometry::TestAux::copy_idx(mesh.get_index_set<1,0>(), v_e0);
  Geometry::TestAux::copy_idx(mesh.get_index_set<2,0>(), v_q0);
  Geometry::TestAux::copy_idx(mesh.get_index_set<2,1>(), e_q0);

}

void fill_cell_set(quad_cell_subset_type_& cell, int face)
{
  switch(face)
  {
  case 0:
    cell.get_target_set<0>()[0] = 0;
    cell.get_target_set<0>()[1] = 1;
    cell.get_target_set<1>()[0] = 0;
    break;

  case 1:
    cell.get_target_set<0>()[0] = 2;
    cell.get_target_set<0>()[1] = 3;
    cell.get_target_set<1>()[0] = 1;
    break;

  case 2:
    cell.get_target_set<0>()[0] = 0;
    cell.get_target_set<0>()[1] = 2;
    cell.get_target_set<1>()[0] = 2;
    break;

  case 3:
    cell.get_target_set<0>()[0] = 1;
    cell.get_target_set<0>()[1] = 3;
    cell.get_target_set<1>()[0] = 3;
    break;
  }
}

template<typename DT1_, typename DT2_, typename DT3_>
TestResult<DT1_> test_check_equal_within_eps(DT1_ l, DT2_ r, DT3_ eps)
{
  return TestResult<DT1_>(l, r, eps);
}

void check_synch_mirror(int rank)
{

  static const Index num_entities[] = {4, 4, 1};
  static const Index num_cellset_entities[] = {2, 1, 0};

  quad_mesh_type_ mesh(num_entities);
  fill_quad_mesh_2d(mesh, double(-1.0));
  quad_cell_subset_type_ cell(num_cellset_entities);
  fill_cell_set(cell, 3);

  quad_trafo_type_ trafo(mesh);
  quad_space_type_ space(trafo);

  Graph dof_adj(Space::DofAdjacency<>::assemble(space));
  Graph dof_mirror(Space::DofMirror::assemble(space, cell));

  DenseVector<Mem::Main, double> target(space.get_num_dofs(), double(rank));
  LAFEM::VectorMirror<Mem::Main, double> target_mirror(dof_mirror);

  //dont use create_buffer(..) from mirror
  DenseVector<Mem::Main, double> sendbuf(target_mirror.size());
  DenseVector<Mem::Main, double> recvbuf(target_mirror.size());

  std::cout << "RANK " << rank << std::endl;
  std::cout << target << std::endl;
#ifndef SERIAL
  Synch<Archs::Parallel, com_exchange>::execute(target, target_mirror, sendbuf, recvbuf, rank == 0 ? 1 : 0, rank == 0 ? 1 : 0);
#endif
  std::cout << "RANK " << rank << std::endl;
  std::cout << target << std::endl;
}

int main(int argc, char* argv[])
{
  int me(0);
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

  check_synch_mirror(me);

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
