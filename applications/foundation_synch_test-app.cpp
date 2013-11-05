//#define SERIAL

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

#include <iostream>
#include <limits>

#ifndef SERIAL
#include <mpi.h>
#endif

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
    //passed = std::abs(l - r) < eps ? true : false;
    passed = (l < r ? r - l : l - r) < eps ? true : false;
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

  Adjacency::Graph dof_adj(Space::DofAdjacency<>::assemble(space));
  Adjacency::Graph dof_mirror(Space::DofMirror::assemble(space, cell));

  DenseVector<Mem::Main, double> target(space.get_num_dofs(), double(rank));
  LAFEM::VectorMirror<Mem::Main, double> target_mirror(dof_mirror);

  //dont use create_buffer(..) from mirror
  DenseVector<Mem::Main, double> sendbuf(target_mirror.size());
  DenseVector<Mem::Main, double> recvbuf(target_mirror.size());

#ifndef SERIAL
  SynchVec<Algo::Generic, Archs::Parallel, com_exchange>::execute(target, target_mirror, sendbuf, recvbuf, rank == 0 ? 1 : 0, rank == 0 ? 1 : 0);
#endif

  TestResult<double> res[4];
#ifndef SERIAL
    res[0] = test_check_equal_within_eps(target(0), rank == 0 ? double(0) : double(1), std::numeric_limits<double>::epsilon());
    res[1] = test_check_equal_within_eps(target(1), rank == 0 ? double(1) : double(0), std::numeric_limits<double>::epsilon());
    res[2] = test_check_equal_within_eps(target(2), rank == 0 ? double(0) : double(1), std::numeric_limits<double>::epsilon());
    res[3] = test_check_equal_within_eps(target(3), rank == 0 ? double(1) : double(0), std::numeric_limits<double>::epsilon());
#else
    res[0] = test_check_equal_within_eps(target(0), double(0), std::numeric_limits<double>::epsilon());
    res[1] = test_check_equal_within_eps(target(1), double(0), std::numeric_limits<double>::epsilon());
    res[2] = test_check_equal_within_eps(target(2), double(0), std::numeric_limits<double>::epsilon());
    res[3] = test_check_equal_within_eps(target(3), double(0), std::numeric_limits<double>::epsilon());
#endif

  bool passed(true);
  for(unsigned long i(0) ; i < 4 ; ++i)
    if(!res[i].passed)
    {
      std::cout << "Failed: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
      passed = false;
      break;
    }

  if(passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_synch-test (Tier-2: vertex-set based exchange (single target, single mirror))" << std::endl;
}

void check_synch_mirrors(int rank)
{

  static const Index num_entities[] = {4, 4, 1};
  static const Index num_cellset_entities[] = {2, 1, 0};

  quad_mesh_type_ mesh(num_entities);
  fill_quad_mesh_2d(mesh, double(-1.0));
  quad_cell_subset_type_ cell(num_cellset_entities);
  fill_cell_set(cell, 3);

  quad_trafo_type_ trafo(mesh);
  quad_space_type_ space(trafo);

  Adjacency::Graph dof_adj(Space::DofAdjacency<>::assemble(space));
  Adjacency::Graph dof_mirror(Space::DofMirror::assemble(space, cell));

  DenseVector<Mem::Main, double> target(space.get_num_dofs(), double(rank));
  LAFEM::VectorMirror<Mem::Main, double> target_mirror(dof_mirror);

  //dont use create_buffer(..) from mirror
  DenseVector<Mem::Main, double> sendbuf(target_mirror.size());
  DenseVector<Mem::Main, double> recvbuf(target_mirror.size());

  std::vector<LAFEM::VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(target_mirror);
  std::vector<LAFEM::DenseVector<Mem::Main, double> > sendbufs;
  sendbufs.push_back(sendbuf);
  std::vector<LAFEM::DenseVector<Mem::Main, double> > recvbufs;
  recvbufs.push_back(recvbuf);
  std::vector<Index> destranks;
  destranks.push_back(rank == 0 ? 1 : 0);
  std::vector<Index> sourceranks;
  sourceranks.push_back(rank == 0 ? 1 : 0);

#ifndef SERIAL
  SynchVec<Algo::Generic, Archs::Parallel, com_exchange>::execute(target, mirrors, sendbufs, recvbufs, destranks, sourceranks);
#endif

  TestResult<double> res[4];
#ifndef SERIAL
    res[0] = test_check_equal_within_eps(target(0), rank == 0 ? double(0) : double(1), std::numeric_limits<double>::epsilon());
    res[1] = test_check_equal_within_eps(target(1), rank == 0 ? double(1) : double(0), std::numeric_limits<double>::epsilon());
    res[2] = test_check_equal_within_eps(target(2), rank == 0 ? double(0) : double(1), std::numeric_limits<double>::epsilon());
    res[3] = test_check_equal_within_eps(target(3), rank == 0 ? double(1) : double(0), std::numeric_limits<double>::epsilon());
#else
    res[0] = test_check_equal_within_eps(target(0), double(0), std::numeric_limits<double>::epsilon());
    res[1] = test_check_equal_within_eps(target(1), double(0), std::numeric_limits<double>::epsilon());
    res[2] = test_check_equal_within_eps(target(2), double(0), std::numeric_limits<double>::epsilon());
    res[3] = test_check_equal_within_eps(target(3), double(0), std::numeric_limits<double>::epsilon());
#endif

  bool passed(true);
  for(unsigned long i(0) ; i < 4 ; ++i)
    if(!res[i].passed)
    {
      std::cout << "Failed: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
      passed = false;
      break;
    }

  if(passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_synch-test (Tier-2: vertex-set based exchange (single target, multiple mirrors))" << std::endl;
}

void check_synch_scal(int rank)
{
#ifndef SERIAL
  int size;
  float value(float(rank + 1));
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  float* send_buffer(new float[1]);
  send_buffer[0] = value;
  float* recv_buffer(new float[1]);

  SynchScal<Parallel, com_allreduce_sqrtsum>::execute(value, *send_buffer, *recv_buffer);

  TestResult<float> res;
  res = test_check_equal_within_eps(value, std::sqrt(3.0f), std::numeric_limits<float>::epsilon());

  if(!res.passed)
  {
    std::cout << "Failed (Tier-2: synch scalar): " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
  }
  else
    std::cout << "PASSED (rank " << rank <<"): foundation_synch_test (Tier-2: synch scalar)" << std::endl;

  delete[] send_buffer;
  delete[] recv_buffer;
#endif
}

int main(int argc, char* argv[])
{
  int me(0);
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
  std::cout<<"CTEST_FULL_OUTPUT"<<std::endl;

  check_synch_mirror(me);
  check_synch_mirrors(me);
  check_synch_scal(me);

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
