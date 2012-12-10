#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/foundation/communication.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/foundation/attribute.hpp>
#include <kernel/foundation/topology.hpp>
#include <kernel/foundation/mesh.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

#include <iostream>
#include <limits>
#ifndef SERIAL
#  include <mpi.h>
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

template<typename DT1_, typename DT2_, typename DT3_>
TestResult<DT1_> test_check_equal_within_eps(DT1_ l, DT2_ r, DT3_ eps)
{
  return TestResult<DT1_>(l, r, eps);
}

void check_sendrecv(int rank)
{
  float* f(new float[100000]);
  float* recvbuffer(new float[100000]);
  for(Index i(0) ; i < 100000 ; ++i)
  {
    f[i] = rank;
    recvbuffer[i] = rank;
  }

  if(rank < 2)
  {
#ifndef SERIAL
    Comm<Archs::Parallel>::send_recv(f,
                                     100000,
                                     rank == 0 ? 1 : 0,
                                     recvbuffer,
                                     100000,
                                     rank == 0 ? 1 : 0);
#else
    Comm<Archs::Serial>::send_recv(f, 100000, 0, recvbuffer, 100000, 0);
#endif

    TestResult<float> res[100000];
    for(unsigned long i(0) ; i < 100000 ; ++i)
#ifndef SERIAL
      res[i] = test_check_equal_within_eps(recvbuffer[i], rank == 0 ? float(1) : float(0), std::numeric_limits<float>::epsilon());
#else
      res[i] = test_check_equal_within_eps(recvbuffer[i], float(0), std::numeric_limits<float>::epsilon());
#endif

    bool passed(true);
    for(unsigned long i(0) ; i < 100000 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-0: sendrecv)" << std::endl;
  }

  delete[] f;
  delete[] recvbuffer;
}

void check_send_and_recv(int rank)
{
#ifndef FEAST_SERIAL_MODE
  float* f(new float[100000]);
  float* recvbuffer(new float[100000]);
  for(Index i(0) ; i < 100000 ; ++i)
  {
    f[i] = i;
    recvbuffer[i] = rank;
  }

  if(rank == 0)
  {
    Comm<Archs::Parallel>::send(f,
                                     100000,
                                     1);

    std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-0: send and recv)" << std::endl;
  }
  if(rank == 1)
  {
    Comm<Archs::Parallel>::recv(recvbuffer,
                                     100000,
                                     0);
    TestResult<float> res[100000];
    for(unsigned long i(0) ; i < 100000 ; ++i)
      res[i] = test_check_equal_within_eps(recvbuffer[i], i, std::numeric_limits<float>::epsilon());

    bool passed(true);
    for(unsigned long i(0) ; i < 100000 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "Failed (Tier-0: send and recv): " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-0: send and recv)" << std::endl;
  }

  delete[] f;
  delete[] recvbuffer;
#endif
}

void check_bcast(int rank)
{
#ifndef FEAST_SERIAL_MODE
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  float* recvbuffer(new float[size]);
  for(Index i(0) ; i < size ; ++i)
  {
    recvbuffer[i] = i + rank;
  }

  Comm<Archs::Parallel>::bcast(recvbuffer, size, 0);

  TestResult<float> res[size];
  for(unsigned long i(0) ; i < size ; ++i)
    res[i] = test_check_equal_within_eps(recvbuffer[i], i, std::numeric_limits<float>::epsilon());

  bool passed(true);
  for(unsigned long i(0) ; i < size ; ++i)
    if(!res[i].passed)
    {
      std::cout << "Failed (Tier-0: bcast): " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
      passed = false;
      break;
    }

  if(passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-0: bcast)" << std::endl;

  delete[] recvbuffer;
#endif
}

void check_scatter_gather(int rank)
{
#ifndef FEAST_SERIAL_MODE
  int size;
  float value(rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  float* buffer(new float[size]);

  Comm<Archs::Parallel>::gather(&value, 1, buffer, 1, 0);
  if (rank == 0)
    for(Index i(0) ; i < size ; ++i)
    {
      buffer[i] += 1;
    }
  Comm<Archs::Parallel>::scatter(buffer, 1, &value, 1, 0);

  TestResult<float> res;
  res = test_check_equal_within_eps(value, rank+1, std::numeric_limits<float>::epsilon());

  bool passed(true);
  if(!res.passed)
  {
    std::cout << "Failed (Tier-0: scatter_gather): " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
    passed = false;
  }

  if(passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-0: scatter_gather)" << std::endl;

  delete[] buffer;
#endif
}

void check_allgather(int rank)
{
#ifndef FEAST_SERIAL_MODE
  int size;
  float value(rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  float* buffer(new float[size]);

  Comm<Archs::Parallel>::allgather(&value, 1, buffer, 1);

  TestResult<float> res[size];
  for(unsigned long i(0) ; i < size ; ++i)
    res[i] = test_check_equal_within_eps(buffer[i], i, std::numeric_limits<float>::epsilon());

  bool passed(true);
  for(unsigned long i(0) ; i < size ; ++i)
    if(!res[i].passed)
    {
      std::cout << "Failed (Tier-0: allgather): " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
      passed = false;
      break;
    }

  if(passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-0: allgather)" << std::endl;

  delete[] buffer;
#endif
}

void check_allreduce(int rank)
{
#ifndef FEAST_SERIAL_MODE
  int size;
  float value(rank + 1);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  float* send_buffer(new float[1]);
  send_buffer[0] = value;
  float* recv_buffer(new float[1]);

  Comm<Archs::Parallel>::allreduce(send_buffer, 1, recv_buffer);

  TestResult<float> res;
  res = test_check_equal_within_eps(recv_buffer[0], float(3), std::numeric_limits<float>::epsilon());

  if(!res.passed)
  {
    std::cout << "Failed (Tier-0: allreduce): " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
  }
  else
    std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-0: allreduce)" << std::endl;

  delete[] send_buffer;
  delete[] recv_buffer;
#endif
}

void check_reduce(int rank)
{
#ifndef FEAST_SERIAL_MODE
  int size;
  float value(rank + 1);
  float result(0);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Comm<Archs::Parallel>::reduce(&value, &result, 1, MPI_SUM, 0);

  if(rank == 0)
  {
    TestResult<float> res;
    res = test_check_equal_within_eps(result, size*(size+1)/2, std::numeric_limits<float>::epsilon());

    bool passed(true);
    if(!res.passed)
    {
      std::cout << "Failed (Tier-0: reduce): " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
      passed = false;
    }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-0: reduce)" << std::endl;
  }

#endif
}

void check_halo_transfer(int rank)
{
  if(rank < 2)
  {
    Mesh<> m(rank);
    Halo<0, pl_face, Mesh<> > h(m, rank == 0 ? 1 : 0);
    h.add_element_pair(rank, rank);

    //halo at rank 1 is bigger
    if(rank == 1)
      h.add_element_pair(rank, rank);


    Halo<0, pl_face, Mesh<> >::buffer_type_ sendbuf(h.buffer());
    Halo<0, pl_face, Mesh<> >::buffer_type_ recvbuf(h.buffer(1));

    h.to_buffer(sendbuf);

    h.send_recv(
        sendbuf,
        rank == 0 ? 1 : 0,
        recvbuf,
        rank == 0 ? 1 : 0);

    h.from_buffer(recvbuf);

    bool passed(true);
#ifndef SERIAL
    //TestResult<Index> res[rank == 0 ? 4 : 2];
    TestResult<Index> res[4];
    res[0] = test_check_equal_within_eps(h.get_element(0), rank == 0 ? Index(1) : Index(0), Index(1));
    res[1] = test_check_equal_within_eps(h.get_element_counterpart(0), rank == 0 ? Index(1) : Index(0), Index(1));
    if(rank == 0) //rank 0 receives more from rank 1
    {
      res[2] = test_check_equal_within_eps(h.get_element(1), Index(1), Index(1));
      res[3] = test_check_equal_within_eps(h.get_element_counterpart(1), Index(1), Index(1));
    }
    for(unsigned long i(0) ; i < (rank == 0 ? 4 : 2) ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: (rank " << rank << "): " <<  res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-1: halo_transfer)" << std::endl;
#else
    TestResult<Index> res[2];
    res[0] = test_check_equal_within_eps(h.get_element(0), rank, Index(1));
    res[1] = test_check_equal_within_eps(h.get_element_counterpart(0), rank, Index(1));
    for(unsigned long i(0) ; i < 2 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: (rank " << rank << "): " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-1: halo_transfer)" << std::endl;
#endif

  }
}

void check_attribute_transfer(int rank)
{
  if(rank < 2)
  {
    Attribute<double> attr;
    attr.push_back(double(rank));
    attr.push_back(double(rank + 42));

    if(rank == 1)
      attr.push_back(double(rank + 10000));

    Attribute<double>::buffer_type_ sendbuf(attr.buffer());
    Attribute<double>::buffer_type_ recvbuf(attr.buffer(10));

    attr.to_buffer(sendbuf);

    attr.send_recv(
        sendbuf,
        rank == 0 ? 1 : 0,
        recvbuf,
        rank == 0 ? 1 : 0);

    attr.from_buffer(recvbuf);

    bool passed(true);
#ifndef SERIAL
    //TestResult<double> res[rank == 0 ? 3 : 2];
    TestResult<double> res[3];
    res[0] = test_check_equal_within_eps(attr.at(0), rank == 0 ? double(1) : double(0), std::numeric_limits<double>::epsilon());
    res[1] = test_check_equal_within_eps(attr.at(1), rank == 0 ? double(43) : double(42), std::numeric_limits<double>::epsilon());
    if(rank == 0)
      res[2] = test_check_equal_within_eps(attr.at(2), double(10001), std::numeric_limits<double>::epsilon());

    for(unsigned long i(0) ; i < (rank == 0 ? 3 : 2) ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: (rank " << rank << "): " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-1: attr_transfer)" << std::endl;
#else
    TestResult<double> res[2];
    res[0] = test_check_equal_within_eps(attr.at(0), double(0), std::numeric_limits<double>::epsilon());
    res[1] = test_check_equal_within_eps(attr.at(1), double(42), std::numeric_limits<double>::epsilon());
    for(unsigned long i(0) ; i < 2 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: (rank " << rank << "): " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-1: attr_transfer)" << std::endl;
#endif
  }
}

void check_topology_transfer(int rank)
{
  if(rank < 2)
  {
    Foundation::Topology<> t;
    t.push_back();
    t.at(0).push_back(rank == 0 ? 42 : 43);
    t.at(0).push_back(rank == 0 ? 47 : 48);
    t.push_back();
    t.at(1).push_back(rank == 0 ? 52 : 53);
    t.at(1).push_back(rank == 0 ? 57 : 58);

    if(rank == 1)
      t.at(1).push_back(100);

    Topology<>::buffer_type_ sendbuf(t.buffer());
    Topology<>::buffer_type_ recvbuf(t.buffer(1));

    t.to_buffer(sendbuf);

    t.send_recv(
        sendbuf,
        rank == 0 ? 1 : 0,
        recvbuf,
        rank == 0 ? 1 : 0);

    //Foundation::Topology<> t2;
    t.from_buffer(recvbuf);

    bool passed(true);
#ifndef SERIAL
    //TestResult<Index> res[rank == 0 ? 5 : 4];
    TestResult<Index> res[5];
    res[0] = test_check_equal_within_eps(t.at(0).at(0), rank == 0 ? Index(43) : Index(42), Index(1));
    res[1] = test_check_equal_within_eps(t.at(0).at(1), rank == 0 ? Index(48) : Index(47), Index(1));
    res[2] = test_check_equal_within_eps(t.at(1).at(0), rank == 0 ? Index(53) : Index(52), Index(1));
    res[3] = test_check_equal_within_eps(t.at(1).at(1), rank == 0 ? Index(58) : Index(57), Index(1));
    if(rank == 0)
      res[4] = test_check_equal_within_eps(t.at(1).at(2), Index(100), Index(1));

    for(Index i(0) ; i < (rank == 0 ? 5 : 4) ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-1: topology_transfer)" << std::endl;
#else
    TestResult<Index> res[4];
    res[0] = test_check_equal_within_eps(t.at(0).at(0), Index(42), Index(1));
    res[1] = test_check_equal_within_eps(t.at(0).at(1), Index(47), Index(1));
    res[2] = test_check_equal_within_eps(t.at(1).at(0), Index(52), Index(1));
    res[3] = test_check_equal_within_eps(t.at(1).at(1), Index(57), Index(1));

    for(Index i(0) ; i < 2 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-1: topology_transfer)" << std::endl;
#endif
  }
}

void check_mesh_transfer(int rank)
{
  if(rank < 2)
  {
    Foundation::Mesh<Foundation::rnt_2D> m(0);

    //add vertices
    m.add_polytope(Foundation::pl_vertex);
    m.add_polytope(Foundation::pl_vertex);
    m.add_polytope(Foundation::pl_vertex);
    m.add_polytope(Foundation::pl_vertex);
    m.add_polytope(Foundation::pl_vertex);
    m.add_polytope(Foundation::pl_vertex);

    //add edges
    m.add_polytope(Foundation::pl_edge);
    m.add_polytope(Foundation::pl_edge);
    m.add_polytope(Foundation::pl_edge);
    m.add_polytope(Foundation::pl_edge);
    m.add_polytope(Foundation::pl_edge);
    m.add_polytope(Foundation::pl_edge);
    m.add_polytope(Foundation::pl_edge);

    //add faces
    m.add_polytope(Foundation::pl_face);
    m.add_polytope(Foundation::pl_face);

    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 0, 0); //v->e is set automagically
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 0, 1);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 1, 1);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 1, 2);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 2, 0);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 2, 3);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 3, 1);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 3, 4);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 4, 2);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 4, 5);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 5, 3);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 5, 4);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 6, 4);
    m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 6, 5);

    m.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 0); //v->f is set automagically
    m.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 1); //v->f is set automagically
    m.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 3); //v->f is set automagically
    m.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 4); //v->f is set automagically
    m.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 1); //v->f is set automagically
    m.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 2); //v->f is set automagically
    m.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 4); //v->f is set automagically
    m.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 5); //v->f is set automagically

    m.send_recv(m.get_topologies(),
                rank == 0 ? 1 : 0,
                rank == 0 ? 1 : 0);

    TestResult<Index> res[13];
    res[0] = test_check_equal_within_eps(m.get_topologies().at(ipi_vertex_edge).size(), Index(6), Index(1));
    res[1] = test_check_equal_within_eps(m.get_topologies().at(ipi_vertex_face).size(), Index(6), Index(1));
    res[2] = test_check_equal_within_eps(m.get_topologies().at(ipi_face_vertex).size(), Index(2), Index(1));
    res[3] = test_check_equal_within_eps(m.get_topologies().at(ipi_edge_vertex).size(), Index(7), Index(1));
    res[4] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_edge_vertex).at(0).size(), Index(2), Index(1));
    res[5] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_edge_vertex).at(1).size(), Index(2), Index(1));
    res[6] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_edge_vertex).at(2).size(), Index(2), Index(1));
    res[7] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_edge_vertex).at(3).size(), Index(2), Index(1));
    res[8] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_edge_vertex).at(4).size(), Index(2), Index(1));
    res[9] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_edge_vertex).at(5).size(), Index(2), Index(1));
    res[10] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_edge_vertex).at(6).size(), Index(2), Index(1));
    res[11] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_face_vertex).at(0).size(), Index(4), Index(1));
    res[12] = test_check_equal_within_eps<Index>(m.get_topologies().at(ipi_face_vertex).at(1).size(), Index(4), Index(1));

    bool passed(true);
    for(Index i(0) ; i < 13 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-1: mesh_transfer)" << std::endl;

  }
}

void check_halobased_attribute_transfer(int rank)
{
  if(rank < 2)
  {
    Attribute<double> attr;
    for(Index i(0) ; i < 100000 ; ++i)
    {
      attr.push_back(double(rank + i));
    }
    attr.at(10) = double(42 + rank);
    attr.at(100) = double(47 + rank);

    Mesh<> m(rank);
    Halo<0, pl_face, Mesh<> > h(m, rank == 0 ? 1 : 0);
    h.add_element_pair(10, 999);
    h.add_element_pair(100, 999);

    InterfacedComm<Mem::Main, com_exchange>::execute(h, attr);

    TestResult<double> res[2];
#ifndef SERIAL
    res[0] = test_check_equal_within_eps(attr.at(10), rank == 0 ? double(43) : double(42), std::numeric_limits<float>::epsilon());
    res[1] = test_check_equal_within_eps(attr.at(100), rank == 0 ? double(48) : double(47), std::numeric_limits<float>::epsilon());
#else
    res[0] = test_check_equal_within_eps(attr.at(10), double(42), std::numeric_limits<float>::epsilon());
    res[1] = test_check_equal_within_eps(attr.at(100), double(47), std::numeric_limits<float>::epsilon());
#endif
    bool passed(true);
    for(Index i(0) ; i < 2 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-2: halo-based attribute transfer)" << std::endl;
  }
}

void check_halobased_dv_transfer(int rank)
{
  if(rank < 2)
  {
    DenseVector<Mem::Main, double> attr(100000);
    for(Index i(0) ; i < 100000 ; ++i)
    {
      attr(i, (double(rank + i)));
    }
    attr(10, double(42 + rank));
    attr(100, double(47 + rank));

    Mesh<> m(rank);
    Halo<0, pl_face, Mesh<> > h(m, rank == 0 ? 1 : 0);
    h.add_element_pair(10, 999);
    h.add_element_pair(100, 999);

    InterfacedComm<Mem::Main, com_exchange>::execute(h, attr);

    TestResult<double> res[2];
#ifndef SERIAL
    res[0] = test_check_equal_within_eps(attr(10), rank == 0 ? double(43) : double(42), std::numeric_limits<float>::epsilon());
    res[1] = test_check_equal_within_eps(attr(100), rank == 0 ? double(48) : double(47), std::numeric_limits<float>::epsilon());
#else
    res[0] = test_check_equal_within_eps(attr(10), double(42), std::numeric_limits<float>::epsilon());
    res[1] = test_check_equal_within_eps(attr(100), double(47), std::numeric_limits<float>::epsilon());
#endif
    bool passed(true);
    for(Index i(0) ; i < 2 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-2: halo-based lafem dv transfer)" << std::endl;
  }
}

void check_halobased_smcsr_transfer(int rank)
{
  if(rank < 2)
  {
    SparseMatrixCOO<Mem::Main, float> smcoo(1000, 1000);
    smcoo(0, 0, float(rank + 42));
    smcoo(0, 1, float(rank + 43));
    smcoo(1, 0, float(rank + 47));
    smcoo(1, 1, float(rank + 48));

    SparseMatrixCSR<Mem::Main, float> smcsr(smcoo);

    Mesh<> m(rank);
    Halo<0, pl_vertex, Mesh<> > h(m, rank == 0 ? 1 : 0);
    h.add_element_pair(0, 999);
    h.add_element_pair(1, 999);

    InterfacedComm<Mem::Main, com_exchange>::execute(h, smcsr);

    TestResult<float> res[4];
#ifndef SERIAL
    res[0] = test_check_equal_within_eps(smcsr(0, 0), rank == 0 ? float(43) : float(42), std::numeric_limits<float>::epsilon());
    res[1] = test_check_equal_within_eps(smcsr(0, 1), rank == 0 ? float(44) : float(43), std::numeric_limits<float>::epsilon());
    res[2] = test_check_equal_within_eps(smcsr(1, 0), rank == 0 ? float(48) : float(47), std::numeric_limits<float>::epsilon());
    res[3] = test_check_equal_within_eps(smcsr(1, 1), rank == 0 ? float(49) : float(48), std::numeric_limits<float>::epsilon());
#else
    res[0] = test_check_equal_within_eps(smcsr(0, 0), float(42), std::numeric_limits<float>::epsilon());
    res[1] = test_check_equal_within_eps(smcsr(0, 1), float(43), std::numeric_limits<float>::epsilon());
    res[2] = test_check_equal_within_eps(smcsr(1, 0), float(47), std::numeric_limits<float>::epsilon());
    res[3] = test_check_equal_within_eps(smcsr(1, 1), float(48), std::numeric_limits<float>::epsilon());
#endif
    bool passed(true);
    for(Index i(0) ; i < 4 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (Tier-2: halo-based lafem smcsr transfer)" << std::endl;
  }
}

int main(int argc, char* argv[])
{
  int me(0);
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

  check_sendrecv(me);
  check_send_and_recv(me);
  check_bcast(me);
  check_scatter_gather(me);
  check_allgather(me);
  check_allreduce(me);
  check_reduce(me);
  check_halo_transfer(me);
  check_attribute_transfer(me);
  check_topology_transfer(me);
  check_mesh_transfer(me);
  check_halobased_attribute_transfer(me);
  check_halobased_dv_transfer(me);
  check_halobased_smcsr_transfer(me);

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
