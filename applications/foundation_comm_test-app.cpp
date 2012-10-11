#ifndef FEAST_SERIAL_MODE
#include <mpi.h>
#endif
#include <iostream>
#include <limits>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/foundation/communication.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/lafem/dense_vector.hpp>

using namespace FEAST;
using namespace Foundation;
using namespace Archs;

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
#ifndef FEAST_SERIAL_MODE
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
#ifndef FEAST_SERIAL_MODE
      res[i] = test_check_equal_within_eps(recvbuffer[i], rank == 0 ? float(1) : float(0), std::numeric_limits<float>::epsilon());
#else
      res[i] = test_check_equal_within_eps(recvbuffer[i], float(0), std::numeric_limits<float>::epsilon());
#endif

    bool passed(true);
    for(unsigned long i(0) ; i < 100000 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "Failed: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (sendrecv)" << std::endl;
  }
}

void check_halo_transfer(int rank)
{
  if(rank < 2)
  {
    Mesh<> m(rank);
    Halo<0, pl_face, Mesh<> > h(m, rank == 0 ? 1 : 0);
    h.add_element_pair(rank, rank);

    Halo<0, pl_face, Mesh<> >::buffer_type_ sendbuf(h.buffer());
    Halo<0, pl_face, Mesh<> >::buffer_type_ recvbuf(h.buffer());
    h.to_buffer(sendbuf);

    h.send_recv(
        sendbuf,
        rank == 0 ? 1 : 0,
        recvbuf,
        rank == 0 ? 1 : 0);

    h.from_buffer(recvbuf);

    TestResult<Index> res[2];
#ifndef FEAST_SERIAL_MODE
    res[0] = test_check_equal_within_eps(h.get_element(0), rank == 0 ? Index(1) : Index(0), Index(1));
    res[1] = test_check_equal_within_eps(h.get_element_counterpart(0), rank == 0 ? Index(1) : Index(0), Index(1));
#else
    res[0] = test_check_equal_within_eps(h.get_element(0), Index(0), Index(1));
    res[1] = test_check_equal_within_eps(h.get_element_counterpart(0), Index(0), Index(1));
#endif

    bool passed(true);
    for(unsigned long i(0) ; i < 2 ; ++i)
      if(!res[i].passed)
      {
        std::cout << "Failed: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "!" << std::endl;
        passed = false;
        break;
      }

    if(passed)
      std::cout << "PASSED (rank " << rank <<"): foundation_comm_test (halo_transfer)" << std::endl;
  }
}

int main(int argc, char* argv[])
{
  int me(0);
#ifndef FEAST_SERIAL_MODE
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

  check_sendrecv(me);
  check_halo_transfer(me);

#ifndef FEAST_SERIAL_MODE
  MPI_Finalize();
#endif

  return 0;
}
