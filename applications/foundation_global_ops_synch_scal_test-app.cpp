#include <kernel/base_header.hpp>
#include <kernel/foundation/comm_base.hpp>
#include <kernel/foundation/global_synch_scal.hpp>
#include <kernel/archs.hpp>

#include <iostream>
#include <limits>

using namespace FEAT;
using namespace Foundation;
using namespace LAFEM;

template<typename DT1_, typename DT2_, typename DT3_>
struct TestResult
{
  TestResult(DT1_ l, DT2_ r, DT3_ eps) :
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

  DT1_ left;
  DT2_ right;
  DT3_ epsilon;
  bool passed;
};

template<typename DT1_, typename DT2_, typename DT3_>
TestResult<DT1_, DT2_, DT3_> test_check_equal_within_eps(DT1_ l, DT2_ r, DT3_ eps)
{
  return TestResult<DT1_, DT2_, DT3_>(l, r, eps);
}

template<typename DT_>
void check_global_synch_scal0(Index rank)
{
  auto nprocs(Comm::size());

  DT_ x(2), r;

  GlobalSynchScal0<Mem::Main>::value(r, x);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(r, x * nprocs, std::numeric_limits<DT_>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_synch_scal0" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

int main(int argc, char* argv[])
{
  int me(0);
#ifdef FEAT_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
  (void)argc;
  (void)argv;
  std::cout<<"CTEST_FULL_OUTPUT"<<std::endl;

//#ifdef FEAT_HAVE_MPI
  check_global_synch_scal0<double>((Index)me);
//#else
//  std::cout << "Parallel tests unavailable on sole process " << me << std::endl;
//#endif

#ifdef FEAT_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
