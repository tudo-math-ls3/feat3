#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the sparse vector blocked class.
*
* \test test description missing
*
* \tparam Mem_
* description missing
*
* \tparam DT_
* description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class SparseVectorBlockedTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  SparseVectorBlockedTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("SparseVectorBlockedTest")
  {
  }

  virtual void run() const
  {
    SparseVectorBlocked<Mem_, DT_, IT_, 2> zero1;
    SparseVectorBlocked<Mem::Main, DT_, IT_, 2> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    SparseVectorBlocked<Mem_, DT_, IT_, 2> a(10);
    TEST_CHECK_EQUAL(a, a);
    Tiny::Vector<DT_, 2> tv1(41);
    Tiny::Vector<DT_, 2> tv2(42);
    a(3, tv1);
    a(3, tv2);
    a(6, tv1);
    a(3, tv1);
    a(1, tv1);
    a(6, tv2);
    TEST_CHECK_EQUAL(a(3).v[0], tv1.v[0]);
    TEST_CHECK_EQUAL(a(1)[0], tv1[1]);
    TEST_CHECK_EQUAL(a(6)[1], tv2[1]);
    TEST_CHECK_EQUAL(a(2)[0], a.zero_element()[0]);

    SparseVectorBlocked<Mem_, DT_, IT_, 2> b(a.clone());
    TEST_CHECK_EQUAL(a, b);
    a(2, tv2);
    TEST_CHECK_NOT_EQUAL(a, b);
  }
};
SparseVectorBlockedTest<Mem::Main, NotSet, float, Index> cpu_sparse_vector_blocked_test_float;
SparseVectorBlockedTest<Mem::Main, NotSet, double, Index> cpu_sparse_vector_blocked_test_double;
#ifdef FEAST_BACKENDS_CUDA
SparseVectorBlockedTest<Mem::CUDA, NotSet, float, Index> cuda_sparse_vector_blocked_test_float;
SparseVectorBlockedTest<Mem::CUDA, NotSet, double, Index> cuda_sparse_vector_blocked_test_double;
#endif
