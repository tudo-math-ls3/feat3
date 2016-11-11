#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/sparse_vector.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse vector class.
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
  typename DT_,
  typename IT_>
class SparseVectorTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  SparseVectorTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseVectorTest")
  {
  }

  virtual void run() const override
  {
    SparseVector<Mem_, DT_, IT_> zero1;
    SparseVector<Mem::Main, DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    SparseVector<Mem_, DT_, IT_> a(10);
    a(3, DT_(7));
    a(3, DT_(3));
    a(6, DT_(1));
    a(5, DT_(6));
    a(6, DT_(8));
    TEST_CHECK_EQUAL(a.used_elements(), Index(3));
    TEST_CHECK_EQUAL(a(3), DT_(3));
    TEST_CHECK_EQUAL(a(2), DT_(0));
    TEST_CHECK_EQUAL(a(5), DT_(6));
    TEST_CHECK_EQUAL(a(6), DT_(8));

    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rng(seed);
    Adjacency::Permutation prm_rnd(a.size(), rng);
    SparseVector<Mem_, DT_, IT_> ap(a.clone());
    ap.permute(prm_rnd);
    prm_rnd = prm_rnd.inverse();
    ap.permute(prm_rnd);
    TEST_CHECK_EQUAL(ap, a);
    TEST_CHECK_EQUAL(ap.used_elements(), Index(3));

    std::stringstream ts;
    a.write_out(FileMode::fm_mtx, ts);
    SparseVector<Mem::Main, DT_, IT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_EQUAL(j, a);

    SparseVector<Mem_, DT_, IT_> b;
    b.convert(a);
    TEST_CHECK_EQUAL(a, b);
    b(6, DT_(1));
    TEST_CHECK_NOT_EQUAL(a, b);
    b.clone(a);
    b(6, DT_(3));
    TEST_CHECK_NOT_EQUAL(a, b);
    TEST_CHECK_NOT_EQUAL((void*)a.elements(), (void*)b.elements());
    TEST_CHECK_NOT_EQUAL((void*)a.indices(), (void*)b.indices());
    b = a.clone();
    TEST_CHECK_NOT_EQUAL((void*)a.elements(), (void*)b.elements());
    TEST_CHECK_NOT_EQUAL((void*)a.indices(), (void*)b.indices());

    SparseVector<Mem::Main, float, unsigned int> c;
    c.convert(a);
    SparseVector<Mem::Main, float, unsigned int> d;
    d.clone(c);
    SparseVector<Mem::Main, float, unsigned int> e;
    e.convert(a);
    TEST_CHECK_EQUAL(d, e);
    c(6, DT_(1));
    TEST_CHECK_NOT_EQUAL(c, e);

    a.format();
    TEST_CHECK_EQUAL(a.used_elements(), Index(3));
    TEST_CHECK_EQUAL(a(2), DT_(0));
    TEST_CHECK_EQUAL(a(3), DT_(0));

    auto op = b.serialise();
    SparseVector<Mem_, DT_, IT_> o(op);
    for (Index i(0) ; i < b.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(o(i), b(i), 1e-5);

    //increase vector size above alloc_increment
    SparseVector<Mem_, DT_, IT_> p(3001);
    for (Index i(1) ; i <= p.size() ; ++i)
    {
      p(p.size() - i, DT_(i));
    }
  }
};
SparseVectorTest<Mem::Main, float, Index> cpu_sparse_vector_test_float;
SparseVectorTest<Mem::Main, double, Index> cpu_sparse_vector_test_double;
//SparseVectorTest<Mem::Main, Index> cpu_sparse_vector_test_index;
#ifdef FEAT_HAVE_CUDA
SparseVectorTest<Mem::CUDA, float, Index> cuda_sparse_vector_test_float;
SparseVectorTest<Mem::CUDA, double, Index> cuda_sparse_vector_test_double;
//SparseVectorTest<Mem::CUDA, Index> cuda_sparse_vector_test_index;
#endif
