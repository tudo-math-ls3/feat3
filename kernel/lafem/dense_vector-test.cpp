#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <list>

#include <cstdio>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the dense vector class.
*
* \test test description missing
*
* \tparam Tag_
* description missing
*
* \tparam DT_
* description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Tag_,
  typename DT_>
class DenseVectorTest
  : public TaggedTest<Tag_, DT_>
{

public:

  DenseVectorTest()
    : TaggedTest<Tag_, DT_>("dense_vector_test")
  {
  }

  virtual void run() const
  {
    DenseVector<Tag_, DT_> a(10, DT_(7));
    DenseVector<Tag_, DT_> b(10, DT_(5));
    b(7, DT_(42));
    DenseVector<Tag_, DT_> c(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c(7), b(7));
    TEST_CHECK_EQUAL(c, b);
    std::list<DenseVector<Tag_, DT_> > list;
    list.push_back(a);
    list.push_back(b);
    list.push_back(c);
    DenseVector<Tag_, DT_> d = a;
    list.push_back(d);
    DenseVector<Tag_, DT_> e(10, DT_(42));
    e = a;
    TEST_CHECK_EQUAL(e(5), a(5));
    TEST_CHECK_EQUAL(e, a);

    DenseVector<Mem::Main, DT_> f(e);
    DenseVector<Mem::Main, DT_> g;
    g = e;
    TEST_CHECK_EQUAL(f, e);
    TEST_CHECK_EQUAL(g, f);
    TEST_CHECK_EQUAL(g, e);

    DenseVector<Mem::Main, DT_> h(g.clone());
    TEST_CHECK_EQUAL(h, g);
    h(1, DT_(5));
    TEST_CHECK_NOT_EQUAL(h, g);
    TEST_CHECK_NOT_EQUAL((unsigned long)h.elements(), (unsigned long)g.elements());

    DenseVector<Tag_, DT_> k(123);
    for (Index i(0) ; i < k.size() ; ++i)
      k(i, DT_(i) / DT_(12));

    k.write_out(fm_exp, "test.exp");
    DenseVector<Tag_, DT_> l(fm_exp, "test.exp");
    for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(l(i), k(i), 1e-5);

    remove("test.exp");
  }
};
DenseVectorTest<Mem::Main, float> cpu_dense_vector_test_float;
DenseVectorTest<Mem::Main, double> cpu_dense_vector_test_double;
DenseVectorTest<Mem::Main, Index> cpu_dense_vector_test_index;
#ifdef FEAST_BACKENDS_CUDA
DenseVectorTest<Mem::CUDA, float> cuda_dense_vector_test_float;
DenseVectorTest<Mem::CUDA, double> cuda_dense_vector_test_double;
DenseVectorTest<Mem::CUDA, Index> cuda_dense_vector_test_index;
#endif
