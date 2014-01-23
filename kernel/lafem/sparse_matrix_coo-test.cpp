#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/util/binary_stream.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the sparse matrix coo class.
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
  typename DT_>
class SparseMatrixCOOTest
  : public TaggedTest<Mem_, DT_>
{

public:

  SparseMatrixCOOTest()
    : TaggedTest<Mem_, DT_>("sparse_matrix_coo_test")
  {
  }

  virtual void run() const
  {
    SparseMatrixCOO<Mem_, DT_> x;
    SparseMatrixCOO<Mem_, DT_> a(10, 10);
    a(5,5,2);
    a(1,2,7);
    a(5,5,2);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);

    a.clear();
    a(1,2,7);
    a(5,5,2);
    a(5,5,2);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);

    a.clear();
    a(1,2,7);
    a(5,5,2);
    a(1,2,7);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);

    SparseMatrixCOO<Mem_, DT_> b(a);
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(a(1,2), b(1,2));
    TEST_CHECK_EQUAL(a(0,2), b(0,2));
    TEST_CHECK_EQUAL(a, b);

    SparseMatrixCOO<Mem_, DT_> c(10, 10);
    c = b;
    TEST_CHECK_EQUAL(c(0,2), b(0,2));
    TEST_CHECK_EQUAL(c(1,2), b(1,2));
    TEST_CHECK_EQUAL(c, b);
    TEST_CHECK_EQUAL(c.used_elements(), b.used_elements());

    c = b.clone();
    TEST_CHECK_EQUAL(c, b);
    c(1,2,3);
    TEST_CHECK_NOT_EQUAL(c, b);

    SparseMatrixCOO<Mem_, DT_> f(10, 10);
    for (Index row(0) ; row < f.rows() ; ++row)
    {
      for (Index col(0) ; col < f.columns() ; ++col)
      {
        if(row == col)
          f(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          f(row, col, DT_(-1));
      }
    }

    BinaryStream bs;
    f.write_out(fm_coo, bs);
    bs.seekg(0);
    SparseMatrixCOO<Mem_, DT_> g(fm_coo, bs);
    TEST_CHECK_EQUAL(g, f);

    std::stringstream ts;
    f.write_out(fm_m, ts);
    SparseMatrixCOO<Mem_, DT_> i(fm_m, ts);
    TEST_CHECK_EQUAL(i, f);

    std::stringstream ms;
    f.write_out(fm_mtx, ms);
    SparseMatrixCOO<Mem_, DT_> j(fm_mtx, ms);
    TEST_CHECK_EQUAL(j, f);
  }
};
SparseMatrixCOOTest<Mem::Main, float> sparse_matrix_coo_test_float;
SparseMatrixCOOTest<Mem::Main, double> sparse_matrix_coo_test_double;
#ifdef FEAST_GMP
SparseMatrixCOOTest<Mem::Main, mpf_class> sparse_matrix_coo_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCOOTest<Mem::CUDA, float> cuda_sparse_matrix_coo_test_float;
SparseMatrixCOOTest<Mem::CUDA, double> cuda_sparse_matrix_coo_test_double;
#endif
