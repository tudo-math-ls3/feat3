#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/util/binary_stream.hpp>

#include <cstdio>
#include <sstream>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the sparse matrix coo class.
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
class SparseMatrixELLTest
  : public TaggedTest<Tag_, DT_>
{

public:

  SparseMatrixELLTest()
    : TaggedTest<Tag_, DT_>("sparse_matrix_ell_test")
  {
  }

  virtual void run() const
  {
    SparseMatrixCOO<Mem::Main, DT_> a(10, 12);
    a(1,2,7);
    a.clear();
    a(1,2,7);
    a(5,5,2);
    SparseMatrixELL<Tag_, DT_> b(a);
    TEST_CHECK_EQUAL(b.used_elements(), 2ul);
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(b(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(b(5, 5), a(5, 5));

    //SparseMatrixCSR<Tag_, DT_> b2(b);
    SparseMatrixCSR<Tag_, DT_> b2(a);
    SparseMatrixELL<Tag_, DT_> b3(b);
    TEST_CHECK_EQUAL(b3, b);

    SparseMatrixELL<Tag_, DT_> z(b);
    TEST_CHECK_EQUAL(z.used_elements(), 2ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z.stride(), b.stride());
    TEST_CHECK_EQUAL(z.num_cols_per_row(), b.num_cols_per_row());
    TEST_CHECK_EQUAL(z(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(z(5, 5), a(5, 5));
    TEST_CHECK_EQUAL(z(1, 3), a(1, 3));

    SparseMatrixELL<Tag_, DT_> c;
    c = b;
    TEST_CHECK_EQUAL(c.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(c(0,2), b(0,2));
    TEST_CHECK_EQUAL(c(1,2), b(1,2));
    TEST_CHECK_EQUAL(c, b);

    SparseMatrixELL<Mem::Main, DT_> e(c);
    TEST_CHECK_EQUAL(e, c);
    e = c;
    TEST_CHECK_EQUAL(e, c);
    e = c.clone();
    TEST_CHECK_EQUAL(e, c);

    SparseMatrixCOO<Mem::Main, DT_> fcoo(10, 10);
    for (Index row(0) ; row < fcoo.rows() ; ++row)
    {
      for (Index col(0) ; col < fcoo.columns() ; ++col)
      {
        if(row == col)
          fcoo(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          fcoo(row, col, DT_(-1));
      }
    }
    SparseMatrixELL<Tag_, DT_> f(fcoo);

    BinaryStream bs;
    f.write_out(fm_ell, bs);
    bs.seekg(0);
    SparseMatrixELL<Tag_, DT_> g(bs);
    TEST_CHECK_EQUAL(g, f);

    std::stringstream ts;
    f.write_out(fm_m, ts);
    SparseMatrixCOO<Mem::Main, DT_> h(f);
    SparseMatrixCOO<Mem::Main, DT_> i(fm_m, ts);
    TEST_CHECK_EQUAL(i, h);
  }
};
SparseMatrixELLTest<Mem::Main, float> cpu_sparse_matrix_ell_test_float;
SparseMatrixELLTest<Mem::Main, double> cpu_sparse_matrix_ell_test_double;
#ifdef FEAST_GMP
SparseMatrixELLTest<Mem::Main, mpf_class> cpu_sparse_matrix_ell_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixELLTest<Mem::CUDA, float> cuda_sparse_matrix_ell_test_float;
SparseMatrixELLTest<Mem::CUDA, double> cuda_sparse_matrix_ell_test_double;
#endif
