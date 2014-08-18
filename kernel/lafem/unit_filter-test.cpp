#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for UnitFilter class template
 *
 * \author Peter Zajac
 */
template<
  typename Algo_,
  typename DT_,
  typename IT_>
class UnitFilterTest
  : public FullTaggedTest<typename Algo_::MemType, Algo_, DT_, IT_>
{
  typedef DenseVector<typename Algo_::MemType, DT_, IT_> VectorType;
  typedef DenseVector<typename Algo_::MemType, IT_, IT_> IVectorType;
  typedef UnitFilter<typename Algo_::MemType, DT_, IT_> FilterType;
public:
  UnitFilterTest()
    : FullTaggedTest<typename Algo_::MemType, Algo_, DT_, IT_>("UnitFilterTest")
  {
  }

  void test_vector() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

    const Index n = 7;
    VectorType a1(n, DT_(1));
    VectorType a2(n, DT_(1));
    VectorType ar(n, DT_(1));
    VectorType b1(n, DT_(1));
    VectorType b2(n, DT_(1));
    VectorType br(n, DT_(1));

    // modify reference results
    ar(0, DT_(0));
    ar(2, DT_(0));
    ar(6, DT_(0));
    br(0, DT_(3));
    br(2, DT_(5));
    br(6, DT_(9));

    // create filter
    FilterType filter(IT_(7));
    filter.add(IT_(0), DT_(3));
    filter.add(IT_(2), DT_(5));
    filter.add(IT_(6), DT_(9));

    // apply the filter
    filter.template filter_def<Algo_>(a1);
    filter.template filter_cor<Algo_>(a2);
    filter.template filter_rhs<Algo_>(b1);
    filter.template filter_sol<Algo_>(b2);

    // subtract reference results
    a1.template axpy<Algo_>(ar, a1, -DT_(1));
    a2.template axpy<Algo_>(ar, a2, -DT_(1));
    b1.template axpy<Algo_>(br, b1, -DT_(1));
    b2.template axpy<Algo_>(br, b2, -DT_(1));

    // check results
    TEST_CHECK_EQUAL_WITHIN_EPS(a1.template norm2<Algo_>(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(a2.template norm2<Algo_>(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b1.template norm2<Algo_>(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b2.template norm2<Algo_>(), DT_(0), tol);
  }

  void test_sparse_matrix_csr() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

    typedef SparseMatrixCSR<typename Algo_::MemType, DT_, IT_> MatrixType;
    IVectorType row_ptr(IT_(8));
    IVectorType col_idx(IT_(18));

    row_ptr(IT_(0), IT_(0));
    row_ptr(IT_(1), IT_(2));
    row_ptr(IT_(2), IT_(4));
    row_ptr(IT_(3), IT_(7));
    row_ptr(IT_(4), IT_(10));
    row_ptr(IT_(5), IT_(13));
    row_ptr(IT_(6), IT_(15));
    row_ptr(IT_(7), IT_(18));
    col_idx(IT_( 0), IT_(0));
    col_idx(IT_( 1), IT_(2));
    col_idx(IT_( 2), IT_(1));
    col_idx(IT_( 3), IT_(5));
    col_idx(IT_( 4), IT_(0));
    col_idx(IT_( 5), IT_(2));
    col_idx(IT_( 6), IT_(4));
    col_idx(IT_( 7), IT_(1));
    col_idx(IT_( 8), IT_(3));
    col_idx(IT_( 9), IT_(4));
    col_idx(IT_(10), IT_(0));
    col_idx(IT_(11), IT_(4));
    col_idx(IT_(12), IT_(6));
    col_idx(IT_(13), IT_(2));
    col_idx(IT_(14), IT_(5));
    col_idx(IT_(15), IT_(1));
    col_idx(IT_(16), IT_(4));
    col_idx(IT_(17), IT_(6));

    VectorType a_data(IT_(18), DT_(7));
    VectorType b_data(IT_(18), DT_(7));
    b_data(IT_( 0), DT_(1));
    b_data(IT_( 1), DT_(0));
    b_data(IT_( 4), DT_(0));
    b_data(IT_( 5), DT_(1));
    b_data(IT_( 6), DT_(0));
    b_data(IT_(15), DT_(0));
    b_data(IT_(16), DT_(0));
    b_data(IT_(17), DT_(1));

    // create a CSR matrix
    MatrixType matrix_a(Index(7), Index(7), col_idx, a_data, row_ptr);
    MatrixType matrix_b(Index(7), Index(7), col_idx, b_data, row_ptr);

    // create filter
    FilterType filter(IT_(7));
    filter.add(IT_(0), DT_(3));
    filter.add(IT_(2), DT_(5));
    filter.add(IT_(6), DT_(9));

    // apply filter onto a
    filter.template filter_mat<Algo_>(matrix_a);

    // subtract reference
    matrix_a.template axpy<Algo_>(matrix_b, matrix_a, -DT_(1));

    // check difference
    TEST_CHECK_EQUAL_WITHIN_EPS(matrix_a.template norm_frobenius<Algo_>(), DT_(0), tol);
  }

  virtual void run() const override
  {
    test_vector();
    test_sparse_matrix_csr();
  }
};

UnitFilterTest<Algo::Generic, double, Index> unit_filter_test_generic_di;
