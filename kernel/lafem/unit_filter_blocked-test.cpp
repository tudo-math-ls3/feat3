#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr_blocked.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for UnitFilterBlocked class template
 *
 * \author Jordi Paul
 *
 **/
template
<
  typename MemType_,
  typename DT_,
  typename IT_,
  int BlockSize_
  >
class UnitFilterBlockedVectorTest
  : public FullTaggedTest<MemType_, DT_, IT_>
{
  typedef Tiny::Vector<DT_, BlockSize_> ValueType;
  typedef DenseVectorBlocked<MemType_, DT_, IT_, BlockSize_> VectorType;
  typedef DenseVectorBlocked<MemType_, IT_, IT_, BlockSize_> IVectorType;
  typedef UnitFilterBlocked<MemType_, DT_, IT_, BlockSize_> FilterType;

public:
  UnitFilterBlockedVectorTest()
    : FullTaggedTest<MemType_, DT_, IT_>("UnitFilterBlockedVectorTest")
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

    const int n = 7;
    VectorType a1(n, DT_(1));
    VectorType a2(n, DT_(1));

    VectorType ar(n, DT_(1));

    VectorType b1(n, DT_(1));
    VectorType b2(n, DT_(1));

    VectorType br(n, DT_(1));

    // modify reference results
    ValueType tmp(DT_(0));

    ar(0, tmp);
    ar(2, tmp);
    ar(6, tmp);

    tmp = DT_(-17);
    tmp(BlockSize_-1) = DT_(711);
    br(0, tmp);
    tmp = DT_(-1);
    br(2, tmp);
    tmp = DT_(-7);
    br(6, tmp);

    // create filter
    FilterType filter(n);

    tmp = DT_(-17);
    tmp(BlockSize_-1) = DT_(711);
    filter.add(IT_(0), tmp);
    tmp = DT_(-1);
    filter.add(IT_(2), tmp);
    tmp = DT_(-7);
    filter.add(IT_(6), tmp);

    // check sizes
    TEST_CHECK_EQUAL(filter.size(), Index(n));
    TEST_CHECK_EQUAL(filter.used_elements(), Index(3));

    FilterType filter2;
    filter2.convert(filter);
    //TEST_CHECK_EQUAL(filter2.get_filter_vector()(2), filter.get_filter_vector()(2));
    filter2.clone(filter);
    //TEST_CHECK_EQUAL(filter2.get_filter_vector()(2), filter.get_filter_vector()(2));
    /// \todo Implement Tiny::Vector operator==

    // apply the filter
    filter.filter_def(a1);
    filter.filter_cor(a2);
    filter.filter_rhs(b1);
    filter.filter_sol(b2);

    // subtract reference results
    a1.axpy(ar, a1, -DT_(1));
    a2.axpy(ar, a2, -DT_(1));
    b1.axpy(br, b1, -DT_(1));
    b2.axpy(br, b2, -DT_(1));

    // check results
    TEST_CHECK_EQUAL_WITHIN_EPS(a1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(a2.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b2.norm2(), DT_(0), tol);
  }
};

UnitFilterBlockedVectorTest<Mem::Main, float, Index, 2> unit_filter_vector_test_generic_fi_2;
UnitFilterBlockedVectorTest<Mem::Main, double, Index, 2> unit_filter_vector_test_generic_di_2;
UnitFilterBlockedVectorTest<Mem::Main, float, Index, 3> unit_filter_vector_test_generic_fi_3;
UnitFilterBlockedVectorTest<Mem::Main, double, Index, 3> unit_filter_vector_test_generic_di_3;
UnitFilterBlockedVectorTest<Mem::Main, float, Index, 4> unit_filter_vector_test_generic_fi_4;
UnitFilterBlockedVectorTest<Mem::Main, double, Index, 4> unit_filter_vector_test_generic_di_4;
#ifdef FEAST_BACKENDS_CUDA
UnitFilterBlockedVectorTest<Mem::CUDA, float, Index, 2> unit_filter_vector_test_cuda_fi_2;
UnitFilterBlockedVectorTest<Mem::CUDA, float, Index, 3> unit_filter_vector_test_cuda_fi_3;
UnitFilterBlockedVectorTest<Mem::CUDA, float, Index, 4> unit_filter_vector_test_cuda_fi_4;
UnitFilterBlockedVectorTest<Mem::CUDA, double, Index, 2> unit_filter_vector_test_cuda_di_2;
UnitFilterBlockedVectorTest<Mem::CUDA, double, Index, 3> unit_filter_vector_test_cuda_di_3;
UnitFilterBlockedVectorTest<Mem::CUDA, double, Index, 4> unit_filter_vector_test_cuda_di_4;
#endif

/**
 * \brief Test class for testing the filter_mat functionality of the UnitFilterBlocked class template
 *
 * \author Jordi Paul
 */
template
<
  typename MemType_,
  typename DT_,
  typename IT_,
  int BlockHeight_,
  int BlockWidth_
  >
class UnitFilterBlockedMatrixTest
  : public FullTaggedTest<MemType_, DT_, IT_>
{
  typedef SparseMatrixCSRBlocked<MemType_, DT_, IT_, BlockHeight_, BlockWidth_> MatrixType;
  typedef typename MatrixType::VectorTypeL VectorTypeR;
  typedef typename MatrixType::ValueType ValueType;
  typedef DenseVector<MemType_, DT_, IT_> VectorType;
  typedef DenseVector<MemType_, IT_, IT_> IVectorType;
  typedef UnitFilterBlocked<MemType_, DT_, IT_, BlockHeight_> FilterType;

  static constexpr int BlockHeight = BlockHeight_;
  static constexpr int BlockWidth = BlockWidth_;

public:
  UnitFilterBlockedMatrixTest()
    : FullTaggedTest<MemType_, DT_, IT_>("UnitFilterBlockedMatrixTest")
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, BlockHeight> VectorValueType;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

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

    VectorType a_data(IT_(BlockHeight)*IT_(BlockWidth)*IT_(18), DT_(7));
    VectorType b_data(IT_(BlockHeight)*IT_(BlockWidth)*IT_(18), DT_(7));

    // create a CSR matrix
    MatrixType matrix_a(Index(7), Index(7), col_idx, a_data, row_ptr);
    MatrixType matrix_b(Index(7), Index(7), col_idx, b_data, row_ptr);

    // Manually correct values in the reference matrix b
    typename MatrixType::ValueType b_tmp(DT_(0));
    // Set off-diagonals to 0
    // Row 1
    matrix_b.val()[1] = b_tmp;
    // Row 3
    matrix_b.val()[4] = b_tmp;
    matrix_b.val()[6] = b_tmp;

    // Row 6
    matrix_b.val()[15] = b_tmp;
    matrix_b.val()[16] = b_tmp;

    // Set the corresponding diagonal blocks to the unit matrix (possibly non-square)
    for(int i(0); i < Math::min(BlockHeight, BlockWidth); ++i)
      b_tmp(i,i) = DT_(1);

    matrix_b.val()[0] = b_tmp;
    matrix_b.val()[5] = b_tmp;
    matrix_b.val()[17] = b_tmp;

    // create filter
    FilterType filter(IT_(7));

    VectorValueType tmp;
    tmp = DT_(3);
    filter.add(IT_(0), tmp);
    tmp = DT_(-17);
    filter.add(IT_(2), tmp);
    tmp = DT_(-1);
    filter.add(IT_(6), tmp);

    // apply filter onto a
    filter.filter_mat(matrix_a);

    // subtract reference
    matrix_a.axpy(matrix_b, matrix_a, -DT_(1));

    // check difference
    TEST_CHECK_EQUAL_WITHIN_EPS(matrix_a.norm_frobenius(), DT_(0), tol);
  }
};

UnitFilterBlockedMatrixTest<Mem::Main, float, Index, 2, 2> unit_filter_matrix_test_generic_fi_22;
UnitFilterBlockedMatrixTest<Mem::Main, float, Index, 3, 3> unit_filter_matrix_test_generic_fi_33;
UnitFilterBlockedMatrixTest<Mem::Main, float, Index, 4, 4> unit_filter_matrix_test_generic_fi_44;
UnitFilterBlockedMatrixTest<Mem::Main, float, Index, 2, 3> unit_filter_matrix_test_generic_fi_23;
UnitFilterBlockedMatrixTest<Mem::Main, float, Index, 4, 3> unit_filter_matrix_test_generic_fi_43;

UnitFilterBlockedMatrixTest<Mem::Main, double, Index, 2, 2> unit_filter_matrix_test_generic_di_22;
UnitFilterBlockedMatrixTest<Mem::Main, double, Index, 3, 3> unit_filter_matrix_test_generic_di_33;
UnitFilterBlockedMatrixTest<Mem::Main, double, Index, 4, 4> unit_filter_matrix_test_generic_di_44;
UnitFilterBlockedMatrixTest<Mem::Main, double, Index, 3, 2> unit_filter_matrix_test_generic_di_32;
UnitFilterBlockedMatrixTest<Mem::Main, double, Index, 3, 4> unit_filter_matrix_test_generic_di_34;
///TODO cuda tests?
