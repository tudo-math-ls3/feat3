// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for UnitFilterBlocked class template
 *
 * \author Jordi Paul
 *
 **/
template
<
  typename DT_,
  typename IT_,
  int BlockSize_
  >
class UnitFilterBlockedVectorTest
  : public UnitTest
{
  typedef Tiny::Vector<DT_, BlockSize_> ValueType;
  typedef DenseVectorBlocked<DT_, IT_, BlockSize_> VectorType;
  typedef DenseVectorBlocked<IT_, IT_, BlockSize_> IVectorType;
  typedef UnitFilterBlocked<DT_, IT_, BlockSize_> FilterType;

public:
  UnitFilterBlockedVectorTest(PreferredBackend backend)
    : UnitTest("UnitFilterBlockedVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~UnitFilterBlockedVectorTest()
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

UnitFilterBlockedVectorTest<float, Index, 2> unit_filter_blocked_vector_test_generic_fi_2(PreferredBackend::generic);
UnitFilterBlockedVectorTest<double, Index, 2> unit_filter_blocked_vector_test_generic_di_2(PreferredBackend::generic);
UnitFilterBlockedVectorTest<float, Index, 3> unit_filter_blocked_vector_test_generic_fi_3(PreferredBackend::generic);
UnitFilterBlockedVectorTest<double, Index, 3> unit_filter_blocked_vector_test_generic_di_3(PreferredBackend::generic);
UnitFilterBlockedVectorTest<float, Index, 4> unit_filter_blocked_vector_test_generic_fi_4(PreferredBackend::generic);
UnitFilterBlockedVectorTest<double, Index, 4> unit_filter_blocked_vector_test_generic_di_4(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
UnitFilterBlockedVectorTest <float, std::uint64_t, 2> mkl_unit_filter_blocked_vector_test_float_uint64_2(PreferredBackend::mkl);
UnitFilterBlockedVectorTest <double, std::uint64_t, 2> mkl_unit_filter_blocked_vector_test_double_uint64_2(PreferredBackend::mkl);
UnitFilterBlockedVectorTest <float, std::uint64_t, 3> mkl_unit_filter_blocked_vector_test_float_uint64_3(PreferredBackend::mkl);
UnitFilterBlockedVectorTest <double, std::uint64_t, 3> mkl_unit_filter_blocked_vector_test_double_uint64_3(PreferredBackend::mkl);
UnitFilterBlockedVectorTest <float, std::uint64_t, 4> mkl_unit_filter_blocked_vector_test_float_uint64_4(PreferredBackend::mkl);
UnitFilterBlockedVectorTest <double, std::uint64_t, 4> mkl_unit_filter_blocked_vector_test_double_uint64_4(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
UnitFilterBlockedVectorTest<__float128, Index, 2> unit_filter_blocked_vector_test_float128_index_2(PreferredBackend::generic);
UnitFilterBlockedVectorTest<__float128, Index, 3> unit_filter_blocked_vector_test_float128_index_3(PreferredBackend::generic);
UnitFilterBlockedVectorTest<__float128, Index, 4> unit_filter_blocked_vector_test_float128_index_4(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
UnitFilterBlockedVectorTest<Half, Index, 2> unit_filter_blocked_vector_test_half_index_2(PreferredBackend::generic);
UnitFilterBlockedVectorTest<Half, Index, 3> unit_filter_blocked_vector_test_half_index_3(PreferredBackend::generic);
UnitFilterBlockedVectorTest<Half, Index, 4> unit_filter_blocked_vector_test_half_index_4(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
UnitFilterBlockedVectorTest<float, Index, 2> unit_filter_blocked_vector_test_cuda_fi_2(PreferredBackend::cuda);
UnitFilterBlockedVectorTest<float, Index, 3> unit_filter_blocked_vector_test_cuda_fi_3(PreferredBackend::cuda);
UnitFilterBlockedVectorTest<float, Index, 4> unit_filter_blocked_vector_test_cuda_fi_4(PreferredBackend::cuda);
UnitFilterBlockedVectorTest<double, Index, 2> unit_filter_blocked_vector_test_cuda_di_2(PreferredBackend::cuda);
UnitFilterBlockedVectorTest<double, Index, 3> unit_filter_blocked_vector_test_cuda_di_3(PreferredBackend::cuda);
UnitFilterBlockedVectorTest<double, Index, 4> unit_filter_blocked_vector_test_cuda_di_4(PreferredBackend::cuda);
#endif

/**
 * \brief Test class for testing the filter_mat functionality of the UnitFilterBlocked class template
 *
 * \author Jordi Paul
 */
template
<
  typename DT_,
  typename IT_,
  int BlockHeight_,
  int BlockWidth_
  >
class UnitFilterBlockedMatrixTest
  : public UnitTest
{
  typedef SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_> MatrixType;
  typedef typename MatrixType::VectorTypeL VectorTypeR;
  typedef typename MatrixType::ValueType ValueType;
  typedef DenseVector<DT_, IT_> VectorType;
  typedef DenseVector<IT_, IT_> IVectorType;
  typedef UnitFilterBlocked<DT_, IT_, BlockHeight_> FilterType;

  static constexpr int BlockHeight = BlockHeight_;
  static constexpr int BlockWidth = BlockWidth_;

public:
  UnitFilterBlockedMatrixTest(PreferredBackend backend)
    : UnitTest("UnitFilterBlockedMatrixTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~UnitFilterBlockedMatrixTest()
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
    VectorType c_data(IT_(BlockHeight)*IT_(BlockWidth)*IT_(18), DT_(7));
    for(IT_ i(0); i < IT_(18); ++i)
    {
      a_data(i, DT_(i+1u));
      b_data(i, DT_(i+1u));
      c_data(i, DT_(i+1u));
    }

    // create a CSR matrix
    VectorType a_data2 = a_data.clone();
    MatrixType matrix_a1(Index(7), Index(7), col_idx, a_data, row_ptr);
    MatrixType matrix_a2(Index(7), Index(7), col_idx, a_data2, row_ptr);
    MatrixType matrix_b(Index(7), Index(7), col_idx, b_data, row_ptr);
    MatrixType matrix_c(Index(7), Index(7), col_idx, c_data, row_ptr);

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
    filter.filter_mat(matrix_a1);

    // subtract reference
    matrix_a1.axpy(matrix_b, matrix_a1, -DT_(1));

    // check difference
    TEST_CHECK_EQUAL_WITHIN_EPS(matrix_a1.norm_frobenius(), DT_(0), tol);


    // Manually correct values in the reference matrix c
    // Row 1
    matrix_c.val()[0] = DT_(3) * matrix_a2.val()[0];
    matrix_c.val()[1] = DT_(3) * matrix_a2.val()[1];
    // Row 3
    matrix_c.val()[4] = DT_(-17) * matrix_a2.val()[4];
    matrix_c.val()[5] = DT_(-17) * matrix_a2.val()[5];
    matrix_c.val()[6] = DT_(-17) * matrix_a2.val()[6];
    // Row 6
    matrix_c.val()[15] = DT_(-1) * matrix_a2.val()[15];
    matrix_c.val()[16] = DT_(-1) * matrix_a2.val()[16];
    matrix_c.val()[17] = DT_(-1) * matrix_a2.val()[17];

    // apply weak filter onto a2
    filter.filter_weak_matrix_rows(matrix_a2, matrix_a2);

    // subtract reference
    matrix_a2.axpy(matrix_c, matrix_a2, -DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(matrix_a2.norm_frobenius(), DT_(0), tol);
  }
};

UnitFilterBlockedMatrixTest<float, Index, 2, 2> unit_filter_matrix_test_generic_fi_22(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<float, Index, 3, 3> unit_filter_matrix_test_generic_fi_33(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<float, Index, 4, 4> unit_filter_matrix_test_generic_fi_44(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<float, Index, 2, 3> unit_filter_matrix_test_generic_fi_23(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<float, Index, 4, 3> unit_filter_matrix_test_generic_fi_43(PreferredBackend::generic);

UnitFilterBlockedMatrixTest<double, Index, 2, 2> unit_filter_matrix_test_generic_di_22(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<double, Index, 3, 3> unit_filter_matrix_test_generic_di_33(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<double, Index, 4, 4> unit_filter_matrix_test_generic_di_44(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<double, Index, 3, 2> unit_filter_matrix_test_generic_di_32(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<double, Index, 3, 4> unit_filter_matrix_test_generic_di_34(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
UnitFilterBlockedMatrixTest <float, std::uint64_t, 2, 2> mkl_unit_filter_blocked_matrix_test_float_uint64_22(PreferredBackend::mkl);
UnitFilterBlockedMatrixTest <float, std::uint64_t, 3, 3> mkl_unit_filter_blocked_matrix_test_float_uint64_33(PreferredBackend::mkl);
UnitFilterBlockedMatrixTest <float, std::uint64_t, 4, 4> mkl_unit_filter_blocked_matrix_test_float_uint64_44(PreferredBackend::mkl);
UnitFilterBlockedMatrixTest <float, std::uint64_t, 2, 3> mkl_unit_filter_blocked_matrix_test_float_uint64_23(PreferredBackend::mkl);
UnitFilterBlockedMatrixTest <float, std::uint64_t, 4, 3> mkl_unit_filter_blocked_matrix_test_float_uint64_43(PreferredBackend::mkl);

UnitFilterBlockedMatrixTest <double, std::uint64_t, 2, 2> mkl_unit_filter_blocked_matrix_test_double_uint64_22(PreferredBackend::mkl);
UnitFilterBlockedMatrixTest <double, std::uint64_t, 3, 3> mkl_unit_filter_blocked_matrix_test_double_uint64_33(PreferredBackend::mkl);
UnitFilterBlockedMatrixTest <double, std::uint64_t, 4, 4> mkl_unit_filter_blocked_matrix_test_double_uint64_44(PreferredBackend::mkl);
UnitFilterBlockedMatrixTest <double, std::uint64_t, 2, 3> mkl_unit_filter_blocked_matrix_test_double_uint64_23(PreferredBackend::mkl);
UnitFilterBlockedMatrixTest <double, std::uint64_t, 4, 3> mkl_unit_filter_blocked_matrix_test_double_uint64_43(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
UnitFilterBlockedMatrixTest<__float128, Index, 2, 2> unit_filter_blocked_matrix_test_float128_index_22(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<__float128, Index, 3, 3> unit_filter_blocked_matrix_test_float128_index_33(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<__float128, Index, 4, 4> unit_filter_blocked_matrix_test_float128_index_44(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<__float128, Index, 2, 3> unit_filter_blocked_matrix_test_float128_index_23(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<__float128, Index, 4, 3> unit_filter_blocked_matrix_test_float128_index_43(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
UnitFilterBlockedMatrixTest<Half, Index, 2, 2> unit_filter_blocked_matrix_test_half_index_22(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<Half, Index, 3, 3> unit_filter_blocked_matrix_test_half_index_33(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<Half, Index, 4, 4> unit_filter_blocked_matrix_test_half_index_44(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<Half, Index, 2, 3> unit_filter_blocked_matrix_test_half_index_23(PreferredBackend::generic);
UnitFilterBlockedMatrixTest<Half, Index, 4, 3> unit_filter_blocked_matrix_test_half_index_43(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
UnitFilterBlockedMatrixTest<float, Index, 2, 2> unit_filter_blocked_matrix_test_cuda_fi_22(PreferredBackend::cuda);
UnitFilterBlockedMatrixTest<float, Index, 3, 3> unit_filter_blocked_matrix_test_cuda_fi_33(PreferredBackend::cuda);
UnitFilterBlockedMatrixTest<float, Index, 4, 4> unit_filter_blocked_matrix_test_cuda_fi_44(PreferredBackend::cuda);
UnitFilterBlockedMatrixTest<float, Index, 2, 3> unit_filter_blocked_matrix_test_cuda_fi_23(PreferredBackend::cuda);
UnitFilterBlockedMatrixTest<float, Index, 4, 3> unit_filter_blocked_matrix_test_cuda_fi_43(PreferredBackend::cuda);

UnitFilterBlockedMatrixTest<double, Index, 2, 2> unit_filter_blocked_matrix_test_cuda_di_22(PreferredBackend::cuda);
UnitFilterBlockedMatrixTest<double, Index, 3, 3> unit_filter_blocked_matrix_test_cuda_di_33(PreferredBackend::cuda);
UnitFilterBlockedMatrixTest<double, Index, 4, 4> unit_filter_blocked_matrix_test_cuda_di_44(PreferredBackend::cuda);
UnitFilterBlockedMatrixTest<double, Index, 2, 3> unit_filter_blocked_matrix_test_cuda_di_23(PreferredBackend::cuda);
UnitFilterBlockedMatrixTest<double, Index, 4, 3> unit_filter_blocked_matrix_test_cuda_di_43(PreferredBackend::cuda);
#endif


/**
 * \brief Test class for UnitFilterBlocked ignore NaNs functionality
 *
 * \author Peter Zajac
 */
template<typename DT_, typename IT_>
class UnitFilterBlockedNansTest
  : public UnitTest
{
  static constexpr int BlockSize_ = 3;
  typedef Tiny::Vector<DT_, BlockSize_> ValueType;
  typedef DenseVectorBlocked<DT_, IT_, BlockSize_> VectorType;
  typedef DenseVectorBlocked<IT_, IT_, BlockSize_> IVectorType;
  typedef UnitFilterBlocked<DT_, IT_, BlockSize_> FilterType;

public:
  UnitFilterBlockedNansTest(PreferredBackend backend)
    : UnitTest("UnitFilterBlockedNansTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~UnitFilterBlockedNansTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));
    const DT_ nan = Math::nan<DT_>();

    const Index n = 7;

    VectorType x(n);
    ValueType* vx = x.elements();
    for(Index i(0); i < n; ++i)
    {
      for(int j(0); j < BlockSize_; ++j)
        vx[i][j] = DT_(10u*i + Index(j));
    }

    // create filter with some NaNs in it
    FilterType filter(n, true);

    ValueType va{DT_(1.0), DT_(2.0), DT_(3.0)};
    ValueType vb(va);
    vb[1] = nan;

    filter.add(2u, va);
    filter.add(5u, vb);

    // filter vector
    filter.filter_rhs(x);

    // check filtered values
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[2u][0], DT_(1.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[2u][1], DT_(2.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[2u][2], DT_(3.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[5u][0], DT_(1.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[5u][1], DT_(51.0), tol); // ignored because of NaN
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[5u][2], DT_(3.0), tol);

    // filter vector for defect
    filter.filter_def(x);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[2u][0], DT_(0.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[2u][1], DT_(0.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[2u][2], DT_(0.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[5u][0], DT_(0.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[5u][1], DT_(51.0), tol); // ignored because of NaN
    TEST_CHECK_EQUAL_WITHIN_EPS(vx[5u][2], DT_(0.0), tol);
  }
};

UnitFilterBlockedNansTest<float, Index> unit_filter_blocked_nans_test_generic_fi(PreferredBackend::generic);
UnitFilterBlockedNansTest<double, Index> unit_filter_blocked_nans_test_generic_di(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
UnitFilterBlockedNansTest <float, std::uint64_t> mkl_unit_filter_blocked_nans_test_float_uint64(PreferredBackend::mkl);
UnitFilterBlockedNansTest <double, std::uint64_t> mkl_unit_filter_blocked_nans_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
UnitFilterBlockedNansTest<__float128, Index> unit_filter_blocked_nans_test_float128_index(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
UnitFilterBlockedNansTest<Half, Index> unit_filter_blocked_nans_test_half_index(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
UnitFilterBlockedNansTest<float, Index> unit_filter_blocked_nans_test_cuda_fi(PreferredBackend::cuda);
UnitFilterBlockedNansTest<double, Index> unit_filter_blocked_nans_test_cuda_di(PreferredBackend::cuda);
#endif
