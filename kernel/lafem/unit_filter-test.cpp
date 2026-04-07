// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for UnitFilter class template
 *
 * \author Peter Zajac
 */
template<
  typename DT_,
  typename IT_>
class UnitFilterVectorTest
  : public UnitTest
{
  typedef DenseVector<DT_, IT_> VectorType;
  typedef DenseVector<IT_, IT_> IVectorType;
  typedef UnitFilter<DT_, IT_> FilterType;

public:
  explicit UnitFilterVectorTest(PreferredBackend backend)
    : UnitTest("UnitFilterVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    const Index n = 7;
    VectorType a1(n, DT_(1));
    VectorType a2(n, DT_(1));
    VectorType ar(n, DT_(1));
    VectorType b1(n, DT_(1));
    VectorType b2(n, DT_(1));
    VectorType br(n, DT_(1));

    // modify reference results
    {
      Memory::TypedView<DT_> var(ar.elements_view_rw());
      Memory::TypedView<DT_> vbr(br.elements_view_rw());
      var[0] = DT_(0);
      var[2] = DT_(0);
      var[6] = DT_(0);
      vbr[0] = DT_(3);
      vbr[2] = DT_(5);
      vbr[6] = DT_(9);
    }

    // create filter
    FilterType filter(n, 3);
    {
      Memory::TypedView<IT_> idx = filter.indices_view_w();
      Memory::TypedView<DT_> val = filter.elements_view_w();
      idx[0] = IT_(0);
      idx[1] = IT_(2);
      idx[2] = IT_(6);
      val[0] = DT_(3);
      val[1] = DT_(5);
      val[2] = DT_(9);
    }

    // check sizes
    TEST_CHECK_EQUAL(filter.size(), n);
    TEST_CHECK_EQUAL(filter.num_nzes(), Index(3));

    FilterType filter2;
    filter2.convert(filter);
    TEST_CHECK(filter2.get_filter_vector().same_layout(filter.get_filter_vector()));
    TEST_CHECK_LESS_THAN(filter2.get_filter_vector().max_rel_diff(filter.get_filter_vector()), tol);

    FilterType filter3;
    filter3.clone(filter);
    TEST_CHECK(filter3.get_filter_vector().same_layout(filter.get_filter_vector()));
    TEST_CHECK_LESS_THAN(filter3.get_filter_vector().max_rel_diff(filter.get_filter_vector()), tol);

    // apply the filter
    filter.filter_def(a1);
    filter.filter_cor(a2);
    filter.filter_rhs(b1);
    filter.filter_sol(b2);

    // subtract reference results
    TEST_CHECK_LESS_THAN(a1.max_rel_diff(ar), tol);
    TEST_CHECK_LESS_THAN(a2.max_rel_diff(ar), tol);
    TEST_CHECK_LESS_THAN(b1.max_rel_diff(br), tol);
    TEST_CHECK_LESS_THAN(b2.max_rel_diff(br), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterVectorTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

/**
 * \brief Test class for UnitFilter class template
 *
 * \author Peter Zajac
 */
template<
  typename DT_,
  typename IT_>
class UnitFilterMatrixTest
  : public UnitTest
{
  typedef DenseVector<DT_, IT_> VectorType;
  typedef DenseVector<IT_, IT_> IVectorType;
  typedef UnitFilter<DT_, IT_> FilterType;
public:
  explicit UnitFilterMatrixTest(PreferredBackend backend)
    : UnitTest("UnitFilterMatrixTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    const Index nrows = 7;
    const Index nnzes = 18u;

    typedef SparseMatrixCSR<DT_, IT_> MatrixType;
    IVectorType row_ptr(nrows + Index(1));
    IVectorType col_idx(nnzes);

    {
      Memory::TypedView<IT_> vrow_ptr(row_ptr.elements_view_w());
      vrow_ptr[0] = IT_(0);
      vrow_ptr[1] = IT_(2);
      vrow_ptr[2] = IT_(4);
      vrow_ptr[3] = IT_(7);
      vrow_ptr[4] = IT_(10);
      vrow_ptr[5] = IT_(13);
      vrow_ptr[6] = IT_(15);
      vrow_ptr[7] = IT_(18);
    }

    {
      Memory::TypedView<IT_> vcol_idx(col_idx.elements_view_w());
      vcol_idx[ 0] = IT_(0);
      vcol_idx[ 1] = IT_(2);
      vcol_idx[ 2] = IT_(1);
      vcol_idx[ 3] = IT_(5);
      vcol_idx[ 4] = IT_(0);
      vcol_idx[ 5] = IT_(2);
      vcol_idx[ 6] = IT_(4);
      vcol_idx[ 7] = IT_(1);
      vcol_idx[ 8] = IT_(3);
      vcol_idx[ 9] = IT_(4);
      vcol_idx[10] = IT_(0);
      vcol_idx[11] = IT_(4);
      vcol_idx[12] = IT_(6);
      vcol_idx[13] = IT_(2);
      vcol_idx[14] = IT_(5);
      vcol_idx[15] = IT_(1);
      vcol_idx[16] = IT_(4);
      vcol_idx[17] = IT_(6);
    }

    VectorType a_data(nnzes, DT_(7));
    VectorType b_data(nnzes, DT_(7));
    VectorType c_data(nnzes, DT_(7));
    {
      Memory::TypedView<DT_> a_view(a_data.elements_view_w());
      Memory::TypedView<DT_> b_view(b_data.elements_view_w());
      Memory::TypedView<DT_> c_view(c_data.elements_view_w());
      for(Index i(0); i < nnzes; ++i)
      {
        a_view[i] = DT_(i+1u);
        b_view[i] = DT_(i+1u);
        c_view[i] = DT_(i+1u);
      }
      b_view[ 0] = DT_(1);
      b_view[ 1] = DT_(0);
      b_view[ 4] = DT_(0);
      b_view[ 5] = DT_(1);
      b_view[ 6] = DT_(0);
      b_view[15] = DT_(0);
      b_view[16] = DT_(0);
      b_view[17] = DT_(1);
      c_view[ 0] = DT_(3*1);
      c_view[ 1] = DT_(3*2);
      c_view[ 4] = DT_(5*5);
      c_view[ 5] = DT_(5*6);
      c_view[ 6] = DT_(5*7);
      c_view[15] = DT_(9*16);
      c_view[16] = DT_(9*17);
      c_view[17] = DT_(9*18);
    }

    // create a pair of CSR matrices
    VectorType a_data2 = a_data.clone();
    MatrixType matrix_a1(nrows, nrows, row_ptr, col_idx, a_data);
    MatrixType matrix_a2(nrows, nrows, row_ptr, col_idx, a_data2);
    MatrixType matrix_b(nrows, nrows, row_ptr, col_idx, b_data);
    MatrixType matrix_c(nrows, nrows, row_ptr, col_idx, c_data);

    // create filter
    FilterType filter(nrows, 3);
    {
      Memory::TypedView<IT_> idx = filter.indices_view_w();
      Memory::TypedView<DT_> val = filter.elements_view_w();
      idx[0] = IT_(0);
      idx[1] = IT_(2);
      idx[2] = IT_(6);
      val[0] = DT_(3);
      val[1] = DT_(5);
      val[2] = DT_(9);
    }

    // apply filter onto a1
    filter.filter_mat(matrix_a1);

    // subtract reference
    TEST_CHECK_LESS_THAN(matrix_a1.max_rel_diff(matrix_b), tol);
    //matrix_a1.axpy(matrix_b, -DT_(1));
    //TEST_CHECK_EQUAL_WITHIN_EPS(matrix_a1.norm_frobenius(), DT_(0), tol);

    // apply weak filter onto a2
    filter.filter_weak_matrix_rows(matrix_a1, matrix_a2);

    // subtract reference
    TEST_CHECK_LESS_THAN(matrix_a1.max_rel_diff(matrix_c), tol);
    //matrix_a2.axpy(matrix_c, -DT_(1));
    //TEST_CHECK_EQUAL_WITHIN_EPS(matrix_a2.norm_frobenius(), DT_(0), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterMatrixTest, double, std::uint32_t, PreferredBackend::cuda);
#endif
