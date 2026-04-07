// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
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
  int block_size_
  >
class UnitFilterBlockedVectorTest
  : public UnitTest
{
  typedef Tiny::Vector<DT_, block_size_> ValueType;
  typedef DenseVectorBlocked<DT_, IT_, block_size_> VectorType;
  typedef DenseVectorBlocked<IT_, IT_, block_size_> IVectorType;
  typedef UnitFilterBlocked<DT_, IT_, block_size_> FilterType;

public:
  explicit UnitFilterBlockedVectorTest(PreferredBackend backend)
    : UnitTest("UnitFilterBlockedVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
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
      Memory::TypedView<ValueType> var(ar.elements_view_rw());
      Memory::TypedView<ValueType> vbr(br.elements_view_rw());
      var[0] = DT_(0);
      var[2] = DT_(0);
      var[6] = DT_(0);
      vbr[0] = DT_(-17);
      vbr[0][block_size_-1] = DT_(711);
      vbr[2] = DT_(-1);
      vbr[6] = DT_(-7);
    }

    // create filter
    FilterType filter(n, 3);
    {
      Memory::TypedView<IT_> idx = filter.indices_view_w();
      Memory::TypedView<ValueType> val = filter.elements_view_w();
      idx[0] = IT_(0);
      idx[1] = IT_(2);
      idx[2] = IT_(6);
      val[0] = DT_(-17);
      val[0][block_size_-1] = DT_(711);
      val[1] = DT_(-1);
      val[2] = DT_(-7);
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

SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, Index, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, Index, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, Index, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, Index, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, Index, 4, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, Index, 4, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, std::uint64_t, 4, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, std::uint64_t, 4, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, __float128, Index, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, __float128, Index, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, __float128, Index, 4, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, Half, Index, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, Half, Index, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, Half, Index, 4, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, Index, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, Index, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, float, Index, 4, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, Index, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, Index, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(UnitFilterBlockedVectorTest, double, Index, 4, PreferredBackend::cuda);
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
  int block_height_,
  int block_width_
  >
class UnitFilterBlockedMatrixTest
  : public UnitTest
{
  typedef SparseMatrixBCSR<DT_, IT_, block_height_, block_width_> MatrixType;
  typedef typename MatrixType::VectorTypeL VectorTypeR;
  typedef typename MatrixType::ValueType ValueType;
  typedef DenseVector<DT_, IT_> VectorType;
  typedef DenseVector<IT_, IT_> IVectorType;
  typedef UnitFilterBlocked<DT_, IT_, block_height_> FilterType;

  static constexpr int block_height = block_height_;
  static constexpr int block_width = block_width_;

public:
  explicit UnitFilterBlockedMatrixTest(PreferredBackend backend)
    : UnitTest("UnitFilterBlockedMatrixTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_height> VectorValueType;
    typedef Tiny::Matrix<DT_, block_height, block_width> MatrixValueType;

    const DT_ tol = TestSystem::tol<DT_>();

    const Index nrows = 7;
    const Index nnzes = 18u;

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

    VectorType a_data(Index(block_height)*Index(block_width)*nnzes, DT_(7));
    VectorType b_data(Index(block_height)*Index(block_width)*nnzes, DT_(7));
    VectorType c_data(Index(block_height)*Index(block_width)*nnzes, DT_(7));
    {
      Memory::TypedView<DT_> a_view(a_data.elements_view_w());
      Memory::TypedView<DT_> b_view(b_data.elements_view_w());
      Memory::TypedView<DT_> c_view(c_data.elements_view_w());
      for(Index i(0); i < Index(block_height)*Index(block_width)*nnzes; ++i)
      {
        a_view[i] = DT_(i+1u);
        b_view[i] = DT_(i+1u);
        c_view[i] = DT_(i+1u);
      }
    }

    // create a CSR matrix
    VectorType a_data2 = a_data.clone();
    MatrixType matrix_a1(nrows, nrows, row_ptr, col_idx, a_data);
    MatrixType matrix_a2(nrows, nrows, row_ptr, col_idx, a_data2);
    MatrixType matrix_b(nrows, nrows, row_ptr, col_idx, b_data);
    MatrixType matrix_c(nrows, nrows, row_ptr, col_idx, c_data);

    // Manually correct values in the reference matrix b
    {
      Memory::TypedView<MatrixValueType> vb = matrix_b.val_view_rw();

      // Row 1
      vb[0].set_identity();
      vb[1].format();
      // Row 3
      vb[4].format();
      vb[5].set_identity();
      vb[6].format();
      // Row 6
      vb[15].format();
      vb[16].format();
      vb[17].set_identity();
    }

    // create filter
    FilterType filter(Index(7), Index(3));
    {
      Memory::TypedView<IT_> idx = filter.indices_view_w();
      Memory::TypedView<VectorValueType> val = filter.elements_view_w();
      idx[0] = IT_(0);
      idx[1] = IT_(2);
      idx[2] = IT_(6);
      val[0] = DT_(3);
      val[1] = DT_(-17);
      val[2] = DT_(-1);
    }

    // apply filter onto a
    filter.filter_mat(matrix_a1);
    TEST_CHECK_LESS_THAN(matrix_a1.max_rel_diff(matrix_b), tol);

    // Manually correct values in the reference matrix c
    {
      const Memory::TypedView<MatrixValueType> va = matrix_a2.val_view_r();
      Memory::TypedView<MatrixValueType> vc = matrix_c.val_view_rw();

      // Row 1
      vc[0] = DT_(3) * va[0];
      vc[1] = DT_(3) * va[1];
      // Row 3
      vc[4] = DT_(-17) * va[4];
      vc[5] = DT_(-17) * va[5];
      vc[6] = DT_(-17) * va[6];
      // Row 6
      vc[15] = DT_(-1) * va[15];
      vc[16] = DT_(-1) * va[16];
      vc[17] = DT_(-1) * va[17];
    }
    auto matrix_a3 = matrix_a2.clone(LAFEM::CloneMode::Weak);

    // apply weak filter onto a2
    filter.filter_weak_matrix_rows(matrix_a2, matrix_a3);

    // subtract reference
    TEST_CHECK_LESS_THAN(matrix_a2.max_rel_diff(matrix_c), tol);
  }
};

SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 2, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 3, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 4, 4, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 2, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 4, 3, PreferredBackend::generic);

SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 2, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 3, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 4, 4, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 3, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 3, 4, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, std::uint64_t, 2, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, std::uint64_t, 3, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, std::uint64_t, 4, 4, PreferredBackend::mkl);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, std::uint64_t, 2, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, std::uint64_t, 4, 3, PreferredBackend::mkl);

SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, std::uint64_t, 2, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, std::uint64_t, 3, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, std::uint64_t, 4, 4, PreferredBackend::mkl);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, std::uint64_t, 2, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, std::uint64_t, 4, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, __float128, Index, 2, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, __float128, Index, 3, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, __float128, Index, 4, 4, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, __float128, Index, 2, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, __float128, Index, 4, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, Half, Index, 2, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, Half, Index, 3, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, Half, Index, 4, 4, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, Half, Index, 2, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, Half, Index, 4, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 2, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 3, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 4, 4, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 2, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, float, Index, 4, 3, PreferredBackend::cuda);

SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 2, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 3, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 4, 4, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 2, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(UnitFilterBlockedMatrixTest, double, Index, 4, 3, PreferredBackend::cuda);
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
  static constexpr int block_size_ = 3;
  typedef Tiny::Vector<DT_, block_size_> ValueType;
  typedef DenseVectorBlocked<DT_, IT_, block_size_> VectorType;
  typedef DenseVectorBlocked<IT_, IT_, block_size_> IVectorType;
  typedef UnitFilterBlocked<DT_, IT_, block_size_> FilterType;

public:
  explicit UnitFilterBlockedNansTest(PreferredBackend backend)
    : UnitTest("UnitFilterBlockedNansTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ nan = Math::nan<DT_>();

    const Index n = 7;

    VectorType x(n), ref1(n), ref2(n);
    {
      Memory::TypedView<ValueType> vx = x.elements_view_w();
      Memory::TypedView<ValueType> vr1 = ref1.elements_view_w();
      Memory::TypedView<ValueType> vr2 = ref2.elements_view_w();
      for(Index i(0); i < n; ++i)
      {
        for(int j(0); j < block_size_; ++j)
        {
          vr1[i][j] = vr2[i][j] = vx[i][j] = DT_(10u*i + Index(j));
        }
      }
      for(int j(0); j < block_size_; ++j)
      {
        vr1[2][j] = vr1[5][j] = DT_(j+1);
        vr2[2][j] = vr2[5][j] = DT_(0);
      }
      vr1[5u][1] = vr2[5u][1] = DT_(51);
    }

    // create filter with some NaNs in it
    FilterType filter(n, 2, true);
    {
      Memory::TypedView<IT_> idx = filter.indices_view_w();
      Memory::TypedView<ValueType> val = filter.elements_view_w();
      idx[0] = IT_(2);
      idx[1] = IT_(5);
      val[0][0] = DT_(1);
      val[0][1] = DT_(2);
      val[0][2] = DT_(3);
      val[1][0] = DT_(1);
      val[1][1] = nan;
      val[1][2] = DT_(3);
    }


    // filter vector
    filter.filter_rhs(x);

    // check filtered values
    TEST_CHECK_LESS_THAN(x.max_rel_diff(ref1), tol);

    // filter vector for defect
    filter.filter_def(x);
    TEST_CHECK_LESS_THAN(x.max_rel_diff(ref2), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedNansTest, float, Index, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedNansTest, double, Index, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedNansTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedNansTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedNansTest, __float128, Index, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedNansTest, Half, Index, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedNansTest, float, Index, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedNansTest, double, Index, PreferredBackend::cuda);
#endif
