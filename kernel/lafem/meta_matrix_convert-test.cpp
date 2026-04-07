// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/util/random.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/power_diag_matrix.hpp>
#include <kernel/lafem/power_col_matrix.hpp>
#include <kernel/lafem/power_row_matrix.hpp>
#include <kernel/lafem/power_full_matrix.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief MetaMatrixConvertTest
 *
 * This class defines the MetaMatrixConvertTest which transforms a meta-matrix to a scalar matrix
 * It depends on the class "MetaMatrixTestBase"
 *
 * \author Christoph Lohmann
 */

template<typename DT_, typename IT_>
class MetaMatrixConvertTest
  : public MetaMatrixTestBase<DT_, IT_>
{
public:
  typedef DT_ DataType;
  typedef IT_ IndexType;
  typedef SparseMatrixCSR<DT_, IT_> ScalarMatrix;
  typedef MetaMatrixTestBase<DataType, IndexType> BaseClass;
  typedef typename BaseClass::SystemDiagMatrix SystemMatrix;
  typedef typename BaseClass::SystemVector SystemVector;

  explicit MetaMatrixConvertTest(PreferredBackend backend)
    : BaseClass("MetaMatrixConvertTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DataType tol = TestSystem::tol<DataType>();

    // generate a test system: A,x,b
    SystemMatrix mat_sys;
    SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_sys, vec_sol, vec_rhs);

    // test t <- A*x; t <- t - b
    vec_sol.scale(vec_sol, DT_(0));

    ScalarMatrix mat_sys_scalar;
    mat_sys_scalar.convert(mat_sys);

    DenseVector<DT_, IT_> vec_rhs_scalar(mat_sys_scalar.num_rows());
    DenseVector<DT_, IT_> vec_rhs_scalar2(mat_sys_scalar.num_rows());
    DenseVector<DT_, IT_> vec_sol_scalar(mat_sys_scalar.num_cols());
    vec_sol_scalar.format();

    for (Index i(0); i < vec_sol.size(); ++i)
    {
      vec_sol_scalar.elements_view_rw()[i] = DT_(1);
      vec_sol_scalar.copy_to(vec_sol);
      mat_sys.apply(vec_rhs, vec_sol);
      vec_sol_scalar.convert(vec_sol);
      vec_rhs_scalar2.copy(vec_rhs);
      mat_sys_scalar.apply(vec_rhs_scalar, vec_sol_scalar);

      //for (Index j(0); j < vec_rhs_scalar.size(); ++j)
      //{
      //  TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j), vec_rhs_scalar2(j), tol);
      //}
      TEST_CHECK_LESS_THAN(vec_rhs_scalar.max_rel_diff(vec_rhs_scalar2), tol);

      vec_sol_scalar.elements_view_rw()[i] = DT_(0);
    }

    auto mat_sys2 = mat_sys.clone(LAFEM::CloneMode::Weak);
    mat_sys2.format();
    mat_sys_scalar.copy_to(mat_sys2);

    ScalarMatrix mat_sys_scalar2;
    mat_sys_scalar2.convert(mat_sys2);
    TEST_CHECK_LESS_THAN(mat_sys_scalar.max_rel_diff(mat_sys_scalar2), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMatrixConvertTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

/**
 * \brief MetaBCSRMatrixConvertTest
 *
 * This class defines the MetaBCSRMatrixConvertTest which transforms a meta-matrix with a block-csr-matrix to a scalar matrix
 *
 * \author Christoph Lohmann
 */
template<typename DT_, typename IT_>
class MetaBCSRMatrixConvertTest
  : public UnitTest
{
public:
  typedef DT_ DataType;
  typedef IT_ IndexType;

  explicit MetaBCSRMatrixConvertTest(PreferredBackend backend)
    : UnitTest("MetaBCSRMatrixConvertTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DataType tol = TestSystem::tol<DataType>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    typedef SparseMatrixBCSR<DT_, IT_, 2, 3> BCSRMatrix;
    typedef PowerColMatrix<BCSRMatrix, 2> ColMatrix;
    typedef PowerRowMatrix<BCSRMatrix, 2> RowMatrix;
    typedef SaddlePointMatrix<BCSRMatrix, RowMatrix, ColMatrix> SaddleMatrix;

    typedef SparseMatrixCSR<DT_, IT_> CSRMatrix;
    typedef PowerColMatrix<CSRMatrix, 2> ColMatrix2;
    typedef PowerRowMatrix<CSRMatrix, 2> RowMatrix2;
    typedef SaddlePointMatrix<CSRMatrix, RowMatrix2, ColMatrix2> SaddleMatrix2;

    BCSRMatrix mat_a(3, 3, 4);
    {
      Memory::TypedView<IT_> row_ptr(mat_a.row_ptr_view_w());
      Memory::TypedView<IT_> col_idx(mat_a.col_idx_view_w());
      Memory::TypedView<DT_> val(mat_a.val_view_raw_w());
      row_ptr[0] = IT_(0);
      row_ptr[1] = IT_(1);
      row_ptr[2] = IT_(2);
      row_ptr[3] = IT_(4);
      col_idx[0] = IT_(0);
      col_idx[1] = IT_(1);
      col_idx[2] = IT_(1);
      col_idx[3] = IT_(2);
      for (Index i(0) ; i < mat_a.num_nzes_raw() ; ++i)
        val[i] = rng(DT_(1), DT_(10));
    }

    BCSRMatrix mat_b(2, 3, 3);
    {
      Memory::TypedView<IT_> row_ptr(mat_b.row_ptr_view_w());
      Memory::TypedView<IT_> col_idx(mat_b.col_idx_view_w());
      Memory::TypedView<DT_> val(mat_b.val_view_raw_w());
      row_ptr[0] = IT_(0);
      row_ptr[1] = IT_(1);
      row_ptr[2] = IT_(3);
      col_idx[0] = IT_(0);
      col_idx[1] = IT_(1);
      col_idx[2] = IT_(2);
      for (Index i(0) ; i < mat_b.num_nzes_raw() ; ++i)
        val[i] = rng(DT_(1), DT_(10));
    }

    // SaddlePoint<BCSR>
    SaddleMatrix sp_bcsr;
    sp_bcsr.template at<0,0>().clone(mat_a);
    mat_a.scale(mat_a, DT_(0.5));
    sp_bcsr.template at<0,1>().template at<0,0>().clone(mat_a);
    mat_a.scale(mat_a, DT_(4.0));
    sp_bcsr.template at<0,1>().template at<0,1>().clone(mat_a);
    sp_bcsr.template at<1,0>().template at<0,0>().clone(mat_b);
    mat_b.scale(mat_b, DT_(0.75));
    sp_bcsr.template at<1,0>().template at<1,0>().clone(mat_b);

    // SaddlePoint<CSR>
    SaddleMatrix2 sp_csr;
    sp_csr.convert(sp_bcsr);

    TEST_CHECK_EQUAL(sp_bcsr.num_rows_raw(), sp_csr.num_rows_raw());
    TEST_CHECK_EQUAL(sp_bcsr.num_cols_raw(), sp_csr.num_cols_raw());
    TEST_CHECK_EQUAL(sp_bcsr.num_nzes_raw(), sp_csr.num_nzes_raw());

    CSRMatrix csr_1;
    csr_1.convert(sp_bcsr);

    TEST_CHECK_EQUAL(sp_bcsr.num_rows_raw(), csr_1.num_rows_raw());
    TEST_CHECK_EQUAL(sp_bcsr.num_cols_raw(), csr_1.num_cols_raw());
    TEST_CHECK_EQUAL(sp_bcsr.num_nzes_raw(), csr_1.num_nzes_raw());

    CSRMatrix csr_2;
    csr_2.convert(sp_csr);

    TEST_CHECK_EQUAL(sp_csr.num_rows_raw(), csr_2.num_rows_raw());
    TEST_CHECK_EQUAL(sp_csr.num_cols_raw(), csr_2.num_cols_raw());
    TEST_CHECK_EQUAL(sp_csr.num_nzes_raw(), csr_2.num_nzes_raw());

    // compare the two CSR matrices derived from SaddlePoint<BCSR> and SaddlePoint<CSR>
    TEST_CHECK(csr_2.same_layout(csr_1));
    TEST_CHECK_LESS_THAN(csr_2.max_rel_diff(csr_1), tol);

    // copy back to SaddlePoint<BCSR> and verify
    SaddleMatrix sp_bcsr2 = sp_bcsr.clone(LAFEM::CloneMode::Layout);
    sp_bcsr2.format();
    csr_2.copy_to(sp_bcsr2);

    TEST_CHECK(sp_bcsr2.same_layout(sp_bcsr));
    TEST_CHECK_LESS_THAN(sp_bcsr2.max_rel_diff(sp_bcsr), tol);

    // copy back to SaddlePoint<CSR> and verify
    SaddleMatrix2 sp_csr2 = sp_csr.clone(LAFEM::CloneMode::Layout);
    sp_csr2.format();
    csr_1.copy_to(sp_csr2);

    TEST_CHECK(sp_csr2.same_layout(sp_csr));
    TEST_CHECK_LESS_THAN(sp_csr2.max_rel_diff(sp_csr), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaBCSRMatrixConvertTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
