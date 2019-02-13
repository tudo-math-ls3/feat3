// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Meta-Matrix apply test class
 *
 * \test The 'apply' operations of the following class templates:
 *  - SparseMatrixCSR
 *  - SparseMatrixCOO
 *  - SparseMatrixELL
 *  - PowerColMatrix
 *  - PowerRowMatrix
 *  - PowerDiagMatrix
 *  - PowerFullMatrix
 *  - SaddlePointMatrix
 *
 * \author Peter Zajac
 */
template<
  typename DataType_,
  typename IndexType_>
class MetaMatrixApplyTest
  : public MetaMatrixTestBase<DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaMatrixTestBase<DataType_, IndexType_> BaseClass;

   MetaMatrixApplyTest(PreferredBackend backend) :
    BaseClass("MetaMatrixApplyTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~MetaMatrixApplyTest()
  {
  }

  virtual void run() const override
  {
    test_diag();
    test_full();
  }

  void test_diag() const
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.7)));

    // generate a test system: A,x,b
    typename BaseClass::SystemDiagMatrix mat_sys;
    typename BaseClass::SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_sys, vec_sol, vec_rhs);

    // test t <- b - A*x
    typename BaseClass::SystemVector vec_tmp(mat_sys.create_vector_l());
    mat_sys.apply(vec_tmp, vec_sol, vec_rhs, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.norm2(), DataType_(0), tol);

    // test t <- A*x; t <- t - b
    mat_sys.apply(vec_tmp, vec_sol);
    vec_tmp.axpy(vec_rhs, vec_tmp, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.norm2(), DataType_(0), tol);

    // generate densevectors
    DenseVector<DataType, IndexType_> vec_sol_dense, vec_rhs_dense;
    vec_sol_dense.convert(vec_sol);
    vec_rhs_dense.convert(vec_rhs);
    DenseVector<DataType, IndexType_> vec_tmp_dense(mat_sys.rows());

    // test t <- b - A*x with densevectors
    mat_sys.apply(vec_tmp_dense, vec_sol_dense, vec_rhs_dense, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp_dense.norm2(), DataType_(0), tol);

    // test t <- A*x; t <- t - b with densevectors
    mat_sys.apply(vec_tmp_dense, vec_sol_dense);
    vec_tmp_dense.axpy(vec_rhs_dense, vec_tmp_dense, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp_dense.norm2(), DataType_(0), tol);
  }

  void test_full() const
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.7)));

    // generate a test system: A,x,b
    typename BaseClass::SystemFullMatrix mat_sys;
    typename BaseClass::SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_sys, vec_sol, vec_rhs);

    // test t <- b - A*x
    typename BaseClass::SystemVector vec_tmp(mat_sys.create_vector_l());
    mat_sys.apply(vec_tmp, vec_sol, vec_rhs, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.norm2(), DataType_(0), tol);

    // test t <- A*x; t <- t - b
    mat_sys.apply(vec_tmp, vec_sol);
    vec_tmp.axpy(vec_rhs, vec_tmp, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.norm2(), DataType_(0), tol);

    // generate densevectors
    DenseVector<DataType, IndexType_> vec_sol_dense, vec_rhs_dense;
    vec_sol_dense.convert(vec_sol);
    vec_rhs_dense.convert(vec_rhs);
    DenseVector<DataType, IndexType_> vec_tmp_dense(mat_sys.rows());

    // test t <- b - A*x with densevectors
    mat_sys.apply(vec_tmp_dense, vec_sol_dense, vec_rhs_dense, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp_dense.norm2(), DataType_(0), tol);

    // test t <- A*x; t <- t - b with densevectors
    mat_sys.apply(vec_tmp_dense, vec_sol_dense);
    vec_tmp_dense.axpy(vec_rhs_dense, vec_tmp_dense, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp_dense.norm2(), DataType_(0), tol);
  }
};

MetaMatrixApplyTest<float, unsigned long> meta_matrix_apply_test_generic_float_ulong(PreferredBackend::generic);
MetaMatrixApplyTest<double, unsigned long> meta_matrix_apply_test_generic_double_ulong(PreferredBackend::generic);
MetaMatrixApplyTest<float, unsigned int> meta_matrix_apply_test_generic_float_uint(PreferredBackend::generic);
MetaMatrixApplyTest<double, unsigned int> meta_matrix_apply_test_generic_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MetaMatrixApplyTest<float, unsigned long> mkl_meta_matrix_apply_test_float_ulong(PreferredBackend::mkl);
MetaMatrixApplyTest<double, unsigned long> mkl_meta_matrix_apply_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MetaMatrixApplyTest<__float128, unsigned long> meta_matrix_apply_test_float128_ulong(PreferredBackend::generic);
MetaMatrixApplyTest<__float128, unsigned int> meta_matrix_apply_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MetaMatrixApplyTest<Half, unsigned int> meta_matrix_apply_test_half_uint(PreferredBackend::generic);
MetaMatrixApplyTest<Half, unsigned long> meta_matrix_apply_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MetaMatrixApplyTest<float, unsigned long> cuda_meta_matrix_apply_test_float_ulong(PreferredBackend::cuda);
MetaMatrixApplyTest<double, unsigned long> cuda_meta_matrix_apply_test_double_ulong(PreferredBackend::cuda);
MetaMatrixApplyTest<float, unsigned int> cuda_meta_matrix_apply_test_float_uint(PreferredBackend::cuda);
MetaMatrixApplyTest<double, unsigned int> cuda_meta_matrix_apply_test_double_uint(PreferredBackend::cuda);
#endif
