#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

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
template<typename Algo_, typename DataType_, typename IndexType_>
class MetaMatrixApplyTest
  : public MetaMatrixTestBase<Algo_, DataType_, IndexType_>
{
public:
  typedef Algo_ AlgoType;
  typedef DataType_ DataType;
  typedef MetaMatrixTestBase<Algo_, DataType_, IndexType_> BaseClass;

  MetaMatrixApplyTest() : BaseClass("MetaMatrixApplyTest") {}

  virtual void run() const
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
    mat_sys.template apply<AlgoType>(vec_tmp, vec_sol, vec_rhs, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.template norm2<AlgoType>(), DataType_(0), tol);

    // test t <- A*x; t <- t - b
    mat_sys.template apply<AlgoType>(vec_tmp, vec_sol);
    vec_tmp.template axpy<AlgoType>(vec_rhs, vec_tmp, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.template norm2<AlgoType>(), DataType_(0), tol);

    // generate densevectors
    DenseVector<typename AlgoType::MemType, DataType, IndexType_> vec_sol_dense, vec_rhs_dense;
    vec_sol_dense.convert(vec_sol);
    vec_rhs_dense.convert(vec_rhs);
    DenseVector<typename AlgoType::MemType, DataType, IndexType_> vec_tmp_dense(mat_sys.rows());

    // test t <- b - A*x with densevectors
    mat_sys.template apply<AlgoType>(vec_tmp_dense, vec_sol_dense, vec_rhs_dense, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp_dense.template norm2<AlgoType>(), DataType_(0), tol);

    // test t <- A*x; t <- t - b with densevectors
    mat_sys.template apply<AlgoType>(vec_tmp_dense, vec_sol_dense);
    vec_tmp_dense.template axpy<AlgoType>(vec_rhs_dense, vec_tmp_dense, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp_dense.template norm2<AlgoType>(), DataType_(0), tol);
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
    mat_sys.template apply<AlgoType>(vec_tmp, vec_sol, vec_rhs, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.template norm2<AlgoType>(), DataType_(0), tol);

    // test t <- A*x; t <- t - b
    mat_sys.template apply<AlgoType>(vec_tmp, vec_sol);
    vec_tmp.template axpy<AlgoType>(vec_rhs, vec_tmp, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.template norm2<AlgoType>(), DataType_(0), tol);

    // generate densevectors
    DenseVector<typename AlgoType::MemType, DataType, IndexType_> vec_sol_dense, vec_rhs_dense;
    vec_sol_dense.convert(vec_sol);
    vec_rhs_dense.convert(vec_rhs);
    DenseVector<typename AlgoType::MemType, DataType, IndexType_> vec_tmp_dense(mat_sys.rows());

    // test t <- b - A*x with densevectors
    mat_sys.template apply<AlgoType>(vec_tmp_dense, vec_sol_dense, vec_rhs_dense, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp_dense.template norm2<AlgoType>(), DataType_(0), tol);

    // test t <- A*x; t <- t - b with densevectors
    mat_sys.template apply<AlgoType>(vec_tmp_dense, vec_sol_dense);
    vec_tmp_dense.template axpy<AlgoType>(vec_rhs_dense, vec_tmp_dense, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp_dense.template norm2<AlgoType>(), DataType_(0), tol);
  }
};

MetaMatrixApplyTest<Algo::Generic, float, Index> meta_matrix_apply_test_generic_float;
MetaMatrixApplyTest<Algo::Generic, double, Index> meta_matrix_apply_test_generic_double;
#ifdef FEAST_BACKENDS_MKL
MetaMatrixApplyTest<Algo::MKL, float, Index> meta_matrix_apply_test_mkl_float;
MetaMatrixApplyTest<Algo::MKL, double, Index> meta_matrix_apply_test_mkl_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
MetaMatrixApplyTest<Algo::CUDA, float, Index> meta_matrix_apply_test_cuda_float;
MetaMatrixApplyTest<Algo::CUDA, double, Index> meta_matrix_apply_test_cuda_double;
#endif
