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
 *  - SaddlePointMatrix
 *
 * \author Peter Zajac
 */
template<typename Algo_, typename DataType_>
class MetaMatrixApplyTest
  : public MetaMatrixTestBase<Algo_, DataType_>
{
public:
  typedef Algo_ AlgoType;
  typedef DataType_ DataType;
  typedef MetaMatrixTestBase<Algo_, DataType_> BaseClass;
  typedef typename BaseClass::SystemMatrix SystemMatrix;
  typedef typename BaseClass::SystemVector SystemVector;

  MetaMatrixApplyTest() : BaseClass("MetaMatrixApplyTest") {}

  virtual void run() const
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.7)));

    // generate a test system: A,x,b
    SystemMatrix mat_sys;
    SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_sys, vec_sol, vec_rhs);

    // test t <- b - A*x
    SystemVector vec_tmp(vec_rhs.clone());
    mat_sys.template apply<AlgoType>(vec_tmp, vec_sol, vec_rhs, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.template norm2<AlgoType>(), DataType_(0), tol);

    // test t <- A*x; t <- t - b
    mat_sys.template apply<AlgoType>(vec_tmp, vec_sol);
    vec_tmp.template axpy<AlgoType>(vec_rhs, vec_tmp, -DataType_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_tmp.template norm2<AlgoType>(), DataType_(0), tol);
  }
};

MetaMatrixApplyTest<Algo::Generic, float> meta_matrix_apply_test_generic_float;
MetaMatrixApplyTest<Algo::Generic, double> meta_matrix_apply_test_generic_double;
#ifdef FEAST_BACKENDS_MKL
MetaMatrixApplyTest<Algo::MKL, float> meta_matrix_apply_test_mkl_float;
MetaMatrixApplyTest<Algo::MKL, double> meta_matrix_apply_test_mkl_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
MetaMatrixApplyTest<Algo::CUDA, float> meta_matrix_apply_test_cuda_float;
MetaMatrixApplyTest<Algo::CUDA, double> meta_matrix_apply_test_cuda_double;
#endif
