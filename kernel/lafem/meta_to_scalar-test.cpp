#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;


/**
 * \brief MetaToScalarTest
 *
 * This class defines the MetaToScalarTest which transforms a meta-matrix to a scalar matrix
 * It depends on the class "MetaMatrixTestBase"
 *
 * \author Christoph Lohmann
 */
template<typename Algo_, typename MT_>
class MetaToScalarTest
  : public MetaMatrixTestBase<Algo_, typename MT_::DataType, typename MT_::IndexType>
{
public:
  typedef Algo_ AlgoType;
  typedef typename MT_::DataType DataType;
  typedef typename MT_::IndexType IndexType;
  typedef MetaMatrixTestBase<Algo_, DataType, IndexType> BaseClass;
  typedef typename BaseClass::SystemDiagMatrix SystemMatrix;
  typedef typename BaseClass::SystemVector SystemVector;

  typedef typename SystemVector::MemType Mem_;
  typedef DataType DT_;
  typedef IndexType IT_;

  MetaToScalarTest()
    : BaseClass("meta_to_scalar_test: " + MT_::name()) {}

  virtual void run() const
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.7)));

    // generate a test system: A,x,b
    SystemMatrix mat_sys;
    SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_sys, vec_sol, vec_rhs);

    // test t <- A*x; t <- t - b
    vec_sol.template scale<Algo_>(vec_sol, DT_(0));

    MT_ mat_sys_scalar;
    mat_sys_scalar.convert(mat_sys);

    DenseVector<Mem_, DT_, IT_> vec_rhs_scalar(mat_sys_scalar.rows());
    DenseVector<Mem_, DT_, IT_> vec_rhs_scalar2(mat_sys_scalar.rows());

    for (Index i(0); i < vec_sol.size(); ++i)
    {
      vec_sol(i, DT_(1));
      mat_sys.template apply<AlgoType>(vec_rhs, vec_sol);
      DenseVector<Mem_, DT_, IT_> vec_sol_scalar;
      vec_sol_scalar.convert(vec_sol);
      vec_rhs_scalar2.copy(vec_rhs);
      mat_sys_scalar.template apply<AlgoType>(vec_rhs_scalar, vec_sol_scalar);

      for (Index j(0); j < vec_rhs_scalar.size(); ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j), vec_rhs_scalar2(j), tol);
      }

      vec_sol(i, DT_(0));
    }

    /**
     * check VecMetaToScalar and "VecScalarToMeta"
     */
    vec_rhs.template scale<Algo_>(vec_rhs, DT_(0));

    for (Index j(0); j < vec_rhs_scalar.size(); ++j)
    {
      vec_rhs_scalar(j, DT_(j));
    }

    vec_rhs_scalar.copy_inv(vec_rhs);

    for (Index j(0); j < vec_rhs_scalar.size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j), vec_rhs(j), tol);
    }

    vec_rhs_scalar.template scale<Algo_>(vec_rhs_scalar, DT_(0));

    vec_rhs_scalar.copy(vec_rhs);

    for (Index j(0); j < vec_rhs_scalar.size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs(j), vec_rhs_scalar(j), tol);
    }
  }
};

MetaToScalarTest<Algo::Generic, SparseMatrixCOO<Mem::Main, float, unsigned int> > cpu_meta_matrix_to_coo_test_generic_float_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixCOO<Mem::Main, double, unsigned int> > cpu_meta_matrix_to_coo_test_generic_double_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixCOO<Mem::Main, float, unsigned long> > cpu_meta_matrix_to_coo_test_generic_float_ulong;
MetaToScalarTest<Algo::Generic, SparseMatrixCOO<Mem::Main, double, unsigned long> > cpu_meta_matrix_to_coo_test_generic_double_ulong;

MetaToScalarTest<Algo::Generic, SparseMatrixCSR<Mem::Main, float, unsigned int> > cpu_meta_matrix_to_csr_test_generic_float_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixCSR<Mem::Main, double, unsigned int> > cpu_meta_matrix_to_csr_test_generic_double_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixCSR<Mem::Main, float, unsigned long> > cpu_meta_matrix_to_csr_test_generic_float_ulong;
MetaToScalarTest<Algo::Generic, SparseMatrixCSR<Mem::Main, double, unsigned long> > cpu_meta_matrix_to_csr_test_generic_double_ulong;

MetaToScalarTest<Algo::Generic, SparseMatrixELL<Mem::Main, float, unsigned int> > cpu_meta_matrix_to_ell_test_generic_float_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixELL<Mem::Main, double, unsigned int> > cpu_meta_matrix_to_ell_test_generic_double_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixELL<Mem::Main, float, unsigned long> > cpu_meta_matrix_to_ell_test_generic_float_ulong;
MetaToScalarTest<Algo::Generic, SparseMatrixELL<Mem::Main, double, unsigned long> > cpu_meta_matrix_to_ell_test_generic_double_ulong;

// TODO: Alle fehlenden Tests hinzufuegen
#ifdef FEAST_HAVE_QUADMATH
MetaToScalarTest<Algo::Generic, SparseMatrixCOO<Mem::Main, __float128, unsigned int> > cpu_meta_matrix_to_coo_test_generic_float128_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixCOO<Mem::Main, __float128, unsigned long> > cpu_meta_matrix_to_coo_test_generic_float128_ulong;

MetaToScalarTest<Algo::Generic, SparseMatrixCSR<Mem::Main, __float128, unsigned int> > cpu_meta_matrix_to_csr_test_generic_float128_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixCSR<Mem::Main, __float128, unsigned long> > cpu_meta_matrix_to_csr_test_generic_float128_ulong;

MetaToScalarTest<Algo::Generic, SparseMatrixELL<Mem::Main, __float128, unsigned int> > cpu_meta_matrix_to_ell_test_generic_float128_uint;
MetaToScalarTest<Algo::Generic, SparseMatrixELL<Mem::Main, __float128, unsigned long> > cpu_meta_matrix_to_ell_test_generic_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
MetaToScalarTest<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, float, unsigned int> > cuda_meta_matrix_to_csr_test_generic_float_uint;
MetaToScalarTest<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, unsigned int> > cuda_meta_matrix_to_csr_test_generic_double_uint;
MetaToScalarTest<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, float, unsigned long> > cuda_meta_matrix_to_csr_test_generic_float_ulong;
MetaToScalarTest<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, unsigned long> > cuda_meta_matrix_to_csr_test_generic_double_ulong;

MetaToScalarTest<Algo::CUDA, SparseMatrixELL<Mem::CUDA, float, unsigned int> > cuda_meta_matrix_to_ell_test_generic_float_uint;
MetaToScalarTest<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, unsigned int> > cuda_meta_matrix_to_ell_test_generic_double_uint;
MetaToScalarTest<Algo::CUDA, SparseMatrixELL<Mem::CUDA, float, unsigned long> > cuda_meta_matrix_to_ell_test_generic_float_ulong;
MetaToScalarTest<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, unsigned long> > cuda_meta_matrix_to_ell_test_generic_double_ulong;
#endif
