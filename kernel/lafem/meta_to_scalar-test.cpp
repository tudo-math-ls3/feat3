#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/meta_to_scalar.hpp>
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
  : public MetaMatrixTestBase<Algo_, typename MT_::DataType, Index>
{
public:
  typedef Algo_ AlgoType;
  typedef typename MT_::DataType DataType;
  typedef MetaMatrixTestBase<Algo_, DataType, Index> BaseClass;
  typedef typename BaseClass::SystemDiagMatrix SystemMatrix;
  typedef typename BaseClass::SystemVector SystemVector;

  typedef typename SystemVector::MemType Mem_;
  typedef DataType DT_;

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

    MT_ mat_sys_scalar (MatMetaToScalar<Algo_>::template value<MT_::layout_id>(mat_sys));

    DenseVector<Mem_, DT_> vec_rhs_scalar(mat_sys_scalar.rows());
    DenseVector<Mem_, DT_> vec_rhs_scalar2(mat_sys_scalar.rows());

    for (Index i(0); i < vec_sol.size(); ++i)
    {
      vec_sol(i, DT_(1));
      mat_sys.template apply<AlgoType>(vec_rhs, vec_sol);
      VecMetaToScalar<Algo_>::copy(vec_rhs_scalar2, vec_rhs);

      mat_sys_scalar.template apply<AlgoType>(vec_rhs_scalar, VecMetaToScalar<Algo_>::value(vec_sol));

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

    VecMetaToScalar<Algo_>::copy(vec_rhs, vec_rhs_scalar);

    for (Index j(0); j < vec_rhs_scalar.size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j), vec_rhs(j), tol);
    }

    vec_rhs_scalar.template scale<Algo_>(vec_rhs_scalar, DT_(0));

    VecMetaToScalar<Algo_>::copy(vec_rhs_scalar, vec_rhs);

    for (Index j(0); j < vec_rhs_scalar.size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs(j), vec_rhs_scalar(j), tol);
    }
  }
};

MetaToScalarTest<Algo::Generic, SparseMatrixCOO<Mem::Main, double> > meta_matrix_to_coo_test_generic_double;
MetaToScalarTest<Algo::Generic, SparseMatrixCSR<Mem::Main, double> > meta_matrix_to_csr_test_generic_double;
MetaToScalarTest<Algo::Generic, SparseMatrixELL<Mem::Main, double> > meta_matrix_to_ell_test_generic_double;
