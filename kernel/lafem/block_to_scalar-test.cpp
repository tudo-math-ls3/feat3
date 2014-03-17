#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/block_to_scalar.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;


/**
 * \brief BlockToScalarTest
 *
 * This class defines the BlockToScalarTest which transforms a block-matrix to a scalar matrix
 * It depends on the class "MetaMatrixTestBase"
 *
 * \author Christoph Lohmann
 */
template<typename Algo_, typename MT_>
class BlockToScalarTest
  : public MetaMatrixTestBase<Algo_, typename MT_::DataType>
{
public:
  typedef Algo_ AlgoType;
  typedef typename MT_::DataType DataType;
  typedef MetaMatrixTestBase<Algo_, DataType> BaseClass;
  typedef typename BaseClass::SystemMatrix SystemMatrix;
  typedef typename BaseClass::SystemVector SystemVector;

  typedef typename SystemVector::MemType Mem_;
  typedef DataType DT_;

  BlockToScalarTest()
    : BaseClass("block_to_scalar_test") {}

  virtual void run() const
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.7)));

    // generate a test system: A,x,b
    SystemMatrix mat_sys;
    SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_sys, vec_sol, vec_rhs);

    // test t <- A*x; t <- t - b
    vec_sol.template scale<Algo_>(vec_sol, DT_(0));

    MT_ mat_sys_scalar (MatBlockToScalar<Algo_>::template value<MT_::LayoutType>(mat_sys));

    DenseVector<Mem_, DT_> vec_rhs_scalar(mat_sys_scalar.rows());
    DenseVector<Mem_, DT_> vec_rhs_scalar2(mat_sys_scalar.rows());

    for (Index i(0); i < vec_sol.size(); ++i)
    {
      vec_sol(i, DT_(1));
      mat_sys.template apply<AlgoType>(vec_rhs, vec_sol);
      VecBlockToScalar<Algo_>::copy(vec_rhs_scalar2, vec_rhs);

      mat_sys_scalar.template apply<AlgoType>(vec_rhs_scalar, VecBlockToScalar<Algo_>::value(vec_sol));

      for (Index j(0); j < vec_rhs_scalar.size(); ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j), vec_rhs_scalar2(j), tol);
      }

      vec_sol(i, DT_(0));
    }
  }
};

BlockToScalarTest<Algo::Generic, SparseMatrixCOO<Mem::Main, double> > meta_matrix_to_coo_test_generic_double;
BlockToScalarTest<Algo::Generic, SparseMatrixCSR<Mem::Main, double> > meta_matrix_to_csr_test_generic_double;
BlockToScalarTest<Algo::Generic, SparseMatrixELL<Mem::Main, double> > meta_matrix_to_ell_test_generic_double;
