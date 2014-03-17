#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/block_to_scalar.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>

#include <kernel/lafem/bicgstab.hpp>
#include <kernel/lafem/preconditioner.hpp>

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
    : BaseClass("block_to_scalar_test: " + MT_::name()) {}

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


/********************************************************************************/


template <SparsePreconType Type_>
struct Precon;

template <>
struct Precon<SparsePreconType::pt_none>
{
  template <typename Algo_, typename MT_, typename VT_>
  static Preconditioner<Algo_, MT_, VT_> * get(const MT_ & /*sys*/, const Index /*opt*/)
  {
    Preconditioner<Algo_, MT_, VT_> * t = new NonePreconditioner<Algo_, MT_, VT_> ();
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_jacobi>
{
  template <typename Algo_, typename MT_, typename VT_>
  static Preconditioner<Algo_, MT_, VT_> * get(const MT_ & sys, const Index /*opt*/)
  {
    Preconditioner<Algo_, MT_, VT_> * t = new JacobiPreconditioner<Algo_, MT_, VT_> (sys, typename VT_::DataType(1));
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_gauss_seidel>
{
  template <typename Algo_, typename MT_, typename VT_>
  static Preconditioner<Algo_, MT_, VT_> * get(const MT_ & sys, const Index /*opt*/)
  {
    Preconditioner<Algo_, MT_, VT_> * t = new GaussSeidelPreconditioner<Algo_, MT_, VT_> (sys, typename VT_::DataType(1));
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_polynomial>
{
  template <typename Algo_, typename MT_, typename VT_>
  static Preconditioner<Algo_, MT_, VT_> * get(const MT_ & sys, const Index opt)
  {
    Preconditioner<Algo_, MT_, VT_> * t = new PolynomialPreconditioner<Algo_, MT_, VT_> (sys, 20, opt == 0);
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_ilu>
{
  template <typename Algo_, typename MT_, typename VT_>
  static Preconditioner<Algo_, MT_, VT_> * get(const MT_ & sys, const Index opt)
  {
    Preconditioner<Algo_, MT_, VT_> * t = new ILUPreconditioner<Algo_, MT_, VT_> (sys, opt);
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_sor>
{
  template <typename Algo_, typename MT_, typename VT_>
  static Preconditioner<Algo_, MT_, VT_> * get(const MT_ & sys, const Index /*opt*/)
  {
    Preconditioner<Algo_, MT_, VT_> * t = new SORPreconditioner<Algo_, MT_, VT_> (sys);
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_ssor>
{
  template <typename Algo_, typename MT_, typename VT_>
  static Preconditioner<Algo_, MT_, VT_> * get(const MT_ & sys, const Index /*opt*/)
  {
    Preconditioner<Algo_, MT_, VT_> * t = new SSORPreconditioner<Algo_, MT_, VT_> (sys);
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_spai>
{
  template <typename Algo_, typename MT_, typename VT_>
  static Preconditioner<Algo_, MT_, VT_> * get(const MT_ & sys, const Index opt)
  {
    Preconditioner<Algo_, MT_, VT_> * t;

    const bool transpose(opt%2 == 0);
    const bool start_layout((opt / 2)%2 == 0);
    const Index max_iter(opt / 4);

    if (start_layout)
    {
      t = new SPAIPreconditioner<Algo_, MT_, VT_> (sys, 2, max_iter, 1e-2, 100, 1e-3, 1e-3, transpose);
    }
    else
    {
      t = new SPAIPreconditioner<Algo_, MT_, VT_> (sys, sys.layout(), max_iter, 1e-2, 10, 1e-3, 1e-3, transpose);
    }
    return t;
  }
};

/**
 * \brief BiCGStabTest
 *
 * This class defines the BiCGStabTest with different preconditioners and matrix-types
 *
 * \author Christoph Lohmann
 */
template<typename PSF_, SparsePreconType PType_, typename Algo_, typename MT_, typename VT_>
class BiCGStabSaddlePointTest
  : public MetaMatrixTestBase<Algo_, typename MT_::DataType>
{
private:
  const Index _opt;

public:
  typedef typename VT_::DataType   DT_;
  typedef typename VT_::MemType    Mem_;
  typedef MetaMatrixTestBase<Algo_, DT_> BaseClass;
  typedef typename BaseClass::SystemMatrix SystemMatrix;
  typedef typename BaseClass::SystemVector SystemVector;

  BiCGStabSaddlePointTest(String pname, Index opt = 0)
    : BaseClass("bicgstab_saddle_point_test: "
                + MT_::name()
                + " "
                + pname
                + " opt = "
                + stringify(opt)), _opt(opt)
  {
  }

  virtual void run() const
  {
    SystemMatrix mat_sys;
    SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_sys, vec_sol, vec_rhs);

    MT_ mat_sys_scalar (MatBlockToScalar<Algo_>::template value<MT_::LayoutType>(mat_sys));

    Index size(mat_sys_scalar.rows());
    VT_ x(size, DT_(1));
    VT_ ref(VecBlockToScalar<Algo_>::value(vec_sol));
    VT_ ref_local(ref);
    VT_ b(size);
    mat_sys_scalar.template apply<Algo_>(b, ref);

    Preconditioner<Algo_, MT_, VT_> * precond(Precon<PType_>::template get<Algo_, MT_, VT_>(mat_sys_scalar, _opt));
    BiCGStab<Algo_>::value(x, mat_sys_scalar, b, *precond, 1000, 1e-12);
    delete precond;

    // check, if the result is correct
    for (Index i(0) ; i < size ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(x(i), ref_local(i), 1e-8);
    }

  }

};


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_none_double("none");

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_none_double("none");

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_none_double("none");


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_spai_double_2("spai", 2);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_spai_double_2("spai", 2);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_spai_double_2("spai", 2);


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_spai_double_3("spai", 3);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_spai_double_3("spai", 3);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_spai_double_3("spai", 3);


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_spai_double_82("spai", 82);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_spai_double_82("spai", 82);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_spai_double_82("spai", 82);


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_spai_double_83("spai", 83);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_spai_double_83("spai", 83);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_spai_double_83("spai", 83);
