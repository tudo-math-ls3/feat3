#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>

#include <kernel/lafem/bicgstab.hpp>
#include <kernel/lafem/preconditioner.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template <SparsePreconType Type_>
struct Precon;

template <>
struct Precon<SparsePreconType::pt_none>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & /*sys*/, const Index /*opt*/)
  {
    Preconditioner<MT_, VT_> * t = new NonePreconditioner<MT_, VT_> ();
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_jacobi>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & sys, const Index /*opt*/)
  {
    Preconditioner<MT_, VT_> * t = new JacobiPreconditioner<MT_, VT_> (sys, typename VT_::DataType(1));
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_gauss_seidel>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & sys, const Index /*opt*/)
  {
    Preconditioner<MT_, VT_> * t = new GaussSeidelPreconditioner<MT_, VT_> (sys, typename VT_::DataType(1));
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_polynomial>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & sys, const Index opt)
  {
    Preconditioner<MT_, VT_> * t = new PolynomialPreconditioner<MT_, VT_> (sys, 20, opt == 0);
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_ilu>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & sys, const Index opt)
  {
    Preconditioner<MT_, VT_> * t = new ILUPreconditioner<MT_, VT_> (sys, opt);
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_sor>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & sys, const Index /*opt*/)
  {
    Preconditioner<MT_, VT_> * t = new SORPreconditioner<MT_, VT_> (sys);
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_ssor>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & sys, const Index /*opt*/)
  {
    Preconditioner<MT_, VT_> * t = new SSORPreconditioner<MT_, VT_> (sys);
    return t;
  }
};

template <>
struct Precon<SparsePreconType::pt_spai>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & sys, const Index opt)
  {
    Preconditioner<MT_, VT_> * t;

    const bool transpose(opt%2 == 0);
    const bool start_layout((opt / 2)%2 == 0);
    const Index max_iter(opt / 4);

    if (start_layout)
    {
      t = new SPAIPreconditioner<MT_, VT_> (sys, 2, max_iter, 1e-2, 100, 1e-3, 1e-3, transpose);
    }
    else
    {
      t = new SPAIPreconditioner<MT_, VT_> (sys, sys.layout(), max_iter, 1e-2, 10, 1e-3, 1e-3, transpose);
    }
    return t;
  }
};


/**
 * \brief BiCGStabSaddlePointTest
 *
 * This class defines the BiCGStabTest with different preconditioners for a Saddle-Point matrix.
 * At the moment first we convert the Saddle-Point matrix to a scalar matrix
 * and solve the scalar linear system with this matrix.
 * Notice: The Saddle-Point matrix has a zero-block at the bottom-right,
 * so we can't use all preconditioners (zeros on the main diagonal).
 *
 * \author Christoph Lohmann
 */
template<typename PSF_, SparsePreconType PType_,  typename MT_, typename VT_>
class BiCGStabSaddlePointTest
  : public MetaMatrixTestBase<typename MT_::MemType, typename MT_::DataType, Index>
{
private:
  const Index _opt;

public:
  typedef typename VT_::DataType   DT_;
  typedef typename VT_::MemType    Mem_;
  typedef MetaMatrixTestBase<typename MT_::MemType, DT_, Index> BaseClass;
  typedef typename BaseClass::SystemDiagMatrix SystemMatrix;
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

    MT_ mat_sys_scalar;
    mat_sys_scalar.convert(mat_sys);

    Index size(mat_sys_scalar.rows());
    VT_ x(size, DT_(1));
    VT_ ref;
    ref.convert(vec_sol);
    VT_ ref_local;
    ref_local.convert(ref);
    VT_ b(size);
    mat_sys_scalar.apply(b, ref);

    Preconditioner<MT_, VT_> * precond(Precon<PType_>::template get<MT_, VT_>(mat_sys_scalar, _opt));
    BiCGStab::value(x, mat_sys_scalar, b, *precond, 1000, 1e-12);
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
                        SparseMatrixCSR<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_none_double("none");

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_none,
                        SparseMatrixCOO<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_none_double("none");

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_none,
                        SparseMatrixELL<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_none_double("none");


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixCSR<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_spai_double_2("spai", 2);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixCOO<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_spai_double_2("spai", 2);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixELL<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_spai_double_2("spai", 2);


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixCSR<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_spai_double_3("spai", 3);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixCOO<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_spai_double_3("spai", 3);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixELL<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_spai_double_3("spai", 3);


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixCSR<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_spai_double_82("spai", 82);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixCOO<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_spai_double_82("spai", 82);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixELL<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_spai_double_82("spai", 82);


BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixCSR<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_csr_spai_double_83("spai", 83);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixCOO<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_coo_spai_double_83("spai", 83);

BiCGStabSaddlePointTest<PointstarFactoryFE<double>,
                        SparsePreconType::pt_spai,
                        SparseMatrixELL<Mem::Main, double>,
                        DenseVector<Mem::Main, double> >
bicgstab_saddle_point_test_cpu_ell_spai_double_83("spai", 83);
