#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/lafem/bicgstab.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

template <SparsePreconType Type_>
struct Precon;

/**
 * none preconditioner
 */
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

/**
 * Jacobi preconditioner
 */
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

/**
 * Gauss-Seidel preconditioner
 */
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

/**
 * ILU preconditioner
 */
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

/**
 * SOR preconditioner
 */
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

/**
 * SSOR preconditioner
 */
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

/**
 * SPAI preconditioner
 */
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
      t = new SPAIPreconditioner<MT_, VT_> (sys, 2, max_iter, 1e-2, 10, 1e-3, 1e-3, transpose);
    }
    else
    {
      t = new SPAIPreconditioner<MT_, VT_> (sys, sys.layout(), max_iter, 1e-2, 10, 1e-3, 1e-3, transpose);
    }
    return t;
  }
};

/**
 * Polynomial preconditioner for Mem::Main
 */
template <>
struct Precon<SparsePreconType::pt_polynomial>
{
  template <typename MT_, typename VT_>
  static Preconditioner<MT_, VT_> * get(const MT_ & sys, const Index opt)
  {
    const Index iter_poly(2);
    Preconditioner<MT_, VT_> * t;

    if (opt%10 == 0)
    {
      auto * temp = Precon<SparsePreconType::pt_none>::template get<MT_, VT_>(sys,(opt - opt%10) / 10);
      t = new PolynomialPreconditioner<MT_, VT_> (sys, iter_poly, temp);
    }
    else if (opt%10 == 1)
    {
      auto * temp = Precon<SparsePreconType::pt_jacobi>::template get<MT_, VT_>(sys,(opt - opt%10) / 10);
      t = new PolynomialPreconditioner<MT_, VT_> (sys, iter_poly, temp);
    }
    else if (opt%10 == 2)
    {
      auto * temp = Precon<SparsePreconType::pt_gauss_seidel>::template get<MT_, VT_>(sys,(opt - opt%10) / 10);
      t = new PolynomialPreconditioner<MT_, VT_> (sys, iter_poly, temp);
    }
    else if (opt%10 == 3)
    {
      auto * temp = Precon<SparsePreconType::pt_sor>::template get<MT_, VT_>(sys,(opt - opt%10) / 10);
      t = new PolynomialPreconditioner<MT_, VT_> (sys, iter_poly, temp);
    }
    else if (opt%10 == 4)
    {
      auto * temp = Precon<SparsePreconType::pt_ssor>::template get<MT_, VT_>(sys,(opt - opt%10) / 10);
      t = new PolynomialPreconditioner<MT_, VT_> (sys, iter_poly, temp);
    }
    else if (opt%10 == 5)
    {
      auto * temp = Precon<SparsePreconType::pt_ilu>::template get<MT_, VT_>(sys,(opt - opt%10) / 10);
      t = new PolynomialPreconditioner<MT_, VT_> (sys, iter_poly, temp);
    }
    else if (opt%10 == 6)
    {
      auto * temp = Precon<SparsePreconType::pt_spai>::template get<MT_, VT_>(sys,(opt - opt%10) / 10);
      t = new PolynomialPreconditioner<MT_, VT_> (sys, iter_poly, temp);
    }
    else
    {
      throw InternalError(__func__, __FILE__, __LINE__, "Parameter " + stringify(opt) + " not allowed!");
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
template<typename PSF_, SparsePreconType PType_, typename MT_>
class BiCGStabTest
  : public FullTaggedTest<typename MT_::MemType,
                          typename MT_::DataType,
                          typename MT_::IndexType>
{
private:
  const Index _opt;

public:
  typedef typename MT_::DataType   DT_;
  typedef typename MT_::IndexType  IT_;
  typedef typename MT_::MemType    Mem_;
  typedef DenseVector<Mem_, DT_, IT_> VT_;

  BiCGStabTest(String pname, Index opt = 0)
    : FullTaggedTest<Mem_, DT_, IT_> ("bicgstab_test: "
                                             + MT_::name()
                                             + " "
                                             + pname
                                             + " opt = "
                                             + stringify(opt)), _opt(opt)
  {
  }

  virtual void run() const override
  {
    PSF_ factory(13);
    MT_ sys;
    sys.convert(factory.matrix_csr());

    Index size(sys.rows());
    VT_ x(size, DT_(1));
    VT_ ref;
    ref.convert(factory.vector_q2_bubble());
    VT_ ref_local(ref.size());
    VT_ b(size);
    sys.apply(b, ref);

    Preconditioner<MT_, VT_> * precond(Precon<PType_>::template get<MT_, VT_>(sys, _opt));
    BiCGStab::value(x, sys, b, *precond, 1000, 1e-12);
    delete precond;

    ref_local.copy(ref);
    // check, if the result is correct
    for (Index i(0) ; i < size ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(x(i), ref_local(i), 1e-10);
    }
  }

};

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_none_double_uint("none");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_none_double_uint("none");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_none_double_uint("none");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_jac_double_uint("jacobi");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_jac_double_uint("jacobi");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_jac_double_uint("jacobi");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_gs_double_uint("gauss-seidel");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_gs_double_uint("gauss-seidel");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_gs_double_uint("gauss-seidel");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_poly_double_uint_0("polynmial", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_poly_double_uint_0("polynmial", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_poly_double_uint_0("polynmial", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_poly_double_uint_1("polynmial", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_poly_double_uint_1("polynmial", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_poly_double_uint_1("polynmial", 1);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_poly_double_uint_2("polynmial", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_poly_double_uint_2("polynmial", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_poly_double_uint_2("polynmial", 2);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_poly_double_uint_3("polynmial", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_poly_double_uint_3("polynmial", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_poly_double_uint_3("polynmial", 3);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_poly_double_uint_4("polynmial", 4);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_poly_double_uint_4("polynmial", 4);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_poly_double_uint_4("polynmial", 4);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_poly_double_uint_5("polynmial", 5);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_poly_double_uint_5("polynmial", 5);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_poly_double_uint_5("polynmial", 5);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_poly_double_uint_6("polynmial", 6);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_poly_double_uint_6("polynmial", 6);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_poly_double_uint_6("polynmial", 6);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_ilu_double_uint_0("ilu", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_ilu_double_uint_0("ilu", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_ilu_double_uint_0("ilu", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_ilu_double_uint_10("ilu", 10);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_ilu_double_uint_10("ilu", 10);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_ilu_double_uint_10("ilu", 10);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_sor_double_uint("sor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_sor_double_uint("sor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_sor_double_uint("sor");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_ssor_double_uint("ssor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_ssor_double_uint("ssor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_ssor_double_uint("ssor");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_spai_double_uint_0("spai", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_spai_double_uint_0("spai", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_spai_double_uint_0("spai", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_spai_double_uint_1("spai", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_spai_double_uint_1("spai", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_spai_double_uint_1("spai", 1);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_spai_double_uint_2("spai", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_spai_double_uint_2("spai", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_spai_double_uint_2("spai", 2);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_spai_double_uint_3("spai", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_spai_double_uint_3("spai", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_spai_double_uint_3("spai", 3);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_spai_double_uint_80("spai", 80);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_spai_double_uint_80("spai", 80);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_spai_double_uint_80("spai", 80);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_spai_double_uint_81("spai", 81);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_spai_double_uint_81("spai", 81);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_spai_double_uint_81("spai", 81);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_spai_double_uint_82("spai", 82);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_spai_double_uint_82("spai", 82);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_spai_double_uint_82("spai", 82);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_csr_spai_double_uint_83("spai", 83);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_coo_spai_double_uint_83("spai", 83);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned int> >
bicgstab_test_cpu_ell_spai_double_uint_83("spai", 83);



BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_none_double_ulong("none");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_none_double_ulong("none");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_none_double_ulong("none");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_jac_double_ulong("jacobi");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_jac_double_ulong("jacobi");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_jac_double_ulong("jacobi");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_gs_double_ulong("gauss-seidel");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_gs_double_ulong("gauss-seidel");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_gs_double_ulong("gauss-seidel");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_poly_double_ulong_0("polynmial", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_poly_double_ulong_0("polynmial", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_poly_double_ulong_0("polynmial", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_poly_double_ulong_1("polynmial", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_poly_double_ulong_1("polynmial", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_poly_double_ulong_1("polynmial", 1);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_poly_double_ulong_2("polynmial", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_poly_double_ulong_2("polynmial", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_poly_double_ulong_2("polynmial", 2);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_poly_double_ulong_3("polynmial", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_poly_double_ulong_3("polynmial", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_poly_double_ulong_3("polynmial", 3);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_poly_double_ulong_4("polynmial", 4);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_poly_double_ulong_4("polynmial", 4);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_poly_double_ulong_4("polynmial", 4);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_poly_double_ulong_5("polynmial", 5);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_poly_double_ulong_5("polynmial", 5);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_poly_double_ulong_5("polynmial", 5);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_poly_double_ulong_6("polynmial", 6);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_poly_double_ulong_6("polynmial", 6);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_poly_double_ulong_6("polynmial", 6);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_ilu_double_ulong_0("ilu", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_ilu_double_ulong_0("ilu", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_ilu_double_ulong_0("ilu", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_ilu_double_ulong_10("ilu", 10);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_ilu_double_ulong_10("ilu", 10);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_ilu_double_ulong_10("ilu", 10);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_sor_double_ulong("sor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_sor_double_ulong("sor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_sor_double_ulong("sor");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_ssor_double_ulong("ssor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_ssor_double_ulong("ssor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_ssor_double_ulong("ssor");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_spai_double_ulong_0("spai", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_spai_double_ulong_0("spai", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_spai_double_ulong_0("spai", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_spai_double_ulong_1("spai", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_spai_double_ulong_1("spai", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_spai_double_ulong_1("spai", 1);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_spai_double_ulong_2("spai", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_spai_double_ulong_2("spai", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_spai_double_ulong_2("spai", 2);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_spai_double_ulong_3("spai", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_spai_double_ulong_3("spai", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_spai_double_ulong_3("spai", 3);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_spai_double_ulong_80("spai", 80);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_spai_double_ulong_80("spai", 80);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_spai_double_ulong_80("spai", 80);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_spai_double_ulong_81("spai", 81);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_spai_double_ulong_81("spai", 81);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_spai_double_ulong_81("spai", 81);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_spai_double_ulong_82("spai", 82);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_spai_double_ulong_82("spai", 82);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_spai_double_ulong_82("spai", 82);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCSR<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_csr_spai_double_ulong_83("spai", 83);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixCOO<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_coo_spai_double_ulong_83("spai", 83);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             SparseMatrixELL<Mem::Main, double, unsigned long> >
bicgstab_test_cpu_ell_spai_double_ulong_83("spai", 83);


/**
 * \brief BiCGStabBandedTest
 *
 * This class defines the BiCGStabBandedTest with different preconditioners and a banded matrix
 *
 * \author Christoph Lohmann
 */
template<SparsePreconType PType_, typename VT_>
class BiCGStabBandedTest
  : public FullTaggedTest<typename VT_::MemType,
                          typename VT_::DataType, typename VT_::IndexType>
{
private:
  const Index _opt;

public:
  typedef typename VT_::DataType   DT_;
  typedef typename VT_::MemType    Mem_;
  typedef typename VT_::IndexType  IT_;
  typedef SparseMatrixBanded<Mem_, DT_, IT_> MatrixType;

  BiCGStabBandedTest(String pname, Index opt = 0)
    : FullTaggedTest<Mem_, DT_, IT_> ("bicgstab_banded_test: "
                                             + pname
                                             + " opt = "
                                             + stringify(opt)), _opt(opt)
  {
  }

  virtual void run() const override
  {
    std::vector<IT_> num_of_subintervalls;
    num_of_subintervalls.push_back(2);
    num_of_subintervalls.push_back(21);
    num_of_subintervalls.push_back(12);
    num_of_subintervalls.push_back(34);

    std::vector<DT_> dimensions;
    dimensions.push_back(DT_(3.0));
    dimensions.push_back(DT_(1.5));
    dimensions.push_back(DT_(48.3));

    // generate FD matrix
    PointstarFactoryFD2<DT_, IT_> factory(num_of_subintervalls, dimensions);
    MatrixType sys;
    sys.convert(factory.matrix_banded());

    Index size(sys.rows());
    VT_ x(size, DT_(0));
    VT_ ref;
    ref.convert(factory.vector_q2_bubble());
    VT_ b(size);
    sys.apply(b, ref);

    Preconditioner<MatrixType, VT_> * precond(Precon<PType_>::template get<MatrixType, VT_>(sys, _opt));
    BiCGStab::value(x, sys, b, *precond, 1000, 1e-12);
    delete precond;

    // check, if the result is correct
    for (Index i(0) ; i < size ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(x(i), ref(i), 1e-5);
    }
  }
};

BiCGStabBandedTest<SparsePreconType::pt_none,
                   DenseVector<Mem::Main, double, unsigned int> >
bicgstabbanded_test_cpu_none_double_uint("none");

BiCGStabBandedTest<SparsePreconType::pt_jacobi,
                   DenseVector<Mem::Main, double, unsigned int> >
bicgstabbanded_test_cpu_jac_double_uint("jacobi");

BiCGStabBandedTest<SparsePreconType::pt_none,
                   DenseVector<Mem::Main, double, unsigned long> >
bicgstabbanded_test_cpu_none_double_ulong("none");

BiCGStabBandedTest<SparsePreconType::pt_jacobi,
                   DenseVector<Mem::Main, double, unsigned long> >
bicgstabbanded_test_cpu_jac_double_ulong("jacobi");


/**
 * \brief ILUTest
 *
 * This class defines an ILUTest in which the LU-decomposition is computed externally
 *
 * \author Christoph Lohmann
 */
template<typename MT_>
class ILUTest
  : public FullTaggedTest<typename MT_::MemType, typename MT_::DataType, typename MT_::IndexType>
{
private:
  typedef typename MT_::DataType  DT_;
  typedef typename MT_::IndexType IT_;
  typedef typename MT_::MemType   Mem_;
  typedef DenseVector<Mem_, DT_, IT_> VectorType;

public:

  ILUTest()
    : FullTaggedTest<Mem_, DT_, IT_>("ilu_test" + MT_::name())
  {
  }

  virtual void run() const override
  {
    Index size(1000);
    VectorType x(size);
    VectorType ref(size);
    VectorType b(size);
    VectorType tmp(size);

    // Define solution vector
    for (Index i(0) ; i < size ; ++i)
    {
      x(i, DT_(42.34 * (i%7 + 1)));
    }

    // Define matrices
    SparseMatrixCOO<Mem_, DT_> cL (size, size);
    SparseMatrixCOO<Mem_, DT_> cU (size, size);
    SparseMatrixCOO<Mem_, DT_> cLU(size, size);

    for (Index i(0) ; i < size ; ++i)
    {
      cLU(i, i, DT_(40 + (i%2)));
      cU(i, i, DT_(40 + (i%2)));
      cL(i, i, DT_(1.0));
    }
    for (Index i(1) ; i < size ; ++i)
    {
      cLU(i, i - 1, DT_(0.01 - 0.5 * (i%2)));
      cL(i, i - 1, DT_(0.01 - 0.5 * (i%2)));
    }
    for (Index i(0) ; i < size - 1; ++i)
    {
      cU(i, i + 1, DT_(-0.7 - 0.2 * (i%3)));
      cLU(i, i + 1, DT_(-0.7 - 0.2 * (i%3)));
    }

    MT_ L;
    L.convert(cL);
    MT_ U;
    U.convert(cU);
    MT_ LU;
    LU.convert(cLU);

    // calculate the reference-solution
    //tmp.product_matvec(U, x);
    U.apply(tmp, x);
    //b.product_matvec(L, tmp);
    L.apply(b, tmp);

    // save reference-solution
    ref.copy(x);

    ILUPreconditioner<MT_, VectorType> precond(LU);
    precond.apply(x, b);

    // check, if the result is coorect
    for (Index i(0) ; i < x.size() ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(ref(i), x(i), 1e-5);
    }
  }
};

ILUTest<SparseMatrixCSR<Mem::Main, double, unsigned int> > ilu_test_cpu_csr_double_uint;
ILUTest<SparseMatrixCSR<Mem::Main, double, unsigned long> > ilu_test_cpu_csr_double_ulong;
ILUTest<SparseMatrixCOO<Mem::Main, double, unsigned int> > ilu_test_cpu_coo_double_uint;
ILUTest<SparseMatrixCOO<Mem::Main, double, unsigned long> > ilu_test_cpu_coo_double_ulong;
ILUTest<SparseMatrixELL<Mem::Main, double, unsigned int> > ilu_test_cpu_ell_double_uint;
ILUTest<SparseMatrixELL<Mem::Main, double, unsigned long> > ilu_test_cpu_ell_double_ulong;
