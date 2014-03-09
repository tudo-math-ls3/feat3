#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/bicgstab.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

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
      t = new SPAIPreconditioner<Algo_, MT_, VT_> (sys, 2, max_iter, 1e-2, 10, 1e-3, 1e-3, transpose);
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
class BiCGStabTest
  : public TaggedTest<typename VT_::MemType,
                      typename VT_::DataType,
                      Algo_>
{
private:
  const Index _opt;

public:
  typedef typename VT_::DataType   DT_;
  typedef typename VT_::MemType    Mem_;

  BiCGStabTest(String pname, Index opt = 0)
    : TaggedTest<Mem_, DT_, Algo_> ("bicgstab_test: "
                                    + MT_::name()
                                    + " "
                                    + pname
                                    + " opt = "
                                    + stringify(opt)), _opt(opt)
  {
  }

  virtual void run() const
  {
    PSF_ factory(13);
    MT_ sys(factory.matrix_csr());

    Index size(sys.rows());
    VT_ x(size, DT_(1));
    VT_ ref(factory.vector_q2_bubble());
    VT_ ref_local(ref);
    VT_ b(size);
    sys.template apply<Algo_>(b, ref);

    Preconditioner<Algo_, MT_, VT_> * precond(Precon<PType_>::template get<Algo_, MT_, VT_>(sys, _opt));
    BiCGStab<Algo_>::value(x, sys, b, *precond, 1000, 1e-12);
    delete precond;

    // check, if the result is correct
    for (Index i(0) ; i < size ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(x(i), ref_local(i), 1e-8);
    }

  }

};


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_none_double("none");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_none_double("none");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_none_double("none");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_jac_double("jacobi");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_jac_double("jacobi");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_jac_double("jacobi");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_gs_double("gauss-seidel");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_gs_double("gauss-seidel");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_gs_double("gauss-seidel");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_poly_double("polynmial");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_poly_double("polynmial");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_poly_double("polynmial");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_ilu_double_0("ilu", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_ilu_double_0("ilu", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_ilu_double_0("ilu", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_ilu_double_10("ilu", 10);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_ilu_double_10("ilu", 10);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_ilu_double_10("ilu", 10);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_sor_double("sor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_sor_double("sor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_sor_double("sor");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_ssor_double("ssor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_ssor_double("ssor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_ssor_double("ssor");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_spai_double_0("spai", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_spai_double_0("spai", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_spai_double_0("spai", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_spai_double_1("spai", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_spai_double_1("spai", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_spai_double_1("spai", 1);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_spai_double_2("spai", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_spai_double_2("spai", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_spai_double_2("spai", 2);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_spai_double_3("spai", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_spai_double_3("spai", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_spai_double_3("spai", 3);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_spai_double_80("spai", 80);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_spai_double_80("spai", 80);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_spai_double_80("spai", 80);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_spai_double_81("spai", 81);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_spai_double_81("spai", 81);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_spai_double_81("spai", 81);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_spai_double_82("spai", 82);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_spai_double_82("spai", 82);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_spai_double_82("spai", 82);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCSR<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_csr_spai_double_83("spai", 83);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixCOO<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_coo_spai_double_83("spai", 83);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::Generic,
             SparseMatrixELL<Mem::Main, double>,
             DenseVector<Mem::Main, double> >
bicgstab_test_cpu_ell_spai_double_83("spai", 83);


#ifdef FEAST_BACKENDS_CUDA
BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_none_double("none");

/*BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_none_double("none");*/

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_none,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_none_double("none");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_jac_double("jacobi");

/*BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_jac_double("jacobi");*/

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_jacobi,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_jac_double("jacobi");


/*BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_gs_double("gauss-seidel");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_gs_double("gauss-seidel");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_gauss_seidel,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_gs_double("gauss-seidel");*/


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_poly_double("polynmial");

/*BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_poly_double("polynmial");*/

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_polynomial,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_poly_double("polynmial");


/*BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_ilu_double_0("ilu", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_ilu_double_0("ilu", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_ilu_double_0("ilu", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_ilu_double_10("ilu", 10);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_ilu_double_10("ilu", 10);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ilu,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_ilu_double_10("ilu", 10);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_sor_double("sor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_sor_double("sor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_sor,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_sor_double("sor");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_ssor_double("ssor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_ssor_double("ssor");

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_ssor,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_ssor_double("ssor");


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_spai_double_0("spai", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_spai_double_0("spai", 0);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_spai_double_0("spai", 0);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_spai_double_1("spai", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_spai_double_1("spai", 1);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_spai_double_1("spai", 1);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_spai_double_2("spai", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_spai_double_2("spai", 2);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_spai_double_2("spai", 2);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_spai_double_3("spai", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_spai_double_3("spai", 3);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_spai_double_3("spai", 3);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_spai_double_80("spai", 80);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_spai_double_80("spai", 80);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_spai_double_80("spai", 80);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_spai_double_81("spai", 81);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_spai_double_81("spai", 81);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_spai_double_81("spai", 81);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_spai_double_82("spai", 82);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_spai_double_82("spai", 82);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_spai_double_82("spai", 82);


BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCSR<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_csr_spai_double_83("spai", 83);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixCOO<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_coo_spai_double_83("spai", 83);

BiCGStabTest<PointstarFactoryFE<double>,
             SparsePreconType::pt_spai,
             Algo::CUDA,
             SparseMatrixELL<Mem::CUDA, double>,
             DenseVector<Mem::CUDA, double> >
bicgstab_test_cuda_ell_spai_double_83("spai", 83);*/
#endif


/**
 * \brief ILUTest
 *
 * This class defines an ILUTest in which the LU-decomposition is computed externally
 *
 * \author Christoph Lohmann
 */
template<typename Algo_, typename MT_>
class ILUTest
  : public TaggedTest<typename MT_::MemType, typename MT_::DataType, Algo_>
{
private:
  typedef typename MT_::DataType DT_;
  typedef typename MT_::MemType Mem_;

public:

  ILUTest()
    : TaggedTest<Mem_, DT_, Algo_>("ilu_test" + MT_::name())
  {
  }

  virtual void run() const
  {
    Index size(1000);
    DenseVector<Mem_, DT_> x(size);
    DenseVector<Mem_, DT_> ref(size);
    DenseVector<Mem_, DT_> b(size);
    DenseVector<Mem_, DT_> tmp(size);

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

    MT_ L(cL);
    MT_ U(cU);
    MT_ LU(cLU);

    // calculate the reference-solution
    //tmp.template product_matvec<Algo_>(U, x);
    U.template apply<Algo_>(tmp, x);
    //b.template product_matvec<Algo_>(L, tmp);
    L.template apply<Algo_>(b, tmp);

    // save reference-solution
    ref.copy(x);

    ILUPreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> > precond(LU);
    precond.apply(x, b);

    // check, if the result is coorect
    for (Index i(0) ; i < x.size() ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(ref(i), x(i), 1e-5);
    }
  }
};

ILUTest<Algo::Generic, SparseMatrixCSR<Mem::Main, double> > ilu__test_cpu_src_double;
ILUTest<Algo::Generic, SparseMatrixCOO<Mem::Main, double> > ilu_test_cpu_coo_double;
ILUTest<Algo::Generic, SparseMatrixELL<Mem::Main, double> > ilu_test_cpu_ell_double;


/*#ifdef FEAST_BACKENDS_CUDA
ILUTest<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double> > ilu_test_cuda_src_double;
ILUTest<Algo::CUDA, SparseMatrixCOO<Mem::CUDA, double> > ilu_test_cuda_coo_double;
ILUTest<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double> > ilu_test_cuda_ell_double;
#endif*/
