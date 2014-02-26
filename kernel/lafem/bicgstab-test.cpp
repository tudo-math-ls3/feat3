#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/bicgstab.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

#include <kernel/lafem/richardson.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;


/**
 * \brief BiCGStabTest
 *
 * This class defines the BiCGStabTest with different preconditioners and matrix-types
 *
 * \author Christoph Lohmann
 */
template< typename PSF_, typename PT_>
class BiCGStabTest
  : public TaggedTest<typename PT_::MemType,
                      typename PT_::DataType,
                      typename PT_::AlgoType>
{
private:
  const int opt;

public:
  typedef typename PT_::AlgoType   Algo_;
  typedef typename PT_::MatrixType MT_;
  typedef typename PT_::DataType   DT_;
  typedef typename PT_::MemType    Mem_;
  typedef DenseVector<Mem_, DT_>   VT_;

  BiCGStabTest(int opt = 0)
    : TaggedTest<Mem_, DT_, Algo_> ("bicgstab_test: "
                                    + MT_::name()
                                    + " "
                                    + PT_::name()
                                    + " opt = "
                                    + stringify(opt)), opt(opt)
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

    // solver-paramters
    Index max_iter = 5000;
    DT_ eps = 1e-12;

    // Define preconditioners for every matrix-type and solve the system
    if (typeid(PT_) == typeid(NonePreconditioner<Algo_, MT_, VT_>))
    {
      NonePreconditioner<Algo_, MT_, VT_> precond;
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(JacobiPreconditioner<Algo_, MT_, VT_>))
    {
      JacobiPreconditioner<Algo_, MT_, VT_> precond(sys, DT_(1));
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(GaussSeidelPreconditioner<Algo_, MT_, VT_>))
    {
      GaussSeidelPreconditioner<Algo_, MT_,
                                VT_> precond(sys, DT_(1));
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(PolynomialPreconditioner<Algo_, MT_, VT_>))
    {
      PolynomialPreconditioner<Algo_, MT_,
                               VT_> precond(sys, 20, opt == 0);
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(ILUPreconditioner<Algo_, MT_, VT_>))
    {
      ILUPreconditioner<Algo_, MT_,
                        VT_> precond(sys, opt);
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(SORPreconditioner<Algo_, MT_, VT_>))
    {
      SORPreconditioner<Algo_, MT_,
                        VT_> precond(sys);
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(SSORPreconditioner<Algo_, MT_, VT_>))
    {
      SSORPreconditioner<Algo_, MT_,
                         VT_> precond(sys);
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(SPAIPreconditioner<Algo_, MT_, VT_>))
    {
      if (opt%2 == 0)
      {
        SPAIPreconditioner<Algo_, MT_, VT_> precond(sys, 2, opt / 2);
        BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
      }
      else
      {
        SPAIPreconditioner<Algo_, MT_, VT_> precond(sys, sys.layout(), (opt - 1) / 2);
        BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
      }
    }
    else
    {
      throw InternalError(__func__, __FILE__, __LINE__, "Preconditioner and Matrix have different matrix-types!");
    }


    // check, if the result is correct
    for (Index i(0) ; i < size ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(x(i), ref_local(i), 1e-8);
    }

  }

};


BiCGStabTest<PointstarFactoryFD<double>,
             NonePreconditioner<Algo::Generic,
                                SparseMatrixCSR<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_csr_none_double;

BiCGStabTest<PointstarFactoryFD<double>,
             NonePreconditioner<Algo::Generic,
                                SparseMatrixCOO<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_coo_none_double;

BiCGStabTest<PointstarFactoryFD<double>,
             NonePreconditioner<Algo::Generic,
                                SparseMatrixELL<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_ell_none_double;


BiCGStabTest<PointstarFactoryFD<double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixCSR<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_csr_jacobi_double;

BiCGStabTest<PointstarFactoryFD<double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixCOO<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_coo_jacobi_double;

BiCGStabTest<PointstarFactoryFD<double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixELL<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_ell_jacobi_double;


BiCGStabTest<PointstarFactoryFD<double>,
             GaussSeidelPreconditioner<Algo::Generic,
                                       SparseMatrixCSR<Mem::Main, double>,
                                       DenseVector<Mem::Main, double> > >
bicgstab_test_csr_gaussSeidel_double;

BiCGStabTest<PointstarFactoryFD<double>,
             GaussSeidelPreconditioner<Algo::Generic,
                                       SparseMatrixCOO<Mem::Main, double>,
                                       DenseVector<Mem::Main, double> > >
bicgstab_test_coo_gaussSeidel_double;

BiCGStabTest<PointstarFactoryFD<double>,
             GaussSeidelPreconditioner<Algo::Generic,
                                       SparseMatrixELL<Mem::Main, double>,
                                       DenseVector<Mem::Main, double> > >
bicgstab_test_ell_gaussSeidel_double;


BiCGStabTest<PointstarFactoryFD<double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixCSR<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_csr_poly_double;

BiCGStabTest<PointstarFactoryFD<double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixCOO<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_coo_poly_double;

BiCGStabTest<PointstarFactoryFD<double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixELL<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_ell_poly_double;


// BiCGStabTest<PointstarFactoryFD<double>,
//              PolynomialPreconditioner<Algo::Generic,
//                                       SparseMatrixCSR<Mem::Main, double>,
//                                       DenseVector<Mem::Main, double> > >
// bicgstab_test_csr_poly_no_double(1);

// BiCGStabTest<PointstarFactoryFD<double>,
//              PolynomialPreconditioner<Algo::Generic,
//                                       SparseMatrixCOO<Mem::Main, double>,
//                                       DenseVector<Mem::Main, double> > >
// bicgstab_test_coo_poly_no_double(1);

// BiCGStabTest<PointstarFactoryFD<double>,
//              PolynomialPreconditioner<Algo::Generic,
//                                       SparseMatrixELL<Mem::Main, double>,
//                                       DenseVector<Mem::Main, double> > >
// bicgstab_test_ell_poly_no_double(1);


BiCGStabTest<PointstarFactoryFD<double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixCSR<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_csr_ilu_0_double;

BiCGStabTest<PointstarFactoryFD<double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixCOO<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_coo_ilu_0_double;

BiCGStabTest<PointstarFactoryFD<double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixELL<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_ell_ilu_0_double;


BiCGStabTest<PointstarFactoryFD<double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixCSR<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_csr_ilu_10_double(10);

BiCGStabTest<PointstarFactoryFD<double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixCOO<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_coo_ilu_10_double(10);

BiCGStabTest<PointstarFactoryFD<double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixELL<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_ell_ilu_10_double(10);


BiCGStabTest<PointstarFactoryFD<double>,
             SORPreconditioner<Algo::Generic,
                               SparseMatrixCSR<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_csr_sor_double;

BiCGStabTest<PointstarFactoryFD<double>,
             SORPreconditioner<Algo::Generic,
                               SparseMatrixCOO<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_coo_sor_double;

BiCGStabTest<PointstarFactoryFD<double>,
             SORPreconditioner<Algo::Generic,
                               SparseMatrixELL<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_ell_sor_double;


BiCGStabTest<PointstarFactoryFD<double>,
             SSORPreconditioner<Algo::Generic,
                                SparseMatrixCSR<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_csr_ssor_double;

BiCGStabTest<PointstarFactoryFD<double>,
             SSORPreconditioner<Algo::Generic,
                                SparseMatrixCOO<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_coo_ssor_double;

BiCGStabTest<PointstarFactoryFD<double>,
             SSORPreconditioner<Algo::Generic,
                                SparseMatrixELL<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_ell_ssor_double;


BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixCSR<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_csr_spai_double_0;

BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixCOO<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_coo_spai_double_0;

BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixELL<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_ell_spai_double_0;


BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixCSR<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_csr_spai_double_1(1);

BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixCOO<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_coo_spai_double_1(1);

BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixELL<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_ell_spai_double_1(1);


BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixCSR<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_csr_spai_double_20(20);

BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixCOO<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_coo_spai_double_20(20);

BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixELL<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_ell_spai_double_20(20);


BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixCSR<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_csr_spai_double_21(21);

BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixCOO<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_coo_spai_double_21(21);

BiCGStabTest<PointstarFactoryFD<double>,
             SPAIPreconditioner<Algo::Generic,
                                SparseMatrixELL<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_ell_spai_double_21(21);

/*
#ifdef FEAST_BACKENDS_CUDA
BiCGStabTest<PointstarFactoryFD<double>,
             NonePreconditioner<Algo::CUDA,
                                SparseMatrixCSR<Mem::CUDA, double>,
                                DenseVector<Mem::CUDA, double> > >
cuda_bicgstab_test_csr_none_double(0);

BiCGStabTest<PointstarFactoryFD<double>,
             NonePreconditioner<Algo::CUDA,
                                SparseMatrixELL<Mem::CUDA, double>,
                                DenseVector<Mem::CUDA, double> > >
cuda_bicgstab_test_ell_none_double(0);
#endif
*/


/**
 * \brief ILUTest
 *
 * This class defines an ILUTest in which the LU-decomposition is computed externally
 *
 * \author Christoph Lohmann
 */
template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename MT_>
class ILUTest
  : public TaggedTest<Mem_, DT_, Algo_>
{

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
    SparseMatrixCOO<Mem::Main, DT_> cL (size, size);
    SparseMatrixCOO<Mem::Main, DT_> cU (size, size);
    SparseMatrixCOO<Mem::Main, DT_> cLU(size, size);

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

ILUTest<Mem::Main, Algo::Generic, double,
        SparseMatrixCSR<Mem::Main, double> > ilu_test_src_double;
ILUTest<Mem::Main, Algo::Generic, double,
        SparseMatrixCOO<Mem::Main, double> > ilu_test_coo_double;
ILUTest<Mem::Main, Algo::Generic, double,
        SparseMatrixELL<Mem::Main, double> > ilu_test_ell_double;
