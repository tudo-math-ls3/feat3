#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/bicgstab.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/algorithm.hpp>

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
template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename MT_,
  typename PT_>
class BiCGStabTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
private:
  const int opt;

public:

  BiCGStabTest(int opt = 0)
    : TaggedTest<Mem_, DT_, Algo_>("bicgstab_test: " + MT_::type_name() + " "
                                    + PT_::type_name() + " opt = "
                                    + stringify(opt)), opt(opt)
  {
  }

  virtual void run() const
  {
    Index size(1000);
    DenseVector<Mem_, DT_> x(size, DT_(1));
    DenseVector<Mem_, DT_> b(size);
    DenseVector<Mem_, DT_> ref(size, DT_(42));

    // Define matrix
    SparseMatrixCOO<Mem::Main, DT_> csys(size, size);

    // opt = 1: alternative matrix for polynomial preconditioner without scaling
    if (typeid(PT_) == typeid(PolynomialPreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> >) && opt == 1)
    {
      for (Index i(0) ; i < size ; ++i)
        csys(i, i, DT_(0.4 + 0.1 * (i%2)));
      for (Index i(1) ; i < size ; ++i)
        csys(i - 1, i, DT_(-0.1 - 0.02 * (i%3)));
      for (Index i(0) ; i < size - 1; ++i)
        csys(i + 1, i, DT_(-0.2 - 0.05 * (i%2)));
    }
    // default case
    else
    {
      for (Index i(0) ; i < size ; ++i)
      {
        for (Index j(0) ; j < size ; ++j)
        {
          if (i == j)
            csys(i, j, DT_(4 + (i%2)));
          else if (i - 1 == j)
            csys(i, j, DT_(-1 - 0.2 * (i%3)));
          else if (i + 1 == j)
            csys(i, j, DT_(-2 - 0.5 * (i%2)));
          else if (i + 5 == j)
            csys(i, j, DT_(0.05 * (i%2)));
          else if (i - 5 == j)
            csys(i, j, DT_(- 0.5 * (i%2)));
          else if (i + 3 == j)
          {
            if (i%2 == 0)
              csys(i, j, DT_(0.234 * (i%7)));
          }
        }
      }
    }
    MT_ sys(csys);

    // calculate the reference-solution
    b.template product_matvec<Algo_>(sys, ref);

    // solver-paramters
    Index max_iter = 1000;
    DT_ eps = 1e-15;


    // Define preconditioners for every matrix-type and solve the system
    if (typeid(PT_) == typeid(NonePreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> >))
    {
      NonePreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> > precond;
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(JacobiPreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> >))
    {
      JacobiPreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> > precond(sys, DT_(1));
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(GaussSeidelPreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> >))
    {
      GaussSeidelPreconditioner<Algo_, MT_,
                                DenseVector<Mem_, DT_> > precond(sys, DT_(1));
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(PolynomialPreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> >))
    {
      if (opt == 0)
      {
        PolynomialPreconditioner<Algo_, MT_,
                                 DenseVector<Mem_, DT_> > precond(sys, 20, true);
        BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
      }
      else if(opt == 1)
      {
        PolynomialPreconditioner<Algo_, MT_,
                                 DenseVector<Mem_, DT_> > precond(sys, 20, false);
        BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
      }
    }
    else if (typeid(PT_) == typeid(ILUPreconditioner<Algo_, MT_, DenseVector<Mem_, DT_> >))
    {
      ILUPreconditioner<Algo_, MT_,
                        DenseVector<Mem_, DT_> > precond(sys, opt);
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else
    {
      throw InternalError(__func__, __FILE__, __LINE__, "Preconditioner and Matrix have different matrix-types!");
    }


    // check, if the result is correct
    for (Index i(0) ; i < x.size() ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(x(i), ref(i), 1e-10);
    }

  }

};


BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCSR<Mem::Main, double>,
             NonePreconditioner<Algo::Generic,
                                SparseMatrixCSR<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_csr_none_double;

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             NonePreconditioner<Algo::Generic,
                                SparseMatrixCOO<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_coo_none_double;

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             NonePreconditioner<Algo::Generic,
                                SparseMatrixELL<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_ell_none_double;


BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCSR<Mem::Main, double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixCSR<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_csr_jacobi_double;

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixCOO<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_coo_jacobi_double;

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixELL<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_ell_jacobi_double;



 BiCGStabTest<Mem::Main, Algo::Generic, double,
              SparseMatrixCSR<Mem::Main, double>,
              GaussSeidelPreconditioner<Algo::Generic,
                                        SparseMatrixCSR<Mem::Main, double>,
                                        DenseVector<Mem::Main, double> > >
 bicgstab_test_csr_gaussSeidel_double;

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             GaussSeidelPreconditioner<Algo::Generic,
                                       SparseMatrixCOO<Mem::Main, double>,
                                       DenseVector<Mem::Main, double> > >
bicgstab_test_coo_gaussSeidel_double;

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             GaussSeidelPreconditioner<Algo::Generic,
                                       SparseMatrixELL<Mem::Main, double>,
                                       DenseVector<Mem::Main, double> > >
bicgstab_test_ell_gaussSeidel_double;



BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCSR<Mem::Main, double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixCSR<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_csr_poly_double(0);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixCOO<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_coo_poly_double(0);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixELL<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_ell_poly_double(0);


BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCSR<Mem::Main, double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixCSR<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_csr_poly_no_double(1);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixCOO<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_coo_poly_no_double(1);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             PolynomialPreconditioner<Algo::Generic,
                                      SparseMatrixELL<Mem::Main, double>,
                                      DenseVector<Mem::Main, double> > >
bicgstab_test_ell_poly_no_double(1);


BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCSR<Mem::Main, double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixCSR<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_csr_ilu_0_double;

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixCOO<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_coo_ilu_0_double;

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             ILUPreconditioner<Algo::Generic,
                               SparseMatrixELL<Mem::Main, double>,
                               DenseVector<Mem::Main, double> > >
bicgstab_test_ell_ilu_0_double;


//BiCGStabTest<Mem::Main, Algo::Generic, double,
//             SparseMatrixCSR<Mem::Main, double>,
//             ILUPreconditioner<Algo::Generic,
//                               SparseMatrixCSR<Mem::Main, double>,
//                               DenseVector<Mem::Main, double> > >
//bicgstab_test_csr_ilu_10_double(10);

//BiCGStabTest<Mem::Main, Algo::Generic, double,
//             SparseMatrixCOO<Mem::Main, double>,
//             ILUPreconditioner<Algo::Generic,
//                               SparseMatrixCOO<Mem::Main, double>,
//                               DenseVector<Mem::Main, double> > >
//bicgstab_test_coo_ilu_10_double(10);

//BiCGStabTest<Mem::Main, Algo::Generic, double,
//             SparseMatrixELL<Mem::Main, double>,
//             ILUPreconditioner<Algo::Generic,
//                               SparseMatrixELL<Mem::Main, double>,
//                               DenseVector<Mem::Main, double> > >
//bicgstab_test_ell_ilu_10_double(10);


#ifdef FEAST_GMP
//BiCGStabTest<Mem::Main, Algo::Generic, mpf_class> bicgstab_test_mpf_class;
#endif
/*#ifdef FEAST_BACKENDS_CUDA
BiCGStabTest<Mem::CUDA, Algo::CUDA, double,
             SparseMatrixCSR<Mem::Main, double>,
             NonePreconditioner<Algo::CUDA,
                                SparseMatrixCSR<Mem::CUDA, double>,
                                DenseVector<Mem::CUDA, double> > >
cuda_bicgstab_test_csr_none_double(0);

BiCGStabTest<Mem::CUDA, Algo::CUDA, double,
             SparseMatrixELL<Mem::Main, double>,
             NonePreconditioner<Algo::CUDA,
                                SparseMatrixELL<Mem::CUDA, double>,
                                DenseVector<Mem::CUDA, double> > >
cuda_bicgstab_test_ell_none_double(0);
#endif*/


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
    : TaggedTest<Mem_, DT_, Algo_>("ilu_test" + MT_::type_name())
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
    tmp.template product_matvec<Algo_>(U, x);
    b.template product_matvec<Algo_>(L, tmp);

    // save reference-solution
    copy(ref, x);

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
        SparseMatrixELL<Mem::Main, double> > ilu_test_ell_double;
