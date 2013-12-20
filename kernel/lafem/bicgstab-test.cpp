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


template<
  typename Arch_,
  typename Algo_,
  typename DT_,
  typename MT_,
  typename PT_>
class BiCGStabTest
  : public TaggedTest<Arch_, DT_, Algo_>
{
private:
  const Index matrix_struct;

public:

  BiCGStabTest(Index matrix_struct = 0)
    : TaggedTest<Arch_, DT_, Algo_>("bicgstab_test: " + MT_::type_name() + " "
                                    + PT_::type_name() + " matrix_struct = "
                                    + std::to_string(matrix_struct)), matrix_struct(matrix_struct)
  {
  }

  virtual void run() const
  {
    Index size(1025);
    DenseVector<Arch_, DT_> x(size, DT_(1));
    DenseVector<Arch_, DT_> b(size);
    DenseVector<Arch_, DT_> ref(size, DT_(42));

    // Define matrix
    SparseMatrixCOO<Mem::Main, DT_> csys(size, size);

    if (matrix_struct == 0)
    {
      for (Index i(0) ; i < size ; ++i)
        csys(i, i, DT_(4 + (i%2)));
      for (Index i(1) ; i < size ; ++i)
        csys(i - 1, i, DT_(-1 - 0.2 * (i%3)));
      for (Index i(0) ; i < size - 1; ++i)
        csys(i + 1, i, DT_(-2 - 0.5 * (i%2)));
    }

    // matix_struct = 1 for polynomial preconditioner without scaling
    else if (matrix_struct == 1)
    {
      for (Index i(0) ; i < size ; ++i)
        csys(i, i, DT_(0.4 + 0.1 * (i%2)));
      for (Index i(1) ; i < size ; ++i)
        csys(i - 1, i, DT_(-0.1 - 0.02 * (i%3)));
      for (Index i(0) ; i < size - 1; ++i)
        csys(i + 1, i, DT_(-0.2 - 0.05 * (i%2)));
    }

    MT_ sys(csys);

    // calculate the reference-solution
    b.template product_matvec<Algo_>(sys, ref);

    // solver-paramters
    Index max_iter = 1000;
    DT_ eps = 1e-15;


    // Define preconditioners for every matrix-type and solve the system
    if (typeid(PT_) == typeid(NonePreconditioner<Algo_, MT_, DenseVector<Arch_,
                              DT_> >))
    {
      NonePreconditioner<Algo_, MT_, DenseVector<Arch_, DT_> > precond;
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(JacobiPreconditioner<Algo_, MT_, DenseVector<Arch_,
                              DT_> >))
    {
      JacobiPreconditioner<Algo_, MT_, DenseVector<Arch_, DT_> > precond(sys, DT_(1));
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(GaussSeidelPreconditioner<Algo_, MT_, DenseVector<Arch_,
                              DT_> >))
    {
      GaussSeidelPreconditioner<Algo_, MT_,
                                DenseVector<Arch_, DT_> > precond(sys, DT_(1));
      BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
    }
    else if (typeid(PT_) == typeid(PolynomialPreconditioner<Algo_, MT_, DenseVector<Arch_,
                              DT_> >))
    {
      if (matrix_struct == 0)
      {
        PolynomialPreconditioner<Algo_, MT_,
                                 DenseVector<Arch_, DT_> > precond(sys, 20, true);
        BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
      }
      else if(matrix_struct == 1)
      {
        PolynomialPreconditioner<Algo_, MT_,
                                 DenseVector<Arch_, DT_> > precond(sys, 20, false);
        BiCGStab<Algo_>::value(x, sys, b, precond, max_iter, eps);
      }
    }


    // check, if the result is coorect
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
bicgstab_test_csr_none_double(0);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             NonePreconditioner<Algo::Generic,
                                SparseMatrixCOO<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_coo_none_double(0);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             NonePreconditioner<Algo::Generic,
                                SparseMatrixELL<Mem::Main, double>,
                                DenseVector<Mem::Main, double> > >
bicgstab_test_ell_none_double(0);


BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCSR<Mem::Main, double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixCSR<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_csr_jacobi_double(0);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixCOO<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_coo_jacobi_double(0);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             JacobiPreconditioner<Algo::Generic,
                                  SparseMatrixELL<Mem::Main, double>,
                                  DenseVector<Mem::Main, double> > >
bicgstab_test_ell_jacobi_double(0);



 BiCGStabTest<Mem::Main, Algo::Generic, double,
              SparseMatrixCSR<Mem::Main, double>,
              GaussSeidelPreconditioner<Algo::Generic,
                                        SparseMatrixCSR<Mem::Main, double>,
                                        DenseVector<Mem::Main, double> > >
 bicgstab_test_csr_gaussSeidel_double(0);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixCOO<Mem::Main, double>,
             GaussSeidelPreconditioner<Algo::Generic,
                                       SparseMatrixCOO<Mem::Main, double>,
                                       DenseVector<Mem::Main, double> > >
bicgstab_test_coo_gaussSeidel_double(0);

BiCGStabTest<Mem::Main, Algo::Generic, double,
             SparseMatrixELL<Mem::Main, double>,
             GaussSeidelPreconditioner<Algo::Generic,
                                       SparseMatrixELL<Mem::Main, double>,
                                       DenseVector<Mem::Main, double> > >
bicgstab_test_ell_gaussSeidel_double(0);



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




#ifdef FEAST_GMP
//BiCGStabTest<Mem::Main, Algo::Generic, mpf_class> bicgstab_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
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
#endif
