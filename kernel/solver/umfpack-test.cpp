// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>

#ifdef FEAT_HAVE_UMFPACK
#include <kernel/solver/umfpack.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

class UmfpackTest :
  public TestSystem::FullTaggedTest<Mem::Main, double, Index>
{
public:
  typedef SparseMatrixCSR<Mem::Main, double, Index> MatrixType;
  typedef DenseVector<Mem::Main, double, Index> VectorType;

  UmfpackTest() :
    TestSystem::FullTaggedTest<Mem::Main, double, Index>("UmfpackTest")
  {
  }

  virtual ~UmfpackTest()
  {
  }

  virtual void run() const override
  {
    const double tol = Math::pow(Math::eps<double>(), 0.6);

    // create a pointstar factory
    PointstarFactoryFD<double> psf(17);

    // create a CSR matrix
    MatrixType mat_sys(psf.matrix_csr());

    // create a reference solution vector
    VectorType vec_ref(psf.eigenvector_min());

    // create an rhs vector
    VectorType vec_rhs(mat_sys.create_vector_r());
    vec_rhs.scale(vec_ref, psf.lambda_min());

    // create an empty solution vector
    VectorType vec_sol(mat_sys.create_vector_r());
    vec_sol.format();

    {
      // create an Umfpack solver
      Umfpack umfpack(mat_sys);
      // initialize
      umfpack.init();
      // solve
      umfpack.apply(vec_sol, vec_rhs);
      // release
      umfpack.done();
    }

    // subtract reference solution
    vec_sol.axpy(vec_ref, vec_sol, -1.0);

    // compute the norm
    double nrm2 = vec_sol.norm2();

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(nrm2, 0.0, tol);
  }
};

UmfpackTest umfpack_test;


class UmfpackMeanTest :
  public TestSystem::FullTaggedTest<Mem::Main, double, Index>
{
public:
  typedef SparseMatrixCSR<Mem::Main, double, Index> MatrixType;
  typedef DenseVector<Mem::Main, double, Index> VectorType;

  UmfpackMeanTest() :
    TestSystem::FullTaggedTest<Mem::Main, double, Index>("UmfpackMeanTest")
  {
  }

  virtual ~UmfpackMeanTest()
  {
  }

  virtual void run() const override
  {
    const double tol = Math::pow(Math::eps<double>(), 0.8);
    const double pi = Math::pi<double>();

    // number of points
    const Index n = 17;

    // create a 1D pointstar factory
    PointstarFactoryFD<double> psf(n, 1);

    // create a CSR matrix
    MatrixType mat_sys(psf.matrix_csr());

    // scale the a_11 and a_nn by 1/2 to emulate Neumann BCs
    {
      double* data = mat_sys.val();
      data[0] *= 0.5;
      data[mat_sys.used_elements()-1] *= 0.5;
    }

    // create a one-vector as Lagrange multiplier
    VectorType vec_one(n, 1.0);

    // create a reference solution vector: x_i := cos(pi*i/(n-1))
    VectorType vec_ref(n, 0.0);
    for(Index i(0); i < n; ++i)
    {
      vec_ref(i, Math::cos(pi * double(i) / double(n-1)));
    }

    // create an rhs vector
    VectorType vec_rhs(mat_sys.create_vector_r());
    mat_sys.apply(vec_rhs, vec_ref);

    // create an empty solution vector
    VectorType vec_sol(mat_sys.create_vector_r());
    vec_sol.format();

    {
      // create an UmfpackMean solver
      UmfpackMean umfpack(mat_sys, vec_one);
      // initialize
      umfpack.init();
      // solve
      umfpack.apply(vec_sol, vec_rhs);
      // release
      umfpack.done();
    }

    // subtract reference solution
    vec_sol.axpy(vec_ref, vec_sol, -1.0);

    // compute the norm
    double nrm2 = vec_sol.norm2();

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(nrm2, 0.0, tol);
  }
};

UmfpackMeanTest umfpack_mean_test;

template<typename DT_>
class GenericUmfpackTest :
  public TestSystem::FullTaggedTest<Mem::Main, DT_, Index>
{
public:
  typedef SparseMatrixCSR<Mem::Main, DT_, Index> MatrixType;
  typedef DenseVector<Mem::Main, DT_, Index> VectorType;

  GenericUmfpackTest() :
    TestSystem::FullTaggedTest<Mem::Main, DT_, Index>("GenericUmfpackTest")
  {
  }

  virtual ~GenericUmfpackTest()
  {
  }

  virtual void run() const override
  {
    // same tolerance for all DT_, since UMFPACK is internally double precision only
    const DT_ tol = DT_(1E-8);

    // create a pointstar factory
    PointstarFactoryFD<DT_> psf(17);

    // create a CSR matrix
    MatrixType mat_sys(psf.matrix_csr());

    // create a reference solution vector
    VectorType vec_ref(psf.eigenvector_min());

    // create an rhs vector
    VectorType vec_rhs(mat_sys.create_vector_r());
    vec_rhs.scale(vec_ref, psf.lambda_min());

    // create an empty solution vector
    VectorType vec_sol(mat_sys.create_vector_r());
    vec_sol.format();

    {
      // create a GenericUmfpack solver
      GenericUmfpack<MatrixType> umfpack(mat_sys);
      // initialize
      umfpack.init();
      // solve
      umfpack.apply(vec_sol, vec_rhs);
      // release
      umfpack.done();
    }

    // subtract reference solution
    vec_sol.axpy(vec_ref, vec_sol, -DT_(1));

    // compute the norm
    DT_ nrm2 = vec_sol.norm2();

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(nrm2, DT_(0), tol);
  }
};

GenericUmfpackTest<double> generic_umfpack_test_double;
#ifdef FEAT_HAVE_QUADMATH
GenericUmfpackTest<__float128> generic_umfpack_test_float128;
#endif
#endif // FEAT_HAVE_UMFPACK
