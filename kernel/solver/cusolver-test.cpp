// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/solver/cusolver.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

#ifdef FEAT_HAVE_CUDA

class CuSolverLUTest :
  public TestSystem::UnitTest
{
public:
  CuSolverLUTest() :
    TestSystem::UnitTest("CuSolverQRTest", "double", "unsigned int", PreferredBackend::cuda)
  {
  }

  virtual ~CuSolverLUTest()
  {
  }

  virtual void run() const override
  {
    typedef SparseMatrixCSR<double, unsigned int> MatrixType;
    typedef DenseVector<double, unsigned int> VectorType;

    const double tol = Math::pow(Math::eps<double>(), 0.6);

    // create a pointstar factory
    PointstarFactoryFD<double, unsigned int> psf(17);

    // create a CSR matrix
    MatrixType mat_sys;
    mat_sys.convert(psf.matrix_csr());

    // create a reference solution vector
    VectorType vec_ref;
    vec_ref.convert(psf.eigenvector_min());

    // create an rhs vector
    VectorType vec_rhs;
    vec_rhs.convert(mat_sys.create_vector_r());
    vec_rhs.scale(vec_ref, psf.lambda_min());

    // create an empty solution vector
    VectorType vec_sol;
    vec_sol.convert(mat_sys.create_vector_r());
    vec_sol.format();

    CuSolverLU lusolver(mat_sys);
    lusolver.apply(vec_sol, vec_rhs);

    // subtract reference solution
    vec_sol.axpy(vec_ref, -1.0);

    // compute the norm
    double nrm2 = vec_sol.norm2();

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(nrm2, 0.0, tol);
  }
};
CuSolverLUTest cusolverlu_test;

class CuSolverQRTest :
  public TestSystem::UnitTest
{
public:
  CuSolverQRTest() :
    TestSystem::UnitTest("CuSolverLUTest", "double", "unsigned int", PreferredBackend::cuda)
  {
  }

  virtual ~CuSolverQRTest()
  {
  }

  virtual void run() const override
  {
    typedef SparseMatrixCSR<double, unsigned int> MatrixType;
    typedef DenseVector<double, unsigned int> VectorType;

    const double tol = Math::pow(Math::eps<double>(), 0.6);

    // create a pointstar factory
    PointstarFactoryFD<double, unsigned int> psf(17);

    // create a CSR matrix
    MatrixType mat_sys;
    mat_sys.convert(psf.matrix_csr());

    // create a reference solution vector
    VectorType vec_ref;
    vec_ref.convert(psf.eigenvector_min());

    // create an rhs vector
    VectorType vec_rhs;
    vec_rhs.convert(mat_sys.create_vector_r());
    vec_rhs.scale(vec_ref, psf.lambda_min());

    // create an empty solution vector
    VectorType vec_sol;
    vec_sol.convert(mat_sys.create_vector_r());
    vec_sol.format();

    CuSolverQR qrsolver(mat_sys);
    qrsolver.apply(vec_sol, vec_rhs);

    // subtract reference solution
    vec_sol.axpy(vec_ref, -1.0);

    // compute the norm
    double nrm2 = vec_sol.norm2();

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(nrm2, 0.0, tol);
  }
};
CuSolverQRTest cusolverqr_test;

#endif //FEAT_HAVE_CUDA
