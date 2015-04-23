#include <test_system/test_system.hpp>
#ifdef FEAST_HAVE_CUSOLVER
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/cusolver.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

class CuSolverLUTest :
  public TestSystem::TaggedTest<Mem::CUDA, double>
{
public:
  CuSolverLUTest() : TestSystem::TaggedTest<Mem::CUDA, double>("CuSolverQRTest") {}

  virtual void run() const
  {
    typedef SparseMatrixCSR<Mem::Main, double, unsigned int> MatrixType;
    typedef DenseVector<Mem::Main, double, unsigned int> VectorType;

    const double tol = Math::pow(Math::eps<double>(), 0.6);

    // create a pointstar factory
    PointstarFactoryFD<double> psf(17);

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

    CuSolverLU::solve(vec_sol, mat_sys, vec_rhs);

    // subtract reference solution
    vec_sol.axpy(vec_ref, vec_sol, -1.0);

    // compute the norm
    double nrm2 = vec_sol.norm2();

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(nrm2, 0.0, tol);
  }
};
CuSolverLUTest cusolverlu_test;

class CuSolverQRTest :
  public TestSystem::TaggedTest<Mem::CUDA, double>
{
public:
  CuSolverQRTest() : TestSystem::TaggedTest<Mem::CUDA, double>("CuSolverLUTest") {}

  virtual void run() const
  {
    typedef SparseMatrixCSR<Mem::CUDA, double, unsigned int> MatrixType;
    typedef DenseVector<Mem::CUDA, double, unsigned int> VectorType;

    const double tol = Math::pow(Math::eps<double>(), 0.6);

    // create a pointstar factory
    PointstarFactoryFD<double> psf(17);

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

    CuSolverQR::solve(vec_sol, mat_sys, vec_rhs);

    // subtract reference solution
    vec_sol.axpy(vec_ref, vec_sol, -1.0);

    // compute the norm
    double nrm2 = vec_sol.norm2();

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(nrm2, 0.0, tol);
  }
};
CuSolverQRTest cusolverqr_test;

#endif // FEAST_HAVE_CUSOLVER
