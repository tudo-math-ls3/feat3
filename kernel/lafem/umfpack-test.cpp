#include <test_system/test_system.hpp>
#ifdef FEAST_HAVE_UMFPACK
#include <kernel/lafem/umfpack.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

class UmfpackTest :
  public TestSystem::TaggedTest<Mem::Main, double>
{
public:
  typedef SparseMatrixCSR<Mem::Main, double> MatrixType;
  typedef DenseVector<Mem::Main, double> VectorType;

  UmfpackTest() : TestSystem::TaggedTest<Mem::Main, double>("UmfpackTest") {}

  virtual void run() const
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
      // create an UMFPACK solver
      Umfpack umfpack(&mat_sys);
      // solve
      umfpack.solve(vec_sol, vec_rhs);
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
#endif // FEAST_HAVE_UMFPACK
