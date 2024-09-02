// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>

#ifdef FEAT_HAVE_MKL
#include <kernel/solver/mkl_dss.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

class MKLDSSTest :
  public TestSystem::UnitTest
{
public:
  typedef SparseMatrixCSR<double, Index> MatrixType;
  typedef DenseVector<double, Index> VectorType;

  MKLDSSTest() :
    TestSystem::UnitTest("MKLDSSTest")
  {
  }

  virtual ~MKLDSSTest()
  {
  }

  virtual void run() const override
  {
    const double tol = Math::pow(Math::eps<double>(), 0.6);

    // set backend to MKL
    Backend::set_preferred_backend(PreferredBackend::mkl);

    // create a pointstar factory
    PointstarFactoryFD<double, Index> psf(3);

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
      // create a solver object
      MKLDSS solver(mat_sys);
      // initialize
      solver.init();
      // solve
      solver.apply(vec_sol, vec_rhs);
      // release
      solver.done();
    }

    // subtract reference solution
    vec_sol.axpy(vec_ref, -1.0);

    // compute the norm
    double nrm2 = vec_sol.norm2();

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(nrm2, 0.0, tol);
  }
};

MKLDSSTest mkl_dss_test;

#endif // FEAT_HAVE_MKL
