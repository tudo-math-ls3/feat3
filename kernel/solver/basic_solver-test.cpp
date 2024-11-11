// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/bicgstabl.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/rgcr.hpp>
#include <kernel/solver/pcr.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/sor_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
//#include <kernel/solver/spai_precond.hpp>
#include <kernel/solver/polynomial_precond.hpp>
#include <kernel/solver/matrix_precond.hpp>
#include <kernel/solver/pcgnr.hpp>
#include <kernel/solver/pcgnrilu.hpp>
#include <kernel/solver/idrs.hpp>
#include <kernel/solver/multigrid.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

template<
  typename DataType_,
  typename IndexType_>
  class BasicSolverTest :
  public UnitTest
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;
  typedef SparseMatrixCSR<DataType, IndexType> MatrixType;
  typedef typename MatrixType::VectorTypeR VectorType;
  typedef NoneFilter<DataType, IndexType> FilterType;

public:
  BasicSolverTest(PreferredBackend backend) :
    UnitTest("BasicSolverTest", Type::Traits<DataType>::name(), Type::Traits<IndexType>::name(), backend)
  {
  }

  virtual ~BasicSolverTest()
  {
  }

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs, const Index ref_iters, const Index iter_tol = 2u) const
  {
    #ifdef FEAT_HAVE_QUADMATH
      // a bit more tolerance for quad precision
      const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.4));
    #else
      const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.5));
    #endif

    // set solver name and plot summary
    std::cout << "\n";
    solver.set_plot_name(name);
    solver.set_plot_mode(PlotMode::summary);

    // initialize solver
    solver.init();

    // solve
    Status status = solver.apply(vec_sol, vec_rhs);
    TEST_CHECK_MSG(status_success(status), name + String(": apply failed with status = ") + stringify(status));

    // release solver
    solver.done();

    // check against reference solution
    vec_sol.axpy(vec_ref, -DataType(1));
    const DataType d = vec_sol.norm2sqr();
    TEST_CHECK_MSG(d <= tol, name + ": failed to reach tolerance\n"
      + "result: " + stringify_fp_sci(d) + "; expected result <= " + stringify(tol));

    // check number of iterations (only in double precision)
    if(sizeof(DataType) == 8u)
    {
      const Index n = solver.get_num_iter();
      TEST_CHECK_MSG((n <= ref_iters + iter_tol) && (n + iter_tol >= ref_iters),
        name + ": performed " + stringify(n) + " iterations; expected "
        + stringify(ref_iters) + " +/- " + stringify(iter_tol));
    }
  }

  virtual void run() const override
  {
    const Index m = 17;
    const Index d = 2;

    // create a pointstar factory
    PointstarFactoryFD<DataType, IndexType> psf(m, d);

    // create 5-point star CSR matrix
    SparseMatrixCSR<DataType, IndexType> csr_mat(psf.matrix_csr());

    // create a Q2 bubble vector
    DenseVector<DataType, IndexType> q2b_vec(psf.vector_q2_bubble());

    // create a NoneFilter
    FilterType filter;

    // convert to system matrix type
    MatrixType matrix;
    matrix.convert(csr_mat);

    // convert bubble vector
    VectorType vec_ref;
    vec_ref.convert(q2b_vec);

    // compute rhs vector
    VectorType vec_rhs(vec_ref.clone(CloneMode::Layout));
    matrix.apply(vec_rhs, vec_ref);

    // initialize sol vector
    VectorType vec_sol(vec_ref.clone(CloneMode::Layout));

    // test plain CG
    {
      auto solver = Solver::new_pcg(matrix, filter);
      test_solver("CG", *solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCG-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_pcg(matrix, filter, precon);
      test_solver("PCG-JAC", *solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCG-POLY(3)
    {
      auto precon = Solver::new_polynomial_precond(matrix, filter, 3);
      auto solver = Solver::new_pcg(matrix, filter, precon);
      test_solver("PCG-POLY(3)", *solver, vec_sol, vec_ref, vec_rhs, 11);
    }

    // test BICGSTAB-SSOR
    {
      auto precon = Solver::new_ssor_precond(this->get_preferred_backend(), matrix, filter);
      auto solver = Solver::new_bicgstab(matrix, filter, precon);
      test_solver("BICGSTAB-SSOR", *solver, vec_sol, vec_ref, vec_rhs, this->get_preferred_backend() == PreferredBackend::cuda ? 20 : 13);
    }

    // test plain CR
    {
      // create a CR solver
      auto solver = Solver::new_pcr(matrix, filter);
      test_solver("CR", *solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCR-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_pcr(matrix, filter, precon);
      solver->set_tol_rel(1e-7);
      test_solver("PCR-JAC", *solver, vec_sol, vec_ref, vec_rhs, 27);
    }

    // test PCR-SSOR
    if (this->get_preferred_backend() != PreferredBackend::cuda)
    {
      auto precon = Solver::new_ssor_precond(this->get_preferred_backend(), matrix, filter);
      auto solver = Solver::new_pcr(matrix, filter, precon);
      test_solver("PCR-SSOR", *solver, vec_sol, vec_ref, vec_rhs, 19);
    }

    // test BICGStab-left-SPAI
//     {
//       auto precon = Solver::new_spai_precond(matrix, filter);
//       auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::left);
//       test_solver("BiCGStab-left-SPAI", *solver, vec_sol, vec_ref, vec_rhs, 12);
//     }

    // test Richardson-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_richardson(matrix, filter, DataType(1.0), precon);
      solver->set_max_iter(1500);
      test_solver("Richardson-JAC", *solver, vec_sol, vec_ref, vec_rhs, 1175);
    }

    // test BiCGStab-right-SOR(1) aka GS
    {
      auto precon = Solver::new_sor_precond(this->get_preferred_backend(), matrix, filter, DataType(1));
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-right-SOR(1)", *solver, vec_sol, vec_ref, vec_rhs, Backend::get_preferred_backend()!=PreferredBackend::cuda ? 29 : 33);
    }

    // test BiCGStab-left-ILU(0)
    if (this->get_preferred_backend() != PreferredBackend::cuda)
    {
      auto precon = Solver::new_ilu_precond(this->get_preferred_backend(), matrix, filter, Index(0));
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-Left-ILU(0)", *solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStabL-ILU(0) L=1 -> BiCGStab
    if (this->get_preferred_backend() != PreferredBackend::cuda)
    {
      auto precon = Solver::new_ilu_precond(this->get_preferred_backend(), matrix, filter, Index(0));
      auto solver = Solver::new_bicgstabl(matrix, filter, 1, precon, BiCGStabLPreconVariant::left);
      test_solver("BiCGStabL(1)-left-ILU", *solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStabL-ILU(0) L=4
    if (this->get_preferred_backend() != PreferredBackend::cuda)
    {
      auto precon = Solver::new_ilu_precond(this->get_preferred_backend(), matrix, filter, Index(0));
      auto solver = Solver::new_bicgstabl(matrix, filter, 4, precon, BiCGStabLPreconVariant::right);
      test_solver("BiCGStabL(4)-right-ILU", *solver, vec_sol, vec_ref, vec_rhs, 3);
    }

    // test IDR(4)-ILU(0)
    if (this->get_preferred_backend() != PreferredBackend::cuda)
    {
      auto precon = Solver::new_ilu_precond(this->get_preferred_backend(), matrix, filter, Index(0));
      auto solver = Solver::new_idrs(matrix, filter, 4, precon);
      solver->reset_shadow_space(false);
      test_solver("IDR(4)-ILU(0)", *solver, vec_sol, vec_ref, vec_rhs, 20);
    }

    // test FGMRES-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_fgmres(matrix, filter, 16, DataType(0), precon);
      test_solver("FGMRES(16)-JAC", *solver, vec_sol, vec_ref, vec_rhs, 48);
    }

    // test GMRES-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_gmres(matrix, filter, 16, DataType(0), precon);
      test_solver("GMRES(16)-JAC", *solver, vec_sol, vec_ref, vec_rhs, 48);
    }

    // test PCG-jac-matrix
    /*{
      SparseMatrixCOO<:Main, DataType, IndexType> coo_jac(csr_mat.rows(), csr_mat.columns());
      for (Index i(0) ; i < csr_mat.rows() ; ++i)
      {
        coo_jac(i, i, DataType(1) / csr_mat(i,i));
      }
      MatrixType jac_matrix;
      jac_matrix.convert(coo_jac);

      auto precon = Solver::new_matrix_precond(jac_matrix, filter);
      auto solver = Solver::new_pcg(matrix, filter, precon);
      test_solver("PCG-MATRIX(JAC)", *solver, vec_sol, vec_ref, vec_rhs, 28);
    }*/

    // test PCGNR-Jac-Jac
    {
      auto precon_l = Solver::new_jacobi_precond(matrix, filter);
      auto precon_r = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_pcgnr(matrix, filter, precon_l, precon_r);
      solver->set_tol_rel(1e-7);
      test_solver("PCGNR-JAC-JAC", *solver, vec_sol, vec_ref, vec_rhs, 42);
    }

    // test PCGNRILU
    // if (this->get_preferred_backend() != PreferredBackend::cuda)
    {
      auto solver = Solver::new_pcgnrilu(matrix, filter, 0);
      test_solver("PCGNRILU", *solver, vec_sol, vec_ref, vec_rhs, 31);
    }

    // test RGCR-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_rgcr(matrix, filter, precon);
      test_solver("RGCR-JAC", *solver, vec_sol, vec_ref, vec_rhs, 28);
    }
    /*
   #if defined(FEAT_CUDA_VERSION_MAJOR) && (FEAT_CUDA_VERSION_MAJOR >= 8)
       // test FGMRES-SPAI
       {
         auto precon = Solver::new_spai_precond(matrix, filter);
         auto solver = Solver::new_fgmres(matrix, filter, 16, 0.0, precon);
         test_solver("FGMRES(16)-SPAI", *solver, vec_sol, vec_ref, vec_rhs, 32);
       }

       // test RGCR-SPAI
       {
         auto precon = Solver::new_spai_precond(matrix, filter);
         auto solver = Solver::new_rgcr(matrix, filter, precon);
         test_solver("RGCR-SPAI", *solver, vec_sol, vec_ref, vec_rhs, 17);
       }
   #endif*/

   // test BiCGStab-Jacobi
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter, DataType(0.5));
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::right);
      solver->set_tol_rel(1e-7);
      test_solver("BiCGStab-right-Jacobi(0.5)", *solver, vec_sol, vec_ref, vec_rhs, 19);
    }

    // test BiCGStab-SOR
    {
      auto precon = Solver::new_sor_precond(this->get_preferred_backend(), matrix, filter, DataType(1));
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-left-SOR", *solver, vec_sol, vec_ref, vec_rhs, Backend::get_preferred_backend()!=PreferredBackend::cuda ? 29 : 37);
    }

    // test BiCGStab-SSOR
    {
      auto precon = Solver::new_ssor_precond(this->get_preferred_backend(), matrix, filter);
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-right-SSOR", *solver, vec_sol, vec_ref, vec_rhs, Backend::get_preferred_backend()!=PreferredBackend::cuda ? 13 : 20);
    }
  }
};

BasicSolverTest <double, std::uint32_t> basic_solver_test_double_uint32(PreferredBackend::generic);
BasicSolverTest <double, std::uint64_t> basic_solver_test_double_uint64(PreferredBackend::generic);
//BasicSolverTest<float, unsigned int> basic_solver_test_float_uint(PreferredBackend::generic);
//BasicSolverTest<float, unsigned long> basic_solver_test_float_ulong(PreferredBackend::generic);

#ifdef FEAT_HAVE_MKL
//BasicSolverTest<float, unsigned long> mkl_basic_solver_test_float_ulong(PreferredBackend::mkl);
BasicSolverTest <double, std::uint64_t> mkl_basic_solver_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
BasicSolverTest <__float128, std::uint32_t> basic_solver_test_float128_uint32(PreferredBackend::generic);
BasicSolverTest <__float128, std::uint64_t> basic_solver_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
//BasicSolverTest<Half, unsigned int> basic_solver_test_half_uint(PreferredBackend::generic);
//BasicSolverTest<Half, unsigned long> basic_solver_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
//BasicSolverTest<float, unsigned int> cuda_basic_solver_test_float_uint(PreferredBackend::cuda);
BasicSolverTest <double, std::uint32_t> cuda_basic_solver_test_double_uint32(PreferredBackend::cuda);
//BasicSolverTest<float, unsigned long> cuda_basic_solver_test_float_ulong(PreferredBackend::cuda);
//BasicSolverTest<double, unsigned long> cuda_basic_solver_test_double_ulong(PreferredBackend::cuda);
#endif

/*
template<
  template<typename,typename,typename> class ScalarMatrix_,
  typename MemType_,
  typename DataType_,
  typename IndexType_>
class CUDASolverTest :
  public FullTaggedTest<MemType_, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;
  typedef ScalarMatrix_<MemType_, DataType, IndexType> MatrixType;
  typedef typename MatrixType::VectorTypeR VectorType;
  typedef NoneFilter<MemType_, DataType, IndexType> FilterType;

public:
  CUDASolverTest() :
    FullTaggedTest<MemType_, DataType, IndexType>("CUDASolverTest-" + MatrixType::name())
  {
  }

  virtual ~CUDASolverTest()
  {
  }

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs, const Index ref_iters, const Index iter_tol = 2u) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.5));

    // set solver name and plot summary
    std::cout << "\n";
    solver.set_plot_name(name);
    solver.set_plot_mode(PlotMode::summary);

    // initialize solver
    solver.init();

    // solve
    Status status = solver.apply(vec_sol, vec_rhs);
    TEST_CHECK_MSG(status_success(status), name + String(": apply failed with status = ") + stringify(status));

    // release solver
    solver.done();

    // check against reference solution
    vec_sol.axpy(vec_ref, vec_sol, -DataType(1));
    const DataType d = vec_sol.norm2sqr();
    TEST_CHECK_MSG(d <= tol, name + ": failed to reach tolerance\n"
      + "result: " + stringify_fp_sci(d) + "; expected result <= " + stringify(tol));

    // check number of iterations
    const Index n = solver.get_num_iter();
    TEST_CHECK_MSG((n <= ref_iters+iter_tol) && (n+iter_tol >= ref_iters),
      name + ": performed " + stringify(n) + " iterations; expected "
      + stringify(ref_iters) + " +/- " + stringify(iter_tol));
  }

  virtual void run() const override
  {
    const Index m = 17;
    const Index d = 2;

    // create a pointstar factory
    PointstarFactoryFD<DataType, IndexType> psf(m, d);

    // create 5-point star CSR matrix
    SparseMatrixCSR<Mem::Main, DataType, IndexType> csr_mat(psf.matrix_csr());

    // create a Q2 bubble vector
    DenseVector<Mem::Main, DataType, IndexType> q2b_vec(psf.vector_q2_bubble());

    // create a NoneFilter
    FilterType filter;

    // convert to system matrix type
    MatrixType matrix;
    matrix.convert(csr_mat);

    // convert bubble vector
    VectorType vec_ref;
    vec_ref.convert(q2b_vec);

    // compute rhs vector
    VectorType vec_rhs(vec_ref.clone(CloneMode::Layout));
    matrix.apply(vec_rhs, vec_ref);

    // initialize sol vector
    VectorType vec_sol(vec_ref.clone(CloneMode::Layout));

    // test plain CG
    {
      auto solver = Solver::new_pcg(matrix, filter);
      test_solver("CG", *solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCG-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_pcg(matrix, filter, precon);
      test_solver("PCG-JAC", *solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCG-POLY
    {
      auto precon = Solver::new_polynomial_precond(matrix, filter, 3);
      auto solver = Solver::new_pcg(matrix, filter, precon);
      test_solver("PCG-POLY", *solver, vec_sol, vec_ref, vec_rhs, 11);
    }

    // test Richardson-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_richardson(matrix, filter, DataType(1.0), precon);
      solver->set_max_iter(1500);
      test_solver("Richardson-JAC", *solver, vec_sol, vec_ref, vec_rhs, 1175);
    }

    // test BiCGStab-ILU(0)
    {
      auto precon = Solver::new_ilu_precond(matrix, filter);
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-ILU", *solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStabL-ILU(0) L=1 -> BiCGStab
    {
      auto precon = Solver::new_ilu_precond(matrix, filter);
      auto solver = Solver::new_bicgstabl(matrix, filter, 1, precon, BiCGStabLPreconVariant::left);
      test_solver("BiCGStabL(1)-ILU", *solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStabL-ILU(0) L=4
    {
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      auto solver = Solver::new_bicgstabl(matrix, filter, 4, precon, BiCGStabLPreconVariant::right);
      test_solver("BiCGStabL(4)-ILU", *solver, vec_sol, vec_ref, vec_rhs, 3);
    }
    // test IDR(4)-ILU(0)

    {
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      auto solver = Solver::new_idrs(matrix, filter, 4, precon);
      solver->reset_shadow_space(false);
      test_solver("IDR(4)-ILU(0)", *solver, vec_sol, vec_ref, vec_rhs, 20);
    }

    // test FGMRES-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_fgmres(matrix, filter, 16, DataType(0), precon);
      test_solver("FGMRES(16)-JAC", *solver, vec_sol, vec_ref, vec_rhs, 48);
    }

    // test PCG-Matrix(jac)
    {
      SparseMatrixCOO<Mem::Main, DataType, IndexType> coo_jac(csr_mat.rows(), csr_mat.columns());
      for (Index i(0) ; i < csr_mat.rows() ; ++i)
      {
        coo_jac(i, i, DataType(1) / csr_mat(i,i));
      }
      MatrixType jac_matrix;
      jac_matrix.convert(coo_jac);

      auto precon = Solver::new_matrix_precond(jac_matrix, filter);
      auto solver = Solver::new_pcg(matrix, filter, precon);
      test_solver("PCG-MATRIX(JAC)", *solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCGNR-Jac-Jac
    {
      auto precon_l = Solver::new_jacobi_precond(matrix, filter);
      auto precon_r = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_pcgnr(matrix, filter, precon_l, precon_r);
      test_solver("PCGNR-JAC-JAC", *solver, vec_sol, vec_ref, vec_rhs, 43);
    }

#if defined(FEAT_CUDA_VERSION_MAJOR) && (FEAT_CUDA_VERSION_MAJOR >= 8)
    // test FGMRES-SPAI
    {
      auto precon = Solver::new_spai_precond(matrix, filter);
      auto solver = Solver::new_fgmres(matrix, filter, 16, 0.0, precon);
      test_solver("FGMRES(16)-SPAI", *solver, vec_sol, vec_ref, vec_rhs, 32);
    }

    // test RGCR-SPAI
    {
      auto precon = Solver::new_spai_precond(matrix, filter);
      auto solver = Solver::new_rgcr(matrix, filter, precon);
      test_solver("RGCR-SPAI", *solver, vec_sol, vec_ref, vec_rhs, 17);
    }
#endif


#ifdef FEAT_HAVE_CUSOLVER
    // test BiCGStab-Jacobi
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter, DataType(0.5));
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-right-Jacobi(0.5)", *solver, vec_sol, vec_ref, vec_rhs, 21);
    }

    // test BiCGStab-SOR
    {
      auto precon = Solver::new_sor_precond(matrix, filter, DataType(1));
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-left-SOR", *solver, vec_sol, vec_ref, vec_rhs, 37);
    }

    // test BiCGStab-SSOR
    {
      auto precon = Solver::new_ssor_precond(matrix, filter);
      auto solver = Solver::new_bicgstab(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-right-SSOR", *solver, vec_sol, vec_ref, vec_rhs, 20);
    }
#endif // FEAT_HAVE_CUSOLVER
  }
};

#ifdef FEAT_HAVE_CUDA
//CUDASolverTest<SparseMatrixCSR, Mem::CUDA, double, std::uint64_t> cuda_solver_csr_generic_double_uint64;
CUDASolverTest<SparseMatrixCSR, Mem::CUDA, double, std::uint32_t> cuda_solver_csr_generic_double_uint32;
#endif

template<
  template<typename,typename,typename> class ScalarMatrix_,
  typename MemType_,
  typename DataType_,
  typename IndexType_>
class BandedSolverTest :
  public FullTaggedTest<MemType_, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;
  typedef ScalarMatrix_<MemType_, DataType, IndexType> MatrixType;
  typedef typename MatrixType::VectorTypeR VectorType;
  typedef NoneFilter<MemType_, DataType, IndexType> FilterType;

public:
  BandedSolverTest() :
    FullTaggedTest<MemType_, DataType, IndexType>("BandedSolverTest-" + MatrixType::name())
  {
  }

  virtual ~BandedSolverTest()
  {
  }

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs, const Index ref_iters) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.5));

    // set solver name and plot summary
    std::cout << "\n";
    solver.set_plot_name(name);
    solver.set_plot_mode(PlotMode::summary);

    // initialize solver
    solver.init();

    // solve
    Status status = solver.apply(vec_sol, vec_rhs);
    TEST_CHECK_MSG(status_success(status), (String("Failed to solve: '") + name + ("'")));

    // release solver
    solver.done();

    // check against reference solution
    vec_sol.axpy(vec_ref, vec_sol, -DataType(1));
    DataType d = vec_sol.norm2sqr();
    TEST_CHECK_EQUAL_WITHIN_EPS(d, DataType(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(solver.get_num_iter(), ref_iters, 3);
  }

  virtual void run() const override
  {
    // create band matrix from pointstar structure fe
    std::vector<IndexType> num_of_nodes;
    num_of_nodes.push_back(13);
    num_of_nodes.push_back(13);
    SparseMatrixBanded<Mem::Main, DataType, IndexType> bm(PointstarStructureFE::template value<DataType>(1, num_of_nodes));
    for (Index i(0) ; i < bm.get_elements_size().at(0) ; ++i)
      bm.val()[i] = DataType_(-1);
    for (Index i(4 * bm.rows()) ; i < 5 * bm.rows() ; ++i)
      bm.val()[i] = DataType_(4);

    // create a Q2 bubble vector
    DenseVector<Mem::Main, DataType, IndexType> q2b_vec(bm.rows());
    for (Index i (0) ; i < q2b_vec.size() ; ++i)
      q2b_vec(i, DataType(i%10) / DataType(100));

    // create a NoneFilter
    FilterType filter;

    // convert to system matrix type
    MatrixType matrix;
    matrix.convert(bm);

    // convert bubble vector
    VectorType vec_ref;
    vec_ref.convert(q2b_vec);

    // compute rhs vector
    VectorType vec_rhs(vec_ref.clone(CloneMode::Layout));
    matrix.apply(vec_rhs, vec_ref);

    // initialize sol vector
    VectorType vec_sol(vec_ref.clone(CloneMode::Layout));

    // test plain CG
    {
      auto solver = Solver::new_pcg(matrix, filter);
      solver->set_max_iter(200);
      test_solver("CG", *solver, vec_sol, vec_ref, vec_rhs, 105);
    }

    // test PCG-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      auto solver = Solver::new_pcg(matrix, filter, precon);
      solver->set_max_iter(200);
      test_solver("PCG-JAC", *solver, vec_sol, vec_ref, vec_rhs, 105);
    }
  }
};

BandedSolverTest<SparseMatrixBanded, Mem::Main, double, Index> banded_solver_csr_generic_double_index;


template<typename MemType_, typename DataType_, typename IndexType_>
class BCSRSolverTest :
  public TestSystem::FullTaggedTest<MemType_, DataType_, IndexType_>
{
public:
  static constexpr int dim = 2;
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;
  typedef LAFEM::SparseMatrixBCSR<MemType_, DataType, IndexType, dim, dim> MatrixType;
  typedef LAFEM::DenseVectorBlocked<MemType_, DataType, IndexType, dim> VectorType;
  typedef LAFEM::NoneFilterBlocked<MemType_, DataType, IndexType, dim> FilterType;

public:
  BCSRSolverTest() :
    TestSystem::FullTaggedTest<MemType_, DataType, IndexType>("BCSRSolverTest")
  {
  }

  virtual ~BCSRSolverTest()
  {
  }

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs, const Index ref_iters, const Index iter_tol = 2u) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.5));

    // set solver name and plot summary
    std::cout << "\n";
    solver.set_plot_name(name);
    solver.set_plot_mode(PlotMode::summary);

    // initialize solver
    solver.init();

    // solve
    Status status = solver.apply(vec_sol, vec_rhs);
    TEST_CHECK_MSG(status_success(status), name + String(": apply failed with status = ") + stringify(status));

    // release solver
    solver.done();

    // check against reference solution
    vec_sol.axpy(vec_ref, vec_sol, -DataType(1));
    const DataType d = vec_sol.norm2sqr();
    TEST_CHECK_MSG(d <= tol, name + ": failed to reach tolerance\n"
      + "result: " + stringify_fp_sci(d) + "; expected result <= " + stringify(tol));

    // check number of iterations
    const Index n = solver.get_num_iter();
    TEST_CHECK_MSG((n <= ref_iters+iter_tol) && (n+iter_tol >= ref_iters),
      name + ": performed " + stringify(n) + " iterations; expected "
      + stringify(ref_iters) + " +/- " + stringify(iter_tol));
  }

  virtual void run() const override
  {
    const Index m = 7;
    const Index d = 2;

    // create a pointstar factory
    LAFEM::PointstarFactoryFD<DataType, IndexType> psf(m, d);

    // create 5-point star CSR matrix
    LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> csr_mat(psf.matrix_csr());

    // create a Q2 bubble vector
    LAFEM::DenseVector<Mem::Main, DataType, IndexType> q2b_vec(psf.vector_q2_bubble());

    // construct blocked matrix in main memory
    LAFEM::SparseMatrixBCSR<Mem::Main, DataType, IndexType, dim, dim> matrix_main(csr_mat.layout());
    {
      const IndexType num_rows = IndexType(csr_mat.rows());
      const IndexType* row_ptr = csr_mat.row_ptr();
      const IndexType* col_idx = csr_mat.col_ind();
      const auto* data_a = csr_mat.val();
      auto* data_b = matrix_main.val();

      // create the following matrix:
      // ( A  0 )
      // (-I  A )
      for(IndexType i(0); i < num_rows; ++i)
      {
        for(IndexType j(row_ptr[i]); j < row_ptr[i+1]; ++j)
        {
          data_b[j][0][0] = data_b[j][1][1] = data_a[j];
          data_b[j][0][1] = DataType(0);
          data_b[j][1][0] = DataType(col_idx[j] == i ? -1 : 0);
        }
      }
    }

    // construct blocked vector in main memory
    LAFEM::DenseVectorBlocked<Mem::Main, DataType, IndexType, dim> vec_ref_main(q2b_vec.size());
    {
      const Index n = q2b_vec.size();
      const auto* data_x = q2b_vec.elements();
      auto* data_y = vec_ref_main.elements();

      for(Index i(0); i < n; ++i)
      {
        data_y[i][0] = DataType(0);
        data_y[i][1] = data_x[i];
      }
    }

    // create filter
    FilterType filter;

    // convert matrix
    MatrixType matrix;
    matrix.convert(matrix_main);

    // convert reference vector
    VectorType vec_ref;
    vec_ref.convert(vec_ref_main);

    // compute rhs vector
    VectorType vec_rhs(vec_ref.clone(LAFEM::CloneMode::Layout));
    matrix.apply(vec_rhs, vec_ref);

    // initialize sol vector
    VectorType vec_sol(vec_ref.clone(LAFEM::CloneMode::Layout));

    // test Richardson-ILU
    {
      auto precon = Solver::new_ilu_precond(matrix, filter);
      auto solver = Solver::new_richardson(matrix, filter, DataType(0.9), precon);
      test_solver("Richardson-ILU", *solver, vec_sol, vec_ref, vec_rhs, 43);
    }

    // test Richardson-SOR
    {
      auto precon = Solver::new_sor_precond(matrix, filter, DataType(1.2));
      auto solver = Solver::new_richardson(matrix, filter, DataType(0.9), precon);
      test_solver("Richardson-SOR", *solver, vec_sol, vec_ref, vec_rhs, 83);
    }

    // test Richardson-SSOR
    {
      auto precon = Solver::new_ssor_precond(matrix, filter, DataType(1));
      auto solver = Solver::new_richardson(matrix, filter, DataType(0.9), precon);
      test_solver("Richardson-SSOR", *solver, vec_sol, vec_ref, vec_rhs, 71);
    }
  }
};
BCSRSolverTest<Mem::Main, double, Index> bcsr_solver_test_main_double_index; */
