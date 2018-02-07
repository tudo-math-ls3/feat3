#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/bicgstabl.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/rgcr.hpp>
#include <kernel/solver/pcr.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/sor_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/spai_precond.hpp>
#include <kernel/solver/polynomial_precond.hpp>
#include <kernel/solver/matrix_precond.hpp>
#include <kernel/solver/pcgnr.hpp>
#include <kernel/solver/pcgnrilu.hpp>
#include <kernel/solver/idrs.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

template<
  template<typename,typename,typename> class ScalarMatrix_,
  typename MemType_,
  typename DataType_,
  typename IndexType_>
class BasicSolverTest :
  public FullTaggedTest<MemType_, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;
  typedef ScalarMatrix_<MemType_, DataType, IndexType> MatrixType;
  typedef typename MatrixType::VectorTypeR VectorType;
  typedef NoneFilter<MemType_, DataType, IndexType> FilterType;

public:
  BasicSolverTest() :
    FullTaggedTest<MemType_, DataType, IndexType>("BasicSolverTest-" + MatrixType::name())
  {
  }

  virtual ~BasicSolverTest()
  {
  }

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs, const Index ref_iters) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.5));

    // initialise solver
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
    TEST_CHECK_EQUAL_WITHIN_EPS(solver.get_num_iter(), ref_iters, 2);
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

    // initialise sol vector
    VectorType vec_sol(vec_ref.clone(CloneMode::Layout));

    // test plain CG
    {
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter);
      test_solver("CG", solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCG-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a PCG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-JAC", solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCG-POLY
    {
      auto precon = Solver::new_polynomial_precond(matrix, filter, 3);
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-POLY", solver, vec_sol, vec_ref, vec_rhs, 11);
    }

    // test PCG-SSOR
    {
      // create a SSOR preconditioner
      auto precon = Solver::new_ssor_precond(matrix, filter);
      // create a PCG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-SSOR", solver, vec_sol, vec_ref, vec_rhs, 19);
    }

    // test plain CR
    {
      // create a CR solver
      PCR<MatrixType, FilterType> solver(matrix, filter);
      test_solver("CR", solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCR-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a PCR solver
      PCR<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCR-JAC", solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCR-SSOR
    {
      // create a SSOR preconditioner
      auto precon = Solver::new_ssor_precond(matrix, filter);
      // create a PCR solver
      PCR<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCR-SSOR", solver, vec_sol, vec_ref, vec_rhs, 19);
    }

    // test FGMRES-SPAI
    {
      // create an SPAI preconditioner
      auto precon = Solver::new_spai_precond(matrix, filter);
      // create a FMGRES solver
      FGMRES<MatrixType, FilterType> solver(matrix, filter, 16, 0.0, precon);
      test_solver("FGMRES(16)-SPAI", solver, vec_sol, vec_ref, vec_rhs, 32);
    }

    // test Richardson-SOR
    {
      // create a SOR preconditioner
      auto precon = Solver::new_sor_precond(matrix, filter, DataType(1.7));
      // create a Richardson solver
      Richardson<MatrixType, FilterType> solver(matrix, filter, DataType(1.0), precon);
      solver.set_max_iter(1000);
      test_solver("Richardson-SOR(1.7)", solver, vec_sol, vec_ref, vec_rhs, 71);
    }

    // test BiCGStab-left-ILU(0)
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-ILU(0)", solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStabL-ILU(0) L=1 -> BiCGStab
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      // create a BiCGStab solver
      BiCGStabL<MatrixType, FilterType> solver(matrix, filter, 1, precon, BiCGStabLPreconVariant::left);
      test_solver("BiCGStab-ILU", solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStabL-ILU(0) L=4
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      // create a BiCGStab solver
      BiCGStabL<MatrixType, FilterType> solver(matrix, filter, 4, precon, BiCGStabLPreconVariant::right);
      test_solver("BiCGStab-ILU", solver, vec_sol, vec_ref, vec_rhs, 3);
    }


    // test BiCGStab-right-SOR(1) aka GS
    {
      auto precon = Solver::new_sor_precond(matrix, filter, DataType(1));
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-right-SOR(1)", solver, vec_sol, vec_ref, vec_rhs, 29);
    }

    // test IDR(4)-ILU(0)
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      // create a IDRS solver
      IDRS<MatrixType, FilterType> solver(matrix, filter, 4, precon);
      test_solver("IDR(4)-ILU(0)", solver, vec_sol, vec_ref, vec_rhs, 20);
    }
    // test PCG-jac-matrix
    {
      SparseMatrixCOO<Mem::Main, DataType, IndexType> coo_jac(csr_mat.rows(), csr_mat.columns());
      for (Index i(0) ; i < csr_mat.rows() ; ++i)
      {
        coo_jac(i, i, DataType(1) / csr_mat(i,i));
      }
      MatrixType jac_matrix;
      jac_matrix.convert(coo_jac);
      // create a Matrix preconditioner
      auto precon = Solver::new_matrix_precond(jac_matrix, filter);
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-JAC", solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCGNR-Jac-Jac
    {
      auto precon_l = Solver::new_jacobi_precond(matrix, filter);
      auto precon_r = Solver::new_jacobi_precond(matrix, filter);
      PCGNR<MatrixType, FilterType> solver(matrix, filter, precon_l, precon_r);
      test_solver("PCGNR-JAC-JAC", solver, vec_sol, vec_ref, vec_rhs, 43);
    }

    // test PCGNRILU
    {
      PCGNRILU<MatrixType, FilterType> solver(matrix, filter, 0);
      test_solver("PCGNRILU", solver, vec_sol, vec_ref, vec_rhs, 31);
    }

    // test RGCR-JAC
    {
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      RGCR<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("RGCR-JAC", solver, vec_sol, vec_ref, vec_rhs, 28);
    }
  }
};

BasicSolverTest<SparseMatrixCSR, Mem::Main, double, unsigned int> basic_solver_csr_generic_double_uint;
BasicSolverTest<SparseMatrixELL, Mem::Main, double, unsigned int> basic_solver_ell_generic_double_uint;
BasicSolverTest<SparseMatrixCSR, Mem::Main, double, unsigned long> basic_solver_csr_generic_double_ulong;
BasicSolverTest<SparseMatrixELL, Mem::Main, double, unsigned long> basic_solver_ell_generic_double_ulong;

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
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs, const Index ref_iters) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.5));

    // initialise solver
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
    TEST_CHECK_EQUAL_WITHIN_EPS(solver.get_num_iter(), ref_iters, 2);
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

    // initialise sol vector
    VectorType vec_sol(vec_ref.clone(CloneMode::Layout));

    // test plain CG
    {
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter);
      test_solver("CG", solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCG-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-JAC", solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCG-POLY
    {
      auto precon = Solver::new_polynomial_precond(matrix, filter, 3);
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-POLY", solver, vec_sol, vec_ref, vec_rhs, 11);
    }

    // test Richardson-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a Richardson solver
      Richardson<MatrixType, FilterType> solver(matrix, filter, DataType(1.0), precon);
      solver.set_max_iter(2000);
      test_solver("Richardson-JAC", solver, vec_sol, vec_ref, vec_rhs, 1175);
    }

    // test BiCGStab-ILU(0)
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter);
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-ILU", solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStabL-ILU(0) L=1 -> BiCGStab
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter);
      // create a BiCGStab solver
      BiCGStabL<MatrixType, FilterType> solver(matrix, filter, 1, precon, BiCGStabLPreconVariant::left);
      test_solver("BiCGStab-ILU", solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStabL-ILU(0) L=4
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      // create a BiCGStab solver
      BiCGStabL<MatrixType, FilterType> solver(matrix, filter, 4, precon, BiCGStabLPreconVariant::right);
      test_solver("BiCGStab-ILU", solver, vec_sol, vec_ref, vec_rhs, 3);
    }
    // test IDR(4)-ILU(0)

    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      // create a IDRS solver
      IDRS<MatrixType, FilterType> solver(matrix, filter, 4, precon);
      test_solver("IDR(4)-ILU(0)", solver, vec_sol, vec_ref, vec_rhs, 20);
    }

    // test FGMRES-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a Richardson solver
      FGMRES<MatrixType, FilterType> solver(matrix, filter, 16, DataType(0), precon);
      solver.set_max_iter(2000);
      test_solver("FGMRES-JAC", solver, vec_sol, vec_ref, vec_rhs, 48);
    }

    // test PCG-jac-matrix
    {
      SparseMatrixCOO<Mem::Main, DataType, IndexType> coo_jac(csr_mat.rows(), csr_mat.columns());
      for (Index i(0) ; i < csr_mat.rows() ; ++i)
      {
        coo_jac(i, i, DataType(1) / csr_mat(i,i));
      }
      MatrixType jac_matrix;
      jac_matrix.convert(coo_jac);
      // create a Matrix preconditioner
      auto precon = Solver::new_matrix_precond(jac_matrix, filter);
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-JAC", solver, vec_sol, vec_ref, vec_rhs, 28);
    }

    // test PCGNR-Jac-Jac
    {
      auto precon_l = Solver::new_jacobi_precond(matrix, filter);
      auto precon_r = Solver::new_jacobi_precond(matrix, filter);
      PCGNR<MatrixType, FilterType> solver(matrix, filter, precon_l, precon_r);
      test_solver("PCGNR-JAC-JAC", solver, vec_sol, vec_ref, vec_rhs, 43);
    }

    // test FGMRES-SPAI
    {
      // create an SPAI preconditioner
      auto precon = Solver::new_spai_precond(matrix, filter);
      // create a FMGRES solver
      FGMRES<MatrixType, FilterType> solver(matrix, filter, 16, 0.0, precon);
      test_solver("FGMRES(16)-SPAI", solver, vec_sol, vec_ref, vec_rhs, 32);
    }

    // test RGCR-SPAI
    {
      auto precon = Solver::new_spai_precond(matrix, filter);
      RGCR<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("RGCR-SPAI", solver, vec_sol, vec_ref, vec_rhs, 17);
    }



#ifdef FEAT_HAVE_CUSOLVER
    // test BiCGStab-SOR
    {
      // create a SOR preconditioner
      auto precon = Solver::new_sor_precond(matrix, filter, DataType(1));
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-left-SOR", solver, vec_sol, vec_ref, vec_rhs, 37);
    }

    // test BiCGStab-SSOR
    {
      // create a SSOR preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter, DataType(0.5));
      // create a CG solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-right-Jacobi(0.5)", solver, vec_sol, vec_ref, vec_rhs, 21);
    }

    // test BiCGStab-SSOR
    {
      // create a SSOR preconditioner
      auto precon = Solver::new_ssor_precond(matrix, filter);
      // create a CG solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-right-SSOR", solver, vec_sol, vec_ref, vec_rhs, 20);
    }
#endif // FEAT_HAVE_CUSOLVER
  }
};

#ifdef FEAT_HAVE_CUDA
//CUDASolverTest<SparseMatrixCSR, Mem::CUDA, double, unsigned long> cuda_solver_csr_generic_double_ulong;
//CUDASolverTest<SparseMatrixELL, Mem::CUDA, double, unsigned long> cuda_solver_ell_generic_double_ulong;
CUDASolverTest<SparseMatrixCSR, Mem::CUDA, double, unsigned int> cuda_solver_csr_generic_double_uint;
//CUDASolverTest<SparseMatrixELL, Mem::CUDA, double, unsigned int> cuda_solver_ell_generic_double_uint;
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

    // initialise solver
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
    TEST_CHECK_EQUAL_WITHIN_EPS(solver.get_num_iter(), ref_iters, 2);
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

    // initialise sol vector
    VectorType vec_sol(vec_ref.clone(CloneMode::Layout));

    // test plain PCG
    {
      // create a PCG solver
      PCG<MatrixType, FilterType> solver(matrix, filter);
      solver.set_max_iter(1000);
      test_solver("PCG", solver, vec_sol, vec_ref, vec_rhs, 105);
    }

    // test PCG-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      solver.set_max_iter(1000);
      test_solver("PCG-JAC", solver, vec_sol, vec_ref, vec_rhs, 105);
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
  void test_solver(String name, std::shared_ptr<Solver_> solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs, Index ref_iters) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.5));

    // initialise solver
    solver->init();

    // solve
    Solver::Status status = solver->apply(vec_sol, vec_rhs);
    TEST_CHECK_MSG(Solver::status_success(status), (String("Failed to solve: '") + name + ("'")));

    // release solver
    solver->done();

    // check against reference solution
    vec_sol.axpy(vec_ref, vec_sol, -DataType(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_sol.norm2sqr(), DataType(0), tol);

    // check for iteration count
    TEST_CHECK_EQUAL_WITHIN_EPS(solver->get_num_iter(), ref_iters, 2);
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

    // initialise sol vector
    VectorType vec_sol(vec_ref.clone(LAFEM::CloneMode::Layout));

    // test Richardson-ILU
    {
      auto precon = Solver::new_ilu_precond(matrix, filter);
      auto solver = Solver::new_richardson(matrix, filter, DataType(0.9), precon);
      solver->set_max_iter(1000);
      test_solver("Richardson-ILU", solver, vec_sol, vec_ref, vec_rhs, 43);
    }

    // test Richardson-SOR
    {
      auto precon = Solver::new_sor_precond(matrix, filter, DataType(1.2));
      auto solver = Solver::new_richardson(matrix, filter, DataType(0.9), precon);
      solver->set_max_iter(1000);
      test_solver("Richardson-SOR", solver, vec_sol, vec_ref, vec_rhs, 83);
    }

    // test Richardson-SSOR
    {
      auto precon = Solver::new_ssor_precond(matrix, filter, DataType(1));
      auto solver = Solver::new_richardson(matrix, filter, DataType(0.9), precon);
      solver->set_max_iter(1000);
      test_solver("Richardson-SSOR", solver, vec_sol, vec_ref, vec_rhs, 71);
    }
  }
};

BCSRSolverTest<Mem::Main, double, Index> bcsr_solver_test_main_double_index;
