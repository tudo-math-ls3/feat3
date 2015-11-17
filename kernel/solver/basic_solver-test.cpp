#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/sor_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/spai_precond.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::Solver;
using namespace FEAST::TestSystem;

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

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs) const
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
  }

  virtual void run() const
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
      test_solver("CG", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test PCG-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-JAC", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test PCG-SSOR
    {
      // create a SSOR preconditioner
      auto precon = Solver::new_ssor_precond(matrix, filter);
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-SSOR", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test FGMRES-SPAI
    {
      // create an SPAI preconditioner
      auto precon = Solver::new_spai_precond(matrix, filter, matrix.layout());
      // create a FMGRES solver
      FGMRES<MatrixType, FilterType> solver(matrix, filter, 16, 0.0, precon);
      test_solver("FGMRES(16)-SPAI", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test Richardson-SOR
    {
      // create a SOR preconditioner
      auto precon = Solver::new_sor_precond(matrix, filter, DataType(1.7));
      // create a Richardson solver
      Richardson<MatrixType, FilterType> solver(matrix, filter, DataType(1.0), precon);
      solver.set_max_iter(1000);
      test_solver("Richardson-SOR(1.7)", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test BiCGStab-ILU(0)
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("BiCGStab-ILU(0)", solver, vec_sol, vec_ref, vec_rhs);
    }
  }
};

BasicSolverTest<SparseMatrixCSR, Mem::Main, double, Index> basic_solver_csr_generic_double_index;
BasicSolverTest<SparseMatrixELL, Mem::Main, double, Index> basic_solver_ell_generic_double_index;

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

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs) const
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
  }

  virtual void run() const
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
      test_solver("CG", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test PCG-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-JAC", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test Richardson-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a Richardson solver
      Richardson<MatrixType, FilterType> solver(matrix, filter, DataType(1.0), precon);
      solver.set_max_iter(1000);
      test_solver("Richardson-JAC", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test BiCGStab-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("BiCGStab-JAC", solver, vec_sol, vec_ref, vec_rhs);
    }
  }
};

#ifdef FEAST_BACKENDS_CUDA
CUDASolverTest<SparseMatrixCSR, Mem::CUDA, double, unsigned long> cuda_solver_csr_generic_double_ulong;
CUDASolverTest<SparseMatrixELL, Mem::CUDA, double, unsigned long> cuda_solver_ell_generic_double_ulong;
CUDASolverTest<SparseMatrixCSR, Mem::CUDA, double, unsigned int> cuda_solver_csr_generic_double_uint;
CUDASolverTest<SparseMatrixELL, Mem::CUDA, double, unsigned int> cuda_solver_ell_generic_double_uint;
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

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs) const
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
  }

  virtual void run() const
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
      test_solver("PCG", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test PCG-JAC
    {
      // create a Jacobi preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      // create a CG solver
      PCG<MatrixType, FilterType> solver(matrix, filter, precon);
      test_solver("PCG-JAC", solver, vec_sol, vec_ref, vec_rhs);
    }
  }
};

BandedSolverTest<SparseMatrixBanded, Mem::Main, double, Index> banded_solver_csr_generic_double_index;
