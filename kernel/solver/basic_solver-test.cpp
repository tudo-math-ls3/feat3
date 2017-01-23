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
#include <kernel/solver/rgcr.hpp>
#include <kernel/solver/pcr.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/sor_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/spai_precond.hpp>
#include <kernel/solver/matrix_precond.hpp>
#include <kernel/solver/pcgnr.hpp>
#include <kernel/solver/pcgnrilu.hpp>

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
    //const Index nrows(289);
    //SparseMatrixCSR<Mem::Main, DataType, IndexType> csr_mat(nrows, nrows, nrows);
    //IndexType* row_ptr(csr_mat.row_ptr());
    //IndexType* col_ind(csr_mat.col_ind());
    //DataType* val(csr_mat.val());
    //for(Index i(0); i < nrows; ++i)
    //{
    //  row_ptr[i] = i;
    //  col_ind[i] = i;
    //  val[i] = DataType(2);
    //}
    //row_ptr[nrows] = Index(nrows);

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

    //// test plain CG
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

    // test BiCGStab-ILU(0) (left)
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_ilu_precond(matrix, filter, Index(0));
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-ILU(0)", solver, vec_sol, vec_ref, vec_rhs, 12);
    }

    // test BiCGStab-Jacobi (right)
    {
      // create a ILU(0) preconditioner
      auto precon = Solver::new_jacobi_precond(matrix, filter, DataType(0.5));
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-Jacobi(0.5)", solver, vec_sol, vec_ref, vec_rhs, 33);
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
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-ILU", solver, vec_sol, vec_ref, vec_rhs, 12);
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
      auto precon = Solver::new_sor_precond(matrix, filter);
      // create a BiCGStab solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-SOR(left)", solver, vec_sol, vec_ref, vec_rhs, 33);
    }

    // test BiCGStab-SSOR
    {
      // create a SSOR preconditioner
      auto precon = Solver::new_ssor_precond(matrix, filter);
      // create a CG solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::left);
      test_solver("BiCGStab-SSOR(left)", solver, vec_sol, vec_ref, vec_rhs, 21);
    }

    // test BiCGStab-SSOR
    {
      // create a SSOR preconditioner
      auto precon = Solver::new_ssor_precond(matrix, filter);
      // create a CG solver
      BiCGStab<MatrixType, FilterType> solver(matrix, filter, precon, BiCGStabPreconVariant::right);
      test_solver("BiCGStab-SSOR(right)", solver, vec_sol, vec_ref, vec_rhs, 38);
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
