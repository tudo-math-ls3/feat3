#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/proto_solver.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  template<typename,typename,typename> class ScalarMatrix_,
  typename Algo_,
  typename DataType_,
  typename IndexType_>
class ProtoSolverTest :
  public FullTaggedTest<typename Algo_::MemType, Algo_, DataType_, IndexType_>
{
public:
  typedef Algo_ AlgoType;
  typedef typename AlgoType::MemType MemType;
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;
  typedef ScalarMatrix_<MemType, DataType, IndexType> MatrixType;
  typedef typename MatrixType::VectorTypeR VectorType;
  typedef NoneFilter<MemType, DataType, IndexType> FilterType;

public:
  ProtoSolverTest() :
    FullTaggedTest<MemType, AlgoType, DataType, IndexType>("ProtoSolverTest")
  {
  }

  template<typename Solver_>
  void test_solver(String name, Solver_& solver, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.5));

    // initialise solver
    bool okay = solver.init();
    TEST_CHECK_MSG(okay, (String("Failed to initialise solver: '") + name + ("'")));

    // solve
    SolverStatus status = solver.solve(vec_sol, vec_rhs);
    TEST_CHECK_MSG(status_success(status), (String("Failed to solve: '") + name + ("'")));

    // release solver
    solver.done();

    // check against reference solution
    vec_sol.template axpy<AlgoType>(vec_ref, vec_sol, -DataType(1));
    DataType d = vec_sol.template norm2sqr<AlgoType>();
    TEST_CHECK_EQUAL_WITHIN_EPS(d, DataType(0), tol);
  }

#ifdef FEAST_HAVE_UMFPACK
  template<typename IT_>
  void test_umfpack(const SparseMatrixCSR<Mem::Main, double, IT_>& matrix, VectorType& vec_sol, const VectorType& vec_ref, const VectorType& vec_rhs) const
  {
    std::cout << "Testing UMFPACK..." << std::endl;
    UmfpackSolver solver(matrix);
    test_solver("Umfpack", solver, vec_sol, vec_ref, vec_rhs);
  }
#endif // FEAST_HAVE_UMFPACK

  template<typename MT_>
  void test_umfpack(const MT_&, VectorType&, const VectorType&, const VectorType&) const
  {
    // nothing to do
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
    matrix.template apply<AlgoType>(vec_rhs, vec_ref);

    // initialise sol vector
    VectorType vec_sol(vec_ref.clone(CloneMode::Layout));

    // test plain CG
    {
      // create a CG solver
      PCGSolver<AlgoType, MatrixType, FilterType> solver(matrix, filter);
      test_solver("CG", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test PCG-SSOR
    {
      // create a SSOR preconditioner
      PreconWrapper<AlgoType, MatrixType, LAFEM::SSORPreconditioner> precon(matrix);
      // create a CG solver
      PCGSolver<AlgoType, MatrixType, FilterType> solver(matrix, filter, &precon);
      test_solver("PCG-SSOR", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test FGMRES-ILU
    {
      // create an SPAI preconditioner
      PreconWrapper<AlgoType, MatrixType, SPAIPreconditioner> precon(matrix, matrix.layout());
      // create a fix-point solver
      FGMRESSolver<AlgoType, MatrixType, FilterType> solver(matrix, filter, 16, 0.0, &precon);
      test_solver("FGMRES(16)-SPAI", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test Fix-Point-SOR
    {
      // create a SOR preconditioner
      PreconWrapper<AlgoType, MatrixType, SORPreconditioner> precon(matrix, DataType(1.7));
      // create a fix-point solver
      FixPointSolver<AlgoType, MatrixType, FilterType> solver(matrix, filter, &precon);
      solver.set_max_iter(1000);
      test_solver("FixPoint-SOR(1.7)", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test BiCGStab-ILU(0)
    {
      // create a ILU(0) preconditioner
      PreconWrapper<AlgoType, MatrixType, ILUPreconditioner> precon(matrix, Index(0));
      // create a BiCGStab solver
      BiCGStabSolver<AlgoType, MatrixType, FilterType> solver(matrix, filter, &precon);
      test_solver("BiCGStab-ILU(0)", solver, vec_sol, vec_ref, vec_rhs);
    }

    // test UMFPACK
    test_umfpack(matrix, vec_sol, vec_ref, vec_rhs);
  }
};

ProtoSolverTest<SparseMatrixCSR, Algo::Generic, double, Index> proto_solver_csr_generic_double_index;
ProtoSolverTest<SparseMatrixELL, Algo::Generic, double, Index> proto_solver_ell_generic_double_index;
