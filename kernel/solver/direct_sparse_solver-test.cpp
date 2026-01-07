// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>

#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/solver/direct_sparse_solver.hpp>

using namespace FEAT;

template<typename IT_>
inline IT_ isqrt(const IT_ n)
{
  IT_ m = IT_(Math::sqrt(double(n)));
  return (m*m == n ? m : IT_(0));
}

/**
 * \brief Simple test for DirectSparseSolver class
 *
 * \author Peter Zajac
 */
template<typename DT_, typename IT_>
class DirectSparseSolverTest :
  public TestSystem::UnitTest
{
  typedef DT_ DataType;
  typedef IT_ IndexType;

  typedef LAFEM::VectorMirror<DataType, IndexType> MirrorType;
  typedef LAFEM::DenseVector<DataType, IndexType> LocalVectorType;
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalMatrixType;
  typedef LAFEM::UnitFilter<DataType, IndexType> LocalFilterType;

  typedef Global::Gate<LocalVectorType, MirrorType> GateType;
  typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;
  typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;
  typedef Global::Filter<LocalFilterType, MirrorType> GlobalFilterType;

  String _allowed_backends;

public:
  DirectSparseSolverTest() :
    TestSystem::UnitTest("DirectSparseSolverTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name()),
    _allowed_backends()
  {
  }

  explicit DirectSparseSolverTest(const String& allowed_backends) :
    TestSystem::UnitTest("DirectSparseSolverTest [" + allowed_backends + "]", Type::Traits<DT_>::name(), Type::Traits<IT_>::name()),
    _allowed_backends(allowed_backends)
  {
  }

  static MirrorType create_mirror_0(const int n, const int k)
  {
    MirrorType mirror((Index)n, 1);
    mirror.indices()[0] = Index(k);
    return mirror;
  }

  static MirrorType create_mirror_1(const int n, const int m, const int o, const int p)
  {
    MirrorType mirror((Index)n, (Index)m);
    IndexType* idx = mirror.indices();
    for(int i(0); i < m; ++i)
      idx[Index(i)] = Index(o + i * p);
    return mirror;
  }

  static void create_gate(const int np, const int m, GateType& gate, std::vector<MirrorType>& bnds)
  {
    const Dist::Comm& comm = *gate.get_comm();

    // get our process (i,j) coords
    const int ii = comm.rank() / np;
    const int jj = comm.rank() % np;

    // add mirrors for our vertex neighbors

    // lower-left neighbor?
    if((ii > 0) && (jj > 0))
      gate.push((ii-1)*np + (jj-1), create_mirror_0(m*m, 0));
    else
      bnds.push_back(create_mirror_0(m*m, 0));

    // lower-right neighbor?
    if((ii > 0) && (jj+1 < np))
      gate.push((ii-1)*np + (jj+1), create_mirror_0(m*m, m-1));
    else
      bnds.push_back(create_mirror_0(m*m, m-1));

    // upper-left neighbor?
    if((ii+1 < np) && (jj > 0))
      gate.push((ii+1)*np + (jj-1), create_mirror_0(m*m, m*(m-1)));
    else
      bnds.push_back(create_mirror_0(m*m, m*(m-1)));

    // upper-right neighbor?
    if((ii+1 < np) && (jj+1 < np))
      gate.push((ii+1)*np + (jj+1), create_mirror_0(m*m, m*m-1));
    else
      bnds.push_back(create_mirror_0(m*m, m*m-1));

    // add mirror for our edge neighbors

    // lower neighbor?
    if(ii > 0)
      gate.push((ii-1)*np + jj, create_mirror_1(m*m, m, 0, 1));
    else
      bnds.push_back(create_mirror_1(m*m, m, 0, 1));

    // upper neighbor?
    if(ii+1 < np)
      gate.push((ii+1)*np + jj, create_mirror_1(m*m, m, m*(m-1), 1));
    else
      bnds.push_back(create_mirror_1(m*m, m, m*(m-1), 1));

    // left neighbor?
    if(jj > 0)
      gate.push(ii*np + (jj-1), create_mirror_1(m*m, m, 0, m));
    else
      bnds.push_back(create_mirror_1(m*m, m, 0, m));

    // right neighbor?
    if(jj+1 < np)
      gate.push(ii*np + (jj+1), create_mirror_1(m*m, m, m-1, m));
    else
      bnds.push_back(create_mirror_1(m*m, m, m-1, m));

    // compile gate
    gate.compile(LocalVectorType(Index(m*m)));
  }

  void run_test_parallel() const
  {
    static const DT_ tol = TestSystem::tol<DT_>();

    const Dist::Comm comm = Dist::Comm::world();

    // get number of processes in each direction
    const int np = isqrt(comm.size());
    if(np*np != comm.size())
    {
      std::cout << "Comm size is " << comm.size() << " which is not a square number; exiting test\n";
      return; // number of procs is not square
    }

    // Check if the selected solver is UMFPACK, because this only works for 1 process
    if((comm.size() > 1) && (_allowed_backends == "umfpack"))
    {
      std::cout << "UMFPACK backend is only available for 1 process; exiting test\n";
      return;
    }

    // distributed cuDSS is not available on Windows (yet) due to missing communication layer
#ifdef _WIN32
    if((comm.size() > 1) && (_allowed_backends == "cudss"))
    {
      std::cout << "cuDSS backend on Windows platform is only available for 1 process; exiting test\n";
      return;
    }
#endif

    // number of vertices in each direction (at least 3)
    const Index m = Index(Math::max((16 / np) + 1, 3));

    // vector of boundary mirrors
    std::vector<MirrorType> bnds;

    // create gate
    GateType gate(comm);
    create_gate(np, int(m), gate, bnds);

    // create global matrix
    GlobalMatrixType matrix(&gate, &gate);

    // assemble local matrix structure and initialize to random values
    LAFEM::PointstarFactoryFE<DataType, IndexType> psf(m);
    matrix.local() = psf.matrix_csr_neumann();

    // create global filter
    GlobalFilterType filter(m*m);
    for(const auto& mir : bnds)
    {
      const Index n = mir.num_indices();
      const IndexType* idx = mir.indices();
      for(Index i(0); i < n; ++i)
        filter.local().add(idx[i], DT_(0));
    }

    // >>>>> DEBUG >>>>>
    /*
    {
      String s;
      int n1 = gate._ranks.size();
      for(int i(0); i < n1; ++i)
      {
        s += stringify(gate._ranks.at(i)) + ":";
        auto* vi = gate._mirrors.at(i).indices();
        int l = gate._mirrors.at(i).num_indices();
        for(int j(0); j < l; ++j)
          (s += " ") += stringify(vi[j]);
        s += "\n";
      }
      s += "F:";
      auto n2 = filter.local().used_elements();
      auto* x = filter.local().get_indices();
      for(Index i(0); i < n2; ++i)
        (s += " ") += stringify(x[i]);
      comm.allprint(s);
    }
    */
    // <<<<< DEBUG <<<<<

    // create two vectors
    GlobalVectorType vec_rhs = matrix.create_vector_l();
    GlobalVectorType vec_sol = matrix.create_vector_l();
    GlobalVectorType vec_ref = matrix.create_vector_l();

    // initialize reference solution vector
    {
      const IT_ ii = IT_(comm.rank()) / IT_(np);
      const IT_ jj = IT_(comm.rank()) % IT_(np);

      const IT_ n = IT_(isqrt(vec_ref.local().size()));
      DataType* v = vec_ref.local().elements();

      // global size in one dimension; remember 1 DOF overlap
      const IT_ gn = IT_(np) * (n - IT_(1)) + IT_(1);

      // scaling factor for sine bubble = sin(pi*x)*sin(pi*y)
      const DT_ sc = Math::pi<DT_>() / DT_(gn - IT_(1));

      for(IT_ i(0); i < n; ++i)
      {
        const DT_ sy = Math::sin(sc * DT_(ii*(n-IT_(1)) + i));
        for(IT_ j(0); j < n; ++j)
        {
          v[i*n+j] = sy * Math::sin(sc * DT_(jj*(n-IT_(1)) + j));
        }
      }
    }

    // filter reference solution for correct Dirichlet BCs
    filter.filter_sol(vec_ref);

    // compute RHS vector from solution
    matrix.apply(vec_rhs, vec_ref);

    // filter RHS for correct Dirichlet BCs
    filter.filter_rhs(vec_rhs);

    // create direct sparse solver
    auto solver = Solver::new_direct_sparse_solver(matrix, filter);

    // push allowed backends, if given
    if(!_allowed_backends.empty())
    {
      std::cout << "[Parallel] Allowed DirectSparseSolver backends: " << _allowed_backends << "\n";
      solver->push_allowed_backend_list(_allowed_backends);
    }

    // initialize
    solver->init();

    // print the selected backend
    std::cout << "[Parallel] Selected DirectSparseSolver backend: " << stringify(solver->get_selected_backend()) << "\n";

    // solve, release and destroy
    solver->apply(vec_sol, vec_rhs);
    solver->done();
    solver.reset();

    // check difference to reference solution
    TEST_CHECK_LESS_THAN(vec_sol.max_rel_diff(vec_ref), tol);
  }

  void run_test_serial() const
  {
    const DataType tol = TestSystem::tol<DataType>();

    // create a pointstar factory
    LAFEM::PointstarFactoryFD<DataType, IndexType> psf(17);

    // create a CSR matrix
    LocalMatrixType matrix(psf.matrix_csr());

    // create a reference solution vector
    LocalVectorType vec_ref(psf.eigenvector_min());

    // create an rhs vector
    LocalVectorType vec_rhs(matrix.create_vector_r());
    vec_rhs.scale(vec_ref, psf.lambda_min());

    // create an empty solution vector
    LocalVectorType vec_sol(matrix.create_vector_r());
    vec_sol.format();

    LAFEM::NoneFilter<DataType, IndexType> filter;

    // create direct sparse solver
    auto solver = Solver::new_direct_sparse_solver(matrix, filter);

    // push allowed backends, if given
    if(!_allowed_backends.empty())
    {
      std::cout << "[Serial] Allowed DirectSparseSolver backends: " << _allowed_backends << "\n";
      solver->push_allowed_backend_list(_allowed_backends);
    }

    // initialize, solve, release
    solver->init();

    // print the selected backend
    std::cout << "[Serial] Selected DirectSparseSolver backend: " << stringify(solver->get_selected_backend()) << "\n";

    // solve, release and destroy
    solver->apply(vec_sol, vec_rhs);
    solver->done();
    solver.reset();

    // check difference to reference solution
    TEST_CHECK_LESS_THAN(vec_sol.max_rel_diff(vec_ref), tol);
  }

  virtual void run() const override
  {
    run_test_serial();
    run_test_parallel();
  }
}; // class DirectSparseSolverTest<...>

#ifdef FEAT_HAVE_UMFPACK
DirectSparseSolverTest<double, Index> direct_sparse_solver_test_double_index_umfpack("umfpack");
#endif

#ifdef FEAT_HAVE_SUPERLU_DIST
DirectSparseSolverTest<double, Index> direct_sparse_solver_test_double_index_superlu("superlu");
#endif

#ifdef FEAT_HAVE_CUDSS
DirectSparseSolverTest<double, Index> direct_sparse_solver_test_double_index_cudss("cudss");
#endif

#ifdef FEAT_HAVE_MKL
DirectSparseSolverTest<double, Index> direct_sparse_solver_test_double_index_mkldss("mkldss");
#endif

#ifdef FEAT_HAVE_MUMPS
DirectSparseSolverTest<double, Index> direct_sparse_solver_test_double_index_mumps("mumps");
#endif

void direct_sparse_solver_test_dummy()
{
  // just to keep linkers from complaining
}
