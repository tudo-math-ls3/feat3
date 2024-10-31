// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/solver/direct_stokes_solver.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

template<typename DT_, typename IT_>
class DirectStokesSolverTest :
  public TestSystem::UnitTest
{
public:
  explicit DirectStokesSolverTest(PreferredBackend backend) :
    TestSystem::UnitTest("DirectStokesSolverTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.6));

    const IT_ m = 8;
    const IT_ num_cells = m*m;
    const IT_ num_edges = 2*m*(m+1);
    const IT_ off_edgev = m*(m+1);
    const DT_ h = DT_(1) / DT_(m);

    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, 2, 2> MatrixA;
    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, 2, 1> MatrixB;
    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, 1, 2> MatrixD;
    //typedef LAFEM::UnitFilterBlocked<DT_, IT_, 2> VeloFilter;
    typedef LAFEM::NoneFilterBlocked<DT_, IT_, 2> VeloFilter;
    typedef LAFEM::MeanFilter<DT_, IT_> PresFilter;
    typedef LAFEM::TupleFilter<VeloFilter, PresFilter> Filter;

    LAFEM::SaddlePointMatrix<MatrixA, MatrixB, MatrixD> matrix;

    MatrixA& matrix_a = matrix.block_a();
    MatrixB& matrix_b = matrix.block_b();
    MatrixD& matrix_d = matrix.block_d();

    // create div(V) matrix
    matrix_d = MatrixD(num_cells, num_edges, 4*num_cells);
    {
      IT_* row_ptr = matrix_d.row_ptr();
      IT_* col_idx = matrix_d.col_ind();
      auto* val   = matrix_d.val();
      for(IT_ i = 0, k = 0; i < num_cells; ++i)
      {
        row_ptr[i] = k;
        // bottom edge
        col_idx[k] = i;
        val[k](0,0) = DT_(0);
        val[k](0,1) = -h;
        ++k;
        // top edge
        col_idx[k] = i + m;
        val[k](0,0) = DT_(0);
        val[k](0,1) = h;
        ++k;
        // left edge
        col_idx[k] = off_edgev + (i/m)*(m+1) + (i%m);
        val[k](0,0) = -h;
        val[k](0,1) = DT_(0);
        ++k;
        // right edge
        col_idx[k] = off_edgev + (i/m)*(m+1) + (i%m) + 1;
        val[k](0,0) = h;
        val[k](0,1) = DT_(0);
        ++k;
      }
      row_ptr[num_cells] = 4*num_cells;
    }

    // create grad(P) matrix
    matrix_b = matrix_d.transpose();

    // create A matrix
    // note: A represents the bilinearform a(v,psi) = dx^2 v1 + dy^2 v2, which is *NOT* the Laplace operator,
    // because it misses the dy^2 v1 and dx^2 v2 terms; this is just to keep the assembly simpler...
    // this is just a unit test, so it will suffice for this cause
    matrix_a = MatrixA(Adjacency::Graph(Adjacency::RenderType::injectify_sorted, matrix.block_b(), matrix.block_d()));
    matrix_a.format(0.0);
    {
      const IT_* row_ptr = matrix_a.row_ptr();
      const IT_* col_idx = matrix_a.col_ind();
      auto* val = matrix_a.val();
      const DT_ ih2 = DT_(1) / (h*h);
      for(IT_ i = 0; i < off_edgev; ++i)
      {
        for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
        {
          if(col_idx[j] < off_edgev)
            val[j].add_scalar_main_diag((i == col_idx[j] ? DT_(2) : DT_(-1))*ih2);
        }
      }
      for(IT_ i = off_edgev; i < num_edges; ++i)
      {
        for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
        {
          if(col_idx[j] >= off_edgev)
            val[j].add_scalar_main_diag((i == col_idx[j] ? DT_(2) : DT_(-1))*ih2);
        }
      }
    }

    // create vectors vector
    auto vec_rhs = matrix.create_vector_r();
    auto vec_def = matrix.create_vector_r();
    auto vec_sol = matrix.create_vector_r();
    auto vec_ref = matrix.create_vector_r();
    vec_ref.format();
    vec_rhs.format();
    vec_sol.format();
    {
      auto* v = vec_ref.template at<0>().elements();
      auto* p = vec_ref.template at<1>().elements();
      const DT_ pi = Math::pi<DT_>();
      // horizontal edges
      for(IT_ i = 0; i < m*(m+1); ++i)
      {
        const DT_ x = (DT_(i%m)+0.5)/DT_(m);
        const DT_ y = DT_(i/m)/DT_(m);
        v[i][0] = Math::sin(pi*x) * Math::cos(pi*y);
        //v[i][1] = -Math::cos(pi*x) * Math::sin(pi*y);
      }
      // vertical edges
      for(IT_ i = 0; i < m*(m+1); ++i)
      {
        const DT_ x = DT_(i%(m+1))/DT_(m);
        const DT_ y = (DT_(i/(m+1))+0.5)/DT_(m);
        //v[i][0] = Math::sin(pi*x) * Math::cos(pi*y);
        v[i][1] = -Math::cos(pi*x) * Math::sin(pi*y);
      }
      // cells
      for(IT_ i = 0; i < num_cells; ++i)
      {
        const DT_ x = (DT_(i%m)+0.5)/DT_(m);
        const DT_ y = (DT_(i/m)+0.5)/DT_(m);
        p[i] = Math::sqr(Math::cos(pi*x)) + Math::sqr(Math::cos(pi*y)) - DT_(1);
      }
    }

    // create a pressure mean filter
    Filter filter;
    filter.template at<1>() = PresFilter(LAFEM::DenseVector<DT_, IT_>(num_cells, 1.0), LAFEM::DenseVector<DT_, IT_>(num_cells, 1.0));

    // compute RHS from reference solution
    matrix.apply(vec_rhs, vec_ref);

    // solve
    auto solver = Solver::new_direct_stokes_solver(matrix, filter);

    // print backend and selected solver
    std::cout << "Selected Backend: " << Backend::get_preferred_backend() << "\n";
    std::cout << "Selected Solver: " << solver->name() << "\n";

    solver->init();
    solver->apply(vec_sol, vec_rhs);
    solver->done();

    // compute defect and error
    matrix.apply(vec_def, vec_sol, vec_rhs, -1.0);
    vec_ref.axpy(vec_sol, -1.0);

    std::cout << "|rhs| = " << stringify_fp_sci(vec_rhs.max_abs_element()) << "\n";
    std::cout << "|sol| = " << stringify_fp_sci(vec_sol.max_abs_element()) << "\n";
    std::cout << "|def| = " << stringify_fp_sci(vec_def.max_abs_element()) << "\n";
    std::cout << "|err| = " << stringify_fp_sci(vec_ref.max_abs_element()) << "\n";

    TEST_CHECK_EQUAL_WITHIN_EPS(vec_def.max_abs_element(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_ref.max_abs_element(), DT_(0), tol);
  }
};

#ifdef FEAT_HAVE_UMFPACK
DirectStokesSolverTest<double, std::uint32_t> direct_stokes_solver_test_umfpack_double_u32(PreferredBackend::generic);
DirectStokesSolverTest<double, std::uint64_t> direct_stokes_solver_test_umfpack_double_u64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DirectStokesSolverTest<double, std::uint32_t> direct_stokes_solver_test_mkldss_double_u32(PreferredBackend::mkl);
DirectStokesSolverTest<double, std::uint64_t> direct_stokes_solver_test_mkldss_double_u64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_CUDSS
DirectStokesSolverTest<double, std::uint32_t> direct_stokes_solver_test_cudss_double_u32(PreferredBackend::cuda);
DirectStokesSolverTest<double, std::uint64_t> direct_stokes_solver_test_cudss_double_u64(PreferredBackend::cuda);
#endif
