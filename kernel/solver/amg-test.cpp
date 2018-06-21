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
#include <kernel/solver/amg.hpp>
#include <kernel/solver/multigrid.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

template<
  template<typename,typename,typename> class ScalarMatrix_,
  typename MemType_,
  typename DataType_,
  typename IndexType_>
class AMGTest :
  public FullTaggedTest<MemType_, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;
  typedef ScalarMatrix_<MemType_, DataType, IndexType> MatrixType;
  typedef typename MatrixType::VectorTypeR VectorType;
  typedef UnitFilter<MemType_, DataType, IndexType> FilterType;
  typedef Transfer<MatrixType> TransferType;

  template <typename MT_, typename FT_, typename TT_>
  struct Level
  {
    MT_ matrix;
    FT_ filter;
    TT_ transfer;

    Level(MT_ & m, FT_ & f, TT_ & t)
    {
      matrix.clone(m);
      filter.clone(f);
      transfer = t.clone();
      //transfer.clone(t);
    }
  };

public:
  AMGTest() :
    FullTaggedTest<MemType_, DataType, IndexType>("AMGTest-" + MatrixType::name())
  {
  }

  virtual ~AMGTest()
  {
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

    // create a UnitFilter
    FilterType filter(q2b_vec.size());

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
    vec_sol.format();

    std::deque<std::shared_ptr<Level<MatrixType, FilterType, TransferType>>> levels;
    {
      TransferType transfer;
      auto l = std::make_shared<Level<MatrixType, FilterType, TransferType>>(matrix, filter, transfer);
      levels.push_back(l);
    }

    while (levels.back()->matrix.rows() > 25)
    {
      MatrixType coarse_matrix;
      FilterType coarse_filter;
      TransferType transfer;
      Dist::Comm comm = Dist::Comm::self();
      AMGFactory<MatrixType, FilterType, TransferType>::new_coarse_level(levels.back()->matrix, levels.back()->filter, 0.8, coarse_matrix, coarse_filter, levels.back()->transfer, &comm);
      auto l = std::make_shared<Level<MatrixType, FilterType, TransferType>>(coarse_matrix, coarse_filter, transfer);
      levels.push_back(l);
    }
    AMGFactory<MatrixType, FilterType, TransferType>::update_coarse_level((*(++levels.rbegin()))->matrix, (*(++levels.rbegin()))->transfer, levels.back()->matrix);

    //std::reverse(levels.begin(), levels.end());

    auto multigrid_hierarchy = std::make_shared<Solver::MultiGridHierarchy<MatrixType, FilterType, TransferType>>(levels.size());

    for(std::size_t i(0); (i+1) < levels.size(); ++i)
    {
      auto& lvl = *levels.at(i);

      auto jacobi = Solver::new_jacobi_precond(lvl.matrix, lvl.filter);

      auto smoother = Solver::new_richardson(lvl.matrix, lvl.filter, DataType(0.8), jacobi);

      smoother->set_max_iter(4);
      smoother->set_min_iter(4);

      multigrid_hierarchy->push_level(
        lvl.matrix,     // the system matrix for this level
        lvl.filter,     // the system filter for this level
        lvl.transfer,   // the transfer operator for this level
        smoother,       // the pre-smoother
        smoother,       // the post-smoother
        smoother        // the peak-smoother
      );
    }

    {
      auto& lvl = *levels.back();

      auto coarse_solver = Solver::new_pcg(lvl.matrix, lvl.filter);

      multigrid_hierarchy->push_level(
        lvl.matrix,       // the coarse-level system matrix
        lvl.filter,       // the coarse-level system filter
        coarse_solver     // the coarse-level solver
      );
    }

    auto multigrid = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

    auto& lvl_fine = *levels.front();
    auto solver = Solver::new_richardson(lvl_fine.matrix, lvl_fine.filter, DataType(1), multigrid);

    multigrid_hierarchy->init();
    solver->init();
    solver->set_plot_mode(Solver::PlotMode::iter);
    //std::cout<<vec_sol<<vec_rhs<<std::endl;
    Solver::solve(*solver, vec_sol, vec_rhs, lvl_fine.matrix, lvl_fine.filter);
    solver->done();
    multigrid_hierarchy->done();

    const Index ref_iters = 5;
    const Index iter_tol = 1;
    const Index n = solver->get_num_iter();
    TEST_CHECK_MSG((n <= ref_iters+iter_tol) && (n+iter_tol >= ref_iters),
      "AMG: performed " + stringify(n) + " iterations; expected "
      + stringify(ref_iters) + " +/- " + stringify(iter_tol));
  }
};
AMGTest<SparseMatrixCSR, Mem::Main, double, unsigned long> amg_csr_generic_double_ulong;
