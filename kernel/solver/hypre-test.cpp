// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>

#if defined(FEAT_HAVE_HYPRE)

#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/solver/hypre.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/pcg.hpp>

using namespace FEAT;

inline int isqrt(const int n)
{
  int m = int(Math::sqrt(double(n)));
  return (m*m == n ? m : 0);
}

/**
 * \brief Simple test for HYPRE solver wrappers
 *
 * This test ensures that the solvers "borrowed" from HYPRE compile and run.
 * Please note that this test does not check the solver convergence, as highly
 * depends on the number of MPI ranks this test is executed with.
 *
 * \author Peter Zajac
 */
template<typename Mem_, typename DT_, typename IT_>
class HypreTest :
  public TestSystem::FullTaggedTest<Mem_, DT_, IT_>
{
  typedef Mem_ MemType;
  typedef DT_ DataType;
  typedef IT_ IndexType;

  typedef LAFEM::VectorMirror<MemType, DataType, IndexType> MirrorType;
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> LocalVectorType;
  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> LocalMatrixType;
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> LocalFilterType;

  typedef Global::Gate<LocalVectorType, MirrorType> GateType;
  typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;
  typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;
  typedef Global::Filter<LocalFilterType, MirrorType> GlobalFilterType;

public:
  HypreTest() :
    TestSystem::FullTaggedTest<Mem_, DT_, IT_>("HypreTest")
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

  static LocalFilterType create_filter(const int m, const std::vector<MirrorType>& bnds)
  {
    LocalFilterType filter((Index)(m*m));

    for(const auto& mir : bnds)
    {
      const Index n = mir.num_indices();
      const IndexType* idx = mir.indices();
      for(Index i(0); i < n; ++i)
        filter.add(idx[i], DT_(0));
    }
    return filter;
  }

  virtual void run() const override
  {
    const Dist::Comm comm = Dist::Comm::world();

    // get number of processes in each direction
    const int np = int(Math::sqrt(double(comm.size())));
    if(np*np != comm.size())
    {
      std::cout << "Comm size is " << comm.size() << " which is not a square number; exiting test" << std::endl;
      return; // number of procs is not square
    }

    // number of vertices in each direction (at least 3)
    const int m = Math::max((16 / np) + 1, 3);

    // vector of boundary mirrors
    std::vector<MirrorType> bnds;

    // create gate
    GateType gate(comm);
    create_gate(np, m, gate, bnds);

    // create global matrix
    GlobalMatrixType matrix(&gate, &gate);

    // assemble local matrix structure and initialize to random values
    LAFEM::PointstarFactoryFE<DataType, IndexType> psf((Index)m);
    matrix.local() = psf.matrix_csr_neumann();

    // create global filter
    GlobalFilterType filter(create_filter(m, bnds));

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

    // "assemble" and filter rhs vector
    vec_rhs.format(1.0);
    filter.filter_rhs(vec_rhs);

    // test solvers
    test_parasails(matrix, filter, vec_rhs, vec_sol);
    test_euclid(matrix, filter, vec_rhs, vec_sol);
    test_boomeramg(matrix, filter, vec_rhs, vec_sol);
  }

  void test_parasails(const GlobalMatrixType& matrix, const GlobalFilterType& filter,
    const GlobalVectorType& vec_rhs, GlobalVectorType& vec_sol) const
  {
    vec_sol.format();

    auto precond = Solver::new_parasails_precond(matrix, filter);
    auto solver = Solver::new_richardson(matrix, filter, 0.5, precond);

    solver->set_max_iter(5);
    solver->set_tol_rel(1E-5);
    solver->set_plot_name("HYPRE:ParaSails");
    solver->set_plot_mode(Solver::PlotMode::all);

    solver->init();
    solver->correct(vec_sol, vec_rhs);
    solver->done();
  }

  void test_euclid(const GlobalMatrixType& matrix, const GlobalFilterType& filter,
    const GlobalVectorType& vec_rhs, GlobalVectorType& vec_sol) const
  {
    vec_sol.format();

    auto precond = Solver::new_euclid_precond(matrix, filter);
    auto solver = Solver::new_richardson(matrix, filter, 0.9, precond);

    solver->set_max_iter(5);
    solver->set_tol_rel(1E-5);
    solver->set_plot_name("HYPRE:Euclid");
    solver->set_plot_mode(Solver::PlotMode::all);

    solver->init();
    solver->correct(vec_sol, vec_rhs);
    solver->done();
  }

  void test_boomeramg(const GlobalMatrixType& matrix, const GlobalFilterType& filter,
    const GlobalVectorType& vec_rhs, GlobalVectorType& vec_sol) const
  {
    vec_sol.format();

    auto precond = Solver::new_boomeramg(matrix, filter);
    auto solver = Solver::new_richardson(matrix, filter, 1.0, precond);

    solver->set_max_iter(5);
    solver->set_tol_rel(1E-5);
    solver->set_plot_name("HYPRE:BoomerAMG");
    solver->set_plot_mode(Solver::PlotMode::all);

    solver->init();
    solver->correct(vec_sol, vec_rhs);
    solver->done();
  }
}; // class HypreTest<...>

HypreTest<Mem::Main, double, Index> hypre_test_main_double_index;

#endif // defined(FEAT_HAVE_HYPRE)
