#pragma once

#include <kernel/solver/direct_stokes_solver.hpp>

namespace SCALEXA
{
  using namespace FEAT;

  template<typename SolverLevel_, typename ScalarizedLevel_>
  class StokesScalarizeHelper
  {
  public:
    typedef typename SolverLevel_::DataType DataType;
    typedef typename SolverLevel_::IndexType IndexType;
    static constexpr int dim = SolverLevel_::dim;

    typedef Solver::DirectStokesCore<DataType, IndexType,
      typename SolverLevel_::LocalMatrixBlockA,
      typename SolverLevel_::LocalMatrixBlockB,
      typename SolverLevel_::LocalMatrixBlockD> DirectStokesCoreType;

    const SolverLevel_* solver_level;

    std::unique_ptr<DirectStokesCoreType> direct_stokes_core;

    ScalarizedLevel_ scalarized_level;

  public:
    StokesScalarizeHelper() :
      solver_level(nullptr),
      direct_stokes_core()
    {
    }

    void create(const SolverLevel_& solver_level_)
    {
      solver_level = &solver_level_;
      direct_stokes_core.reset(new DirectStokesCoreType(
        solver_level->matrix_sys.local().block_a(),
        solver_level->matrix_sys.local().block_b(),
        solver_level->matrix_sys.local().block_d()));

      // create a scalarized system level and set communicator
      scalarized_level.gate_sys.set_comm(solver_level->gate_sys.get_comm());
    }

    void init_symbolic()
    {
      if(solver_level == nullptr)
        return;

      // set filters
      direct_stokes_core->set_filters(solver_level->filter_sys.local());

      // perform symbolic initialization
      direct_stokes_core->init_symbolic();

      // get number of DOFs in scalarized system
      const Index num_sclr_dofs = direct_stokes_core->get_solver_matrix().rows();

      // create the gate for the scalarized system
      {
        typedef typename SolverLevel_::SystemMirror SystemMirrorType;
        typedef typename SolverLevel_::VeloMirror VeloMirrorType;
        typedef typename SolverLevel_::PresMirror PresMirrorType;
        typedef typename ScalarizedLevel_::SystemMirror SclrMirrorType;

        // loop over all mirrors in the gate
        const std::vector<SystemMirrorType>& sys_mirrors = solver_level->gate_sys.get_mirrors();
        const std::vector<int>& neighbor_ranks = solver_level->gate_sys.get_ranks();
        XASSERTM(neighbor_ranks.size() == sys_mirrors.size(), "velocity mirrors rank count mismatch!");

        // loop over all neighbor ranks
        for(std::size_t ineigh(0); ineigh < neighbor_ranks.size(); ++ineigh)
        {
          const int rank = neighbor_ranks.at(ineigh);

          // get velocity mirror
          const VeloMirrorType& velo_mir = sys_mirrors.at(ineigh).template at<0>();
          const Index velo_size = velo_mir.size();
          const Index velo_num_idx = velo_mir.num_indices();
          const IndexType* velo_mir_idx = velo_mir.indices();

          // get pressure mirror
          const PresMirrorType& pres_mir = sys_mirrors.at(ineigh).template at<1>();
          const Index pres_size = pres_mir.size();
          const Index pres_num_idx = pres_mir.num_indices();
          const IndexType* pres_mir_idx = pres_mir.indices();

          // allocate scalarized mirror
          XASSERT(num_sclr_dofs == Index(dim) * velo_size + pres_size);
          Index num_sclr_idx = Index(dim) * velo_num_idx + pres_num_idx;
          SclrMirrorType sclr_mirror(num_sclr_dofs, num_sclr_idx);
          IndexType* sclr_mir_idx = sclr_mirror.indices();

          // add velocity mirror indices
          Index k(0);
          for(Index i(0); i < velo_num_idx; ++i)
            for(int j(0); j < dim; ++j, ++k)
              sclr_mir_idx[k] = Index(dim)*velo_mir_idx[i] + Index(j);

          // add pressure mirror indices
          Index off = Index(dim) * velo_size;
          for(Index i(0); i < pres_num_idx; ++i, ++k)
            sclr_mir_idx[k] = off + pres_mir_idx[i];

          // push the mirror
          scalarized_level.gate_sys.push(rank, std::move(sclr_mirror));
        }

        // compile gate
        scalarized_level.gate_sys.compile(direct_stokes_core->create_solver_vector());
      }

      // create a unit filter for the scalarized system
      {
        // allocate unit filter
        typedef typename ScalarizedLevel_::LocalSystemFilter LocalScalarizedFilter;
        LocalScalarizedFilter& sclr_unit_filter = scalarized_level.filter_sys.local();
        sclr_unit_filter = LocalScalarizedFilter(num_sclr_dofs);


        const auto& vf = solver_level->filter_velo.local();

        // get indices and values
        const Index n = vf.used_elements();
        if(n > 0u)
        {
          const IndexType* idx = vf.get_indices();
          const Tiny::Vector<DataType, dim>* val = vf.get_values();
          for(Index i(0); i < n; ++i)
            for(int j(0); j < dim; ++j)
              sclr_unit_filter.add(Index(dim)*idx[i] + Index(j), val[i][j]);
        }
      }

      // make a shallow copy of the scalarized matrix
      scalarized_level.matrix_sys.local().clone(direct_stokes_core->get_solver_matrix(), LAFEM::CloneMode::Shallow);
    }

    void init_numeric()
    {
      direct_stokes_core->init_numeric();
      scalarized_level.matrix_sys.local().copy(direct_stokes_core->get_solver_matrix());
    }
  };

  template<typename SolverLevel_, typename ScalarizedLevel_>
  class StokesScalarizeWrapper :
    public Solver::SolverBase<typename SolverLevel_::GlobalSystemVector>
  {
  public:
    typedef Solver::SolverBase<typename SolverLevel_::GlobalSystemVector> BaseClass;
    typedef typename SolverLevel_::GlobalSystemVector GlobalSolverVector;
    typedef typename ScalarizedLevel_::GlobalSystemVector GlobalScalarizedVector;

    StokesScalarizeHelper<SolverLevel_, ScalarizedLevel_>& helper;
    std::shared_ptr<Solver::SolverBase<GlobalScalarizedVector>> scalar_solver;

    GlobalScalarizedVector vec_def_sc, vec_cor_sc;

    StopWatch watch_helper_init_symbolic, watch_helper_init_numeric, watch_solver_init_symbolic, watch_solver_init_numeric, watch_solver_apply;

  public:
    explicit StokesScalarizeWrapper(StokesScalarizeHelper<SolverLevel_, ScalarizedLevel_>& helper_,
      std::shared_ptr<Solver::SolverBase<GlobalScalarizedVector>> scalar_solver_) :
      helper(helper_),
      scalar_solver(scalar_solver_)
    {
    }

    virtual String name() const override
    {
      return "StokesScalarizeWrapper[" + scalar_solver->name() + "]";
    }

    virtual void init_symbolic() override
    {
      BaseClass::init_symbolic();
      watch_helper_init_symbolic.start();
      helper.init_symbolic();
      watch_helper_init_symbolic.stop();

      watch_solver_init_symbolic.start();
      scalar_solver->init_symbolic();
      watch_solver_init_symbolic.stop();

      vec_def_sc = helper.scalarized_level.matrix_sys.create_vector_l();
      vec_cor_sc = helper.scalarized_level.matrix_sys.create_vector_l();
    }

    virtual void init_numeric() override
    {
      BaseClass::init_numeric();

      watch_helper_init_numeric.start();
      helper.init_numeric();
      watch_helper_init_numeric.stop();

      watch_solver_init_numeric.start();
      scalar_solver->init_numeric();
      watch_solver_init_numeric.stop();
    }

    virtual void done_numeric() override
    {
      BaseClass::done_numeric();
      //helper.done_numeric();
      scalar_solver->done_numeric();
    }

    virtual void done_symbolic() override
    {
      BaseClass::done_symbolic();
      //helper.done_symbolic();
      scalar_solver->done_symbolic();
    }

    virtual Solver::Status apply(GlobalSolverVector& vec_cor, const GlobalSolverVector& vec_def) override
    {
      helper.direct_stokes_core->upload_vector(vec_def_sc.local(), vec_def.local().template at<0>(), vec_def.local().template at<1>());

      helper.scalarized_level.filter_sys.filter_def(vec_def_sc);

      watch_solver_apply.start();
      Solver::Status status = scalar_solver->apply(vec_cor_sc, vec_def_sc);
      watch_solver_apply.stop();

      helper.direct_stokes_core->download_vector(vec_cor_sc.local(), vec_cor.local().template at<0>(), vec_cor.local().template at<1>());

      return status;
    }
  };
}
