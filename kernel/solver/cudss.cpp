// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>

#ifdef FEAT_HAVE_CUDSS
#include <kernel/util/cuda_util.hpp>
#include <kernel/solver/cudss.hpp>

#include <cudss.h>

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      struct CUDSSCore
      {
        cudssHandle_t handle;
        cudssConfig_t config;
        cudssData_t data;
        cudssMatrix_t matrix;
        cudssMatrix_t vec_sol;
        cudssMatrix_t vec_rhs;
        Memory::Arbiter row_ptr;
        Memory::Arbiter col_idx;
        //Memory::Arbiter mat_val;
      };
    } // namespace Intern

    CUDSS::CUDSS(const MatrixType& system_matrix) :
      BaseClass(),
      _system_matrix(system_matrix),
      _cudss_core(new Intern::CUDSSCore)
    {
      Intern::CUDSSCore& core = *reinterpret_cast<Intern::CUDSSCore*>(_cudss_core);
      core.handle = reinterpret_cast<cudssHandle_t>(Runtime::get_cudss_handle());
      if(CUDSS_STATUS_SUCCESS != cudssConfigCreate(&core.config))
        throw InternalError(__func__, __FILE__, __LINE__, "cudssConfigCreate failed!");
      if(CUDSS_STATUS_SUCCESS != cudssDataCreate(core.handle, &core.data))
        throw InternalError(__func__, __FILE__, __LINE__, "cudssDataCreate failed!");
    }

    CUDSS::~CUDSS()
    {
      Intern::CUDSSCore* core = reinterpret_cast<Intern::CUDSSCore*>(_cudss_core);
      if(core)
      {
        cudssDataDestroy(core->handle, core->data);
        cudssConfigDestroy(core->config);
        delete core;
      }
    }

    void CUDSS::init_symbolic()
    {
      BaseClass::init_symbolic();

      Intern::CUDSSCore& core = *reinterpret_cast<Intern::CUDSSCore*>(_cudss_core);

      cudssStatus_t ret = CUDSS_STATUS_INTERNAL_ERROR;

      // set the basic configuration
      cudssAlgType_t alg_type = CUDSS_ALG_1; // COLAMD based ordering
      ret = cudssConfigSet(core.config, CUDSS_CONFIG_REORDERING_ALG, &alg_type, sizeof(cudssAlgType_t));
      if(ret != CUDSS_STATUS_SUCCESS)
        throw SolverException("CUDSS: cudssConfigSet() for 'CUDSS_CONFIG_REORDERING_ALG' failed!");

      // get the index arrays of the matrix
      const Index neq = _system_matrix.num_rows();
      const Index nze = _system_matrix.num_nzes();

      // convert indices to 32 bit if necessary
      if constexpr (sizeof(Index) == 4u)
      {
        core.row_ptr = _system_matrix.row_ptr_arbiter().attach();
        core.col_idx = _system_matrix.col_idx_arbiter().attach();
      }
      else
      {
        core.row_ptr = Memory::Arbiter(std::size_t(4u*(neq+1u)));
        core.col_idx = Memory::Arbiter(std::size_t(4u*nze));
        core.row_ptr.template convert<std::uint32_t, Index>(_system_matrix.row_ptr_arbiter(), Memory::Location::main);
        core.col_idx.template convert<std::uint32_t, Index>(_system_matrix.col_idx_arbiter(), Memory::Location::main);
      }

      // get device pointers
      const Memory::TypedView<std::uint32_t> row_ptr_view(core.row_ptr.view(Memory::Location::cuda));
      const Memory::TypedView<std::uint32_t> col_idx_view(core.col_idx.view(Memory::Location::cuda));
      const Memory::TypedView<double> val_view(_system_matrix.val_view_r(Memory::Location::cuda));

      // set matrix data
      ret = cudssMatrixCreateCsr(
        &core.matrix,
        std::int64_t(neq),
        std::int64_t(neq),
        std::int64_t(nze),
        const_cast<std::uint32_t*>(row_ptr_view.get_r()),
        nullptr,
        const_cast<std::uint32_t*>(col_idx_view.get_r()),
        const_cast<double*>(val_view.get_r()),
        CUDA_R_32I,
        CUDA_R_64F,
        CUDSS_MTYPE_GENERAL,
        CUDSS_MVIEW_FULL, // is ignored
        CUDSS_BASE_ZERO);
      if(ret != CUDSS_STATUS_SUCCESS)
        throw SolverException("CUDSS: cudssMatrixCreateCsr() for system matrix failed!");

      // allocate solution vector
      ret = cudssMatrixCreateDn(
        &core.vec_sol,
        std::int64_t(neq),
        1,
        std::int64_t(neq),
        nullptr,
        CUDA_R_64F,
        CUDSS_LAYOUT_COL_MAJOR);

      if(ret != CUDSS_STATUS_SUCCESS)
        throw SolverException("CUDSS: cudssMatrixCreateDn() for solution vector failed!");

      // allocate rhs vector
      ret = cudssMatrixCreateDn(
        &core.vec_rhs,
        std::int64_t(neq),
        1,
        std::int64_t(neq),
        nullptr,
        CUDA_R_64F,
        CUDSS_LAYOUT_COL_MAJOR);

      if(ret != CUDSS_STATUS_SUCCESS)
        throw SolverException("CUDSS: cudssMatrixCreateDn() for rhs vector failed!");

      // perform symbolic factorization
      ret = cudssExecute(
        core.handle,
        CUDSS_PHASE_ANALYSIS,
        core.config,
        core.data,
        core.matrix,
        core.vec_sol,
        core.vec_rhs);
      if(ret != CUDSS_STATUS_SUCCESS)
        throw SolverException("CUDSS: cudssExecute() for phase 'CUDSS_PHASE_ANALYSIS' failed!");

      // wait until CUDA is done
      Util::cuda_synchronize();
    }

    void CUDSS::done_symbolic()
    {
      BaseClass::done_symbolic();
      if(_cudss_core == nullptr)
        return;

      Intern::CUDSSCore& core = *reinterpret_cast<Intern::CUDSSCore*>(_cudss_core);

      cudssMatrixDestroy(core.vec_sol);
      cudssMatrixDestroy(core.vec_rhs);
      cudssMatrixDestroy(core.matrix);

      core.col_idx.release();
      core.row_ptr.release();

      // wait until CUDA is done
      Util::cuda_synchronize();
    }

    void CUDSS::init_numeric()
    {
      BaseClass::init_numeric();

      Intern::CUDSSCore& core = *reinterpret_cast<Intern::CUDSSCore*>(_cudss_core);

      // get a cuda view of the numeric data to ensure that its uploaded to device memory
      const Memory::TypedView<double> val_view(_system_matrix.val_view_r(Memory::Location::cuda));

      // perform numeric factorization
      cudssStatus_t ret = cudssExecute(
        core.handle,
        CUDSS_PHASE_FACTORIZATION,
        core.config,
        core.data,
        core.matrix,
        core.vec_sol,
        core.vec_rhs);
      if(ret != CUDSS_STATUS_SUCCESS)
        throw SolverException("CUDSS: cudssExecute() for phase 'CUDSS_PHASE_FACTORIZATION' failed!");

      val_view.release();

      // wait until CUDA is done
      Util::cuda_synchronize();
    }

    void CUDSS::done_numeric()
    {
      // nothing to do here
      BaseClass::done_numeric();
    }

    Status CUDSS::apply(VectorType& vec_sol, const VectorType& vec_rhs)
    {
      Intern::CUDSSCore& core = *reinterpret_cast<Intern::CUDSSCore*>(_cudss_core);

      cudssStatus_t ret = CUDSS_STATUS_INTERNAL_ERROR;

      Memory::TypedView<double> sol_view(vec_sol.elements_view_w(Memory::Location::cuda));
      Memory::TypedView<double> rhs_view(vec_rhs.elements_view_r(Memory::Location::cuda));

      // set solution vector data array
      ret = cudssMatrixSetValues(core.vec_sol, sol_view.get_w());
      if(ret != CUDSS_STATUS_SUCCESS)
        throw SolverException("CUDSS: cudssMatrixSetValues() failed for vec_sol!");

      // set rhs vector data array
      cudssMatrixSetValues(core.vec_rhs, const_cast<double*>(rhs_view.get_r()));
      if(ret != CUDSS_STATUS_SUCCESS)
        throw SolverException("CUDSS: cudssMatrixSetValues() failed for vec_rhs!");

      // solve
      ret = cudssExecute(
        core.handle,
        CUDSS_PHASE_SOLVE,
        core.config,
        core.data,
        core.matrix,
        core.vec_sol,
        core.vec_rhs);

      // wait until CUDA is done
      Util::cuda_synchronize();

      return (ret == CUDSS_STATUS_SUCCESS ? Status::success : Status::aborted);
    }
  } // namespace Solver
} // namespace FEAT
#else
// insert dummy function to suppress linker warnings
void dummy_function_cudss() {}
#endif // FEAT_HAVE_CUDSS
