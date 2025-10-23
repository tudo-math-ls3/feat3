// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/solver/direct_sparse_solver.hpp>

#ifdef FEAT_HAVE_CUDSS

FEAT_DISABLE_WARNINGS
#include <cudss.h>
FEAT_RESTORE_WARNINGS

namespace FEAT
{
  namespace Solver
  {
    namespace DSS
    {
      static_assert(sizeof(CUDSS_IT) == 4, "DirectSparseSolver: cuDSS: index type size mismatch!");

      /**
       * \brief cuDSS core wrapper class
       *
       * \author Peter Zajac
       */
      class CUDSS_Core
      {
      public:
#ifdef FEAT_HAVE_MPI
        MPI_Comm mpi_comm;
#endif
        cudssHandle_t handle;
        cudssConfig_t config;
        cudssData_t data;
        cudssMatrix_t matrix;
        cudssMatrix_t vec_sol;
        cudssMatrix_t vec_rhs;

        // peak memory usage estimate
        std::int64_t memory_estimates[16];

        // system dimensions
        std::int64_t num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes, num_global_nzes, num_rhs;

        // system data in host memory
        std::vector<CUDSS_IT> row_ptr_host, col_idx_host;
        std::vector<CUDSS_DT> mat_val_host, rhs_val_host, sol_val_host;

        // system data in device memory
        void *row_ptr_dev, *col_idx_dev;
        void *mat_val_dev, *rhs_val_dev, *sol_val_dev;

        explicit CUDSS_Core(const Dist::Comm& comm, Index num_global_dofs_, Index dof_offset_,
          Index num_owned_dofs_, Index num_owned_nzes_, Index num_global_nzes_) :
#ifdef FEAT_HAVE_MPI
          mpi_comm(comm.mpi_comm()),
#endif
          handle(reinterpret_cast<cudssHandle_t>(Runtime::get_cudss_handle())),
          config(nullptr),
          data(nullptr),
          matrix(nullptr),
          vec_sol(nullptr),
          vec_rhs(nullptr),
          num_global_dofs(std::uint32_t(num_global_dofs_)),
          dof_offset(std::uint32_t(dof_offset_)),
          num_owned_dofs(std::uint32_t(num_owned_dofs_)),
          num_owned_nzes(std::uint32_t(num_owned_nzes_)),
          num_global_nzes(std::uint32_t(num_global_nzes_)),
          num_rhs(1),
          row_ptr_host(),
          col_idx_host(),
          mat_val_host(),
          rhs_val_host(),
          sol_val_host(),
          row_ptr_dev(nullptr),
          col_idx_dev(nullptr),
          mat_val_dev(nullptr),
          rhs_val_dev(nullptr),
          sol_val_dev(nullptr)
        {
          // prevent unused parameter warnings
          (void)comm;

          XASSERTM(handle != nullptr, "Failed to retrieve cuDSS handle!");

          if(CUDSS_STATUS_SUCCESS != cudssConfigCreate(&config))
            throw InternalError(__func__, __FILE__, __LINE__, "cudssConfigCreate failed!");

          if(CUDSS_STATUS_SUCCESS != cudssDataCreate(handle, &data))
            throw InternalError(__func__, __FILE__, __LINE__, "cudssDataCreate failed!");

#ifdef FEAT_HAVE_MPI
          if(CUDSS_STATUS_SUCCESS != cudssDataSet(handle, data, CUDSS_DATA_COMM, &mpi_comm, sizeof(MPI_Comm*)))
            throw InternalError(__func__, __FILE__, __LINE__, "cudssDataSet for 'CUDSS_DATA_COMM' failed!");
#endif

          // format memory estimates
          memset(memory_estimates, 0, sizeof(memory_estimates));

          // allocate matrix and vector arrays in host memory
          row_ptr_host.resize(num_owned_dofs+1u, 0u);
          col_idx_host.resize(num_owned_nzes, 0u);
          mat_val_host.resize(num_owned_nzes, 0.0);
          rhs_val_host.resize(num_owned_dofs * num_rhs, 0.0);
          sol_val_host.resize(num_owned_dofs * num_rhs, 0.0);

          // allocate matrix and vector arrays in device memory
          row_ptr_dev = Util::cuda_malloc(sizeof(CUDSS_IT) * (num_owned_dofs+1u));
          col_idx_dev = Util::cuda_malloc(sizeof(CUDSS_IT) * num_owned_nzes);
          mat_val_dev = Util::cuda_malloc(sizeof(CUDSS_DT) * num_owned_nzes);
          rhs_val_dev = Util::cuda_malloc(sizeof(CUDSS_DT) * num_owned_dofs * num_rhs);
          sol_val_dev = Util::cuda_malloc(sizeof(CUDSS_DT) * num_owned_dofs * num_rhs);
        }

        ~CUDSS_Core()
        {
          if(vec_sol)
            cudssMatrixDestroy(vec_sol);
          if(vec_rhs)
            cudssMatrixDestroy(vec_rhs);
          if(matrix)
            cudssMatrixDestroy(matrix);

          if(sol_val_dev)
            Util::cuda_free(sol_val_dev);
          if(rhs_val_dev)
            Util::cuda_free(rhs_val_dev);
          if(mat_val_dev)
            Util::cuda_free(mat_val_dev);
          if(col_idx_dev)
            Util::cuda_free(col_idx_dev);
          if(row_ptr_dev)
            Util::cuda_free(row_ptr_dev);
          if(data)
            cudssDataDestroy(handle, data);
          if(config)
            cudssConfigDestroy(config);

          // wait until CUDA is done
          Util::cuda_synchronize();
        }

        void init_symbolic()
        {
          cudssStatus_t ret = CUDSS_STATUS_INTERNAL_ERROR;

          // first and last DOF owned by this process
          const std::int64_t last_owned_dof  = dof_offset + num_owned_dofs - 1;

          // copy structure host arrays to device
          Util::cuda_copy_host_to_device(row_ptr_dev, row_ptr_host.data(), sizeof(CUDSS_IT) * (num_owned_dofs+1u));
          Util::cuda_copy_host_to_device(col_idx_dev, col_idx_host.data(), sizeof(CUDSS_IT) * num_owned_nzes);

          // set the basic configuration
          //cudssAlgType_t alg_type = CUDSS_ALG_DEFAULT;//CUDSS_ALG_1; // COLAMD based ordering
          //ret = cudssConfigSet(config, CUDSS_CONFIG_REORDERING_ALG, &alg_type, sizeof(cudssAlgType_t));
          //if(ret != CUDSS_STATUS_SUCCESS)
            //throw DirectSparseSolverException("cuDSS", "cudssConfigSet() for 'CUDSS_CONFIG_REORDERING_ALG' failed!");

          // set matrix data
          ret = cudssMatrixCreateCsr(
            &matrix,
            num_global_dofs,
            num_global_dofs,
            num_global_nzes,
            row_ptr_dev,
            nullptr,
            col_idx_dev,
            mat_val_dev,
            CUDA_R_32I,
            CUDA_R_64F,
            CUDSS_MTYPE_GENERAL,
            CUDSS_MVIEW_FULL, // is ignored
            CUDSS_BASE_ZERO);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssMatrixCreateCsr() for system matrix failed!");

          // set row distribution of matrix
          ret = cudssMatrixSetDistributionRow1d(matrix, dof_offset, last_owned_dof);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssMatrixSetDistributionRow1d() for solution vector failed!");

          // allocate solution vector
          ret = cudssMatrixCreateDn(
            &vec_sol,
            num_global_dofs,
            num_rhs,
            num_global_dofs,
            nullptr,
            CUDA_R_64F,
            CUDSS_LAYOUT_COL_MAJOR);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssMatrixCreateDn() for solution vector failed!");

          // set row distribution of solution vector
          ret = cudssMatrixSetDistributionRow1d(vec_sol, dof_offset, last_owned_dof);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssMatrixSetDistributionRow1d() for solution vector failed!");

          // allocate rhs vector
          ret = cudssMatrixCreateDn(
            &vec_rhs,
            num_global_dofs,
            num_rhs,
            num_global_dofs,
            nullptr,
            CUDA_R_64F,
            CUDSS_LAYOUT_COL_MAJOR);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssMatrixCreateDn() for rhs vector failed!");

          // set row distribution of solution vector
          ret = cudssMatrixSetDistributionRow1d(vec_rhs, dof_offset, last_owned_dof);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssMatrixSetDistributionRow1d() for rhs vector failed!");

          // perform symbolic factorization
          ret = cudssExecute(
            handle,
            CUDSS_PHASE_ANALYSIS,
            config,
            data,
            matrix,
            vec_sol,
            vec_rhs);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssExecute() for phase 'CUDSS_PHASE_ANALYSIS' failed!");

          // retrieve memory estimates for current device
          std::size_t bytes_written(0u);
          ret = cudssDataGet(
            handle,
            data,
            CUDSS_DATA_MEMORY_ESTIMATES,
            memory_estimates,
            sizeof(memory_estimates),
            &bytes_written);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssDataGet() for 'CUDSS_DATA_MEMORY_ESTIMATES' failed!");

          // wait until CUDA is done
          Util::cuda_synchronize();
        }

        void init_numeric()
        {
          // copy value host arrays to device
          Util::cuda_copy_host_to_device(mat_val_dev, mat_val_host.data(), sizeof(CUDSS_DT) * num_owned_nzes);

          // perform numeric factorization
          cudssStatus_t ret = cudssExecute(
            handle,
            CUDSS_PHASE_FACTORIZATION,
            config,
            data,
            matrix,
            vec_sol,
            vec_rhs);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssExecute() for phase 'CUDSS_PHASE_FACTORIZATION' failed!");

          // wait until CUDA is done
          Util::cuda_synchronize();
        }

        void solve()
        {
          cudssStatus_t ret = CUDSS_STATUS_INTERNAL_ERROR;

          // copy RHS value array to device
          Util::cuda_copy_host_to_device(rhs_val_dev, rhs_val_host.data(), sizeof(CUDSS_DT) * num_owned_dofs * num_rhs);

          // set solution vector data array
          ret = cudssMatrixSetValues(vec_sol, sol_val_dev);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssMatrixSetValues() failed for vec_sol!");

          // set rhs vector data array
          cudssMatrixSetValues(vec_rhs, rhs_val_dev);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssMatrixSetValues() failed for vec_rhs!");

          // solve
          ret = cudssExecute(
            handle,
            CUDSS_PHASE_SOLVE,
            config,
            data,
            matrix,
            vec_sol,
            vec_rhs);
          if(ret != CUDSS_STATUS_SUCCESS)
            throw DirectSparseSolverException("cuDSS", "cudssExecute() for phase 'CUDSS_PHASE_SOLVE' failed!");

          // copy sol values array to host
          Util::cuda_copy_device_to_host(sol_val_host.data(), sol_val_dev, sizeof(CUDSS_DT) * num_owned_dofs * num_rhs);

          // wait until CUDA is done
          Util::cuda_synchronize();
        }

        std::int64_t get_peak_mem_device() const
        {
          return memory_estimates[1];
        }

        std::int64_t get_peak_mem_host() const
        {
          return memory_estimates[3];
        }
      }; // class CUDSS_Core

      void* create_cudss_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes)
      {
        return new CUDSS_Core(*comm, num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes, num_global_nzes);
      }

      void destroy_cudss_core(void* core)
      {
        XASSERT(core != nullptr);
        delete reinterpret_cast<CUDSS_Core*>(core);
      }

      CUDSS_IT* get_cudss_row_ptr(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<CUDSS_Core*>(core)->row_ptr_host.data();
      }

      CUDSS_IT* get_cudss_col_idx(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<CUDSS_Core*>(core)->col_idx_host.data();
      }

      CUDSS_DT* get_cudss_mat_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<CUDSS_Core*>(core)->mat_val_host.data();
      }

      CUDSS_DT* get_cudss_rhs_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<CUDSS_Core*>(core)->rhs_val_host.data();
      }

      CUDSS_DT* get_cudss_sol_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<CUDSS_Core*>(core)->sol_val_host.data();
      }

      void init_cudss_symbolic(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<CUDSS_Core*>(core)->init_symbolic();
      }

      void init_cudss_numeric(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<CUDSS_Core*>(core)->init_numeric();
      }

      void solve_cudss(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<CUDSS_Core*>(core)->solve();
      }

      std::int64_t get_peak_mem_cudss_host(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<CUDSS_Core*>(core)->get_peak_mem_host();
      }

      std::int64_t get_peak_mem_cudss_device(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<CUDSS_Core*>(core)->get_peak_mem_device();
      }
    } // namespace DSS
  } // namespace Solver
} // namespace FEAT

#else // no FEAT_HAVE_CUDSS

void feat_direct_sparse_solver_cudss_dummy()
{
}

#endif // FEAT_HAVE_CUDSS
