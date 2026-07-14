#include <kernel/base_header.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/synch_vec.hpp>
#include <kernel/lafem/arch/batched_vector_apply.hpp>
#include <kernel/lafem/arch/batched_kernel_helper.cuh>
#include <kernel/solver/base.hpp>
#include <kernel/solver/amavanka_gmres.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/solver/backend_aux/csr_apply_helper.hpp>
#include "cusparse_v2.h"
#include <cublas_v2.h>
#include <cuda/atomic>

#include <kernel/util/cuda_helper.cuh>


namespace FEAT::Solver::VankaGMRES
{

  namespace cg = cooperative_groups;

  struct CuMatDesc
  {
    cusparseSpMatDescr_t descr=0;
    void* buffer=nullptr;
  };


  namespace Intern
  {
    // The matrix value/row/column device pointers are resolved by the caller from the matrix
    // arbiters (cuda read views, kept alive for the descriptor's lifetime) and passed in here.
    template<typename DT_, typename IT_>
    CUDA_HOST void alloc_cusparse_interface(cusparseOperation_t op, cudaDataType_t dt, cudaDataType_t ct, cusparseIndexType_t it,
        cusparseDnVecDescr_t dx, cusparseDnVecDescr_t dr,
        const void* row_ptr, const void* col_ptr, const void* val_ptr, Index num_rows, Index num_cols, Index used_elements,
        CuMatDesc& mat_desc, const DT_& alpha, const DT_& beta)
    {
      const auto algo = CUSPARSE_SPMV_CSR_ALG1;

      cusparseStatus_t cu_status = cusparseCreateCsr(&mat_desc.descr, num_rows, num_cols, used_elements, (void*)row_ptr, (void*)col_ptr, (void*)val_ptr,
        it, it, CUSPARSE_INDEX_BASE_ZERO, dt);
      if(cu_status != CUSPARSE_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cusparseCreateCsr failed with status code: " + stringify(cusparseGetErrorString(cu_status)));

      size_t buffer_size(0);
      cu_status = cusparseSpMV_bufferSize(Util::Intern::cusparse_handle, op, &alpha, mat_desc.descr, dx, &beta, dr, ct, algo, &buffer_size);
      if(cu_status != CUSPARSE_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMV_bufferSize failed with status code: " + stringify(cusparseGetErrorString(cu_status)));

      mat_desc.buffer = Memory::alloc_cuda(buffer_size);
      cu_status = cusparseSpMV_preprocess(Util::Intern::cusparse_handle, op, &alpha, mat_desc.descr, dx, &beta, dr, ct, algo, mat_desc.buffer);
      if(cu_status != CUSPARSE_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMV_preprocess failed with status code: " + stringify(cusparseGetErrorString(cu_status)));
    }

    static void cumat_free(CuMatDesc& mat_desc)
    {
      Memory::free_cuda(mat_desc.buffer);
      mat_desc.buffer = nullptr;
      cusparseDestroySpMat(mat_desc.descr);
      mat_desc.descr = 0;
    }

    template<typename DT_>
    static void apply_cu_spmv(cusparseOperation_t op, cudaDataType_t ct, cusparseDnVecDescr_t dx, cusparseDnVecDescr_t dr, const CuMatDesc& mat_desc, DT_ alpha, DT_ beta)
    {
      cusparseStatus_t status = cusparseSpMV(Util::Intern::cusparse_handle, op, &alpha, mat_desc.descr, dx, &beta, dr, ct, CUSPARSE_SPMV_CSR_ALG1, mat_desc.buffer);
      if (status != CUSPARSE_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMV failed with status code: " + stringify(cusparseGetErrorString(status)));
    }

    template<typename DT_>
    static void cublas_dot(DT_& result, const DT_* x, const DT_* y, Index size, cudaDataType_t dt)
    {
      cublasStatus_t status;
      result = DT_(0);

      status = cublasDotEx(Util::Intern::cublas_handle, int(size), x, dt, 1, y, dt, 1, &result, dt, dt);
      if (status != CUBLAS_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
    }

    template<typename DT_>
    static void cublas_scale(DT_* x, DT_ a, Index size, cudaDataType_t dt)
    {
      cublasStatus_t status = cublasScalEx(Util::Intern::cublas_handle, int(size), &a, dt, x, dt, 1, dt);
      if (status != CUBLAS_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
    }

    template<typename DT_>
    static void cublas_axpy(DT_* r, DT_ a, const DT_* x, Index size, cudaDataType_t dt)
    {
      cublasStatus_t status = cublasAxpyEx(Util::Intern::cublas_handle, int(size), &a, dt, x, dt, 1, r, dt, 1, dt);
      if (status != CUBLAS_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
    }

    template<typename DT_>
    static void cuda_sync_dot(DT_* result, DT_* device_helper, const DT_* x, const DT_* y, const DT_* freqs, Index size, const Dist::Comm& comm, cudaStream_t stream = 0)
    {
      cudaCheck(cudaMemsetAsync(device_helper, 0, sizeof(DT_)), "cudaMemSetAsync");
      constexpr Index blocksize = 512u;
      dim3 grid;
      dim3 block;
      block.x = (unsigned)blocksize;
      grid.x = (unsigned)ceil((size)/(double)(block.x*2*10));
      LAFEM::Arch::cuda_triple_dot<DT_, blocksize><<<grid, block, blocksize*sizeof(DT_), stream>>>(device_helper, x, size, y, freqs);
      cudaCheck(cudaMemcpyAsync(result, device_helper, sizeof(DT_), cudaMemcpyDeviceToHost, stream), "cudaMemCpy");
      cudaStreamSynchronize(stream);
      comm.allreduce(result, result, 1u, Dist::op_sum);
    }

    template <typename DT_, typename IT_>
    __global__ void cuda_unit_filter_def(DT_ * v, const IT_ * sv_indices, const Index ue)
    {
      Index idx = threadIdx.x + blockDim.x * blockIdx.x;
      // grid strided for loop
      for(int i = idx; i < ue; i += blockDim.x * gridDim.x)
      {
        v[sv_indices[i]] = DT_(0);
      }
    }

    template <typename DT_>
    __global__ void cuda_ComponentProduct(DT_ * r, const DT_ * x, const DT_ * y, const Index count)
    {
      Index idx = threadIdx.x + blockDim.x * blockIdx.x;
      for(Index i = idx; i < count; i += Index(blockDim.x*gridDim.x))
        r[i] = x[i] * y[i];
    }

    /**
     * \brief Performs m[k] = alpha * (A*rhs)[k_j] + beta * m[k]
     * where k_j is the actual sparse mirror index
     */
    template<typename DT_, typename IT_>
    __global__ void calculate_mirror_entries(DT_* __restrict__ mirror_buffer, const IT_* __restrict__ mirror_indx, const Index mirror_size, const IT_* __restrict__ row_ptr, const IT_* __restrict__ colidx, const DT_* values, const DT_* rhs, DT_ alpha, DT_ beta)
    {
      Index idx = Index(threadIdx.x + blockDim.x * blockIdx.x);
      // grid strided for loop
      for(Index i = idx; i < mirror_size; i += Index(blockDim.x * gridDim.x))
      {
        DT_ val = DT_(0);
        //TODO: retry with fixed kernel, index was wrong before...
        for(Index k = row_ptr[mirror_indx[i]]; k < row_ptr[mirror_indx[i]+1]; ++k)
        {
          val += values[k] * rhs[colidx[k]];
        }
        mirror_buffer[i] = alpha * val + beta * mirror_buffer[i];
      }
    }

    template<typename DT_, typename IT_>
    __global__ void gather_mirror_entries(DT_* __restrict__ mirror_buffer, const IT_* __restrict__ mirror_indx, const Index mirror_size, const DT_* __restrict__ vec)
    {
      Index idx = Index(threadIdx.x + blockDim.x * blockIdx.x);
      // grid strided for loop
      for(Index i = idx; i < mirror_size; i += Index(blockDim.x * gridDim.x))
      {
        mirror_buffer[i] = vec[mirror_indx[i]];
      }
    }

    template<typename DT_, typename IT_>
    __global__ void scatter_mirror_entries(DT_* __restrict__ vec, DT_* __restrict__ mirror_buffer, const IT_* __restrict__ mirror_indx, const Index mirror_size)
    {
      Index idx = Index(threadIdx.x + blockDim.x * blockIdx.x);
      // grid strided for loop
      for(Index i = idx; i < mirror_size; i += Index(blockDim.x * gridDim.x))
      {
        // to allow for asynchronizity by using atomics
        cuda::atomic_ref<DT_, cuda::thread_scope_device> val(vec[mirror_indx[i]]);
        val.fetch_add(mirror_buffer[i]);
      }
    }

    template<typename DT_, typename IT_>
    static void apply_def_filter(DT_* vec, const IT_* filter_idx, Index filter_ue, Index blocksize, Index gridsize)
    {
      // and filter
      dim3 grid, block;
      block.x = (unsigned)blocksize;
      grid.x = (unsigned)gridsize;
      cuda_unit_filter_def<<<grid, block>>>(vec, filter_idx, filter_ue);
      cudaCheckLastErr()
    }

    template<typename DT_>
    static void scale_vec(DT_* r, const DT_* x, const DT_* y, Index size, cudaStream_t stream = 0)
    {
      Index blocksize = Util::cuda_blocksize_axpy;
      dim3 grid;
      dim3 block;
      block.x = (unsigned)blocksize;
      // in theory, this should be num_rows
      grid.x = (unsigned)ceil((size)/(double)(block.x));

      cuda_ComponentProduct<<<grid, block, 0, stream>>>(r, x, y, size);
      cudaCheckLastErr()
    }

    template<typename DT_>
    static void calculate_batched_dots(DT_* dots, const DT_* target, DT_ const * const* device_elements, Index dot_size, Index num_dots)
    {
      constexpr Index blocksize = 256u;
      // constexpr Index blocksize = 32u;
      constexpr int num_batch = 4;
      dim3 grid;
      dim3 block;
      block.x = (unsigned)blocksize;
      block.y = num_batch;
      grid.x = (unsigned)ceil((dot_size)/(double)(block.x*2*10));
      grid.y = (unsigned)ceil((num_dots)/(double)(block.y));
      LAFEM::Arch::cuda_batched_dot<DT_, blocksize><<<grid, block, blocksize*num_batch*sizeof(DT_)>>>(dots, target, dot_size, device_elements, num_dots);
      cudaCheckLastErr()
    }

    template<typename DT_>
    static void calculate_synced_batched_dots(DT_* dot_helper, DT_* dot_helper_device, const DT_* target, DT_ const * const* device_elements, DT_ const * const freqs, Index dot_size, Index num_dots, const Dist::Comm& comm, cudaStream_t stream = 0)
    {
      constexpr Index blocksize = 256u;
      // constexpr Index blocksize = 32u;
      constexpr int num_batch = 4;
      dim3 grid;
      dim3 block;
      block.x = (unsigned)blocksize;
      block.y = num_batch;
      grid.x = (unsigned)ceil((dot_size)/(double)(block.x*2*10));
      grid.y = (unsigned)ceil((num_dots)/(double)(block.y));
      cudaCheck(cudaMemsetAsync(dot_helper_device, 0, size_t(num_dots)*sizeof(DT_), stream), "cudaMemsetAsync");
      LAFEM::Arch::cuda_batched_triple_dot<DT_, blocksize><<<grid, block, blocksize*num_batch*sizeof(DT_), stream>>>(dot_helper_device, target, dot_size, device_elements, freqs, num_dots);
      cudaCheckLastErr()
      cudaCheck(cudaMemcpyAsync(dot_helper, dot_helper_device, size_t(num_dots)*sizeof(DT_), cudaMemcpyDeviceToHost, stream), "cudaMemcpyAsync");
      cudaStreamSynchronize(stream);
      cudaCheckLastErr()
      comm.allreduce(dot_helper, dot_helper, num_dots, Dist::op_sum);
    }

    template<typename DT_>
    static void add_vecs(DT_* vec_v, DT_ const* const * device_elements, const DT_* dots, Index vec_size, Index num_axpy, DT_ alpha)
    {
      if(vec_size == Index(0))
        return;
      constexpr Index blocksize = 128u;
      dim3 grid;
      dim3 block;
      block.x = (unsigned)blocksize;
      // one block always does a warp size chunk of target (due to load coalescing)
      grid.x = static_cast<unsigned int>((vec_size+32u-1u)/32u);
      LAFEM::Arch::cuda_cumulative_axpy<DT_, blocksize><<<grid, block, blocksize*sizeof(DT_)>>>(vec_v, vec_size, device_elements, dots, num_axpy, alpha);
      cudaCheckLastErr()
    }


    template<typename DT_, typename IT_>
    __global__ void cuda_compare(const DT_* mirror_buffer, const IT_* mirror_indx, Index num_mirror, const DT_* vec_sol)
    {
      constexpr DT_ tol = DT_(1E-8);
      Index idx = threadIdx.x + blockDim.x * blockIdx.x;
      if(idx >= num_mirror)
        return;
      std::printf("Working, values %f, %f, %lu, %lu\n", mirror_buffer[idx], vec_sol[mirror_indx[idx]], idx, Index(mirror_indx[idx]));
      if(abs(mirror_buffer[idx] - vec_sol[mirror_indx[idx]]) > tol)
        std::printf("Error, values %f, %f, %lu, %lu larger\n", mirror_buffer[idx], vec_sol[mirror_indx[idx]], idx, Index(mirror_indx[idx]));
    }

    /// global r = a*A*x   product   be careful about setting b correctly, since this factor adds r
    template<typename DT_, typename IT_, typename Gate_, std::size_t st_>
    CUDA_HOST void apply_synced_mat(std::vector<DT_*>& buffer_host, std::vector<DT_*>& buffer_device, DT_* r, cusparseDnVecDescr_t dr, const DT_* x, cusparseDnVecDescr_t dx,
                      const Gate_& gate, cusparseOperation_t op, cudaDataType_t ct, const CuMatDesc& mat_desc, const CSRHelper::CRSMat<DT_, IT_>& mat, DT_ alpha, DT_ beta, cudaStream_t main_stream, const std::array<cudaStream_t, st_>& streams)
    {
      // fill mirror buffers, first gather values
      const auto& mirrors = gate.get_mirrors();
      const auto& ranks = gate.get_ranks();
      const auto& comm = *gate.get_comm();
      const std::size_t n = ranks.size();
      // acquire device views of the mirror index arrays; kept alive for the whole exchange
      std::vector<Memory::TypedView<IT_>> mirror_idx_views;
      mirror_idx_views.reserve(n);
      for(std::size_t k = 0; k < n; ++k)
        mirror_idx_views.emplace_back(mirrors[k].indices_view_r(Memory::Location::cuda));
      // Index s_id = 0;
      // for(Index k = 0; k < n; ++k)
      // {
      //   dim3 block;
      //   dim3 grid;
      //   block.x = (unsigned)64;
      //   grid.x = (unsigned)ceil((mirrors.at(k).num_indices())/(double)((block.x)));
      //   gather_mirror_entries<<<grid, block, 0, streams[s_id]>>>(buffer_device[k], mirrors[k].indices(), mirrors[k].num_indices(), r);
      //   cudaDeviceSynchronize();
      //   cudaCheckLastErr()
      //   calculate_mirror_entries<<<grid, block, 0, streams[s_id]>>>(buffer_device[k], mirrors[k].indices(), mirrors[k].num_indices(), mat.row, mat.col, mat.val, x, alpha, beta);
      //   cudaDeviceSynchronize();
      //   cudaCheckLastErr()
      //   cudaCheck(cudaMemcpyAsync(buffer_host[2*k], buffer_device[k], mirrors[k].num_indices()*sizeof(DT_), cudaMemcpyDeviceToHost, streams[s_id]), "cudaMemcpyAsync");
      //   s_id = (s_id + 1u)%streams.size();
      // }
        // cudaDeviceSynchronize();
        // cudaCheckLastErr()
      // start cusparse
      apply_cu_spmv(op, ct, dx, dr, mat_desc, alpha, beta);

      cudaDeviceSynchronize();
      Index s_id = 0;
      for(Index k = 0; k < n; ++k)
      {
        dim3 block;
        dim3 grid;
        block.x = (unsigned)64;
        grid.x = (unsigned)ceil((mirrors.at(k).num_indices())/(double)((block.x)));
        gather_mirror_entries<<<grid, block, 0, streams[s_id]>>>(buffer_device[k], mirror_idx_views[k].get_r(), mirrors[k].num_indices(), r);
        // cudaDeviceSynchronize();
        // cudaCheckLastErr()
        // calculate_mirror_entries<<<grid, block, 0, streams[s_id]>>>(buffer_device[k], mirrors[k].indices(), mirrors[k].num_indices(), mat.row, mat.col, mat.val, x, alpha, beta);
        // cudaDeviceSynchronize();
        // cudaCheckLastErr()
        cudaCheck(cudaMemcpyAsync(buffer_host[2*k], buffer_device[k], mirrors[k].num_indices()*sizeof(DT_), cudaMemcpyDeviceToHost, streams[s_id]), "cudaMemcpyAsync");
        s_id = (s_id + 1u)%streams.size();
      }
        // cudaDeviceSynchronize();
        // cudaCheckLastErr()
      // compare entries
      // for(Index k = 0; k < n; ++k)
      // {
      //   std::printf("Using mirror of size %lu \n", mirrors.at(k).num_indices());
      //   const auto& mirr_indx = mirrors.at(k).indices();
      //   for(Index l = 0; l < mirrors.at(k).num_indices(); ++l)
      //   {
      //     if(Math::abs(r[mirr_indx[l]] - buffer_host[k][l]) > DT_(1E-3))
      //     {
      //       std::printf("Err between r %f, m %f, at idx %lu, mirr %lu \n", r[mirr_indx[l]], buffer_host[k][l], Index(mirr_indx[l]), Index(l));
      //     }
      //   }

      //   // dim3 block;
      //   // dim3 grid;
      //   // block.x = (unsigned)64;
      //   // grid.x = (unsigned)ceil((mirrors.at(k).num_indices())/(double)((block.x)));
      //   // cuda_compare<<<grid, block, 0>>>(buffer_device[k], mirrors[k].indices(), mirrors[k].num_indices(), r);
      //   // cudaDeviceSynchronize();
      //   // cudaCheckLastErr()
      // }

      Dist::RequestVector _send_reqs, _recv_reqs;
      // post receives
      _recv_reqs.reserve(n);
      for(std::size_t i(0); i < n; ++i)
      {
        // post receive
        _recv_reqs.push_back(comm.irecv(buffer_host.at(2u*i+1), mirrors.at(i).num_indices(), ranks.at(i)));
      }

      // post sends
      _send_reqs.reserve(n);

      // wait for stream to finish and push send
      s_id = 0;
      for(std::size_t i(0); i < n; ++i)
      {
        // wait for stream to finish
        cudaStreamSynchronize(streams[s_id]);
        // post send
        _send_reqs.push_back(comm.isend(buffer_host.at(2u*i), mirrors.at(i).num_indices(), ranks.at(i)));
        s_id = (s_id + 1u)%streams.size();
      }

      // wait for sync to finish and copy back to device
      // process all pending receives
      s_id = 0;
      for(std::size_t idx(0u); _recv_reqs.wait_any(idx); )
      {
        // copy the buffer
        cudaCheck(cudaMemcpyAsync(buffer_device[idx], buffer_host[2*idx+1], mirrors[idx].num_indices()*sizeof(DT_), cudaMemcpyHostToDevice, streams[s_id]), "cudaMemcpyAsync");
        s_id = (s_id + 1u)%streams.size();
      }

      //wait for all streams to finish
      cudaDeviceSynchronize();
      cudaCheckLastErr()

      s_id = 0;
      // and now, add
      for(Index k = 0; k < n; ++k)
      {
        dim3 block;
        dim3 grid;
        block.x = (unsigned)64;
        grid.x = (unsigned)ceil((mirrors.at(k).num_indices())/(double)((block.x)));
        scatter_mirror_entries<<<grid, block, 0, streams[s_id]>>>(r, buffer_device[k], mirror_idx_views[k].get_r(), mirrors[k].num_indices());
        //TODO: stream iter missing...
      }
      cudaCheckLastErr()

      // wait for all sends to finish
      _send_reqs.wait_all();
    }
  }


  template<typename DT_>
  static FEAT::Solver::Status _analyze_defect(CSRHelper::GMRESStatus<DT_>& stat, DT_ def_prev, DT_ def_init, Index max_iter, Index min_iter, DT_ tol_abs, DT_ tol_rel)
  {
    using FEAT::Solver::Status;
    // ensure that the defect is neither NaN nor infinity
    if(!Math::isfinite(stat.def_norm))
      return Status::aborted;

    // // is diverged?
    // if(this->is_diverged(def_cur))
    //   return Status::diverged;

    // maximum number of iterations performed?
    if(stat.num_iter >= max_iter)
      return Status::max_iter;

    // minimum number of iterations performed?
    if(stat.num_iter < min_iter)
      return Status::progress;

    // is converged?
    if((stat.def_norm <= tol_abs) && ((stat.def_norm <= (tol_rel * def_init))))
      return Status::success;

    // // check for stagnation?
    // if(check_stag && (this->_min_stag_iter > Index(0)))
    // {
    //   // did this iteration stagnate?
    //   if(def_cur >= this->_stag_rate * def_prev)
    //   {
    //     // increment stagnation count
    //     if(++this->_num_stag_iter >= this->_min_stag_iter)
    //       return Status::stagnated;
    //   }
    //   else
    //   {
    //     // this iteration did not stagnate
    //     this->_num_stag_iter = Index(0);
    //   }
    // }

    // continue iterating
    return Status::progress;
  }

  template<typename DT_, typename IT_>
  FEAT::Solver::Status cuda_apply_amavanka_gmres_seq(Memory::Arbiter& vec_sol_arb, const Memory::Arbiter& vec_rhs_arb, const CSRHelper::GMRESApplyInfo<DT_, IT_>& apply_info, CSRHelper::GMRESStatus<DT_>& status)
  {
    using namespace FEAT;
    const auto id_setup = Util::cuda_start_range("GMRES Setup");
    // extract info
    const Index rows = apply_info.size;
    const Index columns = apply_info.size;

    // --- acquire device (cuda) views for all arbiter-backed data; these views are kept alive for the
    //     entire solve so the resolved device pointers stay valid, and are released on function return ---
    Memory::TypedView<DT_> vec_sol_view(vec_sol_arb.view(Memory::Location::cuda, Memory::Access::read_write));
    DT_* vec_sol = vec_sol_view.get_w();
    Memory::TypedView<DT_> vec_rhs_view(vec_rhs_arb.view(Memory::Location::cuda, Memory::Access::read));
    const DT_* vec_rhs = vec_rhs_view.get_r();

    // matrix value/row/column device pointers (read views) and the nnz counts (read from a main view)
    Memory::View sys_val_view(apply_info.mat_sys.val->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View sys_row_view(apply_info.mat_sys.row->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View sys_col_view(apply_info.mat_sys.col->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View vanka_val_view(apply_info.mat_vanka.val->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View vanka_row_view(apply_info.mat_vanka.row->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View vanka_col_view(apply_info.mat_vanka.col->view(Memory::Location::cuda, Memory::Access::read));

    // krylov basis vectors: collect device pointers (read-write) into a host array, keep views alive
    const std::vector<Memory::Arbiter*>& vec_v_arbs = *apply_info.vec_v;
    std::vector<Memory::TypedView<DT_>> vec_v_views;
    vec_v_views.reserve(vec_v_arbs.size());
    std::vector<DT_*> vec_v_host(vec_v_arbs.size());
    for(std::size_t k = 0; k < vec_v_arbs.size(); ++k)
    {
      vec_v_views.emplace_back(vec_v_arbs[k]->view(Memory::Location::cuda, Memory::Access::read_write));
      vec_v_host[k] = vec_v_views[k].get_w();
    }
    // const Index used_elements_def = Index(apply_info.mat_sys.row[rows]);
    // const Index used_elements_vanka = Index(apply_info.mat_vanka.row[rows]);
    // create necessary data arrays on host and device
    // setup workspace for cusparse matrices
    cusparseOperation_t trans = CUSPARSE_OPERATION_NON_TRANSPOSE;
    cudaDataType dt;
    cudaDataType ct; //compute type
    if (typeid(DT_) == typeid(double))
    {
        dt = CUDA_R_64F;
        ct = CUDA_R_64F;
    }
    else if (typeid(DT_) == typeid(float))
    {
        dt = CUDA_R_32F;
        ct = CUDA_R_32F;
    }
    cusparseIndexType_t it;
    if(sizeof(IT_) == 4u)
      it = CUSPARSE_INDEX_32I;
    else if(sizeof(IT_) == 8u)
      it = CUSPARSE_INDEX_64I;
    cusparseStatus_t cu_status;
    const DT_ defect_alpha = DT_(-1);
    const DT_ defect_beta = DT_(1);
    const DT_ vanka_alpha = DT_(1);
    const DT_ vanka_beta = DT_(0);
    const DT_ abs_tol = apply_info.abs_tol;
    const DT_ rel_tol = apply_info.rel_tol;
    const DT_ inner_res_scale = apply_info.res_scale;
    const Index krylov_dim = apply_info.krylov_dim;
    const Index max_iter = apply_info.max_iter;
    const Index min_iter = apply_info.min_iter;

    const Index plot_interval = apply_info.plot_interval;
    const Solver::PlotMode plot_mode = apply_info.plot_mode;

    cusparseDnVecDescr_t dx;
    {
      cu_status = cusparseCreateDnVec(&dx, columns, (void*)vec_sol, dt);
      if (cu_status != CUSPARSE_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
    }

    cusparseDnVecDescr_t dr;
    {
      cu_status = cusparseCreateDnVec(&dr, columns, (void*)vec_rhs, dt);
      if (cu_status != CUSPARSE_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
    }

    CuMatDesc cumat_sys;
    Intern::alloc_cusparse_interface<DT_, IT_>(trans, dt, ct, it, dx, dr, sys_row_view.raw_r(), sys_col_view.raw_r(), sys_val_view.raw_r(),
      apply_info.mat_sys.num_rows, apply_info.mat_sys.num_cols, apply_info.mat_sys.num_nnz, cumat_sys, defect_alpha, defect_beta);

    CuMatDesc cumat_vanka;
    Intern::alloc_cusparse_interface<DT_, IT_>(trans, dt, ct, it, dx, dr, vanka_row_view.raw_r(), vanka_col_view.raw_r(), vanka_val_view.raw_r(),
      apply_info.mat_vanka.num_rows, apply_info.mat_vanka.num_cols, apply_info.mat_vanka.num_nnz, cumat_vanka, vanka_alpha, vanka_beta);

    // allocations for helper arrays
    const Index gauss_size = krylov_dim * (krylov_dim+1u) / 2u;
    DT_ *_c, *_s, *_q, *_h;
    cudaCheck(cudaMallocHost(&_c, size_t(krylov_dim)*sizeof(DT_)), "cudaMallocHost");
    cudaCheck(cudaMallocHost(&_s, size_t(krylov_dim)*sizeof(DT_)), "cudaMallocHost");
    cudaCheck(cudaMallocHost(&_q, size_t(krylov_dim+1)*sizeof(DT_)), "cudaMallocHost");
    cudaCheck(cudaMallocHost(&_h, size_t(gauss_size)*sizeof(DT_)), "cudaMallocHost");
    // DT_ *_c_d, *_s_d, *_q_d, *_h_d;
    DT_ *_q_d, *_h_d;
    // cudaCheck(cudaMalloc(&_c_d, size_t(krylov_dim)*sizeof(DT_)), "cudaMalloc");
    // cudaCheck(cudaMalloc(&_s_d, size_t(krylov_dim)*sizeof(DT_)), "cudaMalloc");
    cudaCheck(cudaMalloc(&_q_d, size_t(krylov_dim+1)*sizeof(DT_)), "cudaMalloc");
    cudaCheck(cudaMalloc(&_h_d, size_t(gauss_size)*sizeof(DT_)), "cudaMalloc");
    // scratch for the batched projection dots plus the appended self-dot (krylov_dim projections + 1)
    DT_ *batch_dots_h, *batch_dots_d;
    cudaCheck(cudaMallocHost(&batch_dots_h, size_t(krylov_dim+1)*sizeof(DT_)), "cudaMallocHost");
    cudaCheck(cudaMalloc(&batch_dots_d, size_t(krylov_dim+1)*sizeof(DT_)), "cudaMalloc");
    DT_** vec_v = vec_v_host.data();
    Memory::TypedView<DT_> vec_w_view(apply_info.vec_w->view(Memory::Location::cuda, Memory::Access::read_write));
    DT_* vec_w = vec_w_view.get_w();
    DT_** vec_v_d;
    cudaCheck(cudaMalloc(&vec_v_d, size_t(krylov_dim+1)*sizeof(DT_*)), "cudaMalloc");

    if(apply_info.prefetch_data)
    {
      cudaMemLocation location;
      {location.type = cudaMemLocationTypeDevice; location.id = FEAT::Util::cuda_get_device_id();}
      // hint to prefetch vector data, if possible
      for(Index k = 0; k < krylov_dim+1; ++k)
      {
        cudaMemPrefetchAsync(vec_v[k], rows*sizeof(DT_), location, 0);
      }
    }

    // cudaCheck(cudaMemsetAsync(_c_d, 0, size_t(krylov_dim)*sizeof(DT_)), "memsetAsync");
    // cudaCheck(cudaMemsetAsync(_s_d, 0, size_t(krylov_dim)*sizeof(DT_)), "memsetAsync");
    cudaCheck(cudaMemsetAsync(_q_d, 0, size_t(krylov_dim+1)*sizeof(DT_)), "memsetAsync");

    cudaCheck(cudaMemcpyAsync(vec_v_d, vec_v, size_t(krylov_dim+1)*sizeof(DT_*), cudaMemcpyHostToDevice), "memcopyAsync");

    std::uninitialized_fill(_c, _c+krylov_dim, DT_(0));
    std::uninitialized_fill(_s, _s+krylov_dim, DT_(0));
    std::uninitialized_fill(_q, _q+krylov_dim+1u, DT_(0));
    std::uninitialized_fill(_h, _h+gauss_size, DT_(0));

    cudaDeviceSynchronize();
    Util::cuda_check_last_error();

    // initialize defect
    DT_ def_init = DT_(0);
    DT_ def_prev = DT_(0);
    DT_ def_cur = DT_(0);


    // calculate norm
    Intern::cublas_dot(def_init, vec_v[0], vec_v[0], rows, dt);
    cudaDeviceSynchronize();
    Util::cuda_check_last_error();

    def_init = Math::sqrt(def_init);
    status.def_init = def_init;


    def_prev = def_cur = def_init;

    status.def_norm = def_cur;
    status.num_iter = 0u;
    FEAT::Solver::Status solver_status = FEAT::Solver::Status::progress;

    Memory::TypedView<IT_> filter_idx_view(apply_info.def_filter.ind->view(Memory::Location::cuda, Memory::Access::read));
    const IT_* filter_idx = filter_idx_view.get_r();
    const Index filter_ue = apply_info.def_filter.used_elements;
    const Index blocksize = Util::cuda_blocksize_misc;
    const Index filter_gridsize = Index(ceil((filter_ue)/(double)(blocksize)));


    if(!Math::isfinite(def_init))
      solver_status = FEAT::Solver::Status::aborted;

    Util::cuda_stop_range(id_setup);

    const auto id_outer = Util::cuda_start_range("GMRES Outer Loop");
    while(solver_status == FEAT::Solver::Status::progress)
    {
      // reset h device
      cudaCheck(cudaMemsetAsync(_h_d, 0, size_t(gauss_size)*sizeof(DT_)), "memsetAsync");
      // set first _q
      _q[0] = def_cur;
      // cudaCheck(cudaMemcpyAsync(_q_d, _q, sizeof(DT_), cudaMemcpyHostToDevice), "cudaMemcopyAsync");
      // normalize v[0]
      Intern::cublas_scale(vec_v[0], DT_(1)/_q[0], rows, dt);
      // reset h
      std::fill(_h, _h+gauss_size, DT_(0));

      const auto id_inner = Util::cuda_start_range("Inner GMRES Loop");
      Index i(0);
      Index h_offset(0);
      while(i < krylov_dim)
      {
        const auto first_section = Util::cuda_start_range("First Section");
        // apply vanka matrix
        {
          cu_status = cusparseDnVecSetValues(dx, (void*)vec_v[i]);
          if (cu_status != CUSPARSE_STATUS_SUCCESS)
            throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
          cu_status = cusparseDnVecSetValues(dr, (void*)vec_w);
          if (cu_status != CUSPARSE_STATUS_SUCCESS)
            throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));

          Intern::apply_cu_spmv(trans, ct, dx, dr, cumat_vanka, vanka_alpha, vanka_beta);
          // and filter
          Intern::apply_def_filter(vec_w, filter_idx, filter_ue, blocksize, filter_gridsize);
        }

        // calculate defect
        {
          // reuse dr
          cu_status = cusparseDnVecSetValues(dx, (void*)vec_v[i+1]);
          if (cu_status != CUSPARSE_STATUS_SUCCESS)
            throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
          Intern::apply_cu_spmv(trans, ct, dr, dx, cumat_sys, DT_(1), DT_(0));
          Intern::apply_def_filter(vec_v[i+1], filter_idx, filter_ue, blocksize, filter_gridsize);
        }

        // Gram-Schmidt process (classical / conditional CGS2) with the self-dot fused into the
        // batched reduction: request i+2 dots so slot i+1 (vec_v_d has krylov_dim+1 pointers) holds
        // <v[i+1], v[i+1]>, from which the orthogonalized norm is recovered algebraically.
        DT_ alpha = DT_(0);
        {
          // the batched dot kernel accumulates atomically, so zero the scratch first
          cudaCheck(cudaMemsetAsync(batch_dots_d, 0, size_t(i+2)*sizeof(DT_)), "memsetAsync");
          Intern::calculate_batched_dots(batch_dots_d, vec_v[i+1], vec_v_d, rows, i+2);
          cudaCheck(cudaMemcpyAsync(batch_dots_h, batch_dots_d, size_t(i+2)*sizeof(DT_), cudaMemcpyDeviceToHost), "cudaMemcopyAsync");
          cudaDeviceSynchronize();
          Util::cuda_check_last_error();
          // subtract the projection (reads the leading i+1 scale entries)
          Intern::add_vecs(vec_v[i+1], vec_v_d, batch_dots_d, rows, i+1, DT_(-1));

          // store projections into the Hessenberg column and recover the norm algebraically
          DT_ sum_sqr = DT_(0);
          for(Index k = 0; k < i+1; ++k)
          {
            _h[h_offset + k] = batch_dots_h[k];
            sum_sqr += Math::sqr(batch_dots_h[k]);
          }
          const DT_ self_dot = batch_dots_h[i+1];
          DT_ nrm_sqr = self_dot - sum_sqr;

          // conditional corrective reorthogonalization (cgs2); the second batched reduction supplies
          // both the corrections and the updated self-dot, so no extra reduction is needed
          if((apply_info.ortho_type == FEAT::Solver::KrylovOrthogonalizerType::cgs2)
            && (apply_info.force_reortho || (nrm_sqr < Math::sqr(apply_info.reortho_tol) * self_dot)))
          {
            cudaCheck(cudaMemsetAsync(batch_dots_d, 0, size_t(i+2)*sizeof(DT_)), "memsetAsync");
            Intern::calculate_batched_dots(batch_dots_d, vec_v[i+1], vec_v_d, rows, i+2);
            cudaCheck(cudaMemcpyAsync(batch_dots_h, batch_dots_d, size_t(i+2)*sizeof(DT_), cudaMemcpyDeviceToHost), "cudaMemcopyAsync");
            cudaDeviceSynchronize();
            Util::cuda_check_last_error();
            Intern::add_vecs(vec_v[i+1], vec_v_d, batch_dots_d, rows, i+1, DT_(-1));
            DT_ sum_sqr2 = DT_(0);
            for(Index k = 0; k < i+1; ++k)
            {
              _h[h_offset + k] += batch_dots_h[k];
              sum_sqr2 += Math::sqr(batch_dots_h[k]);
            }
            nrm_sqr = batch_dots_h[i+1] - sum_sqr2;
          }

          // recover the norm; fall back to an explicit dot on heavy cancellation (guarantees a
          // positive result so the subsequent normalization cannot produce a NaN)
          if(nrm_sqr <= Math::Limits<DT_>::epsilon() * DT_(100))
          {
            Intern::cublas_dot(nrm_sqr, vec_v[i+1], vec_v[i+1], rows, dt);
            cudaDeviceSynchronize();
            Util::cuda_check_last_error();
            nrm_sqr = Math::abs(nrm_sqr);
          }
          alpha = Math::sqrt(nrm_sqr);
        }
        Util::cuda_stop_range(first_section);
        { // open "Second Section" scope; RAII so the convergence/max-iter break below still closes it
        FEAT::Util::ScopedNvtxRange _second_section("Second Section");
        Intern::cublas_scale(vec_v[i+1], DT_(1)/alpha, rows, dt);
        {
          const auto id_givens = Util::cuda_start_range("Givens Rotation");
          // apply Givens rotations
          for(Index k(0); k < i; ++k)
          {
            DT_* h_i = _h + h_offset;
            const DT_ t(h_i[k]);
            h_i[k] = _c[k] * t + _s[k] * h_i[k+1];
            h_i[k+1] = _s[k] * t - _c[k] * h_i[k+1];
          }
          Util::cuda_stop_range(id_givens);
        }

        const auto id_def_push = Util::cuda_start_range("Push defect");
        // compute beta
        const DT_ beta = Math::sqrt(Math::sqr(_h[h_offset+i]) + Math::sqr(alpha));

        // compute next plane rotation
        _s[i] = alpha / beta;
        _c[i] = _h[h_offset+i] / beta;

        // important: back sync h -> h device not required
        _h[h_offset+i] = beta;
        _q[i+1] = _s[i] * _q[i];
        _q[i] *= _c[i];
        h_offset += i+1;

        Util::cuda_stop_range(id_def_push);
        // push our new defect
        if(++i < krylov_dim)
        {
          // get the absolute defect
          def_cur = Math::abs(_q[i]);

          // // did we diverge?
          // if((def_cur > div_abs) || (def_cur > (div_rel * def_init)))
          //   break;

          // minimum number of iterations performed?
          if(!(status.num_iter < min_iter))
          {
            // did we converge?
            if((def_cur <= inner_res_scale * abs_tol) && ((def_cur <= inner_res_scale * (rel_tol * def_init))))
              break;

            // maximum number of iterations performed?
            if(status.num_iter >= max_iter)
              break;
          }

          ++status.num_iter;
          status.def_norm = def_cur;
          // plot?
          if((plot_mode == Solver::PlotMode::iter || plot_mode == Solver::PlotMode::all) && status.num_iter%plot_interval == 0u)
          {
            FEAT::String msg = FEAT::String("CudaGMRES")
              +  "* " + stringify(status.num_iter).pad_front(4)
              + " : " + stringify_fp_sci(def_cur);
            std::cout << msg << "\n";
          }
        }

        } // close "Second Section" scope (RAII)
      }
      Util::cuda_stop_range(id_inner);

      const auto id_update = Util::cuda_start_range("Update Solution");
      {
        const Index n = Math::min(i, krylov_dim);
        const auto id_h_solve = Util::cuda_start_range("Solve Hessenberg");
        // solve H*q = q
        for(Index k(n); k > 0;)
        {
          --k;
          const Index h_offset = (k+1)*k/2u;
          _q[k] /= _h[h_offset + k];
          for(Index j(k); j > 0;)
          {
            --j;
            _q[j] -= _h[h_offset + j] * _q[k];
          }
        }
        Util::cuda_stop_range(id_h_solve);
        // copy _q values
        cudaCheck(cudaMemcpyAsync(_q_d, _q, size_t(n)*sizeof(DT_), cudaMemcpyHostToDevice), "cudaMemcpyAsync");

        // Runtime::SyncGuard sync_guard;
        // update solution
        cudaCheck(cudaMemsetAsync(vec_w, 0, rows*sizeof(DT_)), "cudaMemcpyAsync");
        // _add_vecs(this->_vec_w, this->_vec_v, this->_q, n, DataType(1));
        Intern::add_vecs(vec_w, vec_v_d, _q_d, rows, n, DT_(1));

        // apply vanka matrix
        {
          cu_status = cusparseDnVecSetValues(dx, (void*)vec_w);
          if (cu_status != CUSPARSE_STATUS_SUCCESS)
            throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
          cu_status = cusparseDnVecSetValues(dr, (void*)vec_v[0]);
          if (cu_status != CUSPARSE_STATUS_SUCCESS)
            throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));

          Intern::apply_cu_spmv(trans, ct, dx, dr, cumat_vanka, vanka_alpha, vanka_beta);
          // and filter
          Intern::apply_def_filter(vec_v[0], filter_idx, filter_ue, blocksize, filter_gridsize);
        }
        //update solution
        Intern::cublas_axpy(vec_sol, DT_(1), vec_v[0], rows, dt);

        // compute real residual
        {
          cu_status = cusparseDnVecSetValues(dx, (void*)vec_sol);
          if (cu_status != CUSPARSE_STATUS_SUCCESS)
            throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
          // copy rhs to vec_v
          cudaCheck(cudaMemcpyAsync(vec_v[0], vec_rhs, rows*sizeof(DT_), cudaMemcpyDeviceToDevice), "cudaMemcpyAsync");
          // ath this point, sync 1 -> 0 necessary
          Intern::apply_cu_spmv(trans, ct, dx, dr, cumat_sys, defect_alpha, defect_beta);
          // and sync again...
          Intern::apply_def_filter(vec_v[0], filter_idx, filter_ue, blocksize, filter_gridsize);
        }

        // increase iteration count
        ++status.num_iter;

        def_prev = def_cur;

        // finally evaluate our defect
        Intern::cublas_dot(def_cur, vec_v[0], vec_v[0], rows, dt);
        cudaDeviceSynchronize();
        Util::cuda_check_last_error();

        def_cur = Math::sqrt(def_cur);
        status.def_norm = def_cur;

        solver_status = _analyze_defect(status, def_prev, def_init, max_iter, min_iter, abs_tol, rel_tol);
      }
      Util::cuda_stop_range(id_update);
    }
    Util::cuda_stop_range(id_outer);



    cudaCheck(cudaFreeHost(_c), "cudaFreeHost");
    cudaCheck(cudaFreeHost(_s), "cudaFreeHost");
    cudaCheck(cudaFreeHost(_q), "cudaFreeHost");
    cudaCheck(cudaFreeHost(_h), "cudaFreeHost");
    cudaCheck(cudaFreeHost(batch_dots_h), "cudaFreeHost");
    cudaCheck(cudaFree(batch_dots_d), "cudaFree");
    // cudaCheck(cudaFree(_c_d), "cudaFree");
    // cudaCheck(cudaFree(_s_d), "cudaFree");
    cudaCheck(cudaFree(_q_d), "cudaFree");
    cudaCheck(cudaFree(_h_d), "cudaFree");
    cudaCheck(cudaFree(vec_v_d), "cudaFree");




    Intern::cumat_free(cumat_sys);
    Intern::cumat_free(cumat_vanka);
    cusparseDestroyDnVec(dx);
    cusparseDestroyDnVec(dr);

    return solver_status;
  }


  template<typename DT_, typename MT_>
  CUDA_HOST void _setup_comm_buffer(std::vector<DT_*>& host_buffer, std::vector<DT_*>& device_buffer, const std::vector<MT_>& mirrors)
  {
    host_buffer.resize(mirrors.size()*2u);
    device_buffer.resize(mirrors.size());

    for(Index k = 0; k < mirrors.size(); ++k)
    {
      cudaCheck(cudaMallocHost(&host_buffer[2*k], mirrors[k].num_indices() * sizeof(DT_)), "cudaMallocHost");
      cudaCheck(cudaMallocHost(&host_buffer[2*k+1], mirrors[k].num_indices() * sizeof(DT_)), "cudaMallocHost");
      cudaCheck(cudaMalloc(&device_buffer[k], mirrors[k].num_indices() * sizeof(DT_)), "cudaMalloc");
      cudaCheck(cudaMemsetAsync(device_buffer[k], 0, mirrors[k].num_indices() * sizeof(DT_)), "memsetAsync");
    }

    for(Index k = 0; k < mirrors.size(); ++k)
    {
      std::uninitialized_fill(host_buffer[2*k], host_buffer[2*k]+mirrors[k].num_indices(), DT_(0));
      std::uninitialized_fill(host_buffer[2*k+1], host_buffer[2*k+1]+mirrors[k].num_indices(), DT_(0));
    }
  }

  /// global r = a*A*x   product   be careful about setting b correctly, since this factor adds r
  template<typename DT_, typename IT_, typename Gate_, std::size_t st_>
  CUDA_HOST void _apply_mat(std::vector<DT_*>& buffer_host, std::vector<DT_*>& buffer_device, DT_* r, cusparseDnVecDescr_t dr, const DT_* x, cusparseDnVecDescr_t dx,
                    const Gate_& gate, cusparseOperation_t op, cudaDataType_t ct, const CuMatDesc& mat_desc, const CSRHelper::CRSMat<DT_, IT_>& mat, DT_ alpha, DT_ beta, cudaStream_t main_stream, const std::array<cudaStream_t, st_>& streams)
  {
    cusparseStatus_t cu_status = cusparseDnVecSetValues(dx, (void*)x);
    if (cu_status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
    cu_status = cusparseDnVecSetValues(dr, (void*)r);
    if (cu_status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));

    Intern::apply_synced_mat(buffer_host, buffer_device, r, dr, x, dx, gate, op, ct, mat_desc, mat, alpha, beta, main_stream, streams);

  }

  /// global r = a*A*x + b*y   product
  template<typename DT_, typename IT_, typename Gate_, std::size_t st_>
  CUDA_HOST void _apply_mat(std::vector<DT_*>& buffer_host, std::vector<DT_*>& buffer_device, DT_* r, cusparseDnVecDescr_t dr, const DT_* x, cusparseDnVecDescr_t dx, const DT_* y,
                    const Gate_& gate, cusparseOperation_t op, cudaDataType_t ct, const CuMatDesc& mat_desc, const CSRHelper::CRSMat<DT_, IT_>& mat, DT_ alpha, DT_ beta, cudaStream_t main_stream, const std::array<cudaStream_t, st_>& streams)
  {
    // cusparseStatus_t cu_status = cusparseDnVecSetValues(dx, (void*)x);
    // if (cu_status != CUSPARSE_STATUS_SUCCESS)
    //   throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
    // cu_status = cusparseDnVecSetValues(dr, (void*)r);
    // if (cu_status != CUSPARSE_STATUS_SUCCESS)
    //   throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));

    // scale vector with frequency
    Memory::TypedView<DT_> freqs_view(gate.get_freqs().elements_view_r(Memory::Location::cuda));
    const DT_* freqs = freqs_view.get_r();
    Intern::scale_vec(r, y, freqs, mat.num_cols);

    // and now call standard apply
    _apply_mat(buffer_host, buffer_device, r, dr, x, dx, gate, op, ct, mat_desc, mat, alpha, beta, main_stream, streams);
  }

  // for schwarz type synch
  template<typename DT_, typename IT_, typename Gate_, std::size_t st_>
  CUDA_HOST void _apply_mat_mean(std::vector<DT_*>& buffer_host, std::vector<DT_*>& buffer_device, DT_* r, cusparseDnVecDescr_t dr, const DT_* x, cusparseDnVecDescr_t dx,
                    const Gate_& gate, cusparseOperation_t op, cudaDataType_t ct, const CuMatDesc& mat_desc, const CSRHelper::CRSMat<DT_, IT_>& mat, DT_ alpha, DT_ beta, cudaStream_t main_stream, const std::array<cudaStream_t, st_>& streams)
  {
    // cusparseStatus_t cu_status = cusparseDnVecSetValues(dx, (void*)x);
    // if (cu_status != CUSPARSE_STATUS_SUCCESS)
    //   throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
    // cu_status = cusparseDnVecSetValues(dr, (void*)r);
    // if (cu_status != CUSPARSE_STATUS_SUCCESS)
    //   throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
    // call standard apply
    // cudaMemset(r, 0, mat.num_rows * sizeof(DT_));
    _apply_mat(buffer_host, buffer_device, r, dr, x, dx, gate, op, ct, mat_desc, mat, alpha, beta, main_stream, streams);
        // cudaDeviceSynchronize();
        // cudaCheckLastErr()
        //     if(i == 0)
        //     {
        //       String mt;
        //       for(Index k = 0; k < rows; ++k)
        //       {
        //         mt += stringify(r[k]) + "  ";
        //       }
        //       std::cout << "Rank " << comm.rank() << " :" << mt << "\n-------------\n";
        //     }

    // and now scale with vector
    Memory::TypedView<DT_> freqs_view(gate.get_freqs().elements_view_r(Memory::Location::cuda));
    const DT_* freqs = freqs_view.get_r();
    Intern::scale_vec(r, r, freqs, mat.num_rows);
  }



  template<typename DT_, typename IT_>
  FEAT::Solver::Status cuda_apply_amavanka_gmres_parallel(Memory::Arbiter& vec_sol_arb, const Memory::Arbiter& vec_rhs_arb, const CSRHelper::GMRESApplyInfo<DT_, IT_>& apply_info, CSRHelper::GMRESStatus<DT_>& status)
  {
    XASSERTM(apply_info.gate, "Requires allocated gate");
    using namespace FEAT;
    const auto id_setup = Util::cuda_start_range("GMRES Setup");
    // extract info
    const Index rows = apply_info.size;
    const Index columns = apply_info.size;

    // --- acquire device (cuda) views for all arbiter-backed data; kept alive for the whole solve ---
    Memory::TypedView<DT_> vec_sol_view(vec_sol_arb.view(Memory::Location::cuda, Memory::Access::read_write));
    DT_* vec_sol = vec_sol_view.get_w();
    Memory::TypedView<DT_> vec_rhs_view(vec_rhs_arb.view(Memory::Location::cuda, Memory::Access::read));
    const DT_* vec_rhs = vec_rhs_view.get_r();

    Memory::View sys_val_view(apply_info.mat_sys.val->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View sys_row_view(apply_info.mat_sys.row->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View sys_col_view(apply_info.mat_sys.col->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View vanka_val_view(apply_info.mat_vanka.val->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View vanka_row_view(apply_info.mat_vanka.row->view(Memory::Location::cuda, Memory::Access::read));
    Memory::View vanka_col_view(apply_info.mat_vanka.col->view(Memory::Location::cuda, Memory::Access::read));

    const std::vector<Memory::Arbiter*>& vec_v_arbs = *apply_info.vec_v;
    std::vector<Memory::TypedView<DT_>> vec_v_views;
    vec_v_views.reserve(vec_v_arbs.size());
    std::vector<DT_*> vec_v_host(vec_v_arbs.size());
    for(std::size_t k = 0; k < vec_v_arbs.size(); ++k)
    {
      vec_v_views.emplace_back(vec_v_arbs[k]->view(Memory::Location::cuda, Memory::Access::read_write));
      vec_v_host[k] = vec_v_views[k].get_w();
    }
    // const Index used_elements_def = Index(apply_info.mat_sys.row[rows]);
    // const Index used_elements_vanka = Index(apply_info.mat_vanka.row[rows]);
    // create necessary data arrays on host and device
    // setup workspace for cusparse matrices
    cusparseOperation_t trans = CUSPARSE_OPERATION_NON_TRANSPOSE;
    cudaDataType dt;
    cudaDataType ct; //compute type
    if (typeid(DT_) == typeid(double))
    {
        dt = CUDA_R_64F;
        ct = CUDA_R_64F;
    }
    else if (typeid(DT_) == typeid(float))
    {
        dt = CUDA_R_32F;
        ct = CUDA_R_32F;
    }
    cusparseIndexType_t it;
    if(sizeof(IT_) == 4u)
      it = CUSPARSE_INDEX_32I;
    else if(sizeof(IT_) == 8u)
      it = CUSPARSE_INDEX_64I;
    cusparseStatus_t cu_status;
    const DT_ defect_alpha = DT_(-1);
    const DT_ defect_beta = DT_(1);
    const DT_ vanka_alpha = DT_(1);
    const DT_ vanka_beta = DT_(0);
    const DT_ abs_tol = apply_info.abs_tol;
    const DT_ rel_tol = apply_info.rel_tol;
    const DT_ inner_res_scale = apply_info.res_scale;
    const Index krylov_dim = apply_info.krylov_dim;
    const Index max_iter = apply_info.max_iter;
    const Index min_iter = apply_info.min_iter;

    const Index plot_interval = apply_info.plot_interval;
    const Solver::PlotMode plot_mode = apply_info.plot_mode;
    const auto& gate = *apply_info.gate;
    const auto& comm = *gate.get_comm();
    const auto& ranks = gate.get_ranks();
    const auto& mirrors = gate.get_mirrors();
    Memory::TypedView<DT_> freqs_view(gate.get_freqs().elements_view_r(Memory::Location::cuda));
    const DT_* const freqs = freqs_view.get_r();
    const bool use_multiple_streams = true;

    // allocate buffers for our synchoperations
    std::vector<DT_*> mirror_buffer_h;
    std::vector<DT_*> mirror_buffer_d;

    _setup_comm_buffer(mirror_buffer_h, mirror_buffer_d, mirrors);

    // our "main" execution stream
    cudaStream_t mainstream = 0;
    if(use_multiple_streams)
      cudaCheck(cudaStreamCreate(&mainstream), "cudaStreamCreate");
    // else
    //   mainstream = 0;
    // we need some streams (for now lets do 4)
    std::array<cudaStream_t, 4> worker_streams = {mainstream, mainstream, mainstream, mainstream};
    for(auto& s : worker_streams)
    {
      if(use_multiple_streams)
        cudaCheck(cudaStreamCreate(&s), "cudasStreamCreatre");
    }

    // for cusparse, we always use our mainstream
    cusparseSetStream(Util::Intern::cusparse_handle, mainstream);


    cusparseDnVecDescr_t dx;
    {
      cu_status = cusparseCreateDnVec(&dx, columns, (void*)vec_sol, dt);
      if (cu_status != CUSPARSE_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
    }

    cusparseDnVecDescr_t dr;
    {
      cu_status = cusparseCreateDnVec(&dr, columns, (void*)vec_rhs, dt);
      if (cu_status != CUSPARSE_STATUS_SUCCESS)
        throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(cu_status)));
    }

    CuMatDesc cumat_sys;
    Intern::alloc_cusparse_interface<DT_, IT_>(trans, dt, ct, it, dx, dr, sys_row_view.raw_r(), sys_col_view.raw_r(), sys_val_view.raw_r(),
      apply_info.mat_sys.num_rows, apply_info.mat_sys.num_cols, apply_info.mat_sys.num_nnz, cumat_sys, defect_alpha, defect_beta);

    CuMatDesc cumat_vanka;
    Intern::alloc_cusparse_interface<DT_, IT_>(trans, dt, ct, it, dx, dr, vanka_row_view.raw_r(), vanka_col_view.raw_r(), vanka_val_view.raw_r(),
      apply_info.mat_vanka.num_rows, apply_info.mat_vanka.num_cols, apply_info.mat_vanka.num_nnz, cumat_vanka, vanka_alpha, vanka_beta);

    // allocations for helper arrays
    const Index gauss_size = krylov_dim * (krylov_dim+1u) / 2u;
    DT_ *_c, *_s, *_q, *_h;
    _c = new DT_[size_t(krylov_dim)];
    _s = new DT_[size_t(krylov_dim)];
    // cudaCheck(cudaMallocHost(&_c, size_t(krylov_dim)*sizeof(DT_)), "cudaMallocHost");
    // cudaCheck(cudaMallocHost(&_s, size_t(krylov_dim)*sizeof(DT_)), "cudaMallocHost");
    cudaCheck(cudaMallocHost(&_q, size_t(krylov_dim+1)*sizeof(DT_)), "cudaMallocHost");
    cudaCheck(cudaMallocHost(&_h, size_t(gauss_size)*sizeof(DT_)), "cudaMallocHost");
    // DT_ *_c_d, *_s_d, *_q_d, *_h_d;
    DT_ *_q_d, *_h_d;
    // cudaCheck(cudaMalloc(&_c_d, size_t(krylov_dim)*sizeof(DT_)), "cudaMalloc");
    // cudaCheck(cudaMalloc(&_s_d, size_t(krylov_dim)*sizeof(DT_)), "cudaMalloc");
    cudaCheck(cudaMalloc(&_q_d, size_t(krylov_dim+1)*sizeof(DT_)), "cudaMalloc");
    cudaCheck(cudaMalloc(&_h_d, size_t(gauss_size)*sizeof(DT_)), "cudaMalloc");
    DT_** vec_v = vec_v_host.data();
    Memory::TypedView<DT_> vec_w_view(apply_info.vec_w->view(Memory::Location::cuda, Memory::Access::read_write));
    DT_* vec_w = vec_w_view.get_w();
    DT_** vec_v_d;
    cudaCheck(cudaMalloc(&vec_v_d, size_t(krylov_dim+1)*sizeof(DT_*)), "cudaMalloc");

    DT_* dot_helper, *dot_helper_d;
    cudaCheck(cudaMallocHost(&dot_helper, sizeof(DT_)), "cudaMallocHost");
    cudaCheck(cudaMalloc(&dot_helper_d, sizeof(DT_)), "cudaMalloc");

    // scratch for the batched projection dots plus the appended self-dot (krylov_dim projections + 1)
    DT_ *batch_dots_h, *batch_dots_d;
    cudaCheck(cudaMallocHost(&batch_dots_h, size_t(krylov_dim+1)*sizeof(DT_)), "cudaMallocHost");
    cudaCheck(cudaMalloc(&batch_dots_d, size_t(krylov_dim+1)*sizeof(DT_)), "cudaMalloc");

    if(apply_info.prefetch_data)
    {
      cudaMemLocation location;
      {location.type = cudaMemLocationTypeDevice; location.id = FEAT::Util::cuda_get_device_id();}
      // hint to prefetch vector data, if possible
      for(Index k = 0; k < krylov_dim+1; ++k)
      {
        cudaMemPrefetchAsync(vec_v[k], rows*sizeof(DT_), location, 0);
      }
    }

    // cudaCheck(cudaMemsetAsync(_c_d, 0, size_t(krylov_dim)*sizeof(DT_)), "memsetAsync");
    // cudaCheck(cudaMemsetAsync(_s_d, 0, size_t(krylov_dim)*sizeof(DT_)), "memsetAsync");
    cudaCheck(cudaMemsetAsync(_q_d, 0, size_t(krylov_dim+1)*sizeof(DT_)), "memsetAsync");

    cudaCheck(cudaMemcpyAsync(vec_v_d, vec_v, size_t(krylov_dim+1)*sizeof(DT_*), cudaMemcpyHostToDevice), "memcopyAsync");

    std::uninitialized_fill(_c, _c+krylov_dim, DT_(0));
    std::uninitialized_fill(_s, _s+krylov_dim, DT_(0));
    std::uninitialized_fill(_q, _q+krylov_dim+1u, DT_(0));
    std::uninitialized_fill(_h, _h+gauss_size, DT_(0));

    cudaDeviceSynchronize();
    Util::cuda_check_last_error();

    // initialize defect
    DT_ def_init = DT_(0);
    DT_ def_prev = DT_(0);
    DT_ def_cur = DT_(0);


    // calculate norm
    {
      auto id_dot = Util::cuda_start_range("Single Dot");
      // finally evaluate our defect
      Intern::cuda_sync_dot(dot_helper, dot_helper_d, vec_v[0], vec_v[0], freqs, rows, comm, mainstream);
      Util::cuda_stop_range(id_dot);
    }

    def_init = Math::sqrt(*dot_helper);
    status.def_init = def_init;


    def_prev = def_cur = def_init;

    status.def_norm = def_cur;
    status.num_iter = 0u;
    FEAT::Solver::Status solver_status = FEAT::Solver::Status::progress;

    Memory::TypedView<IT_> filter_idx_view(apply_info.def_filter.ind->view(Memory::Location::cuda, Memory::Access::read));
    const IT_* filter_idx = filter_idx_view.get_r();
    const Index filter_ue = apply_info.def_filter.used_elements;
    const Index blocksize = Util::cuda_blocksize_misc;
    const Index filter_gridsize = Index(ceil((filter_ue)/(double)(blocksize)));


    if(!Math::isfinite(def_init))
      solver_status = FEAT::Solver::Status::aborted;

    Util::cuda_stop_range(id_setup);

    const auto id_outer = Util::cuda_start_range("GMRES Outer Loop");
    while(solver_status == FEAT::Solver::Status::progress)
    {
      // reset h device
      cudaCheck(cudaMemsetAsync(_h_d, 0, size_t(gauss_size)*sizeof(DT_)), "memsetAsync");
      // set first _q
      _q[0] = def_cur;
      // cudaCheck(cudaMemcpyAsync(_q_d, _q, sizeof(DT_), cudaMemcpyHostToDevice), "cudaMemcopyAsync");
      // normalize v[0]
      Intern::cublas_scale(vec_v[0], DT_(1)/_q[0], rows, dt);
      // reset h
      std::fill(_h, _h+gauss_size, DT_(0));

      const auto id_inner = Util::cuda_start_range("Inner GMRES Loop");
      Index i(0);
      Index h_offset(0);
      while(i < krylov_dim)
      {
        // cudaDeviceSynchronize();
        // cudaCheckLastErr()
        //     if(i == 0)
        //     {
        //       String mt;
        //       for(Index k = 0; k < rows; ++k)
        //       {
        //         mt += stringify(vec_v[0][k]) + "  ";
        //       }
        //       std::cout << "Rank " << comm.rank() << " :" << mt << "\n-------------\n";
        //     }
        const auto first_section = Util::cuda_start_range("First Section");
        // apply vanka matrix
        {
          auto id_vanka = Util::cuda_start_range("Vanka Apply");
          _apply_mat_mean(mirror_buffer_h, mirror_buffer_d, vec_w, dr, vec_v[i], dx, gate, trans, ct, cumat_vanka, apply_info.mat_vanka, vanka_alpha, vanka_beta, mainstream, worker_streams);
          // and filter
          Intern::apply_def_filter(vec_w, filter_idx, filter_ue, blocksize, filter_gridsize);
          Util::cuda_stop_range(id_vanka);
        }
        // cudaDeviceSynchronize();
        // cudaCheckLastErr()
        //     if(i == 0)
        //     {
        //       String mt;
        //       for(Index k = 0; k < rows; ++k)
        //       {
        //         mt += stringify(vec_w[k]) + "  ";
        //       }
        //       std::cout << "Rank " << comm.rank() << " :" << mt << "\n-------------\n";
        //     }

        // calculate defect
        {
          auto id_defect = Util::cuda_start_range("Defect Calc");
          // reuse dr
          // cudaMemset(vec_v[i+1], 0, rows * sizeof(DT_));
          _apply_mat(mirror_buffer_h, mirror_buffer_d, vec_v[i+1], dr, vec_w, dx, gate, trans, ct, cumat_sys, apply_info.mat_sys, DT_(1), DT_(0), mainstream, worker_streams);
          Intern::apply_def_filter(vec_v[i+1], filter_idx, filter_ue, blocksize, filter_gridsize);
          Util::cuda_stop_range(id_defect);
        }

        // fused batched Gram-Schmidt (classical / conditional CGS2): i+1 frequency-weighted
        // projections plus the appended self-dot, reduced in a SINGLE all-reduce over i+2 entries
        // (previously the norm needed a second all-reduce per inner step)
        DT_ alpha = DT_(0);
        {
          auto id_batched = Util::cuda_start_range("Batched Dot");
          Intern::calculate_synced_batched_dots(batch_dots_h, batch_dots_d, vec_v[i+1], vec_v_d, freqs, rows, i+2, comm, mainstream);
          // copy the leading i+1 (synchronized) projections back to device for the batched axpy
          cudaCheck(cudaMemcpyAsync(batch_dots_d, batch_dots_h, size_t(i+1)*sizeof(DT_), cudaMemcpyHostToDevice, mainstream), "cudaMemcpyAsync");
          cudaStreamSynchronize(mainstream);
          Util::cuda_stop_range(id_batched);
          Intern::add_vecs(vec_v[i+1], vec_v_d, batch_dots_d, rows, i+1, DT_(-1));

          // store projections into the Hessenberg column and recover the norm algebraically
          DT_ sum_sqr = DT_(0);
          for(Index k = 0; k < i+1; ++k)
          {
            _h[h_offset + k] = batch_dots_h[k];
            sum_sqr += Math::sqr(batch_dots_h[k]);
          }
          const DT_ self_dot = batch_dots_h[i+1];
          DT_ nrm_sqr = self_dot - sum_sqr;

          // conditional corrective reorthogonalization (cgs2); the second batched all-reduce supplies
          // both the corrections and the updated self-dot, so no extra reduction is needed
          if((apply_info.ortho_type == FEAT::Solver::KrylovOrthogonalizerType::cgs2)
            && (apply_info.force_reortho || (nrm_sqr < Math::sqr(apply_info.reortho_tol) * self_dot)))
          {
            Intern::calculate_synced_batched_dots(batch_dots_h, batch_dots_d, vec_v[i+1], vec_v_d, freqs, rows, i+2, comm, mainstream);
            cudaCheck(cudaMemcpyAsync(batch_dots_d, batch_dots_h, size_t(i+1)*sizeof(DT_), cudaMemcpyHostToDevice, mainstream), "cudaMemcpyAsync");
            cudaStreamSynchronize(mainstream);
            Intern::add_vecs(vec_v[i+1], vec_v_d, batch_dots_d, rows, i+1, DT_(-1));
            DT_ sum_sqr2 = DT_(0);
            for(Index k = 0; k < i+1; ++k)
            {
              _h[h_offset + k] += batch_dots_h[k];
              sum_sqr2 += Math::sqr(batch_dots_h[k]);
            }
            nrm_sqr = batch_dots_h[i+1] - sum_sqr2;
          }

          // recover the norm; fall back to an explicit synced dot on heavy cancellation (guarantees
          // a positive result so the subsequent normalization cannot produce a NaN)
          if(nrm_sqr <= Math::Limits<DT_>::epsilon() * DT_(100))
          {
            auto id_dot = Util::cuda_start_range("Single Dot");
            Intern::cuda_sync_dot(dot_helper, dot_helper_d, vec_v[i+1], vec_v[i+1], freqs, rows, comm, mainstream);
            Util::cuda_stop_range(id_dot);
            nrm_sqr = Math::abs(*dot_helper);
          }
          alpha = Math::sqrt(nrm_sqr);
        }
        Util::cuda_stop_range(first_section);
        { // open "Second Section" scope; RAII so the convergence/max-iter break below still closes it
        FEAT::Util::ScopedNvtxRange _second_section("Second Section");
        Intern::cublas_scale(vec_v[i+1], DT_(1)/alpha, rows, dt);
        {
          const auto id_givens = Util::cuda_start_range("Givens Rotation");
          // apply Givens rotations
          for(Index k(0); k < i; ++k)
          {
            DT_* h_i = _h + h_offset;
            const DT_ t(h_i[k]);
            h_i[k] = _c[k] * t + _s[k] * h_i[k+1];
            h_i[k+1] = _s[k] * t - _c[k] * h_i[k+1];
          }
          Util::cuda_stop_range(id_givens);
        }

        const auto id_def_push = Util::cuda_start_range("Push defect");
        // compute beta
        const DT_ beta = Math::sqrt(Math::sqr(_h[h_offset+i]) + Math::sqr(alpha));

        // compute next plane rotation
        _s[i] = alpha / beta;
        _c[i] = _h[h_offset+i] / beta;

        // important: back sync h -> h device not required
        _h[h_offset+i] = beta;
        _q[i+1] = _s[i] * _q[i];
        _q[i] *= _c[i];
        h_offset += i+1;

        Util::cuda_stop_range(id_def_push);
        // push our new defect
        if(++i < krylov_dim)
        {
          // get the absolute defect
          def_cur = Math::abs(_q[i]);

          // // did we diverge?
          // if((def_cur > div_abs) || (def_cur > (div_rel * def_init)))
          //   break;

          // minimum number of iterations performed?
          if(!(status.num_iter < min_iter))
          {
            // did we converge?
            if((def_cur <= inner_res_scale * abs_tol) && ((def_cur <= inner_res_scale * (rel_tol * def_init))))
              break;

            // maximum number of iterations performed?
            if(status.num_iter >= max_iter)
              break;
          }

          ++status.num_iter;
          status.def_norm = def_cur;
          // plot?
          if((plot_mode == Solver::PlotMode::iter || plot_mode == Solver::PlotMode::all) && status.num_iter%plot_interval == 0u)
          {
            FEAT::String msg = FEAT::String("CudaGMRES")
              +  "* " + stringify(status.num_iter).pad_front(4)
              + " : " + stringify_fp_sci(def_cur);
            comm.print(msg);
          }
        }

        } // close "Second Section" scope (RAII)
      }
      Util::cuda_stop_range(id_inner);

      const auto id_update = Util::cuda_start_range("Update Solution");
      {
        const Index n = Math::min(i, krylov_dim);
        const auto id_h_solve = Util::cuda_start_range("Solve Hessenberg");
        // solve H*q = q
        for(Index k(n); k > 0;)
        {
          --k;
          const Index h_offset = (k+1)*k/2u;
          _q[k] /= _h[h_offset + k];
          for(Index j(k); j > 0;)
          {
            --j;
            _q[j] -= _h[h_offset + j] * _q[k];
          }
        }
        Util::cuda_stop_range(id_h_solve);
        // copy _q values
        cudaCheck(cudaMemcpyAsync(_q_d, _q, size_t(n)*sizeof(DT_), cudaMemcpyHostToDevice), "cudaMemcpyAsync");

        // Runtime::SyncGuard sync_guard;
        // update solution
        cudaCheck(cudaMemsetAsync(vec_w, 0, rows*sizeof(DT_)), "cudaMemcpyAsync");
        Intern::add_vecs(vec_w, vec_v_d, _q_d, rows, n, DT_(1));

        // apply vanka matrix
        {
          _apply_mat_mean(mirror_buffer_h, mirror_buffer_d, vec_v[0], dr, vec_w, dx, gate, trans, ct, cumat_vanka, apply_info.mat_vanka, vanka_alpha, vanka_beta, mainstream, worker_streams);
          Intern::apply_def_filter(vec_v[0], filter_idx, filter_ue, blocksize, filter_gridsize);
        }
        //update solution
        Intern::cublas_axpy(vec_sol, DT_(1), vec_v[0], rows, dt);

        // compute real residual
        {
          _apply_mat(mirror_buffer_h, mirror_buffer_d, vec_v[0], dr, vec_sol, dx, vec_rhs, gate, trans, ct, cumat_sys, apply_info.mat_sys, defect_alpha, defect_beta, mainstream, worker_streams);
          Intern::apply_def_filter(vec_v[0], filter_idx, filter_ue, blocksize, filter_gridsize);
        }

        // increase iteration count
        ++status.num_iter;

        def_prev = def_cur;

        {
          auto id_dot = Util::cuda_start_range("Single Dot");
          // finally evaluate our defect
          Intern::cuda_sync_dot(dot_helper, dot_helper_d, vec_v[0], vec_v[0], freqs, rows, comm, mainstream);
          Util::cuda_stop_range(id_dot);
        }

        def_cur = Math::sqrt(*dot_helper);
        status.def_norm = def_cur;

        solver_status = _analyze_defect(status, def_prev, def_init, max_iter, min_iter, abs_tol, rel_tol);
      }
      Util::cuda_stop_range(id_update);
    }
    Util::cuda_stop_range(id_outer);



    // cudaCheck(cudaFreeHost(_c), "cudaFreeHost");
    // cudaCheck(cudaFreeHost(_s), "cudaFreeHost");
    delete[](_c);
    delete[](_s);
    cudaCheck(cudaFreeHost(dot_helper), "cudaFreeHost");
    cudaCheck(cudaFreeHost(batch_dots_h), "cudaFreeHost");
    cudaCheck(cudaFree(batch_dots_d), "cudaFree");
    cudaCheck(cudaFreeHost(_q), "cudaFreeHost");
    cudaCheck(cudaFreeHost(_h), "cudaFreeHost");
    // cudaCheck(cudaFree(_c_d), "cudaFree");
    // cudaCheck(cudaFree(_s_d), "cudaFree");
    cudaCheck(cudaFree(_q_d), "cudaFree");
    cudaCheck(cudaFree(_h_d), "cudaFree");
    cudaCheck(cudaFree(vec_v_d), "cudaFree");
    cudaCheck(cudaFree(dot_helper_d), "cudaFree");


    // free buffer
    for(auto k : mirror_buffer_h)
    {
      cudaCheck(cudaFreeHost(k), "cudaFreeHost");
    }
    mirror_buffer_h.clear();
    for(auto k : mirror_buffer_d)
    {
      cudaCheck(cudaFree(k), "cudaFree");
    }
    mirror_buffer_d.clear();





    Intern::cumat_free(cumat_sys);
    Intern::cumat_free(cumat_vanka);
    cusparseDestroyDnVec(dx);
    cusparseDestroyDnVec(dr);

    // reset cusparse context
    cusparseSetStream(Util::Intern::cusparse_handle, 0);

    if(use_multiple_streams)
    {
      // destroy streams
      for(auto& s : worker_streams)
      {
        cudaCheck(cudaStreamDestroy(s), "cudasStreamDestroy");
      }
      cudaCheck(cudaStreamDestroy(mainstream), "cudasStreamDestroy");
    }



    return solver_status;
  }
}

template FEAT::Solver::Status FEAT::Solver::VankaGMRES::cuda_apply_amavanka_gmres_seq<double, std::uint32_t>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, const FEAT::Solver::CSRHelper::GMRESApplyInfo<double, std::uint32_t>&, FEAT::Solver::CSRHelper::GMRESStatus<double>&);
template FEAT::Solver::Status FEAT::Solver::VankaGMRES::cuda_apply_amavanka_gmres_seq<float, std::uint32_t>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, const FEAT::Solver::CSRHelper::GMRESApplyInfo<float, std::uint32_t>&, FEAT::Solver::CSRHelper::GMRESStatus<float>&);
template FEAT::Solver::Status FEAT::Solver::VankaGMRES::cuda_apply_amavanka_gmres_seq<double, std::uint64_t>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, const FEAT::Solver::CSRHelper::GMRESApplyInfo<double, std::uint64_t>&, FEAT::Solver::CSRHelper::GMRESStatus<double>&);
template FEAT::Solver::Status FEAT::Solver::VankaGMRES::cuda_apply_amavanka_gmres_seq<float, std::uint64_t>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, const FEAT::Solver::CSRHelper::GMRESApplyInfo<float, std::uint64_t>&, FEAT::Solver::CSRHelper::GMRESStatus<float>&);

template FEAT::Solver::Status FEAT::Solver::VankaGMRES::cuda_apply_amavanka_gmres_parallel<double, std::uint32_t>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, const FEAT::Solver::CSRHelper::GMRESApplyInfo<double, std::uint32_t>&, FEAT::Solver::CSRHelper::GMRESStatus<double>&);
template FEAT::Solver::Status FEAT::Solver::VankaGMRES::cuda_apply_amavanka_gmres_parallel<float, std::uint32_t>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, const FEAT::Solver::CSRHelper::GMRESApplyInfo<float, std::uint32_t>&, FEAT::Solver::CSRHelper::GMRESStatus<float>&);
template FEAT::Solver::Status FEAT::Solver::VankaGMRES::cuda_apply_amavanka_gmres_parallel<double, std::uint64_t>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, const FEAT::Solver::CSRHelper::GMRESApplyInfo<double, std::uint64_t>&, FEAT::Solver::CSRHelper::GMRESStatus<double>&);
template FEAT::Solver::Status FEAT::Solver::VankaGMRES::cuda_apply_amavanka_gmres_parallel<float, std::uint64_t>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, const FEAT::Solver::CSRHelper::GMRESApplyInfo<float, std::uint64_t>&, FEAT::Solver::CSRHelper::GMRESStatus<float>&);
