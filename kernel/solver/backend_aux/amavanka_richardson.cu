// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/cuda_helper.cuh>
#include <kernel/util/dist.hpp>
#include <kernel/lafem/arch/matvecmult_bcsr_dense.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/backend_aux/richardson_apply_helper.hpp>
#include <kernel/solver/amavanka_richardson.hpp>

#include <cuda/atomic>

#include <array>
#include <vector>

namespace FEAT::Solver::VankaRichardson
{
  namespace Intern
  {
    // ---- small device kernels (grid-stride loops) --------------------------------------------

    template<typename DT_>
    __global__ void cuda_axpy(DT_* r, const DT_ alpha, const DT_* x, const Index n)
    {
      for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < n; i += Index(blockDim.x * gridDim.x))
        r[i] += alpha * x[i];
    }

    template<typename DT_, typename IT_>
    __global__ void cuda_unit_filter_zero(DT_* v, const IT_* idx, const Index ue)
    {
      for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < ue; i += Index(blockDim.x * gridDim.x))
        v[idx[i]] = DT_(0);
    }

    template<typename DT_>
    __global__ void cuda_comp_prod(DT_* v, const DT_* f, const Index n)
    {
      for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < n; i += Index(blockDim.x * gridDim.x))
        v[i] *= f[i];
    }

    /// blocked/scalar mirror gather: buf[i*bs + d] = v[idx[i]*bs + d]
    template<typename DT_, typename IT_>
    __global__ void cuda_gather(DT_* __restrict__ buf, const IT_* __restrict__ idx, const Index n_nodes, const DT_* __restrict__ v, const int bs)
    {
      const Index total = n_nodes * Index(bs);
      for(Index t = threadIdx.x + blockDim.x * blockIdx.x; t < total; t += Index(blockDim.x * gridDim.x))
      {
        const Index i = t / Index(bs);
        const Index d = t % Index(bs);
        buf[t] = v[Index(idx[i]) * Index(bs) + d];
      }
    }

    /// blocked/scalar mirror scatter-add (atomic, allows asynchronous merge): v[idx[i]*bs + d] += buf[i*bs + d]
    template<typename DT_, typename IT_>
    __global__ void cuda_scatter_add(DT_* __restrict__ v, const DT_* __restrict__ buf, const IT_* __restrict__ idx, const Index n_nodes, const int bs)
    {
      const Index total = n_nodes * Index(bs);
      for(Index t = threadIdx.x + blockDim.x * blockIdx.x; t < total; t += Index(blockDim.x * gridDim.x))
      {
        const Index i = t / Index(bs);
        const Index d = t % Index(bs);
        cuda::atomic_ref<DT_, cuda::thread_scope_device> ref(v[Index(idx[i]) * Index(bs) + d]);
        ref.fetch_add(buf[t]);
      }
    }

    // ---- launch helpers ----------------------------------------------------------------------

    static inline dim3 grid_for(Index n, unsigned block)
    {
      return dim3((unsigned)((n + Index(block) - 1u) / Index(block)));
    }

    template<typename DT_>
    static void launch_axpy(DT_* r, DT_ alpha, const DT_* x, Index n, cudaStream_t s = 0)
    {
      if(n == Index(0)) return;
      const unsigned bs = 256u;
      cuda_axpy<DT_><<<grid_for(n, bs), bs, 0, s>>>(r, alpha, x, n);
      cudaCheckLastErr()
    }

    template<typename DT_, typename IT_>
    static void launch_filter(DT_* v, const IT_* idx, Index ue, cudaStream_t s = 0)
    {
      if(ue == Index(0)) return;
      const unsigned bs = 256u;
      cuda_unit_filter_zero<DT_, IT_><<<grid_for(ue, bs), bs, 0, s>>>(v, idx, ue);
      cudaCheckLastErr()
    }

    template<typename DT_>
    static void launch_comp_prod(DT_* v, const DT_* f, Index n, cudaStream_t s = 0)
    {
      if(n == Index(0)) return;
      const unsigned bs = 256u;
      cuda_comp_prod<DT_><<<grid_for(n, bs), bs, 0, s>>>(v, f, n);
      cudaCheckLastErr()
    }

    // ---- resolved (device-pointer) views of a meta matrix ------------------------------------

    template<typename DT_, typename IT_>
    struct ResolvedBlock
    {
      const DT_* val;
      const IT_* row;
      const IT_* col;
      void* workspace;
      int block_height, block_width;
      Index num_rows, num_cols, num_nnz;
      int out_comp, in_comp;
    };

    /// resolve the cuda views of every non-empty block of a meta matrix; keeps the views alive
    template<typename DT_, typename IT_, int NC_>
    static std::vector<ResolvedBlock<DT_, IT_>> resolve_matrix(
      const RichardsonHelper::MetaMatrix<DT_, IT_, NC_, NC_>& mat,
      std::vector<Memory::TypedView<DT_>>& val_views,
      std::vector<Memory::TypedView<IT_>>& idx_views)
    {
      std::vector<ResolvedBlock<DT_, IT_>> out;
      for(int r = 0; r < NC_; ++r)
      {
        for(int c = 0; c < NC_; ++c)
        {
          const auto& b = mat.at(r, c);
          if(b.empty())
            continue;
          val_views.emplace_back(b.val->view(Memory::Location::cuda, Memory::Access::read));
          idx_views.emplace_back(b.row->view(Memory::Location::cuda, Memory::Access::read));
          idx_views.emplace_back(b.col->view(Memory::Location::cuda, Memory::Access::read));
          ResolvedBlock<DT_, IT_> rb;
          rb.val = val_views.back().get_r();
          rb.row = idx_views[idx_views.size() - 2u].get_r();
          rb.col = idx_views[idx_views.size() - 1u].get_r();
          rb.workspace = b.workspace;
          rb.block_height = b.block_height;
          rb.block_width = b.block_width;
          rb.num_rows = b.num_rows;
          rb.num_cols = b.num_cols;
          rb.num_nnz = b.num_nnz;
          rb.out_comp = r;
          rb.in_comp = c;
          out.push_back(rb);
        }
      }
      return out;
    }

    /**
     * \brief Meta SpMV accumulate: out += alpha * Mat * in.
     *
     * \c out must be pre-initialized (to the rhs for a defect, or to zero for a preconditioner
     * apply); every live block accumulates in place via exec_cuda_impl with y == r. This is the sole
     * "how the SpMV runs" seam — a future fused meta-kernel replaces this body only.
     */
    template<typename DT_, typename IT_>
    static void spmv_accumulate(DT_* const* out_ptr, const DT_* const* in_ptr,
      const std::vector<ResolvedBlock<DT_, IT_>>& blocks, DT_ alpha)
    {
      for(const auto& b : blocks)
      {
        DT_* r = out_ptr[b.out_comp];
        const DT_* x = in_ptr[b.in_comp];
        // r = y + alpha * A * x with y == r  =>  r += alpha * A * x (native (bh x bw) BCSR SpMV),
        // using the block's persistent cuSPARSE workspace (required by the BSR SpMV path on CUDA >= 13)
        LAFEM::Arch::MatVecMultBCSRDense::exec_cuda_impl<DT_, IT_>(
          b.workspace, r, alpha, x, /*y*/ r, b.val, b.col, b.row,
          b.num_rows, b.num_cols, b.num_nnz, b.block_height, b.block_width, false);
      }
    }
  } // namespace Intern

  // ---- single-rank (no synchronization) driver ----------------------------------------------

  template<typename DT_, typename IT_, int NC_>
  FEAT::Solver::Status cuda_apply_amavanka_richardson_seq(
    const RichardsonHelper::RichardsonApplyInfo<DT_, IT_, NC_>& info,
    RichardsonHelper::RichardsonStatus<DT_>& status)
  {
    using namespace Intern;

    // resolve the persistent component vectors once (pointers are stable across iterations)
    std::vector<Memory::TypedView<DT_>> keep; // keeps all vector views alive
    keep.reserve(std::size_t(NC_) * 4u);
    std::array<DT_*, NC_> sol{}, def{}, cor{};
    std::array<const DT_*, NC_> rhs{};
    std::array<Index, NC_> rawlen{};
    for(int c = 0; c < NC_; ++c)
    {
      rawlen[c] = info.components[c].size_raw();
      keep.emplace_back(info.vec_sol[c]->view(Memory::Location::cuda, Memory::Access::read_write));
      sol[c] = keep.back().get_w();
      keep.emplace_back(info.vec_def[c]->view(Memory::Location::cuda, Memory::Access::read_write));
      def[c] = keep.back().get_w();
      keep.emplace_back(info.vec_cor[c]->view(Memory::Location::cuda, Memory::Access::read_write));
      cor[c] = keep.back().get_w();
      keep.emplace_back(info.vec_rhs[c]->view(Memory::Location::cuda, Memory::Access::read));
      rhs[c] = keep.back().get_r();
    }

    // resolve the unit-filter index arrays
    std::vector<Memory::TypedView<IT_>> filter_views;
    filter_views.reserve(NC_);
    std::array<const IT_*, NC_> filt_idx{};
    std::array<Index, NC_> filt_ue{};
    for(int c = 0; c < NC_; ++c)
    {
      filt_ue[c] = info.def_filters[c].used_elements;
      if(filt_ue[c] > Index(0) && info.def_filters[c].ind != nullptr)
      {
        filter_views.emplace_back(info.def_filters[c].ind->view(Memory::Location::cuda, Memory::Access::read));
        filt_idx[c] = filter_views.back().get_r();
      }
    }

    // resolve the matrix blocks once
    std::vector<Memory::TypedView<DT_>> sys_val, van_val;
    std::vector<Memory::TypedView<IT_>> sys_idx, van_idx;
    const auto sys_blocks = resolve_matrix<DT_, IT_, NC_>(info.mat_sys, sys_val, sys_idx);
    const auto van_blocks = resolve_matrix<DT_, IT_, NC_>(info.mat_vanka, van_val, van_idx);

    const Index steps = Math::max(info.max_iter, Index(1));
    for(Index step = 0; step < steps; ++step)
    {
      // d = rhs - A*u : start from d = rhs, then accumulate -A*u in place
      for(int c = 0; c < NC_; ++c)
        cudaCheck(cudaMemcpyAsync(def[c], rhs[c], rawlen[c] * sizeof(DT_), cudaMemcpyDeviceToDevice), "def<-rhs");
      spmv_accumulate<DT_, IT_>(def.data(), sol.data(), sys_blocks, -DT_(1));
      for(int c = 0; c < NC_; ++c)
        launch_filter<DT_, IT_>(def[c], filt_idx[c], filt_ue[c]);

      // cor = V*d : start from cor = 0, then accumulate +V*d in place
      for(int c = 0; c < NC_; ++c)
        cudaCheck(cudaMemsetAsync(cor[c], 0, rawlen[c] * sizeof(DT_)), "cor<-0");
      spmv_accumulate<DT_, IT_>(cor.data(), def.data(), van_blocks, DT_(1));
      for(int c = 0; c < NC_; ++c)
        launch_filter<DT_, IT_>(cor[c], filt_idx[c], filt_ue[c]);

      // u += omega * cor
      for(int c = 0; c < NC_; ++c)
        launch_axpy<DT_>(sol[c], info.omega, cor[c], rawlen[c]);
    }
    cudaCheck(cudaDeviceSynchronize(), "richardson seq sync");

    status.num_iter = steps;
    status.def_norm = DT_(0);
    status.def_init = DT_(0);
    return FEAT::Solver::Status::success;
  }

  // ---- distributed driver: per-component hand-rolled sync_0 / sync_1 -------------------------

  namespace Intern
  {
    /// per-component communication buffers for one gate
    template<typename DT_>
    struct CommBuffers
    {
      std::vector<DT_*> dev;       // one device buffer per neighbor
      std::vector<DT_*> host_send; // pinned host send buffer per neighbor
      std::vector<DT_*> host_recv; // pinned host recv buffer per neighbor
    };

    template<typename DT_, typename IT_>
    static CommBuffers<DT_> setup_buffers(const RichardsonHelper::GateWrapper<DT_, IT_>& g, int bs)
    {
      CommBuffers<DT_> b;
      const std::size_t n = g.ranks.size();
      b.dev.resize(n, nullptr);
      b.host_send.resize(n, nullptr);
      b.host_recv.resize(n, nullptr);
      for(std::size_t k = 0; k < n; ++k)
      {
        if(g.mirror_size[k] == Index(0)) // neighbor with no shared DOFs for this component (e.g. discontinuous pressure)
          continue;
        const std::size_t bytes = g.mirror_size[k] * Index(bs) * sizeof(DT_);
        cudaCheck(cudaMalloc(&b.dev[k], bytes), "cudaMalloc");
        cudaCheck(cudaMallocHost(&b.host_send[k], bytes), "cudaMallocHost");
        cudaCheck(cudaMallocHost(&b.host_recv[k], bytes), "cudaMallocHost");
      }
      return b;
    }

    template<typename DT_>
    static void free_buffers(CommBuffers<DT_>& b)
    {
      for(auto* p : b.dev) if(p) cudaFree(p);
      for(auto* p : b.host_send) if(p) cudaFreeHost(p);
      for(auto* p : b.host_recv) if(p) cudaFreeHost(p);
      b.dev.clear(); b.host_send.clear(); b.host_recv.clear();
    }

    /// additive type-0 synchronization of one component vector \p v over its gate
    template<typename DT_, typename IT_>
    static void sync_0_component(DT_* v, int bs, const RichardsonHelper::GateWrapper<DT_, IT_>& g,
      CommBuffers<DT_>& buf, const std::vector<const IT_*>& mirror_ptr)
    {
      const std::size_t n = g.ranks.size();
      if(n == 0u)
        return;
      Dist::Comm& comm = *g.comm;
      const unsigned block = 128u;

      // gather local shared values into the send buffers and copy to host
      for(std::size_t k = 0; k < n; ++k)
      {
        const Index nodes = g.mirror_size[k];
        if(nodes == Index(0)) // no shared DOFs with this neighbor for this component
          continue;
        const Index total = nodes * Index(bs);
        cuda_gather<DT_, IT_><<<grid_for(total, block), block>>>(buf.dev[k], mirror_ptr[k], nodes, v, bs);
        cudaCheckLastErr()
        cudaCheck(cudaMemcpyAsync(buf.host_send[k], buf.dev[k], total * sizeof(DT_), cudaMemcpyDeviceToHost), "gather D2H");
      }
      cudaCheck(cudaDeviceSynchronize(), "sync_0 gather");

      // exchange (only neighbors that actually share DOFs for this component); req_k maps a request
      // slot back to its neighbor index for the wait_any-driven scatter below
      Dist::RequestVector recv_reqs, send_reqs;
      std::vector<std::size_t> req_k;
      recv_reqs.reserve(n); send_reqs.reserve(n); req_k.reserve(n);
      for(std::size_t k = 0; k < n; ++k)
      {
        if(g.mirror_size[k] == Index(0))
          continue;
        recv_reqs.push_back(comm.irecv(buf.host_recv[k], g.mirror_size[k] * Index(bs), g.ranks[k]));
        req_k.push_back(k);
      }
      for(std::size_t k = 0; k < n; ++k)
      {
        if(g.mirror_size[k] == Index(0))
          continue;
        send_reqs.push_back(comm.isend(buf.host_send[k], g.mirror_size[k] * Index(bs), g.ranks[k]));
      }

      // scatter-add received neighbor contributions back
      for(std::size_t ridx = 0; recv_reqs.wait_any(ridx); )
      {
        const std::size_t k = req_k[ridx];
        const Index nodes = g.mirror_size[k];
        const Index total = nodes * Index(bs);
        cudaCheck(cudaMemcpyAsync(buf.dev[k], buf.host_recv[k], total * sizeof(DT_), cudaMemcpyHostToDevice), "recv H2D");
        cudaCheck(cudaDeviceSynchronize(), "recv H2D sync");
        cuda_scatter_add<DT_, IT_><<<grid_for(total, block), block>>>(v, buf.dev[k], mirror_ptr[k], nodes, bs);
        cudaCheckLastErr()
      }
      cudaCheck(cudaDeviceSynchronize(), "sync_0 scatter");
      send_reqs.wait_all();
    }
  } // namespace Intern

  template<typename DT_, typename IT_, int NC_>
  FEAT::Solver::Status cuda_apply_amavanka_richardson_parallel(
    const RichardsonHelper::RichardsonApplyInfo<DT_, IT_, NC_>& info,
    RichardsonHelper::RichardsonStatus<DT_>& status)
  {
    using namespace Intern;

    // resolve component vectors + rawlen + block sizes
    std::vector<Memory::TypedView<DT_>> keep;
    keep.reserve(std::size_t(NC_) * 5u);
    std::array<DT_*, NC_> sol{}, def{}, cor{};
    std::array<const DT_*, NC_> rhs{};
    std::array<const DT_*, NC_> freqs{};
    std::array<Index, NC_> rawlen{};
    std::array<int, NC_> bsz{};
    for(int c = 0; c < NC_; ++c)
    {
      bsz[c] = info.components[c].block_size;
      rawlen[c] = info.components[c].size_raw();
      keep.emplace_back(info.vec_sol[c]->view(Memory::Location::cuda, Memory::Access::read_write));
      sol[c] = keep.back().get_w();
      keep.emplace_back(info.vec_def[c]->view(Memory::Location::cuda, Memory::Access::read_write));
      def[c] = keep.back().get_w();
      keep.emplace_back(info.vec_cor[c]->view(Memory::Location::cuda, Memory::Access::read_write));
      cor[c] = keep.back().get_w();
      keep.emplace_back(info.vec_rhs[c]->view(Memory::Location::cuda, Memory::Access::read));
      rhs[c] = keep.back().get_r();
      if(!info.gates[c].empty() && info.gates[c].freqs != nullptr)
      {
        keep.emplace_back(info.gates[c].freqs->view(Memory::Location::cuda, Memory::Access::read));
        freqs[c] = keep.back().get_r();
      }
    }

    // resolve unit filters
    std::vector<Memory::TypedView<IT_>> filter_views;
    std::array<const IT_*, NC_> filt_idx{};
    std::array<Index, NC_> filt_ue{};
    for(int c = 0; c < NC_; ++c)
    {
      filt_ue[c] = info.def_filters[c].used_elements;
      if(filt_ue[c] > Index(0) && info.def_filters[c].ind != nullptr)
      {
        filter_views.emplace_back(info.def_filters[c].ind->view(Memory::Location::cuda, Memory::Access::read));
        filt_idx[c] = filter_views.back().get_r();
      }
    }

    // resolve matrix blocks
    std::vector<Memory::TypedView<DT_>> sys_val, van_val;
    std::vector<Memory::TypedView<IT_>> sys_idx, van_idx;
    const auto sys_blocks = resolve_matrix<DT_, IT_, NC_>(info.mat_sys, sys_val, sys_idx);
    const auto van_blocks = resolve_matrix<DT_, IT_, NC_>(info.mat_vanka, van_val, van_idx);

    // resolve mirror index arrays + set up per-component comm buffers (once)
    std::vector<Memory::TypedView<IT_>> mirror_views;
    std::array<std::vector<const IT_*>, NC_> mirror_ptr;
    std::array<CommBuffers<DT_>, NC_> buffers;
    for(int c = 0; c < NC_; ++c)
    {
      const auto& g = info.gates[c];
      for(std::size_t k = 0; k < g.ranks.size(); ++k)
      {
        // a hollow mirror (no shared DOFs for this component, e.g. discontinuous pressure) has a
        // null index arbiter — keep the slot aligned with ranks but do not resolve/use it
        if(g.mirror_size[k] > Index(0) && g.mirror_idx[k] != nullptr)
        {
          mirror_views.emplace_back(g.mirror_idx[k]->view(Memory::Location::cuda, Memory::Access::read));
          mirror_ptr[c].push_back(mirror_views.back().get_r());
        }
        else
          mirror_ptr[c].push_back(nullptr);
      }
      if(!g.empty())
        buffers[c] = setup_buffers<DT_, IT_>(g, bsz[c]);
    }

    const Index steps = Math::max(info.max_iter, Index(1));
    for(Index step = 0; step < steps; ++step)
    {
      // d = rhs - A*u (local), then sync_0 the additive (type-0) result, then filter.
      // Note: the defect is reconciled with sync_0 (no frequency weighting); frequency scaling
      // (from_1_to_0) belongs only to the correction's sync_1 below.
      for(int c = 0; c < NC_; ++c)
        cudaCheck(cudaMemcpyAsync(def[c], rhs[c], rawlen[c] * sizeof(DT_), cudaMemcpyDeviceToDevice), "def<-rhs");
      for(int c = 0; c < NC_; ++c)
      {
        if(!info.gates[c].empty())
        {
          launch_comp_prod<DT_>(def[c], freqs[c], rawlen[c]); // from_1_to_0
        }
      }
      spmv_accumulate<DT_, IT_>(def.data(), sol.data(), sys_blocks, -DT_(1));
      cudaCheck(cudaDeviceSynchronize(), "sys spmv");
      for(int c = 0; c < NC_; ++c)
        sync_0_component<DT_, IT_>(def[c], bsz[c], info.gates[c], buffers[c], mirror_ptr[c]);
      for(int c = 0; c < NC_; ++c)
        launch_filter<DT_, IT_>(def[c], filt_idx[c], filt_ue[c]);

      // cor = V*d (local), then sync_1 (freq scale + sync_0), then filter
      for(int c = 0; c < NC_; ++c)
        cudaCheck(cudaMemsetAsync(cor[c], 0, rawlen[c] * sizeof(DT_)), "cor<-0");
      spmv_accumulate<DT_, IT_>(cor.data(), def.data(), van_blocks, DT_(1));
      cudaCheck(cudaDeviceSynchronize(), "vanka spmv");
      for(int c = 0; c < NC_; ++c)
      {
        if(!info.gates[c].empty())
        {
          launch_comp_prod<DT_>(cor[c], freqs[c], rawlen[c]); // from_1_to_0
          sync_0_component<DT_, IT_>(cor[c], bsz[c], info.gates[c], buffers[c], mirror_ptr[c]);
        }
      }
      for(int c = 0; c < NC_; ++c)
        launch_filter<DT_, IT_>(cor[c], filt_idx[c], filt_ue[c]);

      // u += omega * cor
      for(int c = 0; c < NC_; ++c)
        launch_axpy<DT_>(sol[c], info.omega, cor[c], rawlen[c]);
    }
    cudaCheck(cudaDeviceSynchronize(), "richardson parallel sync");

    for(int c = 0; c < NC_; ++c)
      free_buffers<DT_>(buffers[c]);

    status.num_iter = steps;
    status.def_norm = DT_(0);
    status.def_init = DT_(0);
    return FEAT::Solver::Status::success;
  }

  // ---- explicit instantiations (saddle-point => NC_ == 2) -----------------------------------

  template FEAT::Solver::Status cuda_apply_amavanka_richardson_seq<double, std::uint32_t, 2>(
    const RichardsonHelper::RichardsonApplyInfo<double, std::uint32_t, 2>&, RichardsonHelper::RichardsonStatus<double>&);
  template FEAT::Solver::Status cuda_apply_amavanka_richardson_seq<double, std::uint64_t, 2>(
    const RichardsonHelper::RichardsonApplyInfo<double, std::uint64_t, 2>&, RichardsonHelper::RichardsonStatus<double>&);
  template FEAT::Solver::Status cuda_apply_amavanka_richardson_seq<float, std::uint32_t, 2>(
    const RichardsonHelper::RichardsonApplyInfo<float, std::uint32_t, 2>&, RichardsonHelper::RichardsonStatus<float>&);
  template FEAT::Solver::Status cuda_apply_amavanka_richardson_seq<float, std::uint64_t, 2>(
    const RichardsonHelper::RichardsonApplyInfo<float, std::uint64_t, 2>&, RichardsonHelper::RichardsonStatus<float>&);

  template FEAT::Solver::Status cuda_apply_amavanka_richardson_parallel<double, std::uint32_t, 2>(
    const RichardsonHelper::RichardsonApplyInfo<double, std::uint32_t, 2>&, RichardsonHelper::RichardsonStatus<double>&);
  template FEAT::Solver::Status cuda_apply_amavanka_richardson_parallel<double, std::uint64_t, 2>(
    const RichardsonHelper::RichardsonApplyInfo<double, std::uint64_t, 2>&, RichardsonHelper::RichardsonStatus<double>&);
  template FEAT::Solver::Status cuda_apply_amavanka_richardson_parallel<float, std::uint32_t, 2>(
    const RichardsonHelper::RichardsonApplyInfo<float, std::uint32_t, 2>&, RichardsonHelper::RichardsonStatus<float>&);
  template FEAT::Solver::Status cuda_apply_amavanka_richardson_parallel<float, std::uint64_t, 2>(
    const RichardsonHelper::RichardsonApplyInfo<float, std::uint64_t, 2>&, RichardsonHelper::RichardsonStatus<float>&);
} // namespace FEAT::Solver::VankaRichardson
