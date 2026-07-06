// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/batched_kernel_helper.cuh>
#include <kernel/lafem/arch/batched_vector_apply.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/memory_arbiter.hpp>

#include <cooperative_groups.h>

#include <array>

#include <kernel/util/cuda_helper.cuh>

namespace FEAT::LAFEM::Arch
{
  namespace cg = cooperative_groups;

  namespace
  {
    // Persistent, lazily-grown device scratch reused across the batched-kernel calls, so that the
    // per-call cudaMalloc/cudaFree (each of which synchronizes the device) is avoided. There are two
    // arenas: one for the value scratch (device_dots / device_scales) and one for the device pointer
    // table (other_ptrs). Reuse is safe because all batched ops are issued serially on the default
    // stream, so a buffer is only overwritten after the prior kernel that read it has completed
    // (stream ordering); the host orchestration is single-threaded. The buffers grow monotonically and
    // are released at CUDA teardown by the runtime; they are only used by the GMRES driver, not by the
    // leak-checked kernel unit tests.
    void*  g_scratch_val = nullptr;
    size_t g_scratch_val_bytes = 0u;
    void*  g_scratch_ptr = nullptr;
    size_t g_scratch_ptr_bytes = 0u;

    void* grow_scratch(void*& buf, size_t& cap, size_t bytes)
    {
      if(cap < bytes)
      {
        if(buf != nullptr)
          Memory::free_cuda(buf);
        size_t new_bytes = 2u*bytes;
        buf = Memory::alloc_cuda(new_bytes);
        cap = new_bytes;
      }
      return buf;
    }

    void free_scratch_space()
    {
      if(g_scratch_val)
      {
        Memory::free_cuda(g_scratch_val);
        g_scratch_val_bytes = 0u;
      }
      if(g_scratch_ptr)
      {
        Memory::free_cuda(g_scratch_ptr);
        g_scratch_ptr_bytes = 0u;
      }
    }
  }

  void BatchedScratchSpace::free()
  {
    free_scratch_space();
  }



  template<typename DT_, typename DTT_>
  void BatchedDotProduct::exec_cuda(Memory::Arbiter& local_dots, const Memory::Arbiter& target, Index dot_size, Memory::Arbiter const * const * others, Index num_dots)
  {
    constexpr Index max_size_ = 1024u;
    constexpr Index blocksize = 256u;
    // constexpr Index blocksize = 32u;
    constexpr int num_batch = 4;
    // constexpr int num_batch = 16;

    std::array<Memory::TypedView<DTT_>, max_size_> static_views;
    std::array<DTT_ const*, max_size_> static_p;

    dim3 grid;
    dim3 block;
    block.x = (unsigned)blocksize;
    block.y = num_batch;
    // x grid dimension resoluts dot_size, modulo we want to load multiple values per
    // TODO: change with int variants
    grid.x = (unsigned)ceil((dot_size)/(double)(block.x*2*10));
    // grid.x = (unsigned)1;
    grid.y = (unsigned)ceil((num_dots)/(double)(block.y));

    // the cuda kernel operates on device memory, so it acquires cuda read views itself
    Memory::TypedView<DTT_> target_view(target.view(Memory::Location::cuda, Memory::Access::read));
    const DTT_* target_ptr = target_view.get_r();
    Memory::TypedView<DTT_>* tmp_views = static_views.data();
    const DTT_ ** tmp_p = static_p.data();
    if(num_dots > max_size_)
    {
      using HType = const DTT_*;
      tmp_views = new Memory::TypedView<DTT_>[num_dots];
      tmp_p = new HType[num_dots];
    }
    for(Index l = 0; l < num_dots; ++l)
    {
      tmp_views[l] = Memory::TypedView<DTT_>(others[l]->view(Memory::Location::cuda, Memory::Access::read));
      tmp_p[l] = tmp_views[l].get_r();
    }

    // device pointer table for the "others" (transient device scratch)
    DTT_ const ** other_ptrs = static_cast<DTT_ const**>(grow_scratch(g_scratch_ptr, g_scratch_ptr_bytes, num_dots*sizeof(DTT_*)));
    cudaCheck(cudaMemcpyAsync(other_ptrs, tmp_p, num_dots * sizeof(DTT_*), cudaMemcpyHostToDevice), "cudaMemcpyAsync");

    if constexpr(std::is_same_v<DT_, DTT_>)
    {
      Memory::TypedView<DT_> dots_view(local_dots.view(Memory::Location::cuda, Memory::Access::read_write));
      DT_* device_dots = dots_view.get_w();
      cuda_batched_dot<DTT_, blocksize><<<grid, block, blocksize*num_batch*sizeof(DTT_)>>>(device_dots, target_ptr, dot_size, other_ptrs, num_dots);
    }
    else
    {
      // mixed precision: accumulate in the compute type DTT_ on a scratch buffer, then convert-add into
      // the DT_ result arbiter on the host
      DTT_* device_dots = static_cast<DTT_*>(grow_scratch(g_scratch_val, g_scratch_val_bytes, num_dots*sizeof(DTT_)));
      Memory::format_cuda(device_dots, DTT_(0), num_dots);
      cuda_batched_dot<DTT_, blocksize><<<grid, block, blocksize*num_batch*sizeof(DTT_)>>>(device_dots, target_ptr, dot_size, other_ptrs, num_dots);
      DTT_* tmp = static_cast<DTT_*>(Memory::alloc_main(sizeof(DTT_)*num_dots));
      Memory::memcopy_cuda_to_main(tmp, device_dots, num_dots * sizeof(DTT_));
      Memory::TypedView<DT_> dots_view(local_dots.view(Memory::Location::main, Memory::Access::read_write));
      DT_* dots_ptr = dots_view.get_w();
      for(Index l = 0; l < num_dots; ++l)
      {
        dots_ptr[l] += DT_(tmp[l]);
      }
      Memory::free_main(tmp);
    }
    if(num_dots > max_size_)
    {
      delete[](tmp_p);
      delete[](tmp_views);
    }

    #ifdef FEAT_DEBUG_MODE
      Util::cuda_synchronize();
      cudaError_t last_error(cudaGetLastError());
      if (cudaSuccess != last_error)
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
    #endif

  }

  template<typename DT_, typename DTT_>
  void BatchedTripleDotProduct::exec_cuda(Memory::Arbiter& local_dots, const Memory::Arbiter& target, Index dot_size, Memory::Arbiter const * const * others, const Memory::Arbiter& freqs, Index num_dots)
  {
    constexpr Index max_size_ = 1024u;
    constexpr Index blocksize = 256u;
    // constexpr Index blocksize = 32u;
    constexpr int num_batch = 4;
    // constexpr int num_batch = 16;

    std::array<Memory::TypedView<DTT_>, max_size_> static_views;
    std::array<DTT_ const*, max_size_> static_p;



    dim3 grid;
    dim3 block;
    block.x = (unsigned)blocksize;
    block.y = num_batch;
    // x grid dimension resoluts dot_size, modulo we want to load multiple values per
    grid.x = (unsigned)ceil((dot_size)/(double)(block.x*2*10));
    // grid.x = (unsigned)1;
    grid.y = (unsigned)ceil((num_dots)/(double)(block.y));

    // the cuda kernel operates on device memory, so it acquires cuda read views itself
    Memory::TypedView<DTT_> target_view(target.view(Memory::Location::cuda, Memory::Access::read));
    const DTT_* target_ptr = target_view.get_r();
    Memory::TypedView<DTT_> freqs_view(freqs.view(Memory::Location::cuda, Memory::Access::read));
    const DTT_* freqs_ptr = freqs_view.get_r();
    Memory::TypedView<DTT_>* tmp_views = static_views.data();
    const DTT_ ** tmp_p = static_p.data();
    if(num_dots > max_size_)
    {
      using HType = const DTT_*;
      tmp_views = new Memory::TypedView<DTT_>[num_dots];
      tmp_p = new HType[num_dots];
    }
    for(Index l = 0; l < num_dots; ++l)
    {
      tmp_views[l] = Memory::TypedView<DTT_>(others[l]->view(Memory::Location::cuda, Memory::Access::read));
      tmp_p[l] = tmp_views[l].get_r();
    }

    // device pointer table for the "others" (transient device scratch)
    DTT_ const ** other_ptrs = static_cast<DTT_ const**>(grow_scratch(g_scratch_ptr, g_scratch_ptr_bytes, num_dots*sizeof(DTT_*)));
    cudaCheck(cudaMemcpyAsync(other_ptrs, tmp_p, num_dots * sizeof(DTT_*), cudaMemcpyHostToDevice), "cudaMemcpyAsync");

    if constexpr(std::is_same_v<DT_, DTT_>)
    {
      // accumulate straight into the caller's result arbiter on the device (see batched_dot_cuda)
      Memory::TypedView<DT_> dots_view(local_dots.view(Memory::Location::cuda, Memory::Access::read_write));
      DT_* device_dots = dots_view.get_w();
      cuda_batched_triple_dot<DTT_, blocksize><<<grid, block, blocksize*num_batch*sizeof(DTT_)>>>(device_dots, target_ptr, dot_size, other_ptrs, freqs_ptr, num_dots);
    }
    else
    {
      // mixed precision: accumulate in DTT_ on scratch, then convert-add into the DT_ result arbiter
      DTT_* device_dots = static_cast<DTT_*>(grow_scratch(g_scratch_val, g_scratch_val_bytes, num_dots*sizeof(DTT_)));
      Memory::format_cuda(device_dots, DTT_(0), num_dots);
      cuda_batched_triple_dot<DTT_, blocksize><<<grid, block, blocksize*num_batch*sizeof(DTT_)>>>(device_dots, target_ptr, dot_size, other_ptrs, freqs_ptr, num_dots);
      DTT_* tmp = static_cast<DTT_*>(Memory::alloc_main(sizeof(DTT_)*num_dots));
      Memory::memcopy_cuda_to_main(tmp, device_dots, num_dots * sizeof(DTT_));
      Memory::TypedView<DT_> dots_view(local_dots.view(Memory::Location::main, Memory::Access::read_write));
      DT_* dots_ptr = dots_view.get_w();
      for(Index l = 0; l < num_dots; ++l)
      {
        dots_ptr[l] += DT_(tmp[l]);
      }
      Memory::free_main(tmp);
    }

    if(num_dots > max_size_)
    {
      delete[](tmp_p);
      delete[](tmp_views);
    }

    #ifdef FEAT_DEBUG_MODE
      Util::cuda_synchronize();
      cudaError_t last_error(cudaGetLastError());
      if (cudaSuccess != last_error)
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
    #endif

  }


  template<typename DT_>
  void CumulativeAxpy::exec_cuda(Memory::Arbiter& target, Index vec_size, Memory::Arbiter const * const * others, DT_ const * const local_scales, Index num_axpy, DT_ alpha)
  {
    constexpr Index max_size_ = 1024u;
    constexpr Index blocksize = 128u;
    dim3 grid;
    dim3 block;
    block.x = (unsigned)blocksize;
    // one block always does a warp size chunk of target (due to load coalescing)
    grid.x = (unsigned)ceil((vec_size)/(double)((32)));

    std::array<Memory::TypedView<DT_>, max_size_> static_views;
    std::array<DT_ const*, max_size_> static_p;

    // the cuda kernel operates on device memory, so it acquires its cuda views itself
    Memory::TypedView<DT_> target_view(target.view(Memory::Location::cuda, Memory::Access::read_write));
    DT_* target_ptr = target_view.get_w();

    Memory::TypedView<DT_>* tmp_views = static_views.data();
    const DT_ ** tmp_p = static_p.data();
    if(num_axpy > max_size_)
    {
      using HType = const DT_*;
      tmp_views = new Memory::TypedView<DT_>[num_axpy];
      tmp_p = new HType[num_axpy];
    }
    for(Index l = 0; l < num_axpy; ++l)
    {
      tmp_views[l] = Memory::TypedView<DT_>(others[l]->view(Memory::Location::cuda, Memory::Access::read));
      tmp_p[l] = tmp_views[l].get_r();
    }


    // reuse the persistent scratch arenas instead of allocating/freeing per call
    DT_* device_scales = static_cast<DT_*>(grow_scratch(g_scratch_val, g_scratch_val_bytes, num_axpy*sizeof(DT_)));
    DT_ const ** other_ptrs = static_cast<DT_ const**>(grow_scratch(g_scratch_ptr, g_scratch_ptr_bytes, num_axpy*sizeof(DT_*)));
    // stream ordered calls, so no explicit sync required
    cudaCheck(cudaMemcpyAsync(other_ptrs, tmp_p, num_axpy * sizeof(DT_*), cudaMemcpyHostToDevice), "cudaMemcpyAsync");

    cudaCheck(cudaMemcpyAsync(device_scales, local_scales, num_axpy * sizeof(DT_), cudaMemcpyHostToDevice), "cudaMemcpyAsync");

    cuda_cumulative_axpy<DT_, blocksize><<<grid, block, blocksize*sizeof(DT_)>>>(target_ptr, vec_size, other_ptrs, device_scales, num_axpy, DT_(alpha));

    // release and destroy the in-place constructed views (see the placement-new above)
    if(num_axpy > max_size_)
    {
      delete[](tmp_p);
      delete[](tmp_views);
    }

    #ifdef FEAT_DEBUG_MODE
      Util::cuda_synchronize();
      cudaError_t last_error(cudaGetLastError());
      if (cudaSuccess != last_error)
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
    #endif

  }


}

template void FEAT::LAFEM::Arch::BatchedDotProduct::exec_cuda<double, double>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, FEAT::Index);
template void FEAT::LAFEM::Arch::BatchedDotProduct::exec_cuda<double, float>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, FEAT::Index);
template void FEAT::LAFEM::Arch::BatchedDotProduct::exec_cuda<float, double>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, FEAT::Index);
template void FEAT::LAFEM::Arch::BatchedDotProduct::exec_cuda<float, float>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, FEAT::Index);

template void FEAT::LAFEM::Arch::BatchedTripleDotProduct::exec_cuda<double, double>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, const FEAT::Memory::Arbiter&, FEAT::Index);
template void FEAT::LAFEM::Arch::BatchedTripleDotProduct::exec_cuda<double, float>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, const FEAT::Memory::Arbiter&, FEAT::Index);
template void FEAT::LAFEM::Arch::BatchedTripleDotProduct::exec_cuda<float, double>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, const FEAT::Memory::Arbiter&, FEAT::Index);
template void FEAT::LAFEM::Arch::BatchedTripleDotProduct::exec_cuda<float, float>(FEAT::Memory::Arbiter&, const FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, const FEAT::Memory::Arbiter&, FEAT::Index);

template void FEAT::LAFEM::Arch::CumulativeAxpy::exec_cuda<double>(FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, double const *const, FEAT::Index, double);
template void FEAT::LAFEM::Arch::CumulativeAxpy::exec_cuda<float>(FEAT::Memory::Arbiter&, FEAT::Index, FEAT::Memory::Arbiter const * const *, float const *const, FEAT::Index, float);