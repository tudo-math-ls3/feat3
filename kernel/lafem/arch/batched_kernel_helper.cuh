// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/cuda_util.hpp>
#include <cuda/atomic>
#include <cooperative_groups.h>

#include <kernel/util/cuda_helper.cuh>



namespace FEAT::LAFEM::Arch
{
  namespace cg = cooperative_groups;
  extern __shared__ char array[];

  template<typename DT_, int blockSize_>
  __global__ void cuda_batched_dot(DT_* __restrict__ device_dots, const DT_ * target, const Index dot_size, DT_ const* const* others, const Index num_dots)
  {
    // our x index decides on which target we operate, while our y index decides on which other vector we operate
    // load two elements at the same time
    Index idx = threadIdx.x + blockDim.x * blockIdx.x * 2;
    const Index idy = threadIdx.y + blockDim.y * blockIdx.y;

    const Index grid_size = blockDim.x * 2 * gridDim.x;

    // since we perform cyclic reductions, we can not simply jump out, so we will have to perform empty operations
    // first of all, define our shared array
    // link to our actual, independent slice for our group
    DT_* array_dt = (DT_*)&array[0];
    DT_* loc_dot = &array_dt[blockDim.x*threadIdx.y];
    loc_dot[threadIdx.x] = DT_(0);

    // now calculate dot product
    if(idy < num_dots)
    {
      // do as many inner loops as possible before synchronizing
      while(idx < dot_size)
      {
        if((idx+blockSize_) < dot_size)
        {
          loc_dot[threadIdx.x] += target[idx] * others[idy][idx] + target[idx+blockSize_] * others[idy][idx+blockSize_];
        }
        else
        {
          loc_dot[threadIdx.x] += target[idx] * others[idy][idx];
        }
        idx += grid_size;
      }
    }

    reduce_sum<DT_, blockSize_>(&loc_dot[0], threadIdx.x);
    //TODO: actual do this in second kernel
    if(idy >= num_dots)
      return;

    cuda::atomic_ref<DT_, cuda::thread_scope_device> val(device_dots[idy]);

    // finally, we have to add to our global value (and hope that we started this kernel with such few processes that the global lock is not detrimental)
    if(threadIdx.x == 0)
    {
      val.fetch_add(loc_dot[0]);
    }
  }

  template<typename DT_, int blockSize_>
  __global__ void cuda_triple_dot(DT_* device_dot, const DT_ * const target, const Index dot_size, DT_ const* const other, DT_ const* __restrict__ const freqs)
  {
    // our x index decides on which target we operate, while our y index decides on which other vector we operate
    // load two elements at the same time
    Index idx = threadIdx.x + blockDim.x * blockIdx.x * 2;
    const Index grid_size = blockDim.x * 2 * gridDim.x;

    // since we perform cyclic reductions, we can not simply jump out, so we will have to perform empty operations
    // first of all, define our shared array
    // link to our actual, independent slice for our group
    DT_* array_dt = (DT_*)&array[0];
    DT_* loc_dot = &array_dt[0];
    loc_dot[threadIdx.x] = DT_(0);

    // do as many inner loops as possible before synchronizing
    while(idx < dot_size)
    {
      if((idx+blockSize_) < dot_size)
      {
        loc_dot[threadIdx.x] += target[idx] * freqs[idx] * other[idx] + target[idx+blockSize_] * freqs[idx+blockSize_] * other[idx+blockSize_];
      }
      else
      {
        loc_dot[threadIdx.x] += target[idx] * freqs[idx] * other[idx];
      }
      idx += grid_size;
    }

    //TODO: actual do this in second kernel
    reduce_sum<DT_, blockSize_>(&loc_dot[0], threadIdx.x);

    cuda::atomic_ref<DT_, cuda::thread_scope_device> val(*device_dot);

    // finally, we have to add to our global value (and hope that we started this kernel with such few processes that the global lock is not detrimental)
    if(threadIdx.x == 0)
    {
      val.fetch_add(loc_dot[0]);
    }
  }

  template<typename DT_, int blockSize_>
  __global__ void cuda_batched_triple_dot(DT_* __restrict__ device_dots, const DT_ * target, const Index dot_size, DT_ const* const* others, DT_ const* __restrict__ const freqs, const Index num_dots)
  {
    // our x index decides on which target we operate, while our y index decides on which other vector we operate
    // load two elements at the same time
    Index idx = threadIdx.x + blockDim.x * blockIdx.x * 2;
    const Index idy = threadIdx.y + blockDim.y * blockIdx.y;

    const Index grid_size = blockDim.x * 2 * gridDim.x;

    // since we perform cyclic reductions, we can not simply jump out, so we will have to perform empty operations
    // first of all, define our shared array
    // link to our actual, independent slice for our group
    DT_* array_dt = (DT_*)&array[0];
    DT_* loc_dot = &array_dt[blockDim.x*threadIdx.y];
    loc_dot[threadIdx.x] = DT_(0);

    // now calculate dot product
    if(idy < num_dots)
    {
      // do as many inner loops as possible before synchronizing
      while(idx < dot_size)
      {
        if((idx+blockSize_) < dot_size)
        {
          loc_dot[threadIdx.x] += target[idx] * freqs[idx] * others[idy][idx] + target[idx+blockSize_] * freqs[idx+blockSize_] * others[idy][idx+blockSize_];
        }
        else
        {
          loc_dot[threadIdx.x] += target[idx] * freqs[idx] * others[idy][idx];
        }
        idx += grid_size;
      }
    }

    reduce_sum<DT_, blockSize_>(&loc_dot[0], threadIdx.x);
    if(idy >= num_dots)
      return;

      // todo: use second kernel instead
    cuda::atomic_ref<DT_, cuda::thread_scope_device> val(device_dots[idy]);

    // finally, we have to add to our global value (and hope that we started this kernel with such few processes that the global lock is not detrimental)
    if(threadIdx.x == 0)
    {
      val.fetch_add(loc_dot[0]);
    }
  }

  template<typename DT_, int blockSize_>
  __global__ void cuda_cumulative_axpy(DT_* __restrict__ target, const Index vec_size, DT_ const* const* others, DT_ const * __restrict__ device_scales, const Index num_axpy, const DT_ alpha)
  {
    // Obtain the current thread's cooperative group
    cg::thread_block my_group = cg::this_thread_block();

    // Partition the cooperative group into tiles of size 32
    cg::thread_block_tile<32> warp_group = cg::tiled_partition<32>(my_group);
    const Index idx = warp_group.thread_rank() + blockIdx.x * 32;
    Index idy = warp_group.meta_group_rank()*2;

    DT_* shared_array = (DT_*)&array[0];

    // each process sets his shared value to zero
    shared_array[my_group.thread_rank()] = DT_(0);

    if(idx < vec_size)
    {
      // now load (if possible two values) as long as we can
      while(idy < num_axpy)
      {
        if(idy + 1 < num_axpy)
        {
          shared_array[my_group.thread_rank()] += device_scales[idy] * others[idy][idx] + device_scales[idy+1] * others[idy+1][idx];
        }
        else
          shared_array[my_group.thread_rank()] += device_scales[idy] * others[idy][idx];
        idy += warp_group.meta_group_size()*2;
      }
    }

    // now reduce shared_array down to one warp
    if constexpr(blockSize_ >= 1024)
    {
      my_group.sync();
      if(my_group.thread_rank() < 512)
        shared_array[my_group.thread_rank()] += shared_array[my_group.thread_rank()+512];
    }
    if constexpr(blockSize_ >= 512)
    {
      my_group.sync();
      if(my_group.thread_rank() < 256)
        shared_array[my_group.thread_rank()] += shared_array[my_group.thread_rank()+256];
    }
    if constexpr(blockSize_ >= 256)
    {
      my_group.sync();
      if(my_group.thread_rank() < 128)
        shared_array[my_group.thread_rank()] += shared_array[my_group.thread_rank()+128];
    }
    if constexpr(blockSize_ >= 128)
    {
      my_group.sync();
      if(my_group.thread_rank() < 64)
        shared_array[my_group.thread_rank()] += shared_array[my_group.thread_rank()+64];
    }
    if constexpr(blockSize_ >= 64)
    {
      my_group.sync();
      if(my_group.thread_rank() < 32)
        shared_array[my_group.thread_rank()] += shared_array[my_group.thread_rank()+32];
    }
    if((my_group.thread_rank() < 32) && (idx < vec_size))
    {
      target[idx] += alpha * shared_array[my_group.thread_rank()];
    }
  }

}
