// includes, FEAT
#include <kernel/util/memory_pool.hpp>

#include <cstdio>
#include <fstream>
#include <cassert>

#include <cublas_v2.h>
#include "cusparse_v2.h"

namespace FEAT
{
  namespace Util
  {
    namespace Intern
    {
      cublasHandle_t cublas_handle;
      cusparseHandle_t cusparse_handle;

      template <typename DT_>
      __global__ void cuda_set_memory(DT_ * ptr, const DT_ val, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        ptr[idx] = val;
      }

      template <typename DT1_, typename DT2_>
      __global__ void cuda_convert(DT1_ * dest, const DT2_ * src, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        dest[idx] = src[idx];
      }
    }
  }
}


using namespace FEAT;

// static member initialisation
std::map<void*, Util::Intern::MemoryInfo> MemoryPool<Mem::CUDA>::_pool;
Index MemoryPool<Mem::CUDA>::blocksize_misc = 256;
Index MemoryPool<Mem::CUDA>::blocksize_reduction = 256;
Index MemoryPool<Mem::CUDA>::blocksize_spmv = 256;
Index MemoryPool<Mem::CUDA>::blocksize_axpy = 256;

void MemoryPool<Mem::CUDA>::initialise(int rank, int /*ranks_per_node*/, int /*ranks_per_uma*/, int gpus_per_node)
{
  /// \todo enable non cuda ranks and ensure balance of ranks per numa section
  int device = rank % gpus_per_node;
  if (cudaSuccess != cudaSetDevice(device))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaSetDevice failed!");

  if (cudaSuccess != cudaSetDeviceFlags(cudaDeviceMapHost))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaSetDeviceFlags failed!");

  if (CUBLAS_STATUS_SUCCESS != cublasCreate(&Util::Intern::cublas_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cublasCreate failed!");
  if (CUSPARSE_STATUS_SUCCESS != cusparseCreate(&Util::Intern::cusparse_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseCreate failed!");
  cublasSetPointerMode(Util::Intern::cublas_handle, CUBLAS_POINTER_MODE_HOST);
  cusparseSetPointerMode(Util::Intern::cusparse_handle, CUSPARSE_POINTER_MODE_HOST);
}

void MemoryPool<Mem::CUDA>::finalise()
{
  if (_pool.size() > 0)
  {
    std::cout << stderr << " Error: MemoryPool<CUDA> still contains memory chunks on deconstructor call" << std::endl;
    std::exit(1);
  }

  if (CUBLAS_STATUS_SUCCESS != cublasDestroy(Util::Intern::cublas_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cublasDestroy failed!");
  if (CUSPARSE_STATUS_SUCCESS != cusparseDestroy(Util::Intern::cusparse_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseDestroy failed!");

  if (cudaSuccess != cudaDeviceSynchronize())
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSynchronize failed!");
  if (cudaSuccess != cudaDeviceReset())
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSynchronize failed!");

  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "Pending cuda errors occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
}

template <typename DT_>
DT_ * MemoryPool<Mem::CUDA>::allocate_memory(const Index count)
{
  assert(count != 0);

  DT_ * memory(nullptr);
  if (cudaErrorMemoryAllocation == cudaMalloc((void**)&memory, count * sizeof(DT_)))
    throw InternalError("MemoryPool<CUDA> cuda allocation error (cudaErrorMemoryAllocation)");
  if (memory == nullptr)
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA> allocation error (null pointer returned)");
  Util::Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = count * sizeof(DT_);
  _pool.insert(std::pair<void*, Util::Intern::MemoryInfo>(memory, mi));
  return memory;
}

void MemoryPool<Mem::CUDA>::increase_memory(void * address)
{
  if (address == nullptr)
    return;

  std::map<void*, Util::Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::increase_memory: Memory address not found!");
  else
  {
    it->second.counter = it->second.counter + 1;
  }
}

void MemoryPool<Mem::CUDA>::release_memory(void * address)
{
  if (address == nullptr)
    return;

  std::map<void*, Util::Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::relase_memory: Memory address not found!");
  else
  {
    if(it->second.counter == 1)
    {
      if (cudaSuccess != cudaFree(address))
        throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::release_memory: cudaFree failed!");
      _pool.erase(it);
    }
    else
    {
      it->second.counter = it->second.counter - 1;
    }
  }
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::download(DT_ * dest, const DT_ * const src, const Index count)
{
  if (cudaSuccess != cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyDeviceToHost))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::download failed!");
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::upload(DT_ * dest, const DT_ * const src, const Index count)
{
  if (cudaSuccess != cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyHostToDevice))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::upload failed!");
}

template <typename DT_>
DT_ MemoryPool<Mem::CUDA>::get_element(const DT_ * data, const Index index)
{
  const void * src(data + index);
  DT_ value;
  if (cudaSuccess != cudaMemcpy(&value, src, sizeof(DT_), cudaMemcpyDeviceToHost))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::get_element failed!");
  return value;
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::set_memory(DT_ * address, const DT_ val, const Index count)
{
  Index blocksize = blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((count)/(double)(block.x));
  FEAT::Util::Intern::cuda_set_memory<<<grid, block>>>(address, val, count);
#ifdef FEAT_DEBUG
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::set_memory failed!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::copy(DT_ * dest, const DT_ * src, const Index count)
{
  if (dest == src)
    return;

  if (cudaSuccess != cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyDeviceToDevice))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::copy failed!");
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::convert(DT_ * dest, const DT_ * src, const Index count)
{
  if (dest == src)
    return;

  if (cudaSuccess != cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyDeviceToDevice))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::copy failed!");
}

template <typename DT1_, typename DT2_>
void MemoryPool<Mem::CUDA>::convert(DT1_ * dest, const DT2_ * src, const Index count)
{
  Index blocksize = blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((count)/(double)(block.x));
  FEAT::Util::Intern::cuda_convert<<<grid, block>>>(dest, src, count);
#ifdef FEAT_DEBUG
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::convert failed!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

void MemoryPool<Mem::CUDA>::synchronise()
{
  cudaDeviceSynchronize();
}

void MemoryPool<Mem::CUDA>::reset_device()
{
  cudaDeviceReset();
}

void MemoryPool<Mem::CUDA>::set_blocksize(Index misc, Index reduction, Index spmv, Index axpy)
{
  blocksize_misc = misc;

  blocksize_reduction = reduction;

  blocksize_spmv = spmv;

  blocksize_axpy = axpy;
}

Index MemoryPool<Mem::CUDA>::allocated_memory()
{
  Index bytes(0);
  for (auto& i : _pool)
  {
    bytes += i.second.size;
  }
  return bytes;
}

template float * MemoryPool<Mem::CUDA>::allocate_memory<float>(const Index);
template double * MemoryPool<Mem::CUDA>::allocate_memory<double>(const Index);
template unsigned int * MemoryPool<Mem::CUDA>::allocate_memory<unsigned int>(const Index);
template unsigned long * MemoryPool<Mem::CUDA>::allocate_memory<unsigned long>(const Index);
template unsigned long long * MemoryPool<Mem::CUDA>::allocate_memory<unsigned long long>(const Index);

template void MemoryPool<Mem::CUDA>::download<float>(float *, const float * const, const Index);
template void MemoryPool<Mem::CUDA>::download<double>(double *, const double * const, const Index);
template void MemoryPool<Mem::CUDA>::download<unsigned int>(unsigned int *, const unsigned int * const, const Index);
template void MemoryPool<Mem::CUDA>::download<unsigned long>(unsigned long *, const unsigned long * const, const Index);
template void MemoryPool<Mem::CUDA>::download<unsigned long long>(unsigned long long *, const unsigned long long * const, const Index);

template void MemoryPool<Mem::CUDA>::upload<float>(float *, const float * const, const Index);
template void MemoryPool<Mem::CUDA>::upload<double>(double *, const double * const , const Index);
template void MemoryPool<Mem::CUDA>::upload<unsigned int>(unsigned int *, const unsigned int * const, const Index);
template void MemoryPool<Mem::CUDA>::upload<unsigned long>(unsigned long *, const unsigned long * const, const Index);
template void MemoryPool<Mem::CUDA>::upload<unsigned long long>(unsigned long long *, const unsigned long long * const, const Index);

template float MemoryPool<Mem::CUDA>::get_element(const float *, const Index);
template double MemoryPool<Mem::CUDA>::get_element(const double *, const Index);
template unsigned int MemoryPool<Mem::CUDA>::get_element(const unsigned int *, const Index);
template unsigned long MemoryPool<Mem::CUDA>::get_element(const unsigned long *, const Index);
template unsigned long long MemoryPool<Mem::CUDA>::get_element(const unsigned long long *, const Index);

template void MemoryPool<Mem::CUDA>::set_memory(float * , const float, const Index);
template void MemoryPool<Mem::CUDA>::set_memory(double * , const double, const Index);
template void MemoryPool<Mem::CUDA>::set_memory(unsigned int * , const unsigned int, const Index);
template void MemoryPool<Mem::CUDA>::set_memory(unsigned long * , const unsigned long, const Index);
template void MemoryPool<Mem::CUDA>::set_memory(unsigned long long * , const unsigned long long, const Index);

template void MemoryPool<Mem::CUDA>::copy<float>(float *, const float *, const Index);
template void MemoryPool<Mem::CUDA>::copy<double>(double *, const double *, const Index);
template void MemoryPool<Mem::CUDA>::copy<unsigned int>(unsigned int *, const unsigned int *, const Index);
template void MemoryPool<Mem::CUDA>::copy<unsigned long>(unsigned long *, const unsigned long *, const Index);
template void MemoryPool<Mem::CUDA>::copy<unsigned long long>(unsigned long long *, const unsigned long long *, const Index);

template void MemoryPool<Mem::CUDA>::convert<float>(float *, const float *, const Index);
template void MemoryPool<Mem::CUDA>::convert<double>(double *, const double *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned int>(unsigned int *, const unsigned int *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned long>(unsigned long *, const unsigned long *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned long long>(unsigned long long *, const unsigned long long *, const Index);

template void MemoryPool<Mem::CUDA>::convert<float, double>(float *, const double *, const Index);
template void MemoryPool<Mem::CUDA>::convert<double, float>(double *, const float *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned int, unsigned long>(unsigned int *, const unsigned long *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned int, unsigned long long>(unsigned int *, const unsigned long long *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned long, unsigned int>(unsigned long *, const unsigned int *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned long, unsigned long long>(unsigned long *, const unsigned long long *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned long long, unsigned int>(unsigned long long *, const unsigned int *, const Index);
template void MemoryPool<Mem::CUDA>::convert<unsigned long long, unsigned long>(unsigned long long *, const unsigned long *, const Index);
