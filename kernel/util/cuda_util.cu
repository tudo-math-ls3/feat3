// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/assertion.hpp>
#include "cuda_profiler_api.h"

#ifdef FEAT_HAVE_CUDSS
#include <kernel/solver/cudss.hpp>
#endif

using namespace FEAT;

Index FEAT::Util::cuda_blocksize_misc = 256;
Index FEAT::Util::cuda_blocksize_reduction = 256;
Index FEAT::Util::cuda_blocksize_spmv = 256;
Index FEAT::Util::cuda_blocksize_axpy = 256;
Index FEAT::Util::cuda_blocksize_scalar_assembly = 256;
Index FEAT::Util::cuda_blocksize_blocked_assembly = 128;
Index FEAT::Util::cuda_blocksize_vanka_assembly = 64;

int FEAT::Util::cuda_device_number = 0;

cusparseHandle_t FEAT::Util::Intern::cusparse_handle;
cublasHandle_t FEAT::Util::Intern::cublas_handle;
cublasLtMatmulAlgo_t * FEAT::Util::Intern::cublas_lt_algo_matmat;
bool * FEAT::Util::Intern::cublas_lt_algo_matmat_initialized;
size_t FEAT::Util::Intern::cuda_workspace_size;
void * FEAT::Util::Intern::cuda_workspace;

#ifdef FEAT_HAVE_CUDSS
cudssHandle_t FEAT::Util::Intern::cudss_handle = nullptr;
#endif

namespace FEAT
{
  namespace Util
  {
    namespace Intern
    {
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

void FEAT::Util::cuda_set_device(const int device)
{
  if (cudaSuccess != cudaSetDevice(device))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaSetDevice failed!");
}

void FEAT::Util::cuda_check_last_error()
{
  if(Runtime::SyncGuard::enable_synchronize())
    cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
}

void * FEAT::Util::cuda_get_device_pointer(void * host)
{
  void * device(nullptr);
  if (cudaSuccess != cudaHostGetDevicePointer((void**)&device, host, 0))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaHostGetDevicePointer failed!");
  return device;
}

void * FEAT::Util::cuda_malloc_managed(const Index bytes)
{
  void * memory(nullptr);
  if (bytes == 0)
    return memory;

  auto status = cudaMallocManaged((void**)&memory, bytes);
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_malloc_managed allocation error\n" + stringify(cudaGetErrorString(status)));
  if (memory == nullptr)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_malloc_managed allocation error (null pointer returned)");
  return memory;
}

void * FEAT::Util::cuda_malloc(const Index bytes)
{
  void * memory(nullptr);
  if(bytes == 0)
    return memory;

  auto status = cudaMalloc((void**)&memory, bytes);
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_malloc allocation error\n" + stringify(cudaGetErrorString(status)));
  if (memory == nullptr)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_malloc allocation error (null pointer returned)");
  return memory;
}

void * FEAT::Util::cuda_malloc_host(const Index bytes)
{
  void * memory(nullptr);
  if(bytes == 0)
    return memory;

  auto status = cudaMallocHost((void**)&memory, bytes);
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_malloc allocation error\n" + stringify(cudaGetErrorString(status)));
  if (memory == nullptr)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_malloc allocation error (null pointer returned)");
  return memory;
}

void * FEAT::Util::cuda_get_static_memory(const Index bytes)
{
  if(Intern::cuda_workspace_size < bytes)
  {
    cudaFree(Intern::cuda_workspace);
    auto status = cudaMalloc(&Intern::cuda_workspace, bytes);
    if (status != cudaSuccess)
      throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_get_static_memory allocation error\n" + stringify(cudaGetErrorString(status)));
    if (Intern::cuda_workspace == nullptr)
      throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_get_static_memory allocation error (null pointer returned)");
    Intern::cuda_workspace_size = bytes;
  }
  return Intern::cuda_workspace;
}

void FEAT::Util::cuda_free_static_memory()
{
  if (Intern::cuda_workspace == nullptr)
    return;

  auto status = cudaFree(Intern::cuda_workspace);
  if (cudaSuccess != status)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_free_static_memory: cudaFree failed!\n" + stringify(cudaGetErrorString(status)));
  Intern::cuda_workspace_size = size_t(0u);
  Intern::cuda_workspace = nullptr;
}

void FEAT::Util::cuda_free(void * address)
{
  if (address == nullptr)
    return;

  auto status = cudaFree(address);
  if (cudaSuccess != status)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_free: cudaFree failed!\n" + stringify(cudaGetErrorString(status)));
}

void FEAT::Util::cuda_free_host(void * address)
{
  if (address == nullptr)
    return;

  auto status = cudaFreeHost(address);
  if (cudaSuccess != status)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_free_host: cudaFreeHost failed!\n" + stringify(cudaGetErrorString(status)));
}

void FEAT::Util::cuda_initialize(int rank, int /*ranks_per_node*/, int /*ranks_per_uma*/, int gpus_per_node)
{
  /// \todo enable non cuda ranks and ensure balance of ranks per numa section
  FEAT::Util::cuda_device_number = rank % gpus_per_node;
  if (cudaSuccess != cudaSetDevice(cuda_device_number))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaSetDevice failed!");

  int mm_support = 0;
  if (cudaSuccess != cudaDeviceGetAttribute(&mm_support, cudaDevAttrManagedMemory, cuda_device_number))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaGetAttribute failed!");
  XASSERTM(mm_support == 1, "selected cuda device does not support managed memory!");

  if (CUBLAS_STATUS_SUCCESS != cublasCreate(&Util::Intern::cublas_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cublasCreate failed!");
  if (CUSPARSE_STATUS_SUCCESS != cusparseCreate(&Util::Intern::cusparse_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseCreate failed!");
  if (CUBLAS_STATUS_SUCCESS != cublasSetPointerMode(Util::Intern::cublas_handle, CUBLAS_POINTER_MODE_HOST))
    throw InternalError(__func__, __FILE__, __LINE__, "cublasSetPointerMode failed!");
  if (CUSPARSE_STATUS_SUCCESS != cusparseSetPointerMode(Util::Intern::cusparse_handle, CUSPARSE_POINTER_MODE_HOST))
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseSetPointerMode failed!");

  if (CUBLAS_STATUS_SUCCESS != cublasSetMathMode(Util::Intern::cublas_handle, CUBLAS_TF32_TENSOR_OP_MATH))
    throw InternalError(__func__, __FILE__, __LINE__, "cublasSetMathMode failed!");

#ifdef FEAT_HAVE_CUDSS
  if(CUDSS_STATUS_SUCCESS != cudssCreate(&Util::Intern::cudss_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cudssCreate failed!");
#endif

  Util::Intern::cublas_lt_algo_matmat = new cublasLtMatmulAlgo_t[6];
  Util::Intern::cublas_lt_algo_matmat_initialized = new bool[6];
  for (int i(0) ; i < 6 ; ++i)
  {
    Util::Intern::cublas_lt_algo_matmat_initialized[i] = false;
  }

  //Util::Intern::cuda_workspace_size = 1024ul * 1024ul * 1024ul * 2ul;
  Util::Intern::cuda_workspace_size = 0;
  /*auto status = cudaMalloc(&(Util::Intern::cuda_workspace), Util::Intern::cuda_workspace_size);
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaMalloc failed: " + stringify(cudaGetErrorString(status)));*/
}

void FEAT::Util::cuda_finalize()
{
  if (cudaSuccess != cudaDeviceSynchronize())
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSynchronize failed!");

#ifdef FEAT_HAVE_CUDSS
  if(CUDSS_STATUS_SUCCESS != cudssDestroy(Util::Intern::cudss_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cudssDestroy failed!");
#endif
  if (CUBLAS_STATUS_SUCCESS != cublasDestroy(Util::Intern::cublas_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cublasDestroy failed!");
  if (CUSPARSE_STATUS_SUCCESS != cusparseDestroy(Util::Intern::cusparse_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseDestroy failed!");

  delete[] Util::Intern::cublas_lt_algo_matmat;
  delete[] Util::Intern::cublas_lt_algo_matmat_initialized;

  cuda_free_static_memory();

  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "Pending cuda errors occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
}

void FEAT::Util::cuda_synchronize()
{
  if(Runtime::SyncGuard::enable_synchronize())
  {
    auto status = cudaDeviceSynchronize();
    if (status != cudaSuccess)
      throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSynchronize failed: " + stringify(cudaGetErrorString(status)));
  }
}

void FEAT::Util::cuda_force_synchronize()
{
  auto status = cudaDeviceSynchronize();
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSynchronize failed: " + stringify(cudaGetErrorString(status)));
}

void FEAT::Util::cuda_reset_device()
{
  auto status = cudaDeviceReset();
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceReset failed: " + stringify(cudaGetErrorString(status)));
}

void FEAT::Util::cuda_copy(void * dest, const void * src, const Index bytes)
{
  auto status = cudaMemcpy(dest, src, bytes, cudaMemcpyDefault);
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaMemcpy failed: " + stringify(cudaGetErrorString(status)));
}

void FEAT::Util::cuda_copy_host_to_device(void * dest, const void * src, const Index bytes)
{
  auto status = cudaMemcpy(dest, src, bytes, cudaMemcpyHostToDevice);
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaMemcpy failed: " + stringify(cudaGetErrorString(status)));
}

void FEAT::Util::cuda_copy_device_to_host(void * dest, const void * src, const Index bytes)
{
  auto status = cudaMemcpy(dest, src, bytes, cudaMemcpyDeviceToHost);
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaMemcpy failed: " + stringify(cudaGetErrorString(status)));
}

void FEAT::Util::cuda_copy_device_to_device(void * dest, const void * src, const Index bytes)
{
  auto status = cudaMemcpy(dest, src, bytes, cudaMemcpyDeviceToDevice);
  if (status != cudaSuccess)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaMemcpy failed: " + stringify(cudaGetErrorString(status)));
}

void FEAT::Util::cuda_set_blocksize(Index misc, Index reduction, Index spmv, Index axpy, Index scalar_assembly, Index blocked_assembly)
{
  FEAT::Util::cuda_blocksize_misc = misc;

  FEAT::Util::cuda_blocksize_reduction = reduction;

  FEAT::Util::cuda_blocksize_spmv = spmv;

  FEAT::Util::cuda_blocksize_axpy = axpy;

  FEAT::Util::cuda_blocksize_scalar_assembly = scalar_assembly;

  FEAT::Util::cuda_blocksize_blocked_assembly = blocked_assembly;
}

void FEAT::Util::cuda_reset_algos()
{
  for (int i(0) ; i < 6 ; ++i)
  {
    Util::Intern::cublas_lt_algo_matmat_initialized[i] = false;
  }
}

template <typename DT_>
void FEAT::Util::cuda_set_memory(DT_ * address, const DT_ val, const Index count)
{
  Index blocksize = FEAT::Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((count)/(double)(block.x));
  FEAT::Util::Intern::cuda_set_memory<<<grid, block>>>(address, val, count);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_set_memory failed!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
#ifdef FEAT_HAVE_HALFMATH
template void FEAT::Util::cuda_set_memory(Half * , const Half, const Index);
#endif
template void FEAT::Util::cuda_set_memory(float * , const float, const Index);
template void FEAT::Util::cuda_set_memory(double * , const double, const Index);
template void FEAT::Util::cuda_set_memory(unsigned int * , const unsigned int, const Index);
template void FEAT::Util::cuda_set_memory(unsigned long * , const unsigned long, const Index);
template void FEAT::Util::cuda_set_memory(unsigned long long * , const unsigned long long, const Index);
template void FEAT::Util::cuda_set_memory(int * , const int, const Index);
template void FEAT::Util::cuda_set_memory(long * , const long, const Index);
template void FEAT::Util::cuda_set_memory(long long * , const long long, const Index);

template <typename DT1_, typename DT2_>
void FEAT::Util::cuda_convert(DT1_ * dest, const DT2_ * src, const Index count)
{
  Index blocksize = FEAT::Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((count)/(double)(block.x));
  FEAT::Util::Intern::cuda_convert<<<grid, block>>>(dest, src, count);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_convert failed!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
#ifdef FEAT_HAVE_HALFMATH
template void FEAT::Util::cuda_convert<Half, float>(Half *, const float *, const Index);
template void FEAT::Util::cuda_convert<float, Half>(float *, const Half *, const Index);
template void FEAT::Util::cuda_convert<Half, double>(Half *, const double *, const Index);
template void FEAT::Util::cuda_convert<double, Half>(double *, const Half *, const Index);
#endif
template void FEAT::Util::cuda_convert<float, double>(float *, const double *, const Index);
template void FEAT::Util::cuda_convert<double, float>(double *, const float *, const Index);
template void FEAT::Util::cuda_convert<unsigned int, unsigned long>(unsigned int *, const unsigned long *, const Index);
template void FEAT::Util::cuda_convert<unsigned int, unsigned long long>(unsigned int *, const unsigned long long *, const Index);
template void FEAT::Util::cuda_convert<unsigned long, unsigned int>(unsigned long *, const unsigned int *, const Index);
template void FEAT::Util::cuda_convert<unsigned long, unsigned long long>(unsigned long *, const unsigned long long *, const Index);
template void FEAT::Util::cuda_convert<unsigned long long, unsigned int>(unsigned long long *, const unsigned int *, const Index);
template void FEAT::Util::cuda_convert<unsigned long long, unsigned long>(unsigned long long *, const unsigned long *, const Index);
template void FEAT::Util::cuda_convert<unsigned int, double>(unsigned int *, const double *, const Index);
template void FEAT::Util::cuda_convert<unsigned long, double>(unsigned long *, const double *, const Index);
template void FEAT::Util::cuda_convert<unsigned int, float>(unsigned int *, const float *, const Index);
template void FEAT::Util::cuda_convert<unsigned long, float>(unsigned long *, const float *, const Index);

int FEAT::Util::cuda_get_device_count()
{
  int numDevices(-1);
  if (cudaSuccess != cudaGetDeviceCount(&numDevices))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaGetDeviceCount failed!");
  return numDevices;
}

int FEAT::Util::cuda_get_device_id()
{
  int device(-1);
  if (cudaSuccess != cudaGetDevice(&device))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaGetDevice failed!");
  return device;
}

String FEAT::Util::cuda_get_visible_devices()
{
  String result("");
  int numDevices(-1);
  if (cudaSuccess != cudaGetDeviceCount(&numDevices))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaGetDeviceCount failed!");
  result += "Number of visible cuda devices: " + stringify(numDevices) + "\n" ;

  for (int idevice(0); idevice<numDevices; ++idevice)
  {
    // get device properties
    cudaDeviceProp prop;
    if (cudaSuccess != cudaGetDeviceProperties (&prop, idevice))
      throw InternalError(__func__, __FILE__, __LINE__, "cudaGetDeviceProperties failed!");
    // print out device name and compute capabilities
    result += "Device " + stringify(idevice) + ": " + stringify(prop.name) + "\n";
  }
  return result;
}

std::size_t Util::cuda_get_max_cache_thread()
{
  std::size_t value = 0;
  auto status = cudaDeviceGetLimit(&value, cudaLimit::cudaLimitStackSize);
  if(cudaSuccess != status)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSetLimit failed!");
  return value;
}

void FEAT::Util::cuda_set_max_cache_thread(const std::size_t bytes)
{
  std::size_t value = bytes;
  auto status = cudaDeviceSetLimit(cudaLimit::cudaLimitStackSize, value);
  if(cudaSuccess != status)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSetLimit failed!");
  status = cudaDeviceGetLimit(&value, cudaLimit::cudaLimitStackSize);
  if(cudaSuccess != status)
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceGetLimit failed!");
  if(value != bytes)
   XABORTM("Could not set max cache per thread (" +stringify(value)+") to expected amount(" + stringify(bytes) + ")");
}

void FEAT::Util::cuda_start_profiling()
{
  cudaProfilerStart();
}

void FEAT::Util::cuda_stop_profiling()
{
  cudaProfilerStop();
}

std::size_t FEAT::Util::cuda_get_shared_mem_per_sm()
{
  int max_shared_mem_sm = 0;
  if(cudaDeviceGetAttribute(&max_shared_mem_sm, cudaDevAttrMaxSharedMemoryPerMultiprocessor, Util::cuda_device_number) != cudaSuccess)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceGetAttribute failed!");
  }
  return std::size_t(max_shared_mem_sm);

}

std::size_t FEAT::Util::cuda_get_max_blocks_per_sm()
{
  int max_blocks_per_sm = 0;
  if(cudaDeviceGetAttribute(&max_blocks_per_sm, cudaDevAttrMaxBlocksPerMultiprocessor, Util::cuda_device_number) != cudaSuccess)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceGetAttribute failed!");
  }
  return std::size_t(max_blocks_per_sm);

}

std::size_t FEAT::Util::cuda_get_sm_count()
{
  int max_sm_per_device = 0;
  if(cudaDeviceGetAttribute(&max_sm_per_device, cudaDevAttrMultiProcessorCount, Util::cuda_device_number) != cudaSuccess)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceGetAttribute failed!");
  }
  return std::size_t(max_sm_per_device);

}
