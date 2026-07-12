// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/memory_aux.hpp>
#include <kernel/util/cuda_helper.cuh>
#include <kernel/util/cuda_util.hpp>

namespace FEAT
{
  namespace Memory
  {
    void* alloc_main(std::size_t bytes)
    {
      void* ptr = nullptr;
      if(cudaSuccess != ::cudaMallocHost(&ptr, bytes))
        return nullptr;
      return ptr;
    }

    void free_main(void* ptr)
    {
      if(ptr == nullptr)
        return;
      cudaCheck(::cudaFreeHost(ptr), "Memory::free_main: cudaFreeHost");
    }

    void* alloc_cuda(std::size_t bytes)
    {
      void* ptr = nullptr;
      if(bytes == 0u)
        return ptr;
      cudaCheck(::cudaMalloc(&ptr, bytes), "Memory::alloc_cuda: cudaMalloc");
      return ptr;
    }

    void* alloc_cuda_async(std::size_t bytes, [[maybe_unused]] const Util::StreamWrapper& stream)
    {
      void* ptr = nullptr;
      cudaError_t st;
    #ifdef CUDA_ASYNC_SUPPORTED
      const cudaStream_t s_stream = reinterpret_cast<cudaStream_t>(stream.stream);
      st = Util::cuda_get_stream_ordered_alloc_supported() ? ::cudaMallocAsync(&ptr, bytes, s_stream) : ::cudaMalloc(&ptr, bytes);
    #else
      st = ::cudaMalloc(&ptr, bytes);
    #endif
      cudaCheck(st, "Memory::alloc_cuda_async");
      return ptr;
    }

    void free_cuda(void* ptr)
    {
      if(ptr == nullptr)
        return;
      cudaCheck(::cudaFree(ptr), "Memory::free_cuda: cudaFree");
    }

    void free_cuda_async(void* ptr, [[maybe_unused]] const Util::StreamWrapper& stream)
    {
      if(ptr == nullptr)
        return;
      cudaError_t st;
    #ifdef CUDA_ASYNC_SUPPORTED
      const cudaStream_t s_stream = reinterpret_cast<cudaStream_t>(stream.stream);
      st = Util::cuda_get_stream_ordered_alloc_supported() ? cudaFreeAsync(ptr, s_stream) : cudaFree(ptr);
    #else
      st = ::cudaFree(ptr);
    #endif
      cudaCheck(st, "Memory::free_cuda_async");
    }

    void memset_cuda(void* dst, int value, std::size_t bytes, const Util::StreamWrapper& stream)
    {
      if(bytes <= std::size_t(0))
        return;

      XASSERTM(dst != nullptr, "cuda device memset destination is nullptr");

      if constexpr(async_operations)
      {
        const cudaStream_t s_stream = reinterpret_cast<cudaStream_t>(stream.stream);
        cudaCheck(cudaMemsetAsync(dst, value, bytes, s_stream), "Memory::memset_cuda async");
      }
      else
      {
        cudaCheck(cudaMemset(dst, value, bytes), "Memory::memset_cuda");
      }
    }

    void memcopy_cuda(void* dst, const void* src, std::size_t bytes, const Util::StreamWrapper& stream)
    {
      if(bytes <= std::size_t(0))
        return;

      XASSERTM(dst != nullptr, "cuda device-to-device memcpy destination is nullptr");
      XASSERTM(src != nullptr, "cuda device-to-device memcpy source is nullptr");

      if constexpr(async_operations)
      {
        const cudaStream_t s_stream = reinterpret_cast<cudaStream_t>(stream.stream);
        cudaCheck(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToDevice, s_stream), "Memory::memcopy_cuda async");
      }
      else
      {
        cudaCheck(cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToDevice), "Memory::memcopy_cuda");
      }
    }

    void memcopy_cuda_to_main(void* dst, const void* src, std::size_t bytes)
    {
      if(bytes <= std::size_t(0))
        return;

      XASSERTM(dst != nullptr, "cuda device-to-host memcpy destination is nullptr");
      XASSERTM(src != nullptr, "cuda device-to-host memcpy source is nullptr");

      cudaCheck(cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost), "Memory::memcopy_cuda_to_main");
    }

    void memcopy_cuda_to_main_async(void* dst, const void* src, std::size_t bytes, const Util::StreamWrapper& stream)
    {
      if(bytes <= std::size_t(0))
        return;

      XASSERTM(dst != nullptr, "cuda device-to-host memcpy destination is nullptr");
      XASSERTM(src != nullptr, "cuda device-to-host memcpy source is nullptr");

      if constexpr(async_operations)
      {
        const cudaStream_t s_stream = reinterpret_cast<cudaStream_t>(stream.stream);
        cudaCheck(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, s_stream), "Memory::memcopy_cuda_to_main_async");
      }
      else
        cudaCheck(cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost), "Memory::memcopy_cuda_to_main_async");
    }

    void memcopy_main_to_cuda(void* dst, const void* src, std::size_t bytes, const Util::StreamWrapper& stream)
    {
      if(bytes <= std::size_t(0))
        return;

      XASSERTM(dst != nullptr, "cuda host-to-device memcpy destination is nullptr");
      XASSERTM(src != nullptr, "cuda host-to-device memcpy source is nullptr");

      // should be okay to do on main stream, since we can not delete source memory without sync call
      if constexpr(async_operations)
      {
        const cudaStream_t s_stream = reinterpret_cast<cudaStream_t>(stream.stream);
        cudaCheck(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, s_stream), "Memory::memcopy_main_to_cuda async");
      }
      else
        cudaCheck(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice), "Memory::memcopy_main_to_cuda");
    }

    template <typename Tdst_, typename Tsrc_>
    __global__ void convert_cuda_impl(Tdst_* dst, const Tsrc_* src, const Index count)
    {
      // grid strided for loop
      for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < count; i += blockDim.x * gridDim.x)
        dst[i] = Tdst_(src[i]);
    }

    template<typename Tdst_, typename Tsrc_>
    void convert_cuda(Tdst_* dst, const Tsrc_* src, const Index count)
    {
      if(count <= Index(0))
        return;

      dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
      dim3 grid((static_cast<unsigned int>(count) + block.x - 1u) / block.x); // round up to next integer
      convert_cuda_impl<<<grid, block>>>(dst, src, count);

#ifdef DEBUG
      cudaError_t last_error(cudaGetLastError());
      if (cudaSuccess != last_error)
      {
        std::string strerr("Memory::convert_cuda failed!\n");
        strerr.append(cudaGetErrorString(last_error));
        XABORTM(strerr.c_str());
      }
#endif
    }

#ifdef FEAT_HAVE_HALFMATH
    template void convert_cuda<Half, float>(Half *, const float *, const Index);
    template void convert_cuda<Half, double>(Half *, const double *, const Index);
    template void convert_cuda<float, Half>(float *, const Half *, const Index);
    template void convert_cuda<double, Half>(double *, const Half *, const Index);
#endif

    // Python code for generating the mumbo jumbo below:
    // tt = [
    //   "float",
    //   "double",
    //   "signed int",
    //   "signed long",
    //   "signed long long",
    //   "unsigned int",
    //   "unsigned long",
    //   "unsigned long long"
    // ]
    //
    // for t in tt:
    //   for s in tt:
    //     if t != s:
    //       print("template void convert_cuda<" + t + ", " + s + ">(" + t + "*, const " + s + "*, const Index);")

    template void convert_cuda<float, double>(float*, const double*, const Index);
    template void convert_cuda<float, signed int>(float*, const signed int*, const Index);
    template void convert_cuda<float, signed long>(float*, const signed long*, const Index);
    template void convert_cuda<float, signed long long>(float*, const signed long long*, const Index);
    template void convert_cuda<float, unsigned int>(float*, const unsigned int*, const Index);
    template void convert_cuda<float, unsigned long>(float*, const unsigned long*, const Index);
    template void convert_cuda<float, unsigned long long>(float*, const unsigned long long*, const Index);
    template void convert_cuda<double, float>(double*, const float*, const Index);
    template void convert_cuda<double, signed int>(double*, const signed int*, const Index);
    template void convert_cuda<double, signed long>(double*, const signed long*, const Index);
    template void convert_cuda<double, signed long long>(double*, const signed long long*, const Index);
    template void convert_cuda<double, unsigned int>(double*, const unsigned int*, const Index);
    template void convert_cuda<double, unsigned long>(double*, const unsigned long*, const Index);
    template void convert_cuda<double, unsigned long long>(double*, const unsigned long long*, const Index);
    template void convert_cuda<signed int, float>(signed int*, const float*, const Index);
    template void convert_cuda<signed int, double>(signed int*, const double*, const Index);
    template void convert_cuda<signed int, signed long>(signed int*, const signed long*, const Index);
    template void convert_cuda<signed int, signed long long>(signed int*, const signed long long*, const Index);
    template void convert_cuda<signed int, unsigned int>(signed int*, const unsigned int*, const Index);
    template void convert_cuda<signed int, unsigned long>(signed int*, const unsigned long*, const Index);
    template void convert_cuda<signed int, unsigned long long>(signed int*, const unsigned long long*, const Index);
    template void convert_cuda<signed long, float>(signed long*, const float*, const Index);
    template void convert_cuda<signed long, double>(signed long*, const double*, const Index);
    template void convert_cuda<signed long, signed int>(signed long*, const signed int*, const Index);
    template void convert_cuda<signed long, signed long long>(signed long*, const signed long long*, const Index);
    template void convert_cuda<signed long, unsigned int>(signed long*, const unsigned int*, const Index);
    template void convert_cuda<signed long, unsigned long>(signed long*, const unsigned long*, const Index);
    template void convert_cuda<signed long, unsigned long long>(signed long*, const unsigned long long*, const Index);
    template void convert_cuda<signed long long, float>(signed long long*, const float*, const Index);
    template void convert_cuda<signed long long, double>(signed long long*, const double*, const Index);
    template void convert_cuda<signed long long, signed int>(signed long long*, const signed int*, const Index);
    template void convert_cuda<signed long long, signed long>(signed long long*, const signed long*, const Index);
    template void convert_cuda<signed long long, unsigned int>(signed long long*, const unsigned int*, const Index);
    template void convert_cuda<signed long long, unsigned long>(signed long long*, const unsigned long*, const Index);
    template void convert_cuda<signed long long, unsigned long long>(signed long long*, const unsigned long long*, const Index);
    template void convert_cuda<unsigned int, float>(unsigned int*, const float*, const Index);
    template void convert_cuda<unsigned int, double>(unsigned int*, const double*, const Index);
    template void convert_cuda<unsigned int, signed int>(unsigned int*, const signed int*, const Index);
    template void convert_cuda<unsigned int, signed long>(unsigned int*, const signed long*, const Index);
    template void convert_cuda<unsigned int, signed long long>(unsigned int*, const signed long long*, const Index);
    template void convert_cuda<unsigned int, unsigned long>(unsigned int*, const unsigned long*, const Index);
    template void convert_cuda<unsigned int, unsigned long long>(unsigned int*, const unsigned long long*, const Index);
    template void convert_cuda<unsigned long, float>(unsigned long*, const float*, const Index);
    template void convert_cuda<unsigned long, double>(unsigned long*, const double*, const Index);
    template void convert_cuda<unsigned long, signed int>(unsigned long*, const signed int*, const Index);
    template void convert_cuda<unsigned long, signed long>(unsigned long*, const signed long*, const Index);
    template void convert_cuda<unsigned long, signed long long>(unsigned long*, const signed long long*, const Index);
    template void convert_cuda<unsigned long, unsigned int>(unsigned long*, const unsigned int*, const Index);
    template void convert_cuda<unsigned long, unsigned long long>(unsigned long*, const unsigned long long*, const Index);
    template void convert_cuda<unsigned long long, float>(unsigned long long*, const float*, const Index);
    template void convert_cuda<unsigned long long, double>(unsigned long long*, const double*, const Index);
    template void convert_cuda<unsigned long long, signed int>(unsigned long long*, const signed int*, const Index);
    template void convert_cuda<unsigned long long, signed long>(unsigned long long*, const signed long*, const Index);
    template void convert_cuda<unsigned long long, signed long long>(unsigned long long*, const signed long long*, const Index);
    template void convert_cuda<unsigned long long, unsigned int>(unsigned long long*, const unsigned int*, const Index);
    template void convert_cuda<unsigned long long, unsigned long>(unsigned long long*, const unsigned long*, const Index);

    template <typename T_>
    __global__ void format_cuda_impl(T_* dst, T_ value, const Index count)
    {
      // grid strided for loop
      for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < count; i += blockDim.x * gridDim.x)
        dst[i] = value;
    }

    template<typename T_>
    void format_cuda(T_* dst, const T_& value, const Index count)
    {
      if(count <= Index(0))
        return;

      dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
      dim3 grid((static_cast<unsigned int>(count) + block.x - 1u) / block.x); // round up to next integer
      format_cuda_impl<<<grid, block>>>(dst, value, count);

#ifdef DEBUG
      cudaDeviceSynchronize();
      cudaError_t last_error(cudaGetLastError());
      if (cudaSuccess != last_error)
      {
        std::string strerr("Memory::format_cuda failed!\n");
        strerr.append(cudaGetErrorString(last_error));
        XABORTM(strerr.c_str());
      }
#endif
    }

#ifdef FEAT_HAVE_HALFMATH
    template void format_cuda<Half>(Half *, const Half&, const Index);
#endif
    template void format_cuda<float>(float *, const float&, const Index);
    template void format_cuda<double>(double *, const double&, const Index);
    template void format_cuda<signed int>(signed int *, const signed int&, const Index);
    template void format_cuda<signed long>(signed long *, const signed long&, const Index);
    template void format_cuda<signed long long>(signed long long*, const signed long long&, const Index);
    template void format_cuda<unsigned int>(unsigned int *, const unsigned int&, const Index);
    template void format_cuda<unsigned long>(unsigned long *, const unsigned long&, const Index);
    template void format_cuda<unsigned long long>(unsigned long long*, const unsigned long long&, const Index);
  } // namespace Memory
} // namespace FEAT
