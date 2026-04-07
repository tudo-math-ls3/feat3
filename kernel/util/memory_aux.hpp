// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <memory>

namespace FEAT
{
  /**
   * \brief Memory namespace containing functions and classes related to memory management
   */
  namespace Memory
  {
    /**
     * \brief Allocate a memory block in main heap memory
     *
     * This is effectively a wrapper around ::malloc.
     *
     * \param[in] bytes
     * The size of the memory block in bytes.
     *
     * \returns
     * A pointer to the newly allocated main memory block or \c nullptr, if the memory allocation failed.
     */
    void* alloc_main(std::size_t bytes);

    /**
     * \brief Frees a memory block in main heap memory allocated by alloc_main
     *
     * This is effectively a wrapper around ::free.
     *
     * \param[in] ptr
     * A pointer to the main memory block that is to be freed; must not be \c nullptr.
     */
    void free_main(void* ptr);

    /**
     * \brief Formats all bytes of a memory block in main memory to a given value
     *
     * This is effectively a wrapper around ::memset.
     *
     * \param[inout] ptr
     * A pointer to the main memory block that is to be formatted; must not be \c nullptr if bytes > 0.
     *
     * \param[in] value
     * The value that all bytes are to be set to. Set this to 0 to format all bits to 0 and set
     * it to 0xFF (= 255) to format all bits to 1.
     *
     * \param[in] bytes
     * The size of the memory block in bytes.
     */
    void memset_main(void* ptr, int value, std::size_t bytes);

    /**
     * \brief Formats all bytes of a memory block in main memory to random values
     *
     * \param[inout] ptr
     * A pointer to the main memory block that is to be formatted; must not be \c nullptr if bytes > 0.
     *
     * \param[in] bytes
     * The size of the memory block in bytes.
     *
     * \param[in] seed
     * The seed for the RNG. If set to zero, then the address stored in ptr is used as the seed.
     */
    void memset_random_main(void* ptr, std::size_t bytes, std::uint64_t seed = 0u);

    /**
     * \brief Copies a region of a main memory block into another main memory block.
     *
     * \attention
     * The memory regions [dst, ..., dst + bytes] and [src, ..., src + bytes] must not overlap,
     * however, the trivial case of dst = src is allowed and it will not perform any action.
     *
     * This is effectively a wrapper around ::memcpy.
     *
     * \param[inout] dst
     * A pointer to the main memory block that receives the copy; must not be \c nullptr if bytes > 0.
     *
     * \param[in] src
     * A pointer to the main memory block that acts as the source for the copy; must not be \c nullptr if bytes > 0.
     *
     * \param[in] bytes
     * The number of bytes to be copied from *src to *dst
     */
    void memcopy_main(void* dst, const void* src, std::size_t bytes);

    /// overload for Tdst_ = Tsrc_; performs a memcopy internally
    template<typename T_>
    inline void convert_main(T_* dst, const T_* src, const Index count)
    {
      memcopy_main(dst, src, std::size_t(count) * sizeof(T_));
    }

    /**
     * \brief Converts the elements of an array into element of another array of a different data type
     *
     * \attention
     * The memory regions [dst, ..., dst + count*sizeof(Tdst_)] and [src, ..., src + count*sizeof(Tsrc_)]
     * must not overlap.
     *
     * \param[inout] dst
     * A pointer to the main memory block that receives the converted elements; must not be \c nullptr if count > 0.
     *
     * \param[in] src
     * A pointer to the main memory block that contains the source array; must not be \c nullptr if count > 0.
     *
     * \param[in] count
     * The number of arrays to be converted from *src to *dst
     */
    template<typename Tdst_, typename Tsrc_>
    inline void convert_main(Tdst_* dst, const Tsrc_* src, const Index count)
    {
      if(count <= Index(0))
        return;

      XASSERT(dst != nullptr);
      XASSERT(src != nullptr);

      // check for overlapping memory regions
      {
        const char* vdst = reinterpret_cast<const char*>(dst);
        const char* vsrc = reinterpret_cast<const char*>(src);
        XASSERTM((vdst < vsrc) || (vsrc + std::size_t(count)*sizeof(Tsrc_) <= vdst), "source and destination memory regions overlap");
        XASSERTM((vsrc < vdst) || (vdst + std::size_t(count)*sizeof(Tdst_) <= vsrc), "source and destination memory regions overlap");
      }

      for(Index i = 0u; i < count; ++i)
        dst[i] = Tdst_(src[i]);
    }

    /**
     * \brief Sets all entries of a main memory array to a given value
     *
     * \param[out] dst
     * A \transient pointer to the main memory array whose elements are to be formatted
     *
     * \param[in] value
     * A \transient reference to the value that is to be copied to the array elements
     *
     * \param[in] count
     * The number of elements in the array \p dst
     */
    template<typename T_>
    inline void format_main(T_* dst, const T_& value, const Index count)
    {
      if(count <= Index(0))
        return;

      XASSERT(dst != nullptr);

      for(Index i = 0u; i < count; ++i)
        dst[i] = value;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(FEAT_HAVE_CUDA) || defined(DOXYGEN)
    /**
     * \brief Allocate a memory block in cuda heap memory
     *
     * This is effectively a wrapper around ::cudaMalloc.
     *
     * \param[in] bytes
     * The size of the memory block in bytes.
     *
     * \returns
     * A pointer to the newly allocated cuda memory block or \c nullptr, if the memory allocation failed.
     */
    void* alloc_cuda(std::size_t bytes);

    /**
     * \brief Frees a memory block in cuda heap memory allocated by alloc_cuda
     *
     * This is effectively a wrapper around ::cudaFree.
     *
     * \param[in] ptr
     * A pointer to the cuda memory block that is to be freed; must not be \c nullptr.
     */
    void free_cuda(void* ptr);

    /**
     * \brief Formats all bytes of a memory block in cuda memory to a given value
     *
     * This is effectively a wrapper around ::cudaMemset.
     *
     * \param[inout] ptr
     * A device pointer to the cuda memory block that is to be formatted; must not be \c nullptr if bytes > 0.
     *
     * \param[in] value
     * The value that all bytes are to be set to. Set this to 0 to format all bits to 0 and set
     * it to 0xFF (= 255) to format all bits to 1.
     *
     * \param[in] bytes
     * The size of the memory block in bytes.
     */
    void memset_cuda(void* ptr, int value, std::size_t bytes);

    /**
     * \brief Copies a region of a cuda memory block into another cuda memory block.
     *
     * \attention
     * The memory regions [dst, ..., dst + bytes] and [src, ..., src + bytes] must not overlap,
     * however, the trivial case of dst = src is allowed and it will not perform any action.
     *
     * This is effectively a wrapper around ::cudaMemcpy for cudaMemcpyDeviceToDevice.
     *
     * \param[inout] dst
     * A device pointer to the cuda memory block that receives the copy; must not be \c nullptr if bytes > 0.
     *
     * \param[in] src
     * A device pointer to the cuda memory block that acts as the source for the copy; must not be \c nullptr if bytes > 0.
     *
     * \param[in] bytes
     * The number of bytes to be copied from *src to *dst
     */
    void memcopy_cuda(void* dst, const void* src, std::size_t bytes);


    /**
     * \brief Copies a region of a cuda memory block into a main memory block.
     *
     * This is effectively a wrapper around ::cudaMemcpy for cudaMemcpyDeviceToHost.
     *
     * \param[inout] dst
     * A host pointer to the main memory block that receives the copy; must not be \c nullptr if bytes > 0.
     *
     * \param[in] src
     * A device pointer to the cuda memory block that acts as the source for the copy; must not be \c nullptr if bytes > 0.
     *
     * \param[in] bytes
     * The number of bytes to be copied from *src to *dst
     */
    void memcopy_cuda_to_main(void* dst, const void* src, std::size_t bytes);

    /**
     * \brief Copies a region of a main memory block into a cuda memory block.
     *
     * This is effectively a wrapper around ::cudaMemcpy for cudaMemcpyHostToDevice.
     *
     * \param[inout] dst
     * A device pointer to the cuda memory block that receives the copy; must not be \c nullptr if bytes > 0.
     *
     * \param[in] src
     * A host pointer to the main memory block that acts as the source for the copy; must not be \c nullptr if bytes > 0.
     *
     * \param[in] bytes
     * The number of bytes to be copied from *src to *dst
     */
    void memcopy_main_to_cuda(void* dst, const void* src, std::size_t bytes);

    /**
     * \brief Converts the elements of an array into element of another array of a different data type
     *
     * \attention
     * The memory regions [dst, ..., dst + count*sizeof(Tdst_)] and [src, ..., src + count*sizeof(Tsrc_)]
     * must not overlap.
     *
     * \note
     * This template is explicitly instantiated for the cross-product of the following types:
     * float, double, signed int, signed long, signed long long, unsigned int, unsigned long, unsigned long long
     *
     * \tparam Tdst_
     * Specifies the underlying type of the destination array.
     *
     * \tparam Tsrc_
     * Specifies the underlying type of the source array; must be convertible to Tdst_.
     *
     * \param[inout] dst
     * A device pointer to the cuda memory block that receives the converted elements; must not be \c nullptr if count > 0.
     *
     * \param[in] src
     * A device pointer to the cuda memory block that contains the source array; must not be \c nullptr if count > 0.
     *
     * \param[in] count
     * The number of arrays to be converted from *src to *dst
     */
    template<typename Tdst_, typename Tsrc_>
    void convert_cuda(Tdst_* dst, const Tsrc_* src, const Index count);

    /// overload for Tdst_ = Tsrc_; performs a memcopy internally
    template<typename T_>
    inline void convert_cuda(T_* dst, const T_* src, const Index count)
    {
      memcopy_cuda(dst, src, std::size_t(count) * sizeof(T_));
    }

    /**
     * \brief Sets all entries of a cuda memory array to a given value
     *
     * \param[out] dst
     * A \transient pointer to the cuda memory array whose elements are to be formatted
     *
     * \param[in] value
     * A \transient reference to the value that is to be copied to the array elements
     *
     * \param[in] count
     * The number of elements in the array \p dst
     */
    template<typename T_>
    void format_cuda(T_* dst, const T_& value, const Index count);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#else // no CUDA
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline void* alloc_cuda(std::size_t)
    {
      XABORTM("Memory::alloc_cuda: CUDA is not available");
    }

    inline void free_cuda(void*)
    {
      XABORTM("Memory::free_cuda: CUDA is not available");
    }

    inline void memset_cuda(void*, int, std::size_t)
    {
      XABORTM("Memory::memset_cuda: CUDA is not available");
    }

    inline void memcopy_cuda(void*, const void*, std::size_t)
    {
      XABORTM("Memory::memcopy_cuda: CUDA is not available");
    }

    inline void memcopy_cuda_to_main(void*, const void*, std::size_t)
    {
      XABORTM("Memory::memcopy_cuda_to_main: CUDA is not available");
    }

    inline void memcopy_main_to_cuda(void*, const void*, std::size_t)
    {
      XABORTM("Memory::memcopy_main_to_cuda: CUDA is not available");
    }

    template<typename Tdst_, typename Tsrc_>
    inline void convert_cuda(Tdst_*, const Tsrc_*, Index)
    {
      XABORTM("Memory::convert_cuda: CUDA is not available");
    }

    template<typename T_>
    void format_cuda(T_*, const T_&, const Index)
    {
      XABORTM("Memory::format_cuda: CUDA is not available");
    }
#endif // defined(FEAT_HAVE_CUDA) || defined(DOXYGEN)
  } // namespace Memory
} // namespace FEAT
