// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_CUDA_UTIL_HPP
#define KERNEL_UTIL_CUDA_UTIL_HPP 1

#ifdef FEAT_HAVE_CUDA

// includes, FEAT
#include <kernel/base_header.hpp>

#ifdef __CUDACC__
#include <cusparse_v2.h>
#include <cublas_v2.h>
#include <cublasLt.h>
#endif

namespace FEAT
{
  namespace Util
  {
    /// \cond internal
    namespace Intern
    {
#ifdef __CUDACC__
      extern cusparseHandle_t cusparse_handle;
      extern cublasHandle_t cublas_handle;
      extern cublasLtMatmulAlgo_t * cublas_lt_algo_matmat;
      extern bool * cublas_lt_algo_matmat_initialized;
      extern size_t cuda_workspace_size;
      extern void * cuda_workspace;
#endif
    }
    /// \endcond

    /// cuda threading grid blocksize for miscellaneous ops
    extern Index cuda_blocksize_misc;

    /// cuda threading grid blocksize for reduction type ops
    extern Index cuda_blocksize_reduction;

    /// cuda threading grid blocksize for blas-2 type ops
    extern Index cuda_blocksize_spmv;

    /// cuda threading grid blocksize for blas-1 type ops
    extern Index cuda_blocksize_axpy;

    void cuda_set_blocksize(Index misc, Index reduction, Index spmv, Index axpy);
    void cuda_set_device(const int device);
    void cuda_check_last_error();
    void * cuda_get_device_pointer(void * host);
    void * cuda_malloc_managed(const Index bytes);
    void cuda_free(void * address);
    void cuda_initialize(int rank, int ranks_per_node, int ranks_per_uma, int gpus_per_node);
    void cuda_finalize();
    NOINLINE void cuda_synchronize();
    void cuda_reset_device();
    void cuda_copy(void * dest, const void * src, const Index bytes);
    void cuda_reset_algos();
    template <typename DT_>
    void cuda_set_memory(DT_ * address, const DT_ val, const Index count);
    template <typename DT1_, typename DT2_>
    void cuda_convert(DT1_ * dest, const DT2_ * src, const Index count);
  }
}
#endif // FEAT_HAVE_CUDA

#endif //KERNEL_UTIL_CUDA_UTIL_HPP
