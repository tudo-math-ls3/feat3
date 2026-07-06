// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_bcsr_dense.hpp>
#include <kernel/lafem/arch/matvecmult_csr_dense.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/half.hpp>

#include "cusparse_v2.h"
#include <kernel/util/cuda_helper.cuh>
#include "cusparse_helper.cuh"

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_apply_bcsr(DT_ * r, const DT_ alpha, const DT_ * x, const DT_ * val, const IT_ * col_idx,
    const IT_ * row_ptr, const Index count, const int block_height, const int block_width)
  {
    // grid strided for loop
    for(Index idx = threadIdx.x + blockDim.x * blockIdx.x; idx < count; idx += blockDim.x * gridDim.x)
    {
      /// \todo remove hardcoded number
      DT_ bsum[10];
      for (int j(0) ; j < block_height ; ++j)
      {
        bsum[j] = DT_(0);
      }
      const IT_ end(row_ptr[idx + 1]);
      for (IT_ i(row_ptr[idx]) ; i < end ; ++i)
      {
        for (int h(0) ; h < block_height ; ++h)
        {
          for (int w(0) ; w < block_width ; ++w)
          {
            bsum[h] += val[i * block_height * block_width + h * block_width + w] * x[col_idx[i] * block_width + w];
          }
        }
      }
      for (int j(0) ; j < block_height ; ++j)
      {
        r[idx * block_height + j] += alpha * bsum[j];
      }
    }
  }

  template <typename DT_, typename IT_>
  void bcsr_intern_cuda(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index /*nonzeros*/,
    const int block_height, const int block_width, const bool transposed)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    Index rsize = (transposed ? cols*Index(block_width) : rows*Index(block_height));
    if(y == nullptr)
    {
      cudaCheck(cudaMemsetAsync(r, 0, sizeof(DT_) * rsize), "SPMV Apply memset");
    }
    else if (r != y)
    {
      cudaCheck(cudaMemcpyAsync(r, y, sizeof(DT_) * rsize, cudaMemcpyDeviceToDevice), "SPMV Apply memcpy");
    }

    XASSERTM(!transposed, "transposed BCSR matvec for CUDA not implemented");
    XASSERT(block_height < 10); // see cuda_apply_bcsr implementation

    cuda_apply_bcsr<DT_, IT_><<<grid, block>>>(r, alpha, x, val, col_idx, row_ptr, rows, block_height, block_width);

#ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  }

  // maps a scalar data type to its CUDA 2-wide vector type for coalesced 16-byte (double) / 8-byte
  // (float) loads of the two contiguous values of a 2x1 / 1x2 block. Types without a vector partner
  // (e.g. Half) report available == false and fall back to scalar loads.
  template<typename DT_> struct BcsrVec2Trait { static constexpr bool available = false; typedef DT_ type; };
  template<> struct BcsrVec2Trait<float>  { static constexpr bool available = true; typedef float2 type; };
  template<> struct BcsrVec2Trait<double> { static constexpr bool available = true; typedef double2 type; };

  /**
   * \brief Compile-time-specialized rectangular BCSR matvec kernel (one thread per block-row).
   *
   * In contrast to \c cuda_apply_bcsr, the block dimensions are compile-time template parameters, so
   * the inner loops unroll, the index arithmetic is static and the per-thread accumulator lives in
   * registers (no local-memory spill). The y/beta term is merged inside the kernel
   * (r = (y ? y : 0) + alpha * A*x), so the launcher does not need the separate memset/memcpy pre-pass.
   *
   * This is the thread-per-block-row variant, best suited to the "tall" N x 1 blocks (e.g. the 2x1 /
   * 3x1 gradient block B), which have many short block-rows and therefore already saturate the grid.
   * For the 2x1 shape the two contiguous block values are loaded as one vectorized transaction.
   */
  template <typename DT_, typename IT_, int bh_, int bw_>
  __global__ void cuda_apply_bcsr_ct(DT_ * r, const DT_ alpha, const DT_ * x, const DT_ * y, const DT_ * val,
    const IT_ * col_idx, const IT_ * row_ptr, const Index count)
  {
    // grid-strided loop, one thread per block-row
    for(Index idx = threadIdx.x + blockDim.x * blockIdx.x; idx < count; idx += blockDim.x * gridDim.x)
    {
      DT_ bsum[bh_];
      #pragma unroll
      for(int h(0) ; h < bh_ ; ++h)
        bsum[h] = DT_(0);

      const IT_ end(row_ptr[idx + 1]);
      if constexpr (bh_ == 2 && bw_ == 1 && BcsrVec2Trait<DT_>::available)
      {
        // 2x1: the two block values val[2i], val[2i+1] are contiguous and share the single x entry
        typedef typename BcsrVec2Trait<DT_>::type Vec2;
        const Vec2 * const val2 = reinterpret_cast<const Vec2*>(val);
        for(IT_ i(row_ptr[idx]) ; i < end ; ++i)
        {
          const Vec2 v = val2[i];
          const DT_ xv = x[col_idx[i]];
          bsum[0] += v.x * xv;
          bsum[1] += v.y * xv;
        }
      }
      else
      {
        for(IT_ i(row_ptr[idx]) ; i < end ; ++i)
        {
          const IT_ col(col_idx[i] * IT_(bw_));
          const DT_ * const vblk = val + i * (bh_ * bw_);
          #pragma unroll
          for(int h(0) ; h < bh_ ; ++h)
            #pragma unroll
            for(int w(0) ; w < bw_ ; ++w)
              bsum[h] += vblk[h * bw_ + w] * x[col + IT_(w)];
        }
      }

      if(y != nullptr)
      {
        #pragma unroll
        for(int h(0) ; h < bh_ ; ++h)
          r[idx * Index(bh_) + Index(h)] = y[idx * Index(bh_) + Index(h)] + alpha * bsum[h];
      }
      else
      {
        #pragma unroll
        for(int h(0) ; h < bh_ ; ++h)
          r[idx * Index(bh_) + Index(h)] = alpha * bsum[h];
      }
    }
  }

  /// launcher for the compile-time-specialized rectangular BCSR matvec kernel (thread per block-row)
  template <typename DT_, typename IT_, int bh_, int bw_>
  void bcsr_intern_cuda_ct(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const bool transposed)
  {
    XASSERTM(!transposed, "transposed BCSR matvec for CUDA not implemented");
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    cuda_apply_bcsr_ct<DT_, IT_, bh_, bw_><<<grid, block>>>(r, alpha, x, y, val, col_idx, row_ptr, rows);

#ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  }

  /**
   * \brief Routes a real BCSR matvec to the compile-time-specialized custom kernels.
   *
   * Covers the block shapes that occur in 2D/3D saddle-point systems: the square velocity blocks
   * (2x2 / 3x3), the "tall" gradient blocks (2x1 / 3x1, thread-per-block-row) and the "wide"
   * divergence blocks (1x2 / 1x3, warp-per-block-row). Any other shape falls back to the generic
   * runtime-dimension kernel. Used for the rectangular blocks (which cuSPARSE cannot handle) and, when
   * cuSPARSE is disabled for BCSR (FEAT_NO_CUSPARSE_BCSR), for the square blocks as well.
   */
  template <typename DT_, typename IT_>
  void bcsr_intern_cuda_dispatch(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros,
    const int block_height, const int block_width, const bool transposed)
  {
    switch((block_height << 8) | block_width)
    {
    case (2 << 8) | 2: // square velocity block, 2D
      bcsr_intern_cuda_ct<DT_, IT_, 2, 2>(r, alpha, x, y, val, col_idx, row_ptr, rows, transposed);
      break;
    case (3 << 8) | 3: // square velocity block, 3D
      bcsr_intern_cuda_ct<DT_, IT_, 3, 3>(r, alpha, x, y, val, col_idx, row_ptr, rows, transposed);
      break;
    case (2 << 8) | 1: // gradient block B, 2D
      bcsr_intern_cuda_ct<DT_, IT_, 2, 1>(r, alpha, x, y, val, col_idx, row_ptr, rows, transposed);
      break;
    case (3 << 8) | 1: // gradient block B, 3D
      bcsr_intern_cuda_ct<DT_, IT_, 3, 1>(r, alpha, x, y, val, col_idx, row_ptr, rows, transposed);
      break;
    case (1 << 8) | 2: // divergence block D, 2D (warp-per-row)
      bcsr_intern_cuda_ct<DT_, IT_, 1, 2>(r, alpha, x, y, val, col_idx, row_ptr, rows, transposed);
      break;
    case (1 << 8) | 3: // divergence block D, 3D (warp-per-row)
      bcsr_intern_cuda_ct<DT_, IT_, 1, 3>(r, alpha, x, y, val, col_idx, row_ptr, rows, transposed);
      break;
    default:
      bcsr_intern_cuda<DT_, IT_>(r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, block_height, block_width, transposed);
      break;
    }
  }

/// since bsr is only supported with cuda version 13.0.1, we have to branch at this point
#ifdef FEAT_CUSPARSE_BSR_SPMV_SUPPORTED
  template<typename DT_>
  static inline void _apply_spmv_bcsr_cuda(DT_* r, const DT_* y, cusparseSpMatDescr_t descr, cusparseDnVecDescr_t dr, cusparseDnVecDescr_t dx, cusparseOperation_t trans, const DT_ alpha,
                                          const DT_ beta, cudaDataType ct, void* buffer, Index cols, Index rows, int bheight, int bwidth, bool transposed)
  {
    Index rsize = (transposed ? cols*Index(bwidth) : rows*Index(bheight));
    // since async ops are on stream 0 (i.e. default stream), they are always blocking. This might or might not be a feature
    if(y == nullptr)
    {
      cudaCheck(cudaMemsetAsync(r, 0, sizeof(DT_) * rsize), "SPMV Apply memset");
    }
    else if(r != y)
    {
      cudaCheck(cudaMemcpyAsync(r, y, sizeof(DT_) * rsize, cudaMemcpyDeviceToDevice), "SPMV Apply memcpy");
    }
    cusparseCheck(cusparseSpMV(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_BSR_ALG1, buffer), "cuSparseBCSRSPMV");
  }

  template <typename DT_, typename IT_>
  void bcsr_intern_cusparse_v13(void* workspace, DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr,
                            const Index rows, const Index cols, const Index nonzeros, const int bheight, const int bwidth, const bool transposed)
  {
    const cusparseOperation_t trans = transposed ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;

    cudaDataType dt;
    cudaDataType ct; // compute type
    Intern::cusparse_deduce_data_type<DT_>(dt, ct);

    const cusparseIndexType_t it = Intern::cusparse_deduce_index_type<IT_>();
    const cusparseOrder_t order = CUSPARSE_ORDER_ROW;

    cusparseSpMatDescr_t descr = 0;
    cusparseDnVecDescr_t dx = 0;
    cusparseDnVecDescr_t dr = 0;
    void* buffer = nullptr;
    DT_ beta = DT_(1);

    if(workspace)
    {
      CuSparseSPMVWorkSpace* sp_ws = static_cast<CuSparseSPMVWorkSpace*>(workspace);
      descr = sp_ws->csr_mat;
      // set ptr
      {
        cusparseCheck(cusparseSpMatSetValues(descr, (void*)val), "cusparseSetValues");
      }
      dx = sp_ws->vec_x;
      {
        // generally check if the vectors have the correct size
        {
          Index rsize = transposed ? rows*Index(bheight) : cols*Index(bwidth);
          if(!check_size_dnvec(dx, rsize))
          {
            // delete old ptr
            cusparseCheck(cusparseDestroyDnVec(sp_ws->vec_x), "destroyDnVec");
            // create new one
            cusparseCheck(cusparseCreateDnVec(&sp_ws->vec_x, int64_t(rsize), (void*)x, dt), "Cusparse Create dnvec");
            dx = sp_ws->vec_x;
          }
        }
        cusparseCheck(cusparseDnVecSetValues(dx, (void*)x), "Set x Value");
      }
      dr = sp_ws->vec_y;
      {
        if(rows != cols)
        {
          Index rsize = transposed ? cols*Index(bwidth) : rows*Index(bheight);
          if(!check_size_dnvec(dr, rsize))
          {
            // delete old ptr
            cusparseCheck(cusparseDestroyDnVec(sp_ws->vec_y), "destroyDnVec");
            // create new one
            cusparseCheck(cusparseCreateDnVec(&sp_ws->vec_y, int64_t(rsize), (void*)r, dt), "Cusparse Create dnvec");
            dr = sp_ws->vec_y;
          }
        }
        cusparseCheck(cusparseDnVecSetValues(dr, (void*)r), "Set r Value");
      }
      buffer = sp_ws->buffer;
    }
    else
    {
      cusparseCheck(cusparseCreateBsr(&descr, int64_t(rows), int64_t(cols), int64_t(nonzeros), int64_t(bheight), int64_t(bwidth), (void*)row_ptr, (void*)col_idx, (void*)val, it, it, CUSPARSE_INDEX_BASE_ZERO, dt, order), "cusparseCreateBCSR");
      cusparseCheck(cusparseCreateDnVec(&dx, int64_t(transposed?(rows*Index(bheight)):(cols*Index(bwidth))), (void*)x, dt), "cusparseCreateDnVec");
      cusparseCheck(cusparseCreateDnVec(&dr, int64_t(transposed?(cols*Index(bwidth)):(rows*Index(bheight))), (void*)r, dt), "cusparseCreateDnVec");

      size_t buffer_size(0);
      cusparseCheck(cusparseSpMV_bufferSize(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_BSR_ALG1, &buffer_size), "cuSparseBufferSize");
      buffer = Util::cuda_get_static_memory(buffer_size);
    }

    _apply_spmv_bcsr_cuda(r, y, descr, dr, dx, trans, alpha, beta, ct, buffer, cols, rows, bheight, bwidth, transposed);

    if(!workspace)
    {
      cusparseCheck(cusparseDestroySpMat(descr), "cusparseDestroy");
      cusparseDestroyDnVec(dx);
      cusparseDestroyDnVec(dr);
    }

#ifdef DEBUG
    cudaDeviceSynchronize();
    cudaCheckLastErr()
#endif
  }
#else

  void cusparse_apply_bcsr(cusparseDirection_t dir, cusparseOperation_t trans,
    int m, int n, int nnz,
    const float * alpha, const cusparseMatDescr_t descrA,
    const float * csrVal, const int * csrRowPtr, const int *csrColInd,
    int block_dim,
    const float * x, const float * beta, float * y)
  {
    cusparseStatus_t status;
    status = cusparseSbsrmv(Util::Intern::cusparse_handle, dir, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
      csrColInd, block_dim, x, beta, y);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrmv failed with status code: " + stringify(status));
  }

  void cusparse_apply_bcsr(cusparseDirection_t dir, cusparseOperation_t trans,
    int m, int n, int nnz,
    const double * alpha, const cusparseMatDescr_t descrA,
    const double * csrVal, const int * csrRowPtr, const int *csrColInd,
    int block_dim,
    const double * x, const double * beta, double * y)
  {
    cusparseStatus_t status;
    status = cusparseDbsrmv(Util::Intern::cusparse_handle, dir, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
      csrColInd, block_dim, x, beta, y);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrmv failed with status code: " + stringify(status));
  }

  template <typename DT_, typename IT_>
  void bcsr_intern_cusparse_v12(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr,
                            const Index rows, const Index cols, const Index nonzeros, const int block_size, const bool transposed)
  {
    Index rsize = (transposed ? cols : rows) * Index(block_size);
    if(y == nullptr)
    {
      cudaCheck(cudaMemsetAsync(r, 0, sizeof(DT_) * rsize), "SPMV Apply memset");
    }
    else if(r != y)
    {
      cudaCheck(cudaMemcpyAsync(r, y, sizeof(DT_) * rsize, cudaMemcpyDeviceToDevice), "SPMV Apply memcpy");
    }

    cusparseMatDescr_t descr=0;
    cusparseCreateMatDescr(&descr);
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

    cusparseOperation_t trans = (transposed ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);

    DT_ b = DT_(y == nullptr ? 0 : 1);
    cusparse_apply_bcsr(CUSPARSE_DIRECTION_ROW, trans, (int)rows, (int)cols, (int)nonzeros, &alpha, descr, val, (int*)row_ptr, (int*)col_idx,
      block_size, x, &b, r);

    cusparseDestroyMatDescr(descr);
  }
#endif

  template <typename DT_, typename IT_>
  void MatVecMultBCSRDense::exec_cuda_impl(void* workspace, DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const int block_height, const int block_width, const bool transposed)
  {
    if(block_height*block_width == 1)
    {
      // standard CSR always uses cuSPARSE
      MatVecMultCSRDense::exec_cuda_impl<DT_, IT_>(workspace, r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, transposed);
    }
    else
    {
#ifdef FEAT_NO_CUSPARSE_BCSR
      // cuSPARSE is disabled for real (block_dim>1) BCSR at configure time: always use the custom
      // compile-time-specialized kernels (both square and rectangular blocks)
      bcsr_intern_cuda_dispatch<DT_, IT_>(r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, block_height, block_width, transposed);
#else
      // sadly, only quadratic blocks are allowed in cuSPARSE
      if(block_height == block_width)
      {
  #ifdef FEAT_CUSPARSE_BSR_SPMV_SUPPORTED
        bcsr_intern_cusparse_v13<DT_, IT_>(workspace, r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, block_height, block_width, transposed);
  #else
        if(sizeof(IT_) == 4)
          bcsr_intern_cusparse_v12<DT_, IT_>(r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, block_height, transposed);
        else
          bcsr_intern_cuda_dispatch<DT_, IT_>(r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, block_height, block_width, transposed);
  #endif
      }
      else
      {
        // rectangular blocks: cuSPARSE does not support them, so use the custom kernels
        bcsr_intern_cuda_dispatch<DT_, IT_>(r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, block_height, block_width, transposed);
      }
#endif
    }
  }

  template <typename DT_, typename IT_>
  void MatVecMultBCSRDense::init_workspace_cuda(std::shared_ptr<void>& workspace, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const int block_height, const int block_width, const bool transposed)
  {
    if(block_height*block_width == 1)
    {
      MatVecMultCSRDense::init_workspace_cuda<DT_, IT_>(workspace, val, col_idx, row_ptr, rows, cols, nonzeros, transposed);
      return;
    }
#ifdef FEAT_NO_CUSPARSE_BCSR
    // cuSPARSE disabled for real BCSR: the custom kernels need no persistent workspace
    return;
#endif
#ifdef FEAT_CUSPARSE_BSR_SPMV_SUPPORTED
    // nothing to do in this case
    if(block_height != block_width)
      return;
    Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
    Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
    const cusparseOperation_t trans = transposed ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;
    const cusparseOrder_t order = CUSPARSE_ORDER_ROW;

    cudaDataType dt;
    cudaDataType ct; // compute type
    Intern::cusparse_deduce_data_type<DT_>(dt, ct);

    const cusparseIndexType_t it = Intern::cusparse_deduce_index_type<IT_>();

    // The SpMV buffer size does not depend on the actual vector data, only on the vector sizes and
    // types. As the concrete x/r device pointers are only known at apply time, we hand cuSPARSE a
    // pseudo pointer purely to satisfy the descriptor creation.
    // This can not be a nullptr, since there is an assert statement int the cuda interface
    void* const pseudo_ptr = (void*)0x048;

    // create the persistent CSR matrix descriptor from the (fixed) matrix arrays
    cusparseSpMatDescr_t descr = nullptr;
    // cusparseCheck(cusparseCreateBsr(&descr, int64_t(rows), int64_t(cols), int64_t(nonzeros), (void*)rp_view.get_r(), (void*)ci_view.get_r(), pseudo_ptr, it, it, CUSPARSE_INDEX_BASE_ZERO, dt), "cusparseCreateCsr");
    cusparseCheck(cusparseCreateBsr(&descr, int64_t(rows), int64_t(cols), int64_t(nonzeros), int64_t(block_height), int64_t(block_width), (void*)rp_view.get_r(), (void*)ci_view.get_r(), pseudo_ptr, it, it, CUSPARSE_INDEX_BASE_ZERO, dt, order), "cusparseCreateBCSR");

    cusparseDnVecDescr_t dx;
    cusparseCheck(cusparseCreateDnVec(&dx, int64_t(transposed ? (rows*Index(block_height)) : (cols*Index(block_width))), pseudo_ptr, dt), "cusparseCreateDnVec");

    cusparseDnVecDescr_t dr;
    cusparseCheck(cusparseCreateDnVec(&dr, int64_t(transposed ? (cols*Index(block_width)) : (rows*Index(block_height))), pseudo_ptr, dt), "cusparseCreateDnVec");

    // alpha/beta are only inspected for their type here; the actual values used on apply are irrelevant
    DT_ alpha = DT_(1);
    DT_ beta = DT_(1);

    size_t buffer_size(0);
    cusparseCheck(cusparseSpMV_bufferSize(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_BSR_ALG1, &buffer_size), "cusparseBuffersize");

    // allocate the persistent SpMV buffer
    void* buffer = nullptr;
    if(buffer_size > 0u)
    {
      cudaCheck(cudaMalloc(&buffer, buffer_size), "cudaMalloc");
    }
    cusparseCheck(cusparseSpMV_preprocess(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_BSR_ALG1, buffer), "cuSparseBCSRSPMV");

    // hand ownership of the persistent descriptor and buffer to the type-erased workspace; the
    // CuSparseSPMVWorkSpace destructor releases both once the shared_ptr is reset.
    // Note that the Deleter will correctly call the destructor defined by cusparseworkspace
    auto ws = std::make_shared<CuSparseSPMVWorkSpace>();
    ws->csr_mat = descr;
    ws->vec_x = dx;
    ws->vec_y = dr;
    ws->buffer = buffer;
    ws->buffer_size = buffer_size;
    workspace = std::move(ws);
  #endif
  }

  template void MatVecMultBCSRDense::exec_cuda_impl(void*, float *, const float, const float * const, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const int, const int, const bool);
  template void MatVecMultBCSRDense::exec_cuda_impl(void*, double *, const double, const double * const, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const int, const int, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultBCSRDense::exec_cuda_impl(void*, Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const int, const int, const bool);
  #endif

  template void MatVecMultBCSRDense::exec_cuda_impl(void*, float *, const float, const float * const, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const int, const int, const bool);
  template void MatVecMultBCSRDense::exec_cuda_impl(void*, double *, const double, const double * const, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const int, const int, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultBCSRDense::exec_cuda_impl(void*, Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const int, const int, const bool);
  #endif

  template void MatVecMultBCSRDense::init_workspace_cuda<float, std::uint32_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const int, const int, const bool);
  template void MatVecMultBCSRDense::init_workspace_cuda<double, std::uint32_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const int, const int, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultBCSRDense::init_workspace_cuda<Half, std::uint32_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const int, const int, const bool);
  #endif

  template void MatVecMultBCSRDense::init_workspace_cuda<float, std::uint64_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const int, const int, const bool);
  template void MatVecMultBCSRDense::init_workspace_cuda<double, std::uint64_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const int, const int, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultBCSRDense::init_workspace_cuda<Half, std::uint64_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const int, const int, const bool);
  #endif
} // namespace FEAT::LAFEM::Arch
