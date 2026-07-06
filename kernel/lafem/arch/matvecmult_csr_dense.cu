// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
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
  template<typename DT_>
  static inline void _apply_spmv_csr_cuda(DT_* r, const DT_* y, cusparseSpMatDescr_t descr, cusparseDnVecDescr_t dr, cusparseDnVecDescr_t dx, cusparseOperation_t trans, const DT_ alpha,
                                          const DT_ beta, cudaDataType ct, void* buffer, Index cols, Index rows, bool transposed)
  {
    // since async ops are on stream 0 (i.e. default stream), they are always blocking. This might or might not be a feature
    if(y == nullptr)
    {
      cudaCheck(cudaMemsetAsync(r, 0, sizeof(DT_) * (transposed ? cols : rows)), "SPMV Apply memset");
    }
    else if(r != y)
    {
      cudaCheck(cudaMemcpyAsync(r, y, sizeof(DT_) * (transposed ? cols : rows), cudaMemcpyDeviceToDevice), "SPMV Apply memcpy");
    }
    cusparseCheck(cusparseSpMV(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_CSR_ALG1, buffer), "cuSparseCSRSPMV");
  }

  template <typename DT_, typename IT_>
  void MatVecMultCSRDense::exec_cuda_impl(void* workspace, DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows,
                                          const Index cols, const Index nonzeros, const bool transposed)
  {
    const cusparseOperation_t trans = transposed ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;

    cudaDataType dt;
    cudaDataType ct; // compute type
    Intern::cusparse_deduce_data_type<DT_>(dt, ct);

    const cusparseIndexType_t it = Intern::cusparse_deduce_index_type<IT_>();

    cusparseSpMatDescr_t descr = 0;
    cusparseDnVecDescr_t dx = 0;
    cusparseDnVecDescr_t dr = 0;
    void* buffer = nullptr;
    DT_ beta = DT_(1);

    if(workspace)
    {
      CuSparseSPMVWorkSpace* sp_ws = static_cast<CuSparseSPMVWorkSpace*>(workspace);
      descr = sp_ws->csr_mat;
      // only do this in debug
      #ifdef DEBUG
      {
        // check if ptr match
        int64_t t_rows, t_cols, t_nnz;
        void* t_col_ind;
        void* t_row_off;
        void* t_values;
        cusparseIndexType_t t_row_type, t_col_type;
        cusparseIndexBase_t t_idxbase;
        cudaDataType t_valtype;
        cusparseCheck(cusparseCsrGet(descr, &t_rows, &t_cols, &t_nnz, &t_row_off, &t_col_ind, &t_values, &t_row_type, &t_col_type, &t_idxbase, &t_valtype), "cusparseCSRGet");
        // check if structure is the same
        XASSERT(Index(t_rows) == rows);
        XASSERT(Index(t_cols) == cols);
        XASSERT(Index(t_nnz) == nonzeros);
        XASSERT(t_col_ind == (void*)col_idx);
        XASSERT(t_row_off == (void*)row_ptr);
        XASSERT(dt == t_valtype);
      }
      #endif
      // set ptr
      {
        cusparseCheck(cusparseCsrSetPointers(descr, (void*)row_ptr, (void*)col_idx, (void*)val), "cusparseSetPointers");
      }
      dx = sp_ws->vec_x;
      {
        // if the matrix is rectangular, we have to check if the vector descriptor tracks the correct size
        // due to transpose vs. non_transpose
        if(rows != cols)
        {
          Index rsize = transposed ? rows : cols;
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
          Index rsize = transposed ? cols : rows;
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
      cusparseCheck(cusparseCreateCsr(&descr, int64_t(rows), int64_t(cols), int64_t(nonzeros), (void*)row_ptr, (void*)col_idx, (void*)val, it, it, CUSPARSE_INDEX_BASE_ZERO, dt), "cusparseCreateCSR");
      cusparseCheck(cusparseCreateDnVec(&dx, int64_t(transposed?rows:cols), (void*)x, dt), "cusparseCreateDnVec");
      cusparseCheck(cusparseCreateDnVec(&dr, int64_t(transposed?cols:rows), (void*)r, dt), "cusparseCreateDnVec");

      size_t buffer_size(0);
      cusparseCheck(cusparseSpMV_bufferSize(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_CSR_ALG1, &buffer_size), "cuSparseBufferSize");
      buffer = Util::cuda_get_static_memory(buffer_size);
    }

    _apply_spmv_csr_cuda(r, y, descr, dr, dx, trans, alpha, beta, ct, buffer, cols, rows, transposed);

    if(!workspace)
    {
      cusparseDestroySpMat(descr);
      cusparseDestroyDnVec(dx);
      cusparseDestroyDnVec(dr);
    }

  }


  template <typename DT_, typename IT_>
  void MatVecMultCSRDense::init_workspace_cuda(std::shared_ptr<void>& workspace, const Memory::Arbiter& DOXY(val), const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
  {
    // fetch the persistent device pointers of the matrix arrays; these are assumed to stay valid
    // for the whole lifetime of the workspace (see the warning in the header).
    Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
    Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
    // Memory::TypedView<DT_> a_view(val.view(Memory::Location::cuda, Memory::Access::read));

    const cusparseOperation_t trans = transposed ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;

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
    cusparseCheck(cusparseCreateCsr(&descr, int64_t(rows), int64_t(cols), int64_t(nonzeros), (void*)rp_view.get_r(), (void*)ci_view.get_r(), pseudo_ptr, it, it, CUSPARSE_INDEX_BASE_ZERO, dt), "cusparseCreateCsr");

    cusparseDnVecDescr_t dx;
    cusparseCheck(cusparseCreateDnVec(&dx, int64_t(transposed ? rows : cols), pseudo_ptr, dt), "cusparseCreateDnVec");

    cusparseDnVecDescr_t dr;
    cusparseCheck(cusparseCreateDnVec(&dr, int64_t(transposed ? cols : rows), pseudo_ptr, dt), "cusparseCreateDnVec");

    // alpha/beta are only inspected for their type here; the actual values used on apply are irrelevant
    DT_ alpha = DT_(1);
    DT_ beta = DT_(1);

    size_t buffer_size(0);
    cusparseCheck(cusparseSpMV_bufferSize(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_CSR_ALG1, &buffer_size), "cusparseBuffersize");

    // allocate the persistent SpMV buffer
    void* buffer = nullptr;
    if(buffer_size > 0u)
    {
      cudaCheck(cudaMalloc(&buffer, buffer_size), "cudaMalloc");
    }
    cusparseCheck(cusparseSpMV_preprocess(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_CSR_ALG1, buffer), "cuSparseBCSRSPMV");

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
  }

  template void MatVecMultCSRDense::exec_cuda_impl(void*, float *, const float, const float * const, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRDense::exec_cuda_impl(void*, double *, const double, const double * const, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultCSRDense::exec_cuda_impl(void*, Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  #endif

  template void MatVecMultCSRDense::exec_cuda_impl(void*, float *, const float, const float * const, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRDense::exec_cuda_impl(void*,double *, const double, const double * const, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultCSRDense::exec_cuda_impl(void*, Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  #endif

  template void MatVecMultCSRDense::init_workspace_cuda<float, std::uint32_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRDense::init_workspace_cuda<double, std::uint32_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultCSRDense::init_workspace_cuda<Half, std::uint32_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const bool);
  #endif

  template void MatVecMultCSRDense::init_workspace_cuda<float, std::uint64_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRDense::init_workspace_cuda<double, std::uint64_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultCSRDense::init_workspace_cuda<Half, std::uint64_t>(std::shared_ptr<void>&, const Memory::Arbiter&, const Memory::Arbiter&, const Memory::Arbiter&, const Index, const Index, const Index, const bool);
  #endif
} // namespace FEAT::LAFEM::Arch
