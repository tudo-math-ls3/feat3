#pragma once

#ifndef __CUDACC__
  static_assert(false, "Header should only be included from cuda context");
#endif

#include "cusparse_v2.h"
#include <kernel/util/cuda_helper.cuh>

namespace FEAT::LAFEM::Arch
{
  /// \cond internal
  namespace Intern
  {
    /// Maps a FEAT data type onto the cuSPARSE storage and compute types.
    template<typename DT_>
    inline void cusparse_deduce_data_type(cudaDataType& dt, cudaDataType& ct)
    {
      if(typeid(DT_) == typeid(double))
      {
        dt = CUDA_R_64F;
        ct = CUDA_R_64F;
      }
      else if(typeid(DT_) == typeid(float))
      {
        dt = CUDA_R_32F;
        ct = CUDA_R_32F;
      }
    #ifdef FEAT_HAVE_HALFMATH
      else if(typeid(DT_) == typeid(Half))
      {
        dt = CUDA_R_16F;
        ct = CUDA_R_32F; // cusparseSpMV does not support computation in half, yet
      }
    #endif
      else
        XABORTM("unsupported data type!");
    }

    /// Maps a FEAT index type onto the cuSPARSE index type.
    template<typename IT_>
    inline cusparseIndexType_t cusparse_deduce_index_type()
    {
      if(sizeof(IT_) == 4u)
        return CUSPARSE_INDEX_32I;
      else if(sizeof(IT_) == 8u)
        return CUSPARSE_INDEX_64I;
      else
        XABORTM("unsupported index type!");
    }
  } // namespace Intern
  /// \endcond

  #define cusparseCheck(ans, funcStr) if((ans) != CUSPARSE_STATUS_SUCCESS) throw InternalError(__func__, __FILE__, __LINE__, (funcStr) + String(" with error: \n") + stringify(cusparseGetErrorString((ans))))

  /**
  * \brief Handles per matrix cusparse spmv workspace
  *
  * We hold one one hand host allocated handles for standard cuda applications.
  * It is assumed that any calle overwrite the value arrays inside of csr_mat, vec_x and vec_y.
  * On the the other hand, the col and index arrays in csr_mat are assumed to be conistent with the initialised buffer
  * (if initialiazed).
  * This allowws us to reuse allocated resources for the same layout.
  */
  struct CuSparseSPMVWorkSpace
  {
    cusparseSpMatDescr_t csr_mat = nullptr;
    cusparseDnVecDescr_t vec_x = nullptr;
    cusparseDnVecDescr_t vec_y = nullptr;
    void* buffer = nullptr;
    std::size_t buffer_size = 0;

    CuSparseSPMVWorkSpace() = default;

    ~CuSparseSPMVWorkSpace()
    {
      if(this->csr_mat)
        cusparseDestroySpMat(static_cast<cusparseSpMatDescr_t>(this->csr_mat));
      if(this->vec_x)
        cusparseDestroyDnVec(vec_x);
      if(this->vec_y)
        cusparseDestroyDnVec(vec_y);
      if(this->buffer)
        cudaFree(this->buffer);
    }
  };

  /**
  * \brief Handles per matrix cusparse bsrmv workspace
  *
  * Holds the persistent legacy cuSPARSE matrix descriptor used by cusparseXbsrmv() so it does not
  * have to be recreated for every apply() call. bsrmv reads the value/col/row pointers directly on
  * every call, so there is nothing else to cache here.
  */
  struct CuSparseBSRWorkSpace
  {
    cusparseMatDescr_t descr = nullptr;

    CuSparseBSRWorkSpace() = default;

    ~CuSparseBSRWorkSpace()
    {
      if(this->descr)
        cusparseDestroyMatDescr(this->descr);
    }
  };


  static inline bool check_size_dnvec(cusparseDnVecDescr_t dn_vec, Index target_size)
  {
    int64_t t_size;
    void* t_vals;
    cudaDataType t_val_type;
    cusparseCheck(cusparseDnVecGet(dn_vec, &t_size, &t_vals, &t_val_type), "Get x Value");
    return Index(t_size) == target_size;
  }

}