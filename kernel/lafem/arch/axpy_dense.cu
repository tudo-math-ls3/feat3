// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/axpy_dense.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

// includes, CUDA
#include <cublas_v2.h>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void AxpyDense::exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const Index size)
  {
    cudaDataType dt;
    cudaDataType et;
    if (typeid(DT_) == typeid(double))
    {
      dt = CUDA_R_64F;
      et = CUDA_R_64F;
    }
    else if (typeid(DT_) == typeid(float))
    {
      dt = CUDA_R_32F;
      et = CUDA_R_32F;
    }
#ifdef FEAT_HAVE_HALFMATH
    else if (typeid(DT_) == typeid(Half))
    {
      dt = CUDA_R_16F;
      et = CUDA_R_32F;
    }
#endif
    else
    {
      XABORTM("unsupported data type");
    }

    cublasStatus_t status = cublasAxpyEx(Util::Intern::cublas_handle, int(size), &alpha, et, x, dt, 1, r, dt, 1, et);

#ifdef DEBUG
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
#endif
    (void)status;
  }

#ifdef FEAT_HAVE_HALFMATH
  template void AxpyDense::exec_cuda_impl(Half *, const Half, const Half * const, const Index);
#endif
  template void AxpyDense::exec_cuda_impl(float *, const float, const float * const, const Index);
  template void AxpyDense::exec_cuda_impl(double *, const double, const double * const, const Index);
} // namespace FEAT::LAFEM::Arch
