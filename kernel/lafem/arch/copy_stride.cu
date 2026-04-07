// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/copy_stride.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

// includes, CUDA
#include <cublas_v2.h>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void CopyStride::exec_cuda_impl(DT_ * r, const DT_ * const x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count)
  {
    cudaDataType dt;
    if (typeid(DT_) == typeid(double))
    {
      dt = CUDA_R_64F;
    }
    else if (typeid(DT_) == typeid(float))
    {
      dt = CUDA_R_32F;
    }
#ifdef FEAT_HAVE_HALFMATH
    else if (typeid(DT_) == typeid(Half))
    {
      dt = CUDA_R_16F;
    }
#endif
    else
    {
      XABORTM("unsupported data type");
    }

    cublasStatus_t status = cublasCopyEx(Util::Intern::cublas_handle, int(count), &x[comp_x], dt, stride_x, &r[comp_r], dt, stride_r);
#ifdef DEBUG
    if (CUBLAS_STATUS_SUCCESS != status)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
#endif
    (void)status;
  }

#ifdef FEAT_HAVE_HALFMATH
  template void CopyStride::exec_cuda_impl(Half *, const Half * const, const int, const int, const int, const int, const Index);
#endif
  template void CopyStride::exec_cuda_impl(float *, const float * const, const int, const int, const int, const int, const Index);
  template void CopyStride::exec_cuda_impl(double *, const double * const, const int, const int, const int, const int, const Index);

} // namespace FEAT::LAFEM::Arch
