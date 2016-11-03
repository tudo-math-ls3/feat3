// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/scatter_axpy_prim.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_scatter_axpy_prim_dv_csr(DT_ * v, const DT_* b, const IT_* col_ind, const DT_* val, const IT_* row_ptr, const DT_ alpha, const Index size)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= size)
          return;

          // skip empty rows
          if(row_ptr[idx] == row_ptr[idx + 1])
            return;

          DT_ sum(0);
          const IT_ end(row_ptr[idx + 1]);
          for (Index i(row_ptr[idx]) ; i < end ; ++i)
          {
            sum += val[i] * (b[col_ind[i]]);
          }
          v[idx] += alpha * sum;
      }

      template <typename DT_, typename IT_, int BlockSize_>
      __global__ void cuda_scatter_axpy_prim_dvb_csr(DT_ * v, const DT_* b, const IT_* col_ind, const DT_* val, const IT_* row_ptr, const DT_ alpha, const Index size)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= size)
          return;

        // skip empty rows
        if(row_ptr[idx] == row_ptr[idx + 1])
          return;

        DT_ bsum[BlockSize_];
        for (int j(0) ; j < BlockSize_ ; ++j)
        {
          bsum[j] = DT_(0);
        }

        const IT_ end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          const DT_ vali(val[i]);
          const IT_ coli(col_ind[i] * BlockSize_);
          for (int j(0) ; j < BlockSize_ ; ++j)
          {
            bsum[j] += vali * b[coli + j];
          }
        }
        for (int j(0) ; j < BlockSize_ ; ++j)
        {
          v[idx * BlockSize_ + j] += alpha * bsum[j];
        }
      }
    }
  }
}

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_, typename IT_>
void ScatterAxpyPrim<Mem::CUDA>::dv_csr(DT_ * v, const DT_* b, const IT_* col_ind, const DT_* val, const IT_* row_ptr, const DT_ alpha, const Index size)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_scatter_axpy_prim_dv_csr<<<grid, block>>>(v, b, col_ind, val, row_ptr, alpha, size);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void ScatterAxpyPrim<Mem::CUDA>::dv_csr(float *, const float*, const unsigned long*, const float*, const unsigned long*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dv_csr(double *, const double*, const unsigned long*, const double*, const unsigned long*, const double alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dv_csr(float *, const float*, const unsigned int*, const float*, const unsigned int*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dv_csr(double *, const double*, const unsigned int*, const double*, const unsigned int*, const double alpha, const Index);

template <typename DT_, typename IT_, int BlockSize_>
void ScatterAxpyPrim<Mem::CUDA>::dvb_csr(DT_ * v, const DT_* b, const IT_* col_ind, const DT_* val, const IT_* row_ptr, const DT_ alpha, const Index size)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_scatter_axpy_prim_dvb_csr<DT_, IT_, BlockSize_><<<grid, block>>>(v, b, col_ind, val, row_ptr, alpha, size);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<float, unsigned long, 1>(float *, const float*, const unsigned long*, const float*, const unsigned long*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<double, unsigned long, 1>(double *, const double*, const unsigned long*, const double*, const unsigned long*, const double alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<float, unsigned int, 1>(float *, const float*, const unsigned int*, const float*, const unsigned int*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<double, unsigned int, 1>(double *, const double*, const unsigned int*, const double*, const unsigned int*, const double alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<float, unsigned long, 2>(float *, const float*, const unsigned long*, const float*, const unsigned long*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<double, unsigned long, 2>(double *, const double*, const unsigned long*, const double*, const unsigned long*, const double alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<float, unsigned int, 2>(float *, const float*, const unsigned int*, const float*, const unsigned int*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<double, unsigned int, 2>(double *, const double*, const unsigned int*, const double*, const unsigned int*, const double alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<float, unsigned long, 3>(float *, const float*, const unsigned long*, const float*, const unsigned long*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<double, unsigned long, 3>(double *, const double*, const unsigned long*, const double*, const unsigned long*, const double alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<float, unsigned int, 3>(float *, const float*, const unsigned int*, const float*, const unsigned int*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::CUDA>::dvb_csr<double, unsigned int, 3>(double *, const double*, const unsigned int*, const double*, const unsigned int*, const double alpha, const Index);
