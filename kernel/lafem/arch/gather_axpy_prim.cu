// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/gather_axpy_prim.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_gather_axpy_prim_dv_csr(DT_ * b, const DT_* v, const Index* col_ind, const DT_* val, const Index* row_ptr, const DT_ alpha, const Index size, const Index offset)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= size)
          return;

        // loop over all gather-matrix rows
        DT_ sum(DT_(0));
        for (Index i(row_ptr[idx]) ; i < row_ptr[idx + 1] ; ++i)
        {
          sum += DT_(val[i]) * DT_(v[col_ind[i]]);
        }
        b[offset + idx] += alpha * sum;
      }
    }
  }
}

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void GatherAxpyPrim<Mem::CUDA>::dv_csr(DT_ * v, const DT_* b, const Index* col_ind, const DT_* val, const Index* row_ptr, const DT_ alpha, const Index size, const Index offset)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_gather_axpy_prim_dv_csr<<<grid, block>>>(v, b, col_ind, val, row_ptr, alpha, size, offset);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void GatherAxpyPrim<Mem::CUDA>::dv_csr(float *, const float*, const Index*, const float*, const Index*, const float alpha, const Index, const Index);
template void GatherAxpyPrim<Mem::CUDA>::dv_csr(double *, const double*, const Index*, const double*, const Index*, const double alpha, const Index, const Index);
