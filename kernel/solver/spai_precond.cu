// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

#include <cublas_v2.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>

#include <vector>

//#include <thrust/copy.h>
//#include <thrust/execution_policy.h>
//#include <thrust/reduce.h>

using namespace FEAT;

namespace FEAT
{
    namespace Solver
    {
        namespace Intern
        {

            /***********************************
             *        ERROR HANDLING           *
             ***********************************/

#define SPAI_API_CALL(err) cuda_spai_check_for_api_errors(err,__func__,__FILE__,__LINE__)
            inline void cuda_spai_check_for_api_errors(cudaError_t err, const char* func, const char* file, const long line)
            {
                if (err != cudaSuccess)
                    throw InternalError(func, file, line, "CUDA error occurred in API call: " + stringify(cudaGetErrorString(err)));
            }

#define SPAI_CHECK_KERNEL_ERROR() cuda_spai_check_for_kernel_errors(__func__,__FILE__,__LINE__)
            inline void cuda_spai_check_for_kernel_errors(const char* func, const char* file, const long line)
            {
                // catch launch errors
                cudaError_t last_error(cudaGetLastError());
                if (last_error != cudaSuccess)
                    throw InternalError(func, file, line, "CUDA error occurred in kernel launch: " + stringify(cudaGetErrorString(last_error)));
#ifdef FEAT_DEBUG_MODE
                // catch execution errors -> may affect performance!
                if (cudaDeviceSynchronize() != cudaSuccess)
                    throw InternalError(func, file, line, "CUDA error occurred in kernel execution: " + stringify(cudaGetErrorString(last_error)));
#endif
            }

#define SPAI_CUBLAS_CALL(err) cuda_spai_check_for_cublas_errors(err,__func__,__FILE__,__LINE__)
            inline void cuda_spai_check_for_cublas_errors(cublasStatus_t err, const char* func, const char* file, const long line)
            {
                if (err != CUBLAS_STATUS_SUCCESS)
                    throw InternalError(func, file, line, "CUDA error occurred in cublas call.");
            }


            /***************************
             *        KERNEL           *
             ***************************/

            // calculate lengths of rows from A csr matrix layout
            __global__ void
                cuda_spai_csr_rowlen_kernel(const unsigned int m, const unsigned int *rowptr, unsigned int *rowlens)
            {
                unsigned int myrow = blockIdx.x * blockDim.x + threadIdx.x;
                if (myrow < m) rowlens[myrow] = rowptr[myrow + 1] - rowptr[myrow];
            }

            // write Adags in column major order
            template<class DT_>
            __global__ void
                cuda_spai_construct_adags_kernel(
                    const unsigned int this_chunk_size, const unsigned int first_row_in_chunk, const unsigned int mrl,
                    const unsigned int *rowptr, const unsigned int *colind, const DT_ *Aval, const unsigned int *rowlens,
                    DT_ *Adags, unsigned int *Jks, unsigned int *lensJk)
            {
                const unsigned int kchunk = blockIdx.x * blockDim.x + threadIdx.x;
                if (kchunk < this_chunk_size)
                {
                  const unsigned int k = kchunk + first_row_in_chunk;
                  DT_ *thisadag = Adags + kchunk*mrl*mrl*mrl;
                  unsigned int *thisjk = Jks + kchunk*mrl*mrl;

                  // all entries in row k
                  int nj = 0; // count the size of set Jk
                  for (int i = 0; i < rowlens[k]; i++)
                  {
                        const int ii = colind[rowptr[k] + i];
                        for (int j = 0; j < rowlens[ii]; j++)
                        {
                            const int jj = colind[rowptr[ii] + j];
                            // already in matrix?
                            for (int l = 0; l < nj; l++)
                            {
                                if (thisjk[l] == jj)
                                {
                                    // insert at this position
                                    thisadag[i*mrl*mrl + l] = Aval[rowptr[ii] + j];
                                    /// \todo remove goto
                                    goto nextj; // break inner loop and skip rest of outer loop
                                }
                            }
                            // not yet in matrix -> new line
                            thisjk[nj] = jj;
                            thisadag[i*mrl*mrl + nj] = Aval[rowptr[ii] + j];
                            nj++;
                        nextj:
                            continue;
                        }
                    }
                    // save actual number of rows of Adag_k
                    lensJk[kchunk] = nj;
                }
            }

            // add unit matrix blocks at the end of each Adag
            template<class DT_>
            __global__ void
                cuda_spai_pad_adags_kernel(const unsigned int this_chunk_size, const unsigned int first_row_in_chunk,
                    const unsigned int mrl, const unsigned int *lensJk, const unsigned int *rowlens, DT_ *Adags)
            {
                const unsigned int kchunk = blockIdx.x * blockDim.x + threadIdx.x;
                if (kchunk < this_chunk_size)
                {
                  const unsigned int k = kchunk + first_row_in_chunk;
                  DT_ *thisadag = Adags + kchunk*mrl*mrl*mrl;

                  int blockdim = mrl - rowlens[k];
                  for (int offset = 0; offset < blockdim; offset++)
                    thisadag[(rowlens[k] + offset)*mrl*mrl + (lensJk[kchunk] + offset)] = 1.0;
                }
            }

            // construct right hand sides of the minimization problems
            template<class DT_>
            __global__ void
                cuda_spai_fill_rhs(
                    const unsigned int this_chunk_size, const unsigned int first_row_in_chunk, const unsigned int mrl,
                    DT_ *Erhs, const unsigned int *Jks, const unsigned int *lenJks)
            {
                const unsigned int kchunk = blockIdx.x * blockDim.x + threadIdx.x;
                if (kchunk < this_chunk_size)
                {
                  const unsigned int k = kchunk + first_row_in_chunk;

                  const unsigned int *thisJk = Jks + kchunk*mrl*mrl;
                  DT_ *thisrhs = Erhs + kchunk*mrl*mrl;
                  const unsigned int thislen = *(lenJks + kchunk);

                  for (int i = 0; i < thislen; i++)
                  {
                    if (thisJk[i] == k)
                    {
                      thisrhs[i] = 1.0;
                      // Jks are unique by construction
                      break;
                    }
                  }
                }
            }

            // set batch pointers to the parts of Adags/Erhs
            template<class DT_>
            __global__ void
                cuda_spai_fill_batch_pointers(
                    const unsigned int this_chunk_size, const unsigned int maxrowlen,
                    DT_ **Aptrs, DT_ **Cptrs, DT_ *Adags, DT_ *Erhs)
            {
                const unsigned int kchunk = blockIdx.x * blockDim.x + threadIdx.x;

                if (kchunk < this_chunk_size)
                {
                    Aptrs[kchunk] = Adags + kchunk*maxrowlen*maxrowlen*maxrowlen;
                    Cptrs[kchunk] = Erhs + kchunk*maxrowlen*maxrowlen;
                }
            }

            // copy minimzation results to the CSR matrix M
            template<class DT_>
            __global__ void
                cuda_spai_copy_to_m(
                    const unsigned int this_chunk_size, const unsigned int first_row_in_chunk, const unsigned int maxrowlen,
                    DT_ *M_val, const unsigned int *M_rowptr, const unsigned int *rowlens, const DT_ *Erhs)
            {
                const unsigned int kchunk = blockIdx.x * blockDim.x + threadIdx.x;
                const unsigned int k = kchunk + first_row_in_chunk;

                if (kchunk < this_chunk_size)
                {
                    DT_ *thisrow = M_val + M_rowptr[k]; // dst
                    const DT_ *thisres = Erhs + kchunk*maxrowlen*maxrowlen; // src

                    for (int i = 0; i < rowlens[k]; i++)
                        thisrow[i] = thisres[i];
                }
            }

            /************************************
            *        HELPER FUNCTIONS           *
            *************************************/

            /*
            // helper function to show content of a gpu vector on screen
            template <typename T>
            void cuda_spai_print_gpu_array(const T *gpu_data, const size_t n)
            {
                std::vector<T> host_data(n);
                SPAI_API_CALL(cudaMemcpy(host_data.data(), gpu_data, n*sizeof(T), cudaMemcpyDeviceToHost));
                std::cout << "GPU array of length " << n << ":" << std::endl;
                for (size_t i = 0; i < n; i++)
                    std::cout << "(" << i << ")\t:" << host_data[i] << std::endl;
            }
            */

            // wrap cublas calls for single / double precision
            template <typename DT_>
            void cuda_spai_gelsBatched_wrapper(cublasHandle_t handle,
                cublasOperation_t trans, unsigned int m, unsigned int n, unsigned int nrhs,
                DT_ *Aarray[], unsigned int lda, DT_ *Carray[], unsigned int ldc,
                int *info, int *devInfoArray, unsigned int batchSize);

            template<>
            void cuda_spai_gelsBatched_wrapper(cublasHandle_t handle,
                cublasOperation_t trans, unsigned int m, unsigned int n, unsigned int nrhs,
                double *Aarray[], unsigned int lda, double *Carray[], unsigned int ldc,
                int *info, int *devInfoArray, unsigned int batchSize)
            {
                SPAI_CUBLAS_CALL(cublasDgelsBatched(handle,
                    trans, m, n, nrhs,
                    Aarray, lda, Carray, ldc,
                    info, devInfoArray, batchSize));
                SPAI_CHECK_KERNEL_ERROR();
            }

            template<>
            void cuda_spai_gelsBatched_wrapper(cublasHandle_t handle,
                cublasOperation_t trans, unsigned int m, unsigned int n, unsigned int nrhs,
                float *Aarray[], unsigned int lda, float *Carray[], unsigned int ldc,
                int *info, int *devInfoArray, unsigned int batchSize)
            {
                SPAI_CUBLAS_CALL(cublasSgelsBatched(handle,
                    trans, m, n, nrhs,
                    Aarray, lda, Carray, ldc,
                    info, devInfoArray, batchSize));
                SPAI_CHECK_KERNEL_ERROR();
            }

            /************************************
            *        EXPORTED FUNCTIONS         *
            *************************************/

            unsigned int cuda_spai_maxrowlen(const unsigned int m, const unsigned int *rowptr, unsigned int *rowlens)
            {
                dim3 rwblock(256);
                dim3 rwgrid((unsigned)ceil(m / (double)rwblock.x));

                cuda_spai_csr_rowlen_kernel <<<rwgrid, rwblock >>> (m, rowptr, rowlens);
                SPAI_CHECK_KERNEL_ERROR();
                if (cudaDeviceSynchronize() != cudaSuccess)
                    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred");
                thrust::device_ptr<unsigned int> result_dev = thrust::max_element(thrust::device_ptr<unsigned int>(rowlens),
                    thrust::device_ptr<unsigned int>(rowlens) + m);
                if (cudaDeviceSynchronize() != cudaSuccess)
                    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred");
                unsigned int maxrowlen;
                SPAI_API_CALL(cudaMemcpy(&maxrowlen, result_dev.get(), sizeof(unsigned int), cudaMemcpyDeviceToHost));
                return maxrowlen;
            }

            template <typename DT_>
            void cuda_spai_construct_rows(const unsigned int this_chunk_size, const unsigned int first_row_in_chunk,
                const DT_ *A_val, const unsigned int *A_colind, const unsigned int *A_rowptr, DT_ *M_val,
                const unsigned int *rowlens, const unsigned int maxrowlen)
            {
                dim3 rwblock(256);
                dim3 rwgrid((unsigned)ceil(this_chunk_size / (double)rwblock.x));

                // allocate memory on the GPU
                DT_ *Adags, *Erhs;   // parts of the min problem
                unsigned int *Jks, *lenJks;      // indices corresponding to the rows of Adag
                                        // and the cardinality of these sets
                SPAI_API_CALL(cudaMalloc(&Adags,
                    sizeof(*Adags) * this_chunk_size * maxrowlen * maxrowlen * maxrowlen));
                SPAI_API_CALL(cudaMalloc(&Erhs,
                    sizeof(*Erhs) * this_chunk_size * maxrowlen * maxrowlen));
                SPAI_API_CALL(cudaMalloc(&Jks,
                    sizeof(*Jks) * this_chunk_size * maxrowlen * maxrowlen));
                SPAI_API_CALL(cudaMalloc(&lenJks,
                    sizeof(*lenJks) * this_chunk_size));

                // fill the floating point memory with zeros
                SPAI_API_CALL(cudaMemset(Adags, 0,
                    sizeof(*Adags) * this_chunk_size * maxrowlen * maxrowlen * maxrowlen));
                SPAI_API_CALL(cudaMemset(Erhs, 0,
                    sizeof(*Erhs) * this_chunk_size * maxrowlen * maxrowlen));

                // construct small matrices A^dagger and pad with unit matrices to the same size
                cuda_spai_construct_adags_kernel <<<rwgrid, rwblock>>>(
                    this_chunk_size, first_row_in_chunk, maxrowlen,
                    A_rowptr, A_colind, A_val, rowlens,
                    Adags, Jks, lenJks);
                SPAI_CHECK_KERNEL_ERROR();
                cuda_spai_pad_adags_kernel <<<rwgrid, rwblock>>>(
                    this_chunk_size, first_row_in_chunk, maxrowlen,
                    lenJks, rowlens, Adags);
                SPAI_CHECK_KERNEL_ERROR();

                // create rhs of the minimization problems
                cuda_spai_fill_rhs <<<rwgrid, rwblock>>> (
                    this_chunk_size, first_row_in_chunk, maxrowlen,
                    Erhs, Jks, lenJks);
                SPAI_CHECK_KERNEL_ERROR();

                // create arrays of pointers to the Adags / rhss
                DT_ **Aptrs, **Cptrs;
                SPAI_API_CALL(cudaMalloc(&Aptrs, sizeof(*Aptrs) * this_chunk_size));
                SPAI_API_CALL(cudaMalloc(&Cptrs, sizeof(*Aptrs) * this_chunk_size));
                cuda_spai_fill_batch_pointers <<<rwgrid, rwblock>>> (
                    this_chunk_size, maxrowlen, Aptrs, Cptrs, Adags, Erhs);
                SPAI_CHECK_KERNEL_ERROR();

                // allocate memory for the status vector
                int *devInfoArray;
                SPAI_API_CALL(cudaMalloc(&devInfoArray, sizeof(*devInfoArray)*this_chunk_size));

                // perform the actual batched minimization
                int info;
                cuda_spai_gelsBatched_wrapper(
                    Util::Intern::cublas_handle, CUBLAS_OP_N,
                    maxrowlen*maxrowlen, maxrowlen, 1,
                    Aptrs, maxrowlen*maxrowlen, Cptrs, maxrowlen*maxrowlen,
                    &info, devInfoArray, this_chunk_size);

                // check status flags of batched minimization
                if (info != 0)
                {
                    throw InternalError(__func__, __FILE__, __LINE__,
                        "CUBLAS gelsBatched error: invalid argument no. " + stringify(-info));
                }
                std::vector<int> infoArray(this_chunk_size);
                SPAI_API_CALL(cudaMemcpy(infoArray.data(), devInfoArray, this_chunk_size*sizeof(int), cudaMemcpyDeviceToHost));
                for (unsigned int i = 0; i < this_chunk_size; i++)
                {
                    if(infoArray[i] != 0)
                        throw InternalError(__func__, __FILE__, __LINE__,
                            "CUBLAS gelsBatched error: problem no. " + stringify(i) + " could not be solved.");
                }

                // copy results to M
                cuda_spai_copy_to_m <<<rwgrid, rwblock >>> (
                    this_chunk_size, first_row_in_chunk, maxrowlen, M_val, A_rowptr, rowlens, Erhs);
                SPAI_CHECK_KERNEL_ERROR();

                // free GPU memory
                SPAI_API_CALL(cudaFree(devInfoArray));
                SPAI_API_CALL(cudaFree(Cptrs));
                SPAI_API_CALL(cudaFree(Aptrs));
                SPAI_API_CALL(cudaFree(lenJks));
                SPAI_API_CALL(cudaFree(Jks));
                SPAI_API_CALL(cudaFree(Erhs));
                SPAI_API_CALL(cudaFree(Adags));
            }

            // export explicit instantiations for float and double
            template void cuda_spai_construct_rows<double>(const unsigned int this_chunk_size, const unsigned int first_row_in_chunk,
                const double *A_val, const unsigned int *A_colind, const unsigned int *A_rowptr, double *M_val,
                const unsigned int *rowlens, const unsigned int maxrowlen);

            template void cuda_spai_construct_rows<float>(const unsigned int this_chunk_size, const unsigned int first_row_in_chunk,
                const float *A_val, const unsigned int *A_colind, const unsigned int *A_rowptr, float *M_val,
                const unsigned int *rowlens, const unsigned int maxrowlen);

        } // namespace Intern
    } // namespace Solver
} // namespace FEAST
