#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/solver/amavanka_base.hpp>
#include <kernel/solver/voxel_amavanka.hpp>

#include "cusparse_v2.h"
#include <cuda/std/utility>

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      // Wrapper arrays
      template<int n_>
      struct MacroDofs
      {
        Index* macro_dofs[n_];
      };

      template<int n_>
      struct MDegreeDofs
      {
        Index max_degree_dofs[n_];
      };

      template<int n_>
      struct DofMacros
      {
        Index* dof_macros[n_];
      };

      template<int n_>
      struct MDegreeMacros
      {
        Index max_degree_macros[n_];
      };
    }
    namespace Kernel
    {
      template<typename DT_, typename IT_, int n_>
      __global__ void gather_local_matrices(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_> mat_wrap,
        Intern::MacroDofs<n_> macro_dofs_wrapper, Intern::MDegreeDofs<n_> max_degree_wrapper,
        Index num_macros, Index stride, DT_* __restrict__ _local)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if(idx >= num_macros)
          return;
        typedef DT_ DataType;
        DT_* local = _local + idx*stride*stride;
        // note: actually, blas routines work in column-major storage, but with quadratic matrices and
        // as long as gather and scatter both work on the same storage format (i.e. row-major or column-major), we can simply gather in row-major format
        Intern::VoxelAmaVankaCore::template gather<DT_, IT_, n_, true>(mat_wrap, local, stride, idx, macro_dofs_wrapper.macro_dofs,
                                                                      max_degree_wrapper.max_degree_dofs, Index(0), Index(0),
                                                                      Index(0), Index(0));
      }

      template<typename DT_, typename IT_, int n_, bool skip_singular_>
      __global__ void scatter_local_matrices(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_> vanka_wrap,
        Intern::MacroDofs<n_> macro_dofs_wrapper, Intern::MDegreeDofs<n_> max_degree_wrapper,
        Index stride, const DT_* __restrict__ _local, int* __restrict__ macro_mask,
        const int* __restrict__ coloring_map, Index color_size, int* info_array)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if(idx >= color_size)
          return;
        typedef DT_ DataType;
        const int imacro = coloring_map[idx];
        const DT_* local = _local + imacro*stride*stride;
        // typedef IT_ IndexType;
        if constexpr(skip_singular_)
        {
          const bool singular = info_array[imacro] > 0;
          // set macro regularity mask
          macro_mask[imacro] = (singular ? 0 : 1);

          // scatter local matrix
          if(!singular)
          {
            Intern::VoxelAmaVankaCore::template scatter_add<DT_, IT_, n_, true>(vanka_wrap, local, stride, imacro, macro_dofs_wrapper.macro_dofs,
              max_degree_wrapper.max_degree_dofs, Index(0), Index(0), Index(0), Index(0));
          }
        }
        else
        {
          Intern::VoxelAmaVankaCore::template scatter_add<DT_, IT_, n_, true>(vanka_wrap, local, stride, imacro, macro_dofs_wrapper.macro_dofs,
            max_degree_wrapper.max_degree_dofs, Index(0), Index(0), Index(0), Index(0));
        }
      }

      template<typename DT_>
      __global__ void set_batched_matrix_ptrs(DT_** local_ptrs, DT_* local, Index num_macros, Index stride)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if(idx >= num_macros)
          return;
        local_ptrs[idx] = local + idx*stride*stride;
      }

      template<typename DT_>
      void cublas_getrfBatched(DT_** __restrict__ a_array, int* __restrict__ pivot, int* __restrict__ info_array, Index batch_size, Index leading_dimension, Index n);

      template<>
      void cublas_getrfBatched(double** __restrict__ a_array, int* __restrict__ pivot, int* __restrict__ info_array, Index batch_size, Index leading_dimension, Index n)
      {
        cublasStatus_t status =  cublasDgetrfBatched(Util::Intern::cublas_handle, int(n), a_array, int(leading_dimension), pivot, info_array, int(batch_size));
        if(status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasDgetrfBatched failed with status code: " + stringify(status));
      }

      template<>
      void cublas_getrfBatched(float** __restrict__ a_array, int* __restrict__ pivot, int* __restrict__ info_array, Index batch_size, Index leading_dimension, Index n)
      {
        cublasStatus_t status =  cublasSgetrfBatched(Util::Intern::cublas_handle, int(n), a_array, int(leading_dimension), pivot, info_array, int(batch_size));
        if(status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasSgetrfBatched failed with status code: " + stringify(status));
      }

      template<typename DT_>
      void cublas_getriBatched(DT_** __restrict__ a_array, DT_** __restrict__ c_array, int* __restrict__ pivot, int* __restrict__ info_array, Index batch_size, Index leading_dimension_a, Index n, Index leading_dimension_c);

      template<>
      void cublas_getriBatched(double** __restrict__ a_array, double** __restrict__ c_array, int* __restrict__ pivot, int* __restrict__ info_array, Index batch_size, Index leading_dimension_a, Index n, Index leading_dimension_c)
      {
        cublasStatus_t status =  cublasDgetriBatched(Util::Intern::cublas_handle, int(n), a_array, int(leading_dimension_a), pivot, c_array, int(leading_dimension_c), info_array, int(batch_size));
        if(status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasDgetriBatched failed with status code: " + stringify(status));
      }

      template<>
      void cublas_getriBatched(float** __restrict__ a_array, float** __restrict__ c_array, int* __restrict__ pivot, int* __restrict__ info_array, Index batch_size, Index leading_dimension_a, Index n, Index leading_dimension_c)
      {
        cublasStatus_t status =  cublasSgetriBatched(Util::Intern::cublas_handle, int(n), a_array, int(leading_dimension_a), pivot, c_array, int(leading_dimension_c), info_array, int(batch_size));
        if(status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasSgetriBatched failed with status code: " + stringify(status));
      }

      template<typename DT_>
      void batch_invert_matrices(dim3 grid, dim3 block, Index num_macros, Index stride, Index uniform_mat_size, DT_* __restrict__ local, DT_* __restrict__ local_t, int* __restrict__ pivot, int* __restrict__ info_array)
      {
        DT_** local_ptrs = (DT_**)Util::cuda_malloc(sizeof(DT_*)*num_macros);
        DT_** local_t_ptrs = (DT_**)Util::cuda_malloc(sizeof(DT_*)*num_macros);

        //set matrix ptr arrays
        Kernel::template set_batched_matrix_ptrs<DT_><<<grid, block>>>(local_ptrs, local, num_macros, stride);
        Kernel::template set_batched_matrix_ptrs<DT_><<<grid, block>>>(local_t_ptrs, local_t, num_macros, stride);

        Util::cuda_check_last_error();

        //and now call cublas kernels to invert our matrices
        cublas_getrfBatched(local_ptrs, pivot, info_array, num_macros, uniform_mat_size, uniform_mat_size);
        cublas_getriBatched(local_ptrs, local_t_ptrs, pivot, info_array, num_macros, uniform_mat_size, uniform_mat_size, uniform_mat_size);

        Util::cuda_free((void*)local_t_ptrs);
        Util::cuda_free((void*)local_ptrs);
      }


      template<typename DT_, typename IT_, int n_, bool skip_singular_, bool uniform_macros_>
      __global__ void assemble_unscaled_vanka_device(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_> mat_wrap,
        Intern::CSRTupleMatrixWrapper<DT_, IT_, n_> vanka_wrap,
        Intern::MacroDofs<n_> macro_dofs_wrapper, int* __restrict__ macro_mask,
        Intern::MDegreeDofs<n_> max_degree_wrapper, const int* __restrict__ coloring_map, Index color_size,
        Index num_macros, Index stride, DT_ eps, DT_* __restrict__ _local, DT_* __restrict__ _local_t, Index* __restrict__ _pivot)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if(idx >= color_size)
          return;
        typedef DT_ DataType;
        // typedef IT_ IndexType;
        if constexpr(skip_singular_)
        {
          const int imacro = coloring_map[idx];
          DT_* local = _local + idx*stride*stride;
          DT_* local_t = _local_t + idx*stride*stride;
          Index* pivot = _pivot + idx*stride;
          Index** ptr_dofs = macro_dofs_wrapper.macro_dofs;
          Index* ptr_deg = max_degree_wrapper.max_degree_dofs;
          const cuda::std::pair<Index, Index> nrc = Intern::VoxelAmaVankaCore::template gather<DT_, IT_, n_, uniform_macros_>(mat_wrap, local, stride, Index(imacro), ptr_dofs,
                                                                                                  ptr_deg, Index(0), Index(0),
                                                                                                 Index(0), Index(0));
          // the approach used for checking the regularity of the local matrix is to check whether
          //
          //     || I - A*A^{-1} ||_F^2 < eps
          //
          // we could try to analyse the pivots returned by invert_matrix function instead, but
          // unfortunately this approach sometimes leads to false positives

          // make a backup if checking for singularity
          for(Index i(0); i < nrc.first; ++i)
            for(Index j(0); j < nrc.second; ++j)
              local_t[i*stride+j] = local[i*stride+j];

          // invert matrix
          CudaMath::invert_matrix(nrc.first, stride, local, pivot);

          // compute (squared) Frobenius norm of (I - A*A^{-1})
          DataType norm = DataType(0);
          for(Index i(0); i < nrc.first; ++i)
          {
            for(Index j(0); j < nrc.first; ++j)
            {
              DataType xij = DataType(i == j ? 1 : 0);
              for(Index l(0); l < nrc.first; ++l)
                xij -= local_t[i*stride+l] * local[l*stride+j]; // A_il * (A^{-1})_lj
              norm += xij * xij;
            }
          }

          // is the matrix block singular?
          // Note: we check for !(norm < eps) instead of (norm >= eps),
          // because the latter one evaluates to false if norm is NaN,
          // which would result in a false negative
          const bool singular = !(norm < eps);

          // set macro regularity mask
          macro_mask[imacro] = (singular ? 0 : 1);

          // scatter local matrix
          if(!singular)
          {
            Intern::VoxelAmaVankaCore::template scatter_add<DT_, IT_, n_, uniform_macros_>(vanka_wrap, local, stride, imacro, macro_dofs_wrapper.macro_dofs,
              max_degree_wrapper.max_degree_dofs, Index(0), Index(0), Index(0), Index(0));
          }

          // reformat local matrix
          for(Index i(0); i < nrc.first; ++i)
            for(Index j(0); j < nrc.second; ++j)
              local[i*stride+j] = DataType(0);
        }
        else
        {
          const int imacro = coloring_map[idx];
          DT_* local = _local + idx*stride*stride;
          Index* pivot = _pivot + idx*stride;
          const cuda::std::pair<Index, Index> nrc = Intern::VoxelAmaVankaCore::template gather<DT_, IT_, n_, uniform_macros_>(mat_wrap, local, stride, imacro, macro_dofs_wrapper.macro_dofs,
                                                                                                  max_degree_wrapper.max_degree_dofs, Index(0), Index(0),
                                                                                                 Index(0), Index(0));
          // invert matrix
          CudaMath::invert_matrix(nrc.first, stride, local, pivot);

          Intern::VoxelAmaVankaCore::template scatter_add<DT_, IT_, n_, uniform_macros_>(vanka_wrap, local, stride, imacro, macro_dofs_wrapper.macro_dofs,
                                                 max_degree_wrapper.max_degree_dofs, Index(0), Index(0),
                                                Index(0), Index(0));

          // reformat local matrix
          for(Index i(0); i < nrc.first; ++i)
            for(Index j(0); j < nrc.second; ++j)
              local[i*stride+j] = DataType(0);
        }
      }

      template<typename DT_, typename  IT_, int n_, bool skip_singular_>
      __global__ void scale_vanka_rows_device(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_> vanka_wrap, const DT_ omega,
                                  Intern::DofMacros<n_> wrapper_dof_macros, Intern::MDegreeMacros<n_> wrapper_max_degree_macros, const int* __restrict__ m_mask, const Index all_row_size)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if(idx >= all_row_size * n_)
          return;

        // find meta row and actual row index for the current row
        Index xrow = idx % all_row_size;
        Index meta_col = idx / all_row_size;
        Index rowc = 0;
        Index next_row = vanka_wrap.tensor_counts[0];
        while(next_row <= xrow)
        {
          xrow -= next_row;
          ++rowc;
          next_row = vanka_wrap.tensor_counts[1+rowc];
        }
        const Index max_row_degree_size = wrapper_max_degree_macros.max_degree_macros[rowc]+1;
        const Index* row_img_idx = wrapper_dof_macros.dof_macros[rowc] + xrow*max_row_degree_size;
        const Index real_degree = *row_img_idx;
        const int hw = vanka_wrap.blocksizes[Index(1) + CudaMath::cuda_min(Index(rowc),Index(1))];
        // const int num_rows = int(next_row);
        DT_* vals = vanka_wrap.data_arrays[rowc*n_ + meta_col];
        const IT_* row_ptr = vanka_wrap.row_arrays[rowc*n_+meta_col];
        const IT_* col_idx = vanka_wrap.col_arrays[rowc*n_+meta_col];
        const int hb = vanka_wrap.blocksizes[meta_col +1];
        // careful, num rows is counted against native elements, not raw elements
        Intern::VoxelAmaVankaCore::scale_row<DT_, IT_, skip_singular_>(vals, omega, row_ptr, col_idx, row_img_idx+1, real_degree, hw, hb, xrow, m_mask);

      }
    }












    namespace Arch
    {
      template<typename DT_, typename IT_, int n_>
      void assemble_vanka_device(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
        Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& vanka_wrap, const std::vector<Index*>& d_macro_dofs,
        const std::vector<Index*>& d_dof_macros, int* d_macro_mask, const std::vector<Index>& max_degree_dofs,
        const std::vector<Index>& max_degree_macros, const Adjacency::ColoringDataHandler& coloring_data,
        Index num_macros, Index stride, DT_ omega, DT_ eps, bool skip_singular, bool uniform_macros)
        {
          typedef DT_ DataType;
          const Index blocksize = Util::cuda_blocksize_vanka_assembly;
          dim3 grid;
          dim3 block;
          block.x = (unsigned int)blocksize;
          //this is probably to much data
          grid.x = (unsigned int)(ceil(double(coloring_data.get_max_size())/double(block.x)));
          //further extract internal arrays
          Intern::DofMacros<n_> graphs_dof_macros;
          Intern::MacroDofs<n_> graphs_macro_dofs;
          Intern::MDegreeDofs<n_> degree_dofs;
          Intern::MDegreeMacros<n_> degree_macros;
          for(int i = 0; i < n_; ++i)
          {
            graphs_dof_macros.dof_macros[i] = d_dof_macros.at(i);
            graphs_macro_dofs.macro_dofs[i] = d_macro_dofs.at(i);
            degree_dofs.max_degree_dofs[i] = max_degree_dofs.at(i);
            degree_macros.max_degree_macros[i] = max_degree_macros.at(i);
          }
          auto& coloring_map = coloring_data.get_coloring_maps();
          const Index max_color_size = coloring_data.get_max_size();
          // allocate arrays for local matrix
          DataType* local = (DataType*)Util::cuda_malloc(max_color_size*stride*stride*sizeof(DataType));
          Util::cuda_set_memory(local, DataType(0), max_color_size*stride*stride);
          DataType* local_t = nullptr;
          if(skip_singular)
          {
            local_t = (DataType*)Util::cuda_malloc(max_color_size*stride*stride*sizeof(DataType));
            Util::cuda_set_memory(local_t, DataType(0), max_color_size*stride*stride);
          }
          Index* pivot = (Index*)Util::cuda_malloc(max_color_size*stride*sizeof(Index));
          Util::cuda_set_memory(pivot, Index(0), max_color_size*stride);
          for(Index k = 0; k < coloring_map.size(); ++k)
          {
            if(uniform_macros)
            {
              if(skip_singular)
                Solver::Kernel::template assemble_unscaled_vanka_device<DT_, IT_, n_, true, true><<<grid, block>>>
                                                (mat_wrap, vanka_wrap, graphs_macro_dofs,
                                                d_macro_mask, degree_dofs,
                                                (int*)coloring_map[k], coloring_data.get_color_size(k), num_macros, stride, eps,
                                                local, local_t, pivot);
              else
                Solver::Kernel::template assemble_unscaled_vanka_device<DT_, IT_, n_, false, true><<<grid, block>>>
                                                (mat_wrap, vanka_wrap, graphs_macro_dofs,
                                                d_macro_mask, degree_dofs,
                                                (int*)coloring_map[k], coloring_data.get_color_size(k), num_macros, stride, eps,
                                                local, local_t, pivot);
            }
            else
            {
              if(skip_singular)
                Solver::Kernel::template assemble_unscaled_vanka_device<DT_, IT_, n_, true, false><<<grid, block>>>
                                                (mat_wrap, vanka_wrap, graphs_macro_dofs,
                                                d_macro_mask, degree_dofs,
                                                (int*)coloring_map[k], coloring_data.get_color_size(k), num_macros, stride, eps,
                                                local, local_t, pivot);
              else
                Solver::Kernel::template assemble_unscaled_vanka_device<DT_, IT_, n_, false, false><<<grid, block>>>
                                                (mat_wrap, vanka_wrap, graphs_macro_dofs,
                                                d_macro_mask, degree_dofs,
                                                (int*)coloring_map[k], coloring_data.get_color_size(k), num_macros, stride, eps,
                                                local, local_t, pivot);
            }
          }
          //check for cuda error in our kernel
          Util::cuda_check_last_error();
          Util::cuda_free((void*)pivot);
          Util::cuda_free((void*)local_t);
          Util::cuda_free((void*)local);
          // get max row_size
          const Index all_row_size = vanka_wrap.get_all_rows_size();
          block.x = (unsigned int)Util::cuda_blocksize_spmv;
          grid.x = (unsigned int)(ceil(double(all_row_size*n_)/double(block.x)));
          if(skip_singular)
            Solver::Kernel::template scale_vanka_rows_device<DT_, IT_, n_, true><<<grid, block>>>(vanka_wrap, omega, graphs_dof_macros, degree_macros, d_macro_mask, all_row_size);
          else
            Solver::Kernel::template scale_vanka_rows_device<DT_, IT_, n_, false><<<grid, block>>>(vanka_wrap, omega, graphs_dof_macros, degree_macros, d_macro_mask, all_row_size);

          //check for cuda error in our kernel
          Util::cuda_check_last_error();
        }

      // only works with uniform macros
      template<typename DT_, typename IT_, int n_>
      void assemble_vanka_device_batched(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
        Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& vanka_wrap, const std::vector<Index*>& d_macro_dofs,
        const std::vector<Index*>& d_dof_macros, int* d_macro_mask, const std::vector<Index>& max_degree_dofs,
        const std::vector<Index>& max_degree_macros, const Adjacency::ColoringDataHandler& coloring_data,
        Index num_macros, Index stride, Index actual_matrix_size, DT_ omega, bool skip_singular)
        {
          typedef DT_ DataType;
          const Index blocksize = Util::cuda_blocksize_vanka_assembly;
          dim3 grid;
          dim3 block;
          //further extract internal arrays
          Intern::DofMacros<n_> graphs_dof_macros;
          Intern::MacroDofs<n_> graphs_macro_dofs;
          Intern::MDegreeDofs<n_> degree_dofs;
          Intern::MDegreeMacros<n_> degree_macros;
          for(int i = 0; i < n_; ++i)
          {
            graphs_dof_macros.dof_macros[i] = d_dof_macros.at(i);
            graphs_macro_dofs.macro_dofs[i] = d_macro_dofs.at(i);
            degree_dofs.max_degree_dofs[i] = max_degree_dofs.at(i);
            degree_macros.max_degree_macros[i] = max_degree_macros.at(i);
          }
          auto& coloring_map = coloring_data.get_coloring_maps();
          // allocate arrays for local matrix
          DataType* local = (DataType*)Util::cuda_malloc(num_macros*stride*stride*sizeof(DataType));
          Util::cuda_set_memory(local, DataType(0), num_macros*stride*stride);
          DataType* local_t = (DataType*)Util::cuda_malloc(num_macros*stride*stride*sizeof(DataType));
          Util::cuda_set_memory(local_t, DataType(0), num_macros*stride*stride);
          int* pivot = (int*)Util::cuda_malloc(num_macros*actual_matrix_size*sizeof(int));
          Util::cuda_set_memory(pivot, int(0), num_macros*stride);
          int* info_array = (int*)Util::cuda_malloc(num_macros*sizeof(int));
          Util::cuda_set_memory(info_array, int(0), num_macros);
          // gather and invert
          block.x = (unsigned int)blocksize;
          //this is probably to much data
          grid.x = (unsigned int)(ceil(double(num_macros)/double(block.x)));
          Solver::Kernel::template gather_local_matrices<DT_, IT_, n_><<<grid, block>>>(mat_wrap, graphs_macro_dofs, degree_dofs, num_macros, stride, local);
          Solver::Kernel::batch_invert_matrices(grid, block, num_macros, stride, actual_matrix_size, local, local_t, pivot, info_array);
          // now scatter local matrices... this requires coloring
          block.x = (unsigned int)blocksize;
          //this is probably to much data
          grid.x = (unsigned int)(ceil(double(coloring_data.get_max_size())/double(block.x)));
          for(Index k = 0; k < coloring_map.size(); ++k)
          {
            if(skip_singular)
              Solver::Kernel::template scatter_local_matrices<DT_, IT_, n_, true><<<grid, block>>>
                                              (vanka_wrap, graphs_macro_dofs,
                                              degree_dofs, stride, local_t, d_macro_mask,
                                              (int*)coloring_map[k], coloring_data.get_color_size(k), info_array);
            else
              Solver::Kernel::template scatter_local_matrices<DT_, IT_, n_, false><<<grid, block>>>
                                              (vanka_wrap, graphs_macro_dofs,
                                              degree_dofs, stride, local_t, d_macro_mask,
                                              (int*)coloring_map[k], coloring_data.get_color_size(k), info_array);
          }

          //check for cuda error in our kernel
          Util::cuda_check_last_error();
          Util::cuda_free((void*)info_array);
          Util::cuda_free((void*)pivot);
          Util::cuda_free((void*)local_t);
          Util::cuda_free((void*)local);
          // get max row_size
          const Index all_row_size = vanka_wrap.get_all_rows_size();
          block.x = (unsigned int)Util::cuda_blocksize_spmv;
          grid.x = (unsigned int)(ceil(double(all_row_size*n_)/double(block.x)));
          if(skip_singular)
            Solver::Kernel::template scale_vanka_rows_device<DT_, IT_, n_, true><<<grid, block>>>(vanka_wrap, omega, graphs_dof_macros, degree_macros, d_macro_mask, all_row_size);
          else
            Solver::Kernel::template scale_vanka_rows_device<DT_, IT_, n_, false><<<grid, block>>>(vanka_wrap, omega, graphs_dof_macros, degree_macros, d_macro_mask, all_row_size);

          //check for cuda error in our kernel
          Util::cuda_check_last_error();



        }
    }
  }
}


using namespace FEAT;
using namespace FEAT::Solver;

//########################## oneThreadperMacro kernel ################################################
//#########################1x1 kernels###############################################################
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 1>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 1>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, double, double, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 1>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 1>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, float, float, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 1>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 1>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, double, double, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 1>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 1>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, float, float, bool, bool);
//#########################2x2 kernels###############################################################
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 2>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 2>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, double, double, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 2>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 2>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, float, float, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 2>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 2>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, double, double, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 2>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 2>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, float, float, bool, bool);


/** For copy pasting new nxn kernels...
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, n_>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, n_>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, double, double, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, n_>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, n_>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, float, float, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, n_>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, n_>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, double, double, bool, bool);
template void Arch::assemble_vanka_device(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, n_>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, n_>&, const std::vector<Index*>&,
        const std::vector<Index*>&, int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, float, float, bool, bool);
//*/

//######################### cuBlasbasedKernels ##############################################################
//#########################1x1 kernels###############################################################
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 1>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 1>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, double, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 1>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 1>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, float, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 1>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 1>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, double, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 1>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 1>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, float, bool);


//#########################2x2 kernels###############################################################
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 2>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 2>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, double, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 2>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 2>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, float, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 2>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 2>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, double, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 2>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 2>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, float, bool);




/** For copy pasting new nxn kernels...
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, n_>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint64_t, n_>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, double, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, n_>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint64_t, n_>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, float, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, n_>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<double, std::uint32_t, n_>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, double, bool);
template void Arch::assemble_vanka_device_batched(const FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, n_>&, FEAT::Solver::Intern::CSRTupleMatrixWrapper<float, std::uint32_t, n_>&, const std::vector<Index*>&, const std::vector<Index*>&,
        int*, const std::vector<Index>&, const std::vector<Index>&, const Adjacency::ColoringDataHandler&, Index, Index, Index, float, bool);
//*/