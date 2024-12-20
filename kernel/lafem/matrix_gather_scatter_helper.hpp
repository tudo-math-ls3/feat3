// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/tiny_algebra.hpp>

#ifdef __CUDACC__
#include <cuda/std/type_traits>
#endif



namespace FEAT
{

  namespace Intern
  {
    /**
     * \brief Policy Enum for Gather and Scatter operations
     *
     * \TODO: Documentation
     *
     * \author Maximilian Esser
     */
    enum MatrixGatherScatterPolicy
    {
      useLocalOps = 0,
      useLocalSortHelper = 1,
      useColPtr = 2
    };
  }

  namespace LAFEM
  {
    /**
      * \brief Standalone Matrix Gather and Scatter Axpy Interface
      *
      * This class provides *standalone* implementations for Gather- and Scatter-Axpy
      * operations for different underlying matrix containers.
      * You have to extract the data arrays and choose the correct functions yourself
      * (for now), so be sure these are correct.
      * The functions are also provided as (cuda) device functions, if compiled by
      * a nvcc compile unit.
      *
      * \tparam Space_ The underlying space.
      * \tparam DT_ The datatype to be used.
      * \tparam IT_ The indextype to be used.
      * \tparam policy Matrix Gather Scatter policy. For details see Intern::MatrixGatherScatterPolicy
      *
      * \author Maximilian Esser
      */
    template<typename Space_, typename  DT_, typename IT_, FEAT::Intern::MatrixGatherScatterPolicy policy_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
    struct MatrixGatherScatterHelper DOXY({});

    template<typename Space_, typename  DT_, typename IT_>
    struct MatrixGatherScatterHelper<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
    {
      /// The spacetype
      typedef Space_ SpaceType;
      /// The datatype
      typedef DT_ DataType;
      /// The indextype
      typedef IT_ IndexType;


      /**
       * \brief CSR scatter axpy function
       *
       * \tparam numr_ The number of rows of the local matrix.
       * \tparam numc_ The number of columns of the local matrix.
       *
       * \param[in] loc_mat The local matrix which values are to be scattered.
       * \param[out] matrix_data A pointer to the csr data array to be scattered in. Has to be initialized to the correct size beforehand. Previous values are added on.
       * \param[in] row_map A pointer to the dof row mapping array.
       * \param[in] col_map A pointer to the dof col mapping array.
       * \param[in] matrix_num_rows Number of rows of the csr matrix.
       * \param[in] matrix_num_cols Number of columns of the csr matrix.
       * \param[in] matrix_row_ptr CSR row pointer.
       * \param[in] matrix_col_idx CSR column index ptr.
       * \param[in] alpha Optional scaling parameter.
       * \param[in] dummy_ptr Dummy ptr used for other policies.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Matrix .
       */
      template<typename InnerType_, int numr_, int numc_ = numr_>
      CUDA_HOST_DEVICE static void scatter_matrix_csr(const Tiny::Matrix<InnerType_, numr_, numc_>& loc_mat, InnerType_* matrix_data, const IndexType* row_map, const IndexType* col_map,
                                      [[maybe_unused]] Index matrix_num_rows, [[maybe_unused]] Index matrix_num_cols, const IndexType* matrix_row_ptr,
                                      const IndexType* matrix_col_idx, DataType alpha = DataType(1), [[maybe_unused]] IndexType* dummy_ptr = nullptr)
      {
        #ifndef __CUDACC__
        static_assert(std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #else
        static_assert(::cuda::std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #endif
        IndexType loc_idx_map[numc_];
        for(int i = 0; i < numr_; ++i)
        {
          const Index ix = row_map[i];
          for(IndexType k = matrix_row_ptr[ix]; k < matrix_row_ptr[ix+1]; ++k)
          {
            for(int k_ptr = 0; k_ptr < numc_; ++k_ptr)
            {
              if(matrix_col_idx[k] == col_map[k_ptr])
              {
                loc_idx_map[k_ptr] = k;
                break;
              }
            }
          }

          //now loop over all local columns
          for(int j = 0; j < numc_; ++j)
          {
            Tiny::axpy(matrix_data[loc_idx_map[j]], loc_mat[i][j], alpha);
          }
        }
      }

      #if defined(__CUDACC__) || defined(DOXYGEN)
      /**
       * \brief CSR grouped scatter axpy function
       *
       * \tparam ThreadGroup_ The thread group type, that work together.
       * \tparam numr_, numc_ The row and column size of the internal data, i.e. each object map is mapping to consists of numr_*numc_ entries.
       *
       * \param[in] tg A reference to the thread group object on which to cooperate.
       * \param[in] scatter_size The number of blocked size elements to be scattered
       * \param[in] scatter_offset The offset into the global matrix in which to scatter
       * \param[in] loc_mat The local matrix which values are to be scattered.
       * \param[out] matrix_data A pointer to the csr data array to be scattered in. Has to be initialized to the correct size beforehand. Previous values are added on.
       * \param[in] row_map A pointer to the dof row mapping array.
       * \param[in] col_map A pointer to the dof col mapping array.
       * \param[in] matrix_num_rows Number of rows of the csr matrix.
       * \param[in] matrix_num_cols Number of columns of the csr matrix.
       * \param[in] matrix_row_ptr CSR row pointer.
       * \param[in] matrix_col_idx CSR column index ptr.
       * \param[in] num_data_row Number of InnerTypes to be scattered per local row.
       * \param[in] num_data_col Number of InnerTypes to be scattered per local col.
       * \param[in] alpha Optional scaling parameter.
       * \param[in] dummy_ptr Dummy ptr used for other policies.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Matrix .
       */
      template<typename ThreadGroup_, int numr_, int numc_=numr_>
      CUDA_DEVICE __forceinline__ static void grouped_scatter_matrix_csr(const ThreadGroup_& tg, const int scatter_size, const int scatter_offset, const DataType* loc_mat, DataType* matrix_data, const IndexType* row_map, const IndexType* col_map,
                                      [[maybe_unused]] Index matrix_num_rows, [[maybe_unused]] Index matrix_num_cols, const IndexType* matrix_row_ptr,
                                      const IndexType* matrix_col_idx, int num_data_row, int num_data_col, DataType alpha = DataType(1), [[maybe_unused]] IndexType* dummy_ptr = nullptr)
      {
        for(int idx = tg.thread_rank(); (idx < scatter_size*numr_*numc_)/* && ((idx + scatter_offset*numr_*numc_) < num_data_row*num_data_col*numr_*numc_)*/; idx += tg.num_threads())
        {
          IndexType loc_idx_map = ~IndexType(0);
          const int i = ((idx/(numr_*numc_)+scatter_offset))/num_data_row;
          const int j = ((idx/(numr_*numc_)+scatter_offset))%num_data_row;
          const Index ix = row_map[i];
          // brute force search for the correct value
          for(IndexType k = matrix_row_ptr[ix]; k < matrix_row_ptr[ix+1]; ++k)
          {
            loc_idx_map = matrix_col_idx[k] == col_map[j] ? k : loc_idx_map;
          }

          // ASSERT(loc_idx_map != ~IndexType(0));


          // and now add our value strided to the correct length
          matrix_data[loc_idx_map * numr_ * numc_ + idx%(numr_*numc_)] += alpha * loc_mat[idx];
        }
      }
      #endif

      /**
       * \brief CSR gather axpy function
       *
       * \tparam numr_ The number of rows of the local matrix.
       * \tparam numc_ The number of columns of the local matrix.
       *
       * \param[out] loc_mat The local matrix in which values are to be gathered in. Values are added on.
       * \param[in] matrix_data A pointer to the csr data array to be gathered from.
       * \param[in] row_map A pointer to the dof row mapping array.
       * \param[in] col_map A pointer to the dof col mapping array.
       * \param[in] matrix_num_rows Number of rows of the csr matrix.
       * \param[in] matrix_num_cols Number of columns of the csr matrix.
       * \param[in] matrix_row_ptr CSR row pointer.
       * \param[in] matrix_col_idx CSR column index ptr.
       * \param[in] alpha Optional scaling parameter.
       * \param[in] dummy_ptr Dummy ptr used for other policies.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Matrix .
       */
      template<typename InnerType_, int numr_, int numc_ = numr_>
      CUDA_HOST_DEVICE static void gather_matrix_csr(Tiny::Matrix<InnerType_, numr_, numc_>& loc_mat, const InnerType_* matrix_data, const IndexType* row_map, const IndexType* col_map,
                                      [[maybe_unused]] Index matrix_num_rows, [[maybe_unused]] Index matrix_num_cols, const IndexType* matrix_row_ptr,
                                      const IndexType* matrix_col_idx, DataType alpha = DataType(1), [[maybe_unused]] const IndexType* dummy_ptr = nullptr)
      {
        #ifndef __CUDACC__
        static_assert(std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #else
        static_assert(::cuda::std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #endif
        IndexType loc_idx_map[numc_];
        // loop over all local row entries
        for(int i(0); i < numr_; ++i)
        {
          // fetch row index
          const Index ix = row_map[i];

          // build column pointer for this row entry contribution
          for(IndexType k = matrix_row_ptr[ix]; k < matrix_row_ptr[ix + 1]; ++k)
          {
            for(int k_ptr = 0; k_ptr < numc_; ++k_ptr)
            {
              if(matrix_col_idx[k] == col_map[k_ptr])
              {
                loc_idx_map[k_ptr] = k;
                break;
              }
            }
          }

          // loop over all local column entries
          for(int j(0); j < numc_; ++j)
          {
            Tiny::axpy(loc_mat[i][j], matrix_data[loc_idx_map[j]], alpha);
          }
          // continue with next row entry
        }
      }

      #if defined(__CUDACC__) || defined(DOXYGEN)
      /**
       * \brief CSR grouped gather axpy function
       *
       * \tparam ThreadGroup_ The thread group type, that work together.
       * \tparam numr_, numc_ The row and column size of the internal data, i.e. each object map is mapping to consists of numr_*numc_ entries.
       *
       * \param[in] tg A reference to the thread group object on which to cooperate.
       * \param[in] scatter_size The number of blocked size elements to be scattered
       * \param[in] scatter_offset The offset into the global matrix in which to scatter
       * \param[out] loc_mat The local matrix which values are to be gathered in, has to be correctly initilized.
       * \param[in] matrix_data A pointer to the csr data array to be gathered from.
       * \param[in] row_map A pointer to the dof row mapping array.
       * \param[in] col_map A pointer to the dof col mapping array.
       * \param[in] matrix_num_rows Number of rows of the csr matrix.
       * \param[in] matrix_num_cols Number of columns of the csr matrix.
       * \param[in] matrix_row_ptr CSR row pointer.
       * \param[in] matrix_col_idx CSR column index ptr.
       * \param[in] num_data_row Number of InnerTypes to be scattered per local row.
       * \param[in] num_data_col Number of InnerTypes to be scattered per local col.
       * \param[in] alpha Optional scaling parameter.
       * \param[in] dummy_ptr Dummy ptr used for other policies.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Matrix .
       */
      template<typename ThreadGroup_, int numr_, int numc_=numr_>
      CUDA_HOST_DEVICE static void grouped_gather_matrix_csr(const ThreadGroup_& tg, const int scatter_size, const int scatter_offset, DataType* loc_mat, const DataType* matrix_data, const IndexType* row_map, const IndexType* col_map,
                                      [[maybe_unused]] Index matrix_num_rows, [[maybe_unused]] Index matrix_num_cols, const IndexType* matrix_row_ptr,
                                      const IndexType* matrix_col_idx, int num_data_row, int num_data_col, DataType alpha = DataType(1), [[maybe_unused]] IndexType* dummy_ptr = nullptr)
      {
        for(int idx = tg.thread_rank(); (idx < scatter_size*numr_*numc_) && ((idx + scatter_offset*numr_*numc_) < num_data_row*num_data_col*numr_*numc_); idx += tg.num_threads())
        {
          IndexType loc_idx_map = ~IndexType(0);
          const int i = ((idx/(numr_*numc_)+scatter_offset))/num_data_row;
          const int j = ((idx/(numr_*numc_)+scatter_offset))%num_data_row;
          const Index ix = row_map[i];
          // brute force search for the correct value
          for(IndexType k = matrix_row_ptr[ix]; k < matrix_row_ptr[ix+1]; ++k)
          {
            loc_idx_map = matrix_col_idx[k] == col_map[j] ? k : loc_idx_map;
          }

          // and now add our value strided to the correct length
          loc_mat[idx] += alpha * matrix_data[loc_idx_map * numr_ * numc_ + idx%(numr_*numc_)];
        }
      }
      #endif
    }; // struct MatrixGatherScatterHelper<localOps>

   template<typename Space_, typename  DT_, typename IT_>
   struct MatrixGatherScatterHelper<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>
   {
      typedef DT_ DataType;
      typedef IT_ IndexType;
      typedef Space_ SpaceType;


      /**
       * \brief CSR scatter axpy function
       *
       * \tparam numr_ The number of rows of the local matrix.
       * \tparam numc_ The number of columns of the local matrix.
       *
       * \param[in] loc_mat The local matrix which values are to be scattered.
       * \param[out] matrix_data A pointer to the csr data array to be scattered in. Has to be initialized to the correct size beforehand. Previous values are added on.
       * \param[in] row_map A pointer to the dof row mapping array.
       * \param[in] col_map A pointer to the dof col mapping array.
       * \param[in] matrix_num_rows Number of rows of the csr matrix.
       * \param[in] matrix_num_cols Number of columns of the csr matrix.
       * \param[in] matrix_row_ptr CSR row pointer.
       * \param[in] matrix_col_idx CSR column index ptr.
       * \param[in] alpha Optional scaling parameter.
       * \param[in] col_map_sorter Helper array with locally sorted matrix mapping.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Matrix .
       */
      template<typename InnerType_, int numr_, int numc_ = numr_>
      CUDA_HOST_DEVICE static void scatter_matrix_csr(const Tiny::Matrix<InnerType_, numr_, numc_>& loc_mat, InnerType_* matrix_data, const IndexType* row_map, const IndexType* col_map,
                                      [[maybe_unused]] Index matrix_num_rows, [[maybe_unused]] Index matrix_num_cols, const IndexType* matrix_row_ptr,
                                      const IndexType* matrix_col_idx, DataType alpha, const IndexType* col_map_sorter)
      {
        #ifndef __CUDACC__
        static_assert(std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #else
        static_assert(::cuda::std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #endif
        IndexType loc_idx_map[numc_];

        // loop over all local row entries
        for(int i(0); i < numr_; ++i)
        {
          // fetch row index
          Index k = matrix_row_ptr[row_map[i]];
          for(Index k_ptr = 0; k_ptr < numc_; ++k_ptr)
          {
            const Index real_dof = col_map_sorter[k_ptr];
            //search for our column value, no boundary checks, so be damn sure the value is inside matrix_col_idx
            while(matrix_col_idx[k] < col_map[real_dof])
            {
              ++k;
            }
            loc_idx_map[real_dof] = IndexType(k++);
          }

          // loop over all local column entries
          for(int j(0); j < numc_; ++j)
          {
            Tiny::axpy(matrix_data[loc_idx_map[j]], loc_mat[i][j], alpha);
          }
          // continue with next row entry
        }
      }

      #if defined(__CUDACC__) || defined(DOXYGEN)
      /**
       * \brief CSR grouped scatter axpy function
       *        Does not use the local_sorter array, since useless...
       *
       * \tparam ThreadGroup_ The thread group type, that work together.
       * \tparam numr_, numc_ The row and column size of the internal data, i.e. each object map is mapping to consists of numr_*numc_ entries.
       *
       * \param[in] tg A reference to the thread group object on which to cooperate.
       * \param[in] scatter_size The number of blocked size elements to be scattered
       * \param[in] scatter_offset The offset into the global matrix in which to scatter
       * \param[in] loc_mat The local matrix which values are to be scattered.
       * \param[out] matrix_data A pointer to the csr data array to be scattered in. Has to be initialized to the correct size beforehand. Previous values are added on.
       * \param[in] row_map A pointer to the dof row mapping array.
       * \param[in] col_map A pointer to the dof col mapping array.
       * \param[in] matrix_num_rows Number of rows of the csr matrix.
       * \param[in] matrix_num_cols Number of columns of the csr matrix.
       * \param[in] matrix_row_ptr CSR row pointer.
       * \param[in] matrix_col_idx CSR column index ptr.
       * \param[in] num_data_row Number of InnerTypes to be scattered per local row.
       * \param[in] num_data_col Number of InnerTypes to be scattered per local col.
       * \param[in] alpha Optional scaling parameter.
       * \param[in] dummy_ptr Dummy ptr used for other policies.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Matrix .
       */
      template<typename ThreadGroup_, int numr_, int numc_=numr_>
      CUDA_HOST_DEVICE static __forceinline__ void grouped_scatter_matrix_csr(const ThreadGroup_& tg, const int scatter_size, const int scatter_offset, const DataType* loc_mat, DataType* matrix_data, const IndexType* row_map, const IndexType* col_map,
                                      [[maybe_unused]] Index matrix_num_rows, [[maybe_unused]] Index matrix_num_cols, const IndexType* matrix_row_ptr,
                                      const IndexType* matrix_col_idx, int num_data_row, int num_data_col, DataType alpha = DataType(1), [[maybe_unused]] IndexType* dummy_ptr = nullptr)
      {
        MatrixGatherScatterHelper<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>::template grouped_scatter_matrix_csr<ThreadGroup_, numr_, numc_>(tg, scatter_size, scatter_offset, loc_mat,
              matrix_data, row_map, col_map, matrix_num_rows, matrix_num_cols, matrix_row_ptr, matrix_col_idx, num_data_row, num_data_col, alpha, nullptr);
      }
      #endif

      /**
       * \brief CSR gather axpy function
       *
       * \tparam numr_ The number of rows of the local matrix.
       * \tparam numc_ The number of columns of the local matrix.
       *
       * \param[out] loc_mat The local matrix in which values are to be gathered in. Values are added on.
       * \param[in] matrix_data A pointer to the csr data array to be gathered from.
       * \param[in] row_map A pointer to the dof row mapping array.
       * \param[in] col_map A pointer to the dof col mapping array.
       * \param[in] matrix_num_rows Number of rows of the csr matrix.
       * \param[in] matrix_num_cols Number of columns of the csr matrix.
       * \param[in] matrix_row_ptr CSR row pointer.
       * \param[in] matrix_col_idx CSR column index ptr.
       * \param[in] alpha Optional scaling parameter.
       * \param[in] dummy_ptr Dummy ptr used for other policies.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Matrix .
       */
      template<typename InnerType_, int numr_, int numc_ = numr_>
      CUDA_HOST_DEVICE static void gather_matrix_csr(Tiny::Matrix<InnerType_, numr_, numc_>& loc_mat, const InnerType_* matrix_data, const IndexType* row_map, const IndexType* col_map,
                                      [[maybe_unused]] Index matrix_num_rows, [[maybe_unused]] Index matrix_num_cols, const IndexType* matrix_row_ptr,
                                      const IndexType* matrix_col_idx, DataType alpha, const IndexType* col_map_sorter)
      {
        #ifndef __CUDACC__
        static_assert(std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #else
        static_assert(::cuda::std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #endif

        IndexType loc_idx_map[numc_];

        // loop over all local row entries
        for(int i(0); i < numr_; ++i)
        {
          // fetch row index
          Index k = matrix_row_ptr[row_map[i]];
          for(Index k_ptr = 0; k_ptr < numc_; ++k_ptr)
          {
            const Index real_dof = col_map_sorter[k_ptr];
            //search for our column value, no boundary checks, so be damn sure the value is inside matrix_col_idx
            while(matrix_col_idx[k] < col_map[real_dof])
            {
              ++k;
            }
            loc_idx_map[real_dof] = IndexType(k++);
          }

          // loop over all local column entries
          for(int j(0); j < numc_; ++j)
          {
            Tiny::axpy(loc_mat[i][j], matrix_data[loc_idx_map[j]], alpha);
          }
          // continue with next row entry
        }
      }

      #if defined(__CUDACC__) || defined(DOXYGEN)
      /**
       * \brief CSR grouped gather axpy function
       *
       * \tparam ThreadGroup_ The thread group type, that work together.
       * \tparam numr_, numc_ The row and column size of the internal data, i.e. each object map is mapping to consists of numr_*numc_ entries.
       *
       * \param[in] tg A reference to the thread group object on which to cooperate.
       * \param[in] scatter_size The number of blocked size elements to be scattered
       * \param[in] scatter_offset The offset into the global matrix in which to scatter
       * \param[out] loc_mat The local matrix which values are to be gathered in, has to be correctly initilized.
       * \param[in] matrix_data A pointer to the csr data array to be gathered from.
       * \param[in] row_map A pointer to the dof row mapping array.
       * \param[in] col_map A pointer to the dof col mapping array.
       * \param[in] matrix_num_rows Number of rows of the csr matrix.
       * \param[in] matrix_num_cols Number of columns of the csr matrix.
       * \param[in] matrix_row_ptr CSR row pointer.
       * \param[in] matrix_col_idx CSR column index ptr.
       * \param[in] num_data_row Number of InnerTypes to be scattered per local row.
       * \param[in] num_data_col Number of InnerTypes to be scattered per local col.
       * \param[in] alpha Optional scaling parameter.
       * \param[in] dummy_ptr Dummy ptr used for other policies.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Matrix .
       */
      template<typename ThreadGroup_, int numr_, int numc_=numr_>
      CUDA_HOST_DEVICE static void grouped_gather_matrix_csr(const ThreadGroup_& tg, const int scatter_size, const int scatter_offset, DataType* loc_mat, const DataType* matrix_data, const IndexType* row_map, const IndexType* col_map,
                                      [[maybe_unused]] Index matrix_num_rows, [[maybe_unused]] Index matrix_num_cols, const IndexType* matrix_row_ptr,
                                      const IndexType* matrix_col_idx, int num_data_row, int num_data_col, DataType alpha = DataType(1), [[maybe_unused]] IndexType* dummy_ptr = nullptr)
      {
        MatrixGatherScatterHelper<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>::template grouped_gather_matrix_csr<ThreadGroup_, numr_, numc_>(tg, scatter_size, scatter_offset, loc_mat,
              matrix_data, row_map, col_map, matrix_num_rows, matrix_num_cols, matrix_row_ptr, matrix_col_idx, num_data_row, num_data_col, alpha, nullptr);
      }
      #endif
   }; // struct MatrixGatherScatterHelper
  }
}
