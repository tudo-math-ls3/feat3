#pragma once
#ifndef KERNEL_LAFEM_VECGASC_HELPER_HPP
#define KERNEL_LAFEM_VECGASC_HELPER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/tiny_algebra.hpp>

#ifdef __CUDACC__
#include <cuda/std/type_traits>
#endif

namespace FEAT
{
  namespace LAFEM
  {
    /**
      * \brief Standalone Vector Gather and Scatter Axpy Interface
      *
      * This class provides *standalone* implementations for Gather- and Scatter-Axpy
      * operations for different underlying vector containers.
      * You have to extract the data arrays and choose the correct functions yourself
      * (for now), so be sure these are correct.
      * The functions are also provided as (cuda) device functions, if compiled by
      * a nvcc compile unit.
      * Furthermore, grouped operations, see cuda cooperative groups, are also available.
      *
      * \tparam Space_ The underlying space.
      * \tparam DT_ The datatype to be used.
      * \tparam IT_ The indextype to be used.
      *
      * \author Maximilian Esser
      */
    template<typename Space_, typename DT_, typename IT_>
    struct VectorGatherScatterHelper;

    template<typename Space_, typename DT_, typename IT_>
    struct VectorGatherScatterHelper
    {
      typedef Space_ SpaceType;
      typedef DT_ DataType;
      typedef IT_ IndexType;

      /**
       * \brief Dense Vector scatter axpy function
       *
       * \tparam InnerType_ The inner dataobject type.
       * \tparam numr_ The number of rows of the local matrix.
       *
       * \param[in] loc_vec The local vector which values are to be scattered.
       * \param[out] data A pointer to the vector data array to be scattered in. Has to be initialized to the correct size beforehand. Previous values are added on.
       * \param[in] num_entries Number of entries of the dense vector.
       * \param[in] map The dof mapping used to map local to global values.
       * \param[in] alpha Optional scaling parameter.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Vector .
       */
      template<typename InnerType_, int numr_>
      CUDA_HOST_DEVICE static void scatter_vector_dense(const Tiny::Vector<InnerType_, numr_>& loc_vec, InnerType_* data, [[maybe_unused]] IndexType num_entries, const IndexType* map, DataType alpha = DataType(1))
      {
        #ifndef __CUDACC__
        static_assert(std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #else
        static_assert(::cuda::std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #endif
        // loop over all local entries
        for (int i(0); i < numr_; ++i)
        {
          // get dof index
          Index dof_idx = map[i];
          // ASSERT(dof_idx < num_entries);

          // update vector entry
          Tiny::axpy(data[dof_idx], loc_vec[i], alpha);
        }
      }

      #if defined(__CUDACC__) || defined(DOXYGEN)
      /**
       * \brief Dense Vector grouped scatter axpy function async version
       *
       * \tparam ThreadGroup_ The thread group type, that work together.
       * \tparam inner_size_ The size of the internal data, i.e. each object map is mapping to consists of inner_dim_ entries
       *
       * \param[in] tg A reference to the thread group object on which to cooperate.
       * \param[in] scatter_size The number of blocked size elements to be scattered
       * \param[in] scatter_offset The offset into the global vector in which to scatter
       * \param[in] loc_vec The local vector which values are to be scattered.
       * \param[out] data A pointer to the vector data array to be scattered in. Has to be initialized to the correct size beforehand. Previous values are added on.
       * \param[in] num_entries Number of entries of the dense vector.
       * \param[in] map The dof mapping used to map local to global values.
       * \param[in] data_size Number of blocked elements to be scattered.
       * \param[in] alpha Optional scaling parameter.
       *
       * \note Contrary to the standard method, this methods always uses the blank data array of the underlying containers.
       */
      template<typename ThreadGroup_, int inner_size_>
      CUDA_DEVICE static void __forceinline__ grouped_scatter_vector_dense(const ThreadGroup_& tg, const int scatter_size, const int scatter_offset, const DataType* loc_vec,
                                                            DataType* data, [[maybe_unused]] IndexType num_entries, const IndexType* map,
                                                            int data_size, DataType alpha = DataType(1))
      {
        // stride based for loop
        for(int idx = tg.thread_rank(); (idx < scatter_size * inner_size_) && ((idx + scatter_offset*inner_size_) < data_size*inner_size_); idx += tg.num_threads())
        {
          // get dof index
          Index dof_idx = map[idx/inner_size_];
          // ASSERT(dof_idx < num_entries);

          // update vector entry
          data[dof_idx*inner_size_+(idx%inner_size_)] +=  alpha * loc_vec[idx];
        }
      }

      #endif

      /**
       * \brief Dense Vector gather axpy function
       *
       * \tparam numr_ The number of rows of the local matrix.
       *
       * \param[out] loc_vec The local vector which values are to be gathered in. Previous values are added onto.
       * \param[in] data A pointer to the vector data array to be gathered from.
       * \param[in] num_entries Number of entries of the dense vector.
       * \param[in] map The dof mapping used to map local to global values.
       * \param[in] alpha Optional scaling parameter.
       *
       * \note This method also works for blocked variants. The interface is always used from the *native* variant,
       *         i.e. InnerType_ refers to Tiny::Vector .
       */
      template<typename InnerType_, int numr_>
      CUDA_HOST_DEVICE static void gather_vector_dense(Tiny::Vector<InnerType_, numr_>& loc_vec, const InnerType_* data, [[maybe_unused]] IndexType num_entries, const IndexType* map, DataType alpha = DataType(1))
      {
        #ifndef __CUDACC__
        static_assert(std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #else
        static_assert(::cuda::std::is_same<typename Tiny::Intern::DataTypeExtractor<InnerType_>::MyDataType, DataType>(), "Inner Datatype does not match!");
        #endif
        // loop over all local entries
        for (int i(0); i < numr_; ++i)
        {
          // get dof index
          Index dof_idx = map[i];
          // ASSERT(dof_idx < num_entries);

          // update local vector data
          Tiny::axpy(loc_vec[i], data[dof_idx], alpha);
        }
      }

      #if defined(__CUDACC__) || defined(DOXYGEN)
      /**
       * \brief Dense Vector grouped scatter axpy function
       *
       * \tparam ThreadGroup_ The thread group type, that work together.
       * \tparam inner_size_ The size of the internal data, i.e. each object map is mapping to consists of inner_dim_ entries
       *
       * \param[in] tg A reference to the thread group object on which to cooperate.
       * \param[out] loc_vec The local vector in which the values are to be gathered. Has to be a shared array
       * \param[in] data A pointer to the vector data array to be gather from.
       * \param[in] num_entries Number of entries of the dense vector.
       * \param[in] map The dof mapping used to map local to global values.
       * \param[in] num_data Number of blocked Elements to be gathered.
       * \param[in] alpha Optional scaling parameter.
       *
       * \note Contrary to the standard method, this methods always uses the blank data array of the underlying containers.
       * \note generally requires a synchronization after the call if threadgroup is not a warp
       */
      template<typename ThreadGroup_, int inner_size_>
      CUDA_DEVICE static void __forceinline__ grouped_gather_vector_dense(const ThreadGroup_& tg, DataType* loc_vec,
                                                            const DataType* data, [[maybe_unused]] IndexType num_entries, const IndexType* map,
                                                            int num_data, DataType alpha = DataType(1))
      {
        // stride based for loop
        for(int i = tg.thread_rank(); i < num_data*inner_size_; i += tg.num_threads())
        {
          // get dof index
          Index dof_idx = map[i/inner_size_];
          // ASSERT(dof_idx < num_entries);

          // update vector entry
          loc_vec[i] += alpha * data[dof_idx*inner_size_+(i%inner_size_)];
        }
      }

      #endif
    }; //struct GPUVectorGatherScatterHelper
  }
}

#endif
