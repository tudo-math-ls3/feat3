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
    }; //struct GPUVectorGatherScatterHelper
  }
}

#endif
