#pragma once
#ifndef KERNEL_LAFEM_POINTSTAR_STRUCTURE_HPP
#define KERNEL_LAFEM_POINTSTAR_STRUCTURE_HPP 1

#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>

#include <stack>

namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct PointstarStructureFE
    {
    };

    /**
     * \brief empty Finite-Elements pointstar matrix creator.
     *
     * This class generates the matrix-structure for Finite-Elements on a structured mesh.
     *
     * \author Christoph Lohmann
     */
    template <>
    struct PointstarStructureFE<Algo::Generic>
    {
      /**
       * \brief Generates an empty FE-style pointstar banded matrix
       *
       * \param[in] FE_order
       * order of the Finite Element discretisation
       * \param[in] num_of_subintervalls
       * The vector with number of subintervalls per dimension.
       *
       * \returns
       * The m^d x m^d FE-stye pointstar matrix.
       */
      template<typename DataType_, typename IndexType_>
      static SparseMatrixBanded<Mem::Main, DataType_, IndexType_> value(const Index fe_order,
                                                                        const std::vector<IndexType_> & num_of_subintervalls)
      {
        const IndexType_ * const pnos(num_of_subintervalls.data());

        // get number of dimensions
        const Index d(num_of_subintervalls.size());

        // output of errors if wrong input
        ASSERT(d >= IndexType_(1), "You need at least 1 dimension");

        for (IndexType_ i(0); i < d; ++i)
        {
          ASSERT(pnos[i] >= IndexType_(3), "You need at least 3 subintervalls per dimension");
        }

        // calculate size of matrix and number of offsets
        IndexType_ size(1);
        Index noo(1);
        for (Index i(0); i < d; ++i)
        {
          size *= pnos[i] * fe_order - 1;
          noo *=  2 * fe_order + 1;
        }

        // allocate memory for vectors of matrix
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_offsets(noo);
        DenseVector<Mem::Main, DataType_, IndexType_> vec_val(noo * size);

        // fill offsets-vector
        IndexType_ * const poffsets(vec_offsets.elements());
        const Index h_off((noo - 1) / 2);

        // save position of main-diagonal
        poffsets[h_off] = IndexType_(size - 1);

        for (Index i(0), k(1), m(1); i < d; ++i, k *= 2 * fe_order + 1, m *= pnos[i - 1] * fe_order - 1)
        {
          IndexType_ k1((k - 1) / 2);

          for (IndexType_ j(1); j <= fe_order; ++j)
          {
            for (IndexType_ l(0); l < k; ++l)
            {
              poffsets[h_off - k1 + l + k * j] = poffsets[h_off - k1 + l] + IndexType_(j * m);
              poffsets[h_off - k1 + l - k * j] = poffsets[h_off - k1 + l] - IndexType_(j * m);
            }
          }
        }

        // return the matrix
        return SparseMatrixBanded<Mem::Main, DataType_, IndexType_>(size, size, vec_val, vec_offsets);
      }
    }; // struct PointstarStructureFE

    template <typename Algo_>
    struct PointstarStructureFD
    {
    };

    /**
     * \brief empty Finite-Differences pointstar matrix creator.
     *
     * This class generates the matrix-structure for Finite-Differences on a structured mesh.
     *
     * \author Christoph Lohmann
     */
    template <>
    struct PointstarStructureFD<Algo::Generic>
    {
      /**
       * \brief Generates an empty FD-style pointstar banded matrix
       *
       * \param[in] num_of_subintervalls
       * The vector with number of subintervalls per dimension plus a leading 2.
       *
       * \returns
       * The m^d x m^d FD-stye pointstar matrix.
       */
      template<typename DataType_, typename IndexType_>
      static SparseMatrixBanded<Mem::Main, DataType_, IndexType_> value(const std::vector<IndexType_> & num_of_subintervalls)
      {
        const IndexType_ * const pnos(num_of_subintervalls.data());

        const Index d(num_of_subintervalls.size() - 1);

        // calculate dimension of the matrix
        IndexType_ size(1);
        for (Index i(1); i <= d; ++i)
        {
          size *= pnos[i] - 1;
        }

        // calculate number of offsets
        const Index num_of_offsets(2 * d + 1);

        // allocate memory for vectors of matrix
        DenseVector<Mem::Main, DataType_, IndexType_> vec_val(size * num_of_offsets);
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_offsets(num_of_offsets);

        // fill vec_offsets
        IndexType_ * const poffsets(vec_offsets.elements());

        poffsets[d] = IndexType_(size - 1);

        Index tmp(1);
        for (Index i(0); i < d; ++i)
        {
          tmp *= pnos[i] - 1;
          poffsets[d - 1 - i] = size - IndexType_(1 + tmp);
          poffsets[d + 1 + i] = size + IndexType_(tmp - 1);
        }

        // return the matrix
        return SparseMatrixBanded<Mem::Main, DataType_, IndexType_>(size, size, vec_val, vec_offsets);
      }
    }; // struct PointstarStructureFD

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POINTSTAR_STRUCTURE_HPP
