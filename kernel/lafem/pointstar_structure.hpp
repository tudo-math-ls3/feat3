// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>

#include <stack>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief empty Finite-Elements pointstar matrix creator.
     *
     * This class generates the matrix-structure for Finite-Elements on a structured mesh.
     *
     * \author Christoph Lohmann
     */
    struct PointstarStructureFE
    {
      /**
       * \brief Generates an empty FE-style pointstar banded matrix
       *
       * \param[in] fe_degree
       * Polynomial degree of the Finite Element discretization.
       * I.e. FE_order=1 yields the so called Q1 Matrix with three tridiagonal bands.
       *
       * \param[in] num_of_subintervalls
       * The vector with number of subintervalls per dimension.
       *
       * \returns
       * The m^d x m^d FE-stye pointstar matrix.
       */
      template<typename DataType_, typename IndexType_>
      static SparseMatrixBanded<DataType_, IndexType_> value(const Index fe_degree,
                                                                        const std::vector<IndexType_> & num_of_subintervalls)
      {
        const IndexType_ * const pnos(num_of_subintervalls.data());

        // get number of dimensions
        const Index d(Index(num_of_subintervalls.size()));

        // output of errors if wrong input
        XASSERTM(d >= IndexType_(1), "You need at least 1 dimension");

        for (IndexType_ i(0); i < d; ++i)
        {
          XASSERTM(pnos[i] >= IndexType_(3), "You need at least 3 subintervalls per dimension");
        }

        // calculate size of matrix and number of offsets
        IndexType_ size(1);
        Index noo(1);
        if (fe_degree > Index(0))
        {
          for (Index i(0); i < d; ++i)
          {
            size *= pnos[i] * IndexType_(fe_degree) + 1;
            noo *=  2 * fe_degree + 1;
          }
        }
        else
        {
          for (Index i(0); i < d; ++i)
          {
            size *= pnos[i];
          }
        }

        // allocate memory for vectors of matrix
        DenseVector<IndexType_, IndexType_> vec_offsets(noo);
        DenseVector<DataType_, IndexType_> vec_val(noo * Index(size));

        // fill offsets-vector
        IndexType_ * const poffsets(vec_offsets.elements());
        const Index h_off((noo - 1) / 2);

        // save position of main-diagonal
        poffsets[h_off] = size - IndexType_(1);

        for (Index i(0), k(1), m(1); i < d; ++i, k *= 2 * fe_degree + 1, m *= pnos[i - 1] * fe_degree + 1)
        {
          Index k1((k - 1) / 2);

          for (IndexType_ j(1); j <= fe_degree; ++j)
          {
            for (IndexType_ l(0); l < k; ++l)
            {
              poffsets[h_off - k1 + l + k * j] = poffsets[h_off - k1 + l] + IndexType_(j * m);
              poffsets[h_off - k1 + l - k * j] = poffsets[h_off - k1 + l] - IndexType_(j * m);
            }
          }
        }

        // return the matrix
        return SparseMatrixBanded<DataType_, IndexType_>(Index(size), Index(size), vec_val, vec_offsets);
      }
    }; // struct PointstarStructureFE

    /**
     * \brief empty Finite-Differences pointstar matrix creator.
     *
     * This class generates the matrix-structure for Finite-Differences on a structured mesh.
     *
     * \author Christoph Lohmann
     */
    struct PointstarStructureFD
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
      static SparseMatrixBanded<DataType_, IndexType_> value(const std::vector<IndexType_> & num_of_subintervalls)
      {
        const IndexType_ * const pnos(num_of_subintervalls.data());

        const Index d(Index(num_of_subintervalls.size()) - 1);

        // calculate dimension of the matrix
        IndexType_ size(1);
        for (Index i(1); i <= d; ++i)
        {
          size *= pnos[i] - 1;
        }

        // calculate number of offsets
        const Index num_of_offsets(2 * d + 1);

        // allocate memory for vectors of matrix
        DenseVector<DataType_, IndexType_> vec_val(Index(size) * num_of_offsets);
        DenseVector<IndexType_, IndexType_> vec_offsets(num_of_offsets);

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
        return SparseMatrixBanded<DataType_, IndexType_>(Index(size), Index(size), vec_val, vec_offsets);
      }
    }; // struct PointstarStructureFD

  } // namespace LAFEM
} // namespace FEAT
