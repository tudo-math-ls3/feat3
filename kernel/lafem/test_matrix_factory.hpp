// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <iostream>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/adjacency/dynamic_graph.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/util/random.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Flags for desired test matrix
     */
    enum class TestMatrixFlags : int
    {
      /// generic matrix without special properties, values in range [-1,+1]
      generic           = 0x00,
      /// matrix with symmetric structure
      symmetric_struct  = 0x01,
      /// matrix with symmetric structure and values
      symmetric         = 0x02,
      /// all values non-negative, i.e. values in range [0,1]
      non_negative      = 0x04,
      /// all rows non-empty
      non_empty_rows    = 0x08,
      /// all columns non-empty
      non_empty_cols    = 0x10,
      /// all main diagonal entries non-empty
      non_empty_diag    = 0x20,
      /// all main diagonal entries non-zero
      non_zero_diag     = 0x20,
      /// diagonal dominant matrix
      diagonal_dominant = 0x40
    };

    /// bit-wise AND operator for TestMatrixFlags
    inline TestMatrixFlags operator&(TestMatrixFlags a, TestMatrixFlags b)
    {
      return static_cast<TestMatrixFlags>(static_cast<int>(a) & static_cast<int>(b));
    }

    /// bit-wise OR operator for TestMatrixFlags
    inline TestMatrixFlags operator|(TestMatrixFlags a, TestMatrixFlags b)
    {
      return static_cast<TestMatrixFlags>(static_cast<int>(a) | static_cast<int>(b));
    }

    /// checks whether a & b != 0
    inline bool operator*(TestMatrixFlags a, TestMatrixFlags b)
    {
      return static_cast<int>(a & b) != 0;
    }

    class TestMatrixFactory
    {
    public:
      Random rng;

      TestMatrixFactory() = default;

      Random::SeedType get_rng_seed() const
      {
        return rng.get_seed();
      }

      /**
      * \brief Creates a sparse matrix in CSR format with various properties based on the specified flags.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) to be created and filled.
      *
      * \param[in] nrows
      * The number of rows in the matrix.
      *
      * \param[in] ncols
      * The number of columns in the matrix.
      *
      * \param[in] nnze
      * The number of non-zero entries in the matrix.
      *
      * \param[in] flags
      */

      template<typename DT_, typename IT_>
      void create(SparseMatrixCSR<DT_, IT_>& matrix, Index nrows, Index ncols, Index nnze, TestMatrixFlags flags)
      {
        if(flags * (TestMatrixFlags::non_negative) && flags * (TestMatrixFlags::diagonal_dominant))
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_non_negative_values(matrix);
          _make_diagonal_dominant(matrix);
        }
        else if(flags * (TestMatrixFlags::non_negative) && flags * (TestMatrixFlags::symmetric))
        {
          _create_symmetric_struct(matrix, nrows, ncols, nnze);
          _fill_non_negative_values(matrix);
          _make_symmetric_values(matrix);
        }
        else if(flags * TestMatrixFlags::symmetric_struct)
        {
          _create_symmetric_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
        }
        else if(flags * TestMatrixFlags::symmetric)
        {
          _create_symmetric_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
          _make_symmetric_values(matrix);
        }
        else if(flags * TestMatrixFlags::non_negative)
        {
          _create_generic_struct(matrix, nrows, ncols, nnze);
          _fill_non_negative_values(matrix);
        }
        else if(flags * TestMatrixFlags::non_empty_diag)
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
        }
        else if(flags * TestMatrixFlags::non_zero_diag)
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
          _make_non_zero_diag(matrix);
        }
        else if(flags * TestMatrixFlags::diagonal_dominant)
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
          _make_diagonal_dominant(matrix);
        }
        else
        {
          _create_generic_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
        }
      }

      /**
      * \brief Creates a sparse matrix in BCSR format with various properties based on the specified flags.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) to be created and filled.
      *
      * \param[in] nrows
      * The number of rows in the matrix.
      *
      * \param[in] ncols
      * The number of columns in the matrix.
      *
      * \param[in] nnze
      * The number of non-zero entries in the matrix.
      *
      * \param[in] flags
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void create(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix, Index nrows, Index ncols, Index nnze, TestMatrixFlags flags)
      {
        if(flags * (TestMatrixFlags::non_negative) && flags * (TestMatrixFlags::diagonal_dominant))
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_non_negative_values(matrix);
          _make_diagonal_dominant(matrix);
        }
        else if(flags * (TestMatrixFlags::non_negative) && flags * (TestMatrixFlags::symmetric))
        {
          _create_symmetric_struct(matrix, nrows, ncols, nnze);
          _fill_non_negative_values(matrix);
          _make_symmetric_values(matrix);
        }
        else if(flags * TestMatrixFlags::symmetric_struct)
        {
          _create_symmetric_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
        }
        else if(flags * TestMatrixFlags::symmetric)
        {
          _create_symmetric_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
          _make_symmetric_values(matrix);
        }
        else if(flags * TestMatrixFlags::non_negative)
        {
          _create_generic_struct(matrix, nrows, ncols, nnze);
          _fill_non_negative_values(matrix);
        }
        else if(flags * TestMatrixFlags::non_empty_diag)
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
        }
        else if(flags * TestMatrixFlags::non_zero_diag)
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
          _make_non_zero_diag(matrix);
        }
        else if(flags * TestMatrixFlags::diagonal_dominant)
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
          _make_diagonal_dominant(matrix);
        }
        else
        {
          _create_generic_struct(matrix, nrows, ncols, nnze);
          _fill_generic_values(matrix);
        }
      }

    protected:
      /**
      * \brief Creates a generic sparse matrix structure in CSR format based on the specified dimensions and non-zero entries.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) to be created and initialized with a generic structure.
      *
      * \param[in] nrows
      * The number of rows in the matrix, must be greater than zero.
      *
      * \param[in] ncols
      * The number of columns in the matrix, must be greater than zero.
      *
      * \param[in] nnze
      * The desired number of non-zero entries in the matrix.
      */
      Adjacency::Graph _create_generic_struct(Index nrows, Index ncols, Index nnze)
      {
        XASSERT(nrows > 0u);
        XASSERT(ncols > 0u);
        XASSERT(nnze <= nrows * ncols);
        Adjacency::DynamicGraph dygraph(nrows, ncols);
        Index nm = nrows * ncols - Index(1);
        while(dygraph.get_num_indices() < nnze)
        {
          Index k = rng(Index(0), nm);
          dygraph.insert(k / ncols, k % ncols);
        }
        Adjacency::Graph graph(Adjacency::RenderType::as_is, dygraph);
        return graph;
      }

      template<typename DT_, typename IT_>
      void _create_generic_struct(SparseMatrixCSR<DT_, IT_>& matrix, Index nrows, Index ncols, Index nnze)
      {
        matrix = SparseMatrixCSR<DT_, IT_>(_create_generic_struct(nrows, ncols, nnze));
      }

      /**
      * \brief Creates a generic sparse matrix structure in BCSR format based on the specified dimensions and non-zero entries.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) to be created and initialized with a generic structure.
      *
      * \param[in] nrows
      * The number of rows in the matrix, must be greater than zero.
      *
      * \param[in] ncols
      * The number of columns in the matrix, must be greater than zero.
      *
      * \param[in] nnze
      * The desired number of non-zero entries in the matrix.
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void _create_generic_struct(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix, Index nrows, Index ncols, Index nnze)
      {
        matrix = SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>(_create_generic_struct(nrows, ncols, nnze));
      }

      /**
      * \brief Creates a symmetric sparse matrix structure in CSR format based on the specified dimensions and non-zero entries.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) to be created and initialized with a generic structure.
      *
      * \param[in] nrows
      * The number of rows in the matrix, must be greater than zero.
      *
      * \param[in] ncols
      * The number of columns in the matrix, must be greater than zero.
      *
      * \param[in] nnze
      * The desired number of non-zero entries in the matrix.
      *
      * \note The matrix has to be square
      * There might be nnze +- 1 entries
      */

      Adjacency::Graph _create_symmetric_struct(Index nrows, Index ncols, Index nnze) {
        XASSERT(nrows > 0u);
        XASSERT(ncols > 0u);
        XASSERT(nrows == ncols);
        XASSERT(nnze <= nrows * ncols);
        Adjacency::DynamicGraph dygraph(nrows, ncols);
        Index nm = nrows * ncols - Index(1);
        while(dygraph.get_num_indices() < nnze)
        {
          Index k = rng(Index(0), nm);
          dygraph.insert(k / ncols, k % ncols);
          dygraph.insert(k % ncols, k / ncols);
        }
        Adjacency::Graph graph(Adjacency::RenderType::as_is, dygraph);
        return graph;
      }

      template<typename DT_, typename IT_>
      void _create_symmetric_struct(SparseMatrixCSR<DT_, IT_>& matrix, Index nrows, Index ncols, Index nnze)
      {
        matrix = SparseMatrixCSR<DT_, IT_>(_create_symmetric_struct(nrows, ncols, nnze));
      }

      /**
      * \brief Creates a symmetric sparse matrix structure in BCSR format based on the specified dimensions and non-zero entries.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) to be created and initialized with a generic structure.
      *
      * \param[in] nrows
      * The number of rows in the matrix, must be greater than zero.
      *
      * \param[in] ncols
      * The number of columns in the matrix, must be greater than zero.
      *
      * \param[in] nnze
      * The desired number of non-zero entries in the matrix.
      *
      * \note The matrix as well as the blocks have to be square
      * There might be nnze +-1 entries
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void _create_symmetric_struct(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix, Index nrows, Index ncols, Index nnze)
      {
        XASSERT(BlockWidth_ == BlockHeight_);
        matrix = SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>(_create_symmetric_struct(nrows, ncols, nnze));
      }

      /**
      * \brief Creates a non empty diagonal sparse matrix structure in CSR format based on the specified dimensions and non-zero entries.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) to be created and initialized with a generic structure.
      *
      * \param[in] nrows
      * The number of rows in the matrix, must be greater than zero.
      *
      * \param[in] ncols
      * The number of columns in the matrix, must be greater than zero.
      *
      * \param[in] nnze
      * The desired number of non-zero entries in the matrix.
      *
      * \note the matrix has to be square
      */
      Adjacency::Graph _create_non_empty_diag_struct(Index nrows, Index ncols, Index nnze) {
        XASSERT(nrows > 0u);
        XASSERT(ncols > 0u);
        XASSERT(nrows == ncols);
        XASSERT(nnze <= nrows * ncols);
        XASSERT(nnze >= ncols);
        Adjacency::DynamicGraph dygraph(nrows, ncols);
        Index nm = nrows * ncols - Index(1);
        for(Index i(0); i < ncols; ++i)
        {
          dygraph.insert(i, i);
        }

        while(dygraph.get_num_indices() < nnze)
        {
          Index k = rng(Index(0), nm);
          dygraph.insert(k / ncols, k % ncols);
        }
        Adjacency::Graph graph(Adjacency::RenderType::as_is, dygraph);
        return graph;
      }

      template<typename DT_, typename IT_>
      void _create_non_empty_diag_struct(SparseMatrixCSR<DT_, IT_>& matrix, Index nrows, Index ncols, Index nnze)
      {
        matrix = SparseMatrixCSR<DT_, IT_>(_create_non_empty_diag_struct(nrows, ncols, nnze));
      }

      /**
      * \brief Creates a non empty diagonal sparse matrix structure in BCSR format based on the specified dimensions and non-zero entries.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) to be created and initialized with a generic structure.
      *
      * \param[in] nrows
      * The number of rows in the matrix, must be greater than zero.
      *
      * \param[in] ncols
      * The number of columns in the matrix, must be greater than zero.
      *
      * \param[in] nnze
      * The desired number of non-zero entries in the matrix.
      *
      * \note the matrix as well as the blocks have to be square
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void _create_non_empty_diag_struct(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix, Index nrows, Index ncols, Index nnze)
      {
        XASSERT(BlockWidth_ == BlockHeight_);
        matrix = SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>(_create_non_empty_diag_struct(nrows, ncols, nnze));
      }

      /**
      * \brief Fills the values of a sparse matrix in CSR format with random values between -1 and 1.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) whose non-zero values will be filled
      * with randomly generated values.
      */
      template<typename DT_, typename IT_>
      void _fill_generic_values(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        const Index nnze = matrix.used_elements();
        for(Index i(0); i < nnze; ++i)
          values[i] = rng(-DT_(1), DT_(1));
      }

      /**
      * \brief Fills the values of a sparse matrix in BCSR format with random values between -1 and 1.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) whose non-zero values will be filled
      * with randomly generated values.
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void _fill_generic_values(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix)
      {
        auto * values = matrix.val();
        const Index nnze = matrix.used_elements();
        for(Index i(0); i < nnze; ++i)
          for(int j(0); j < BlockHeight_; ++j)
            for(int k(0); k< BlockWidth_; ++k)
              values[i][j][k] = rng(-DT_(1), DT_(1));
      }

      /**
      * \brief Fills the values of a sparse matrix in CSR format with random non negative values between 0 and 1.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) whose non-zero values will be filled
      * with randomly generated values.
      */
      template<typename DT_, typename IT_>
      void _fill_non_negative_values(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        const Index nnze = matrix.used_elements();
        for(Index i(0); i < nnze; ++i)
          values[i] = rng(DT_(0), DT_(1));
      }

      /**
      * \brief Fills the values of a sparse matrix in BCSR format with random non negative values between 0 and 1.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) whose non-zero values will be filled
      * with randomly generated values.
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void _fill_non_negative_values(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix)
      {
        auto* values = matrix.val();
        const Index nnze = matrix.used_elements();
        for(Index i(0); i < nnze; ++i)
          for(int j(0); j < BlockHeight_; ++j)
            for(int k(0); k< BlockWidth_; ++k)
              values[i][j][k] = rng(DT_(0), DT_(1));
      }

      /**
      * \brief Ensures that the values of a sparse matrix in CSR format are symmetric.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) whose non-zero values will be modified
      * to ensure symmetry.
      *
      * \note
      * The matrix has to be square
      */
      template<typename DT_, typename IT_>
      void _make_symmetric_values(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        const Index nrows = matrix.rows();
        const IT_ * col_idx = matrix.col_ind();
        const IT_ * row_ptr = matrix.row_ptr();

        for (Index i = 0; i < nrows; ++i)
        {
          for (IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
          {
            Index col_j = col_idx[j];
            if(col_j > i)
              break;
            else
            {
              for (IT_ k = row_ptr[col_j]; k < row_ptr[col_j + 1]; ++k)
              {
                if (col_idx[k] == i)
                {
                  values[k] = values[j];
                  break;
                }
              }
            }
          }
        }
      }

      /**
      * \brief Ensures that the values of a sparse matrix in BCSR format are symmetric.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) whose non-zero values will be modified
      * to ensure symmetry.
      *
      * \note
      * The matrix as well as the blocks have to be square
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void _make_symmetric_values(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix)
      {
        auto* values = matrix.val();
        const IT_* col_idx = matrix.col_ind();
        const IT_* row_ptr = matrix.row_ptr();
        const Index nrows = matrix.rows();

        for (Index i = 0; i < nrows; ++i) {
          for (IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            Index col_j = col_idx[j];
            if(col_j > i)
              break;
            else
            {
              for (IT_ k = row_ptr[col_j]; k < row_ptr[col_j + 1]; ++k)
              {
                if (col_idx[k] == i)
                {
                  for (int l = 0; l < BlockHeight_; ++l)
                  {
                    for (int m = 0; m < BlockWidth_; ++m)
                    {
                      values[k][m][l] = values[j][l][m];
                    }
                  }
                  break;
                }
              }
            }
          }
        }
      }

      /**
      * \brief Ensures that the diagonal entries of a sparse matrix in CSR format are non-zero.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) whose diagonal entries will be modified
      * to ensure they are non-zero.
      *
      * \note
      * The matrix has to be square
      */
      template<typename DT_, typename IT_>
      void _make_non_zero_diag(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        const IT_* col_idx = matrix.col_ind();
        const IT_* row_ptr = matrix.row_ptr();
        const Index nrows = matrix.rows();

        for(Index i = 0; i < nrows; ++i)
        {
          for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
          {
            if(col_idx[j] == i)  // Check if this is a diagonal element
            {
              values[j] = rng(DT_(0.1), DT_(1));
              break;
            }
          }
        }
      }

      /**
      * \brief Ensures that the diagonal entries of a sparse matrix in BCSR format are non-zero.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) whose diagonal entries will be modified
      * to ensure they are non-zero.
      *
      * \note
      * The matrix as well as the blocks have to be square
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void _make_non_zero_diag(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix)
      {
        auto* values = matrix.val();
        const IT_* col_idx = matrix.col_ind();
        const IT_* row_ptr = matrix.row_ptr();
        const Index nrows = matrix.rows();

        for(Index i = 0; i < nrows; ++i)
        {
          for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
          {
            if(col_idx[j] == i)  // Check if this is a diagonal element
            {
              for(int l(0); l < BlockHeight_; ++l)
                values[j][l][l] = rng(DT_(0.1), DT_(1));
              break;
            }
          }
        }
      }

      /**
      * \brief Ensures that a sparse matrix in CSR format is diagonal dominant.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in CSR format) whose diagonal entries will be modified
      * to ensure that the matrix is diagonal dominant.
      *
      * \note
      * The matrix has to be square
      */
      template<typename DT_, typename IT_>
      void _make_diagonal_dominant(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        const IT_* col_idx = matrix.col_ind();
        const IT_* row_ptr = matrix.row_ptr();
        const Index nrows = matrix.rows();

        for(Index i = 0; i < nrows; ++i)
        {
          DT_ off_diagonal = 0;
          for(IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
          {
            if(col_idx[j] != i) // Skip the diagonal element
            {
              off_diagonal += std::abs(values[j]);
            }
          }

          for(IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
          {
            if(col_idx[j] == i) // Check if this is the diagonal element
            {
              values[j] = off_diagonal + rng(DT_(0.1), DT_(1)); // Set diagonal value and ensure it's greater than the sum of off-diagonal elements
              break;
            }
          }
        }
      }

      /**
      * \brief Ensures that a sparse matrix in BCSR format is diagonal dominant.
      *
      * \param[in,out] matrix
      * The sparse matrix object (in BCSR format) whose diagonal entries will be modified
      * to ensure that the matrix is diagonal dominant.
      *
      * \note
      * The matrix as well as the blocks have to be square
      */
      template<typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void _make_diagonal_dominant(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix)
      {
        auto* values = matrix.val();
        const IT_* col_idx = matrix.col_ind();
        const IT_* row_ptr = matrix.row_ptr();
        const Index nrows = matrix.rows();

        for(Index i = 0; i < nrows; ++i)
        {
          DT_ off_diagonal = 0;
          for(IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
          {
            if(col_idx[j] != i) // Skip the diagonal element
            {
              for(int k(0); k < BlockWidth_; ++k)
              {
                for(int l(0); l < BlockHeight_; ++l)
                {
                  off_diagonal += std::abs(values[j][k][l]);
                }
              }
            }
            else
            {
              for(int k(0); k < BlockWidth_; ++k)
              {
                for(int l(0); l < BlockHeight_; ++l)
                {
                  if(k != l)
                    off_diagonal += std::abs(values[j][k][l]);
                }
              }
            }
          }

          for(IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
          {
            if(col_idx[j] == i) // Check if this is the diagonal element
            {
              for(int k(0); k < BlockWidth_; ++k)
                values[j][k][k] = off_diagonal + rng(DT_(0.1), DT_(1)); // Set diagonal value and ensure it's greater than the sum of off-diagonal elements
              break;
            }
          }
        }
      }
    }; // class TestMatrixFactory
  } // namespace LAFEM
} // namespace FEAT
