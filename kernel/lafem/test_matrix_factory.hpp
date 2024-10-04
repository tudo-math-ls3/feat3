// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
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

      void seed_rng_by_timer(int rank = 0)
      {
        rng = Random(Random::get_seed(rank));
      }

      template<typename DT_, typename IT_>
      void create(SparseMatrixCSR<DT_, IT_>& matrix, Index nrows, Index ncols, Index nnze, TestMatrixFlags flags)
      {
        _create_generic_struct(matrix, nrows, ncols, nnze);
        _fill_generic_values(matrix);

        if(flags * (TestMatrixFlags::non_negative & TestMatrixFlags::diagonal_dominant))
        {
          _create_non_empty_diag_struct(matrix, nrows, ncols, nnze);
          _fill_non_negative_values(matrix);
          _make_diagonal_dominant(matrix);
        }
        else if(flags * (TestMatrixFlags::non_negative & TestMatrixFlags::symmetric))
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
      }

    protected:
      template<typename DT_, typename IT_>
      void _create_generic_struct(SparseMatrixCSR<DT_, IT_>& matrix, Index nrows, Index ncols, Index nnze)
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
        matrix = SparseMatrixCSR<DT_, IT_>(graph);
      }

      template<typename DT_, typename IT_>
      void _create_symmetric_struct(SparseMatrixCSR<DT_, IT_>& matrix, Index nrows, Index ncols, Index nnze)
      {
        XASSERT(nrows > 0u);
        XASSERT(ncols > 0u);
        XASSERT(nrows == ncols);
        XASSERT(nnze <= nrows * ncols);
        Adjacency::DynamicGraph dygraph(nrows, ncols);
        Index nm = nrows * ncols - Index(1);
        while(dygraph.get_num_indices() < nnze)
        {
          Index k = rng(Index(0), nm);
          if(dygraph.get_num_indices() == nnze - 1)
          {
            dygraph.insert(k / ncols, k / ncols);
          }
          else
          {
            dygraph.insert(k / ncols, k % ncols);
            dygraph.insert(k % ncols, k / ncols);
          }
        }
        Adjacency::Graph graph(Adjacency::RenderType::as_is, dygraph);
        matrix = SparseMatrixCSR<DT_, IT_>(graph);
      }

      template<typename DT_, typename IT_>
      void _create_non_empty_diag_struct(SparseMatrixCSR<DT_, IT_>& matrix, Index nrows, Index ncols, Index nnze)
      {
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
        matrix = SparseMatrixCSR<DT_, IT_>(graph);
      }

      template<typename DT_, typename IT_>
      void _fill_generic_values(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        Index nnze = matrix.used_elements();
        for(Index i(0); i < nnze; ++i)
          values[i] = rng(-DT_(1), DT_(1));
      }

      template<typename DT_, typename IT_>
      void _fill_non_negative_values(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        Index nnze = matrix.used_elements();
        for(Index i(0); i < nnze; ++i)
          values[i] = rng(DT_(0), DT_(1));
      }

      template<typename DT_, typename IT_>
      void _make_symmetric_values(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        Index nnze = matrix.used_elements();
        IT_ * col_idx = matrix.col_ind();
        IT_ * row_ptr = matrix.row_ptr();
        for (Index i(0); i < nnze; ++i) {
          Index row_i = std::upper_bound(row_ptr, row_ptr + matrix.rows(), i) - row_ptr - 1;
          Index col_i = col_idx[i];

          DT_ value = values[i];

          Index start = row_ptr[col_i];
          Index end = row_ptr[col_i + 1];

          for (Index j = start; j < end; ++j) {
            if (col_idx[j] == row_i) {
              values[j] = value;
              break;
            }
          }
        }
      }
      template<typename DT_, typename IT_>
      void _make_non_zero_diag(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        IT_* col_idx = matrix.col_ind();
        IT_* row_ptr = matrix.row_ptr();
        const Index nrows = matrix.rows();
        bool diagonal_entry = false;

        for(Index i = 0; i < nrows; ++i)
        {
          diagonal_entry = false;
          for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
          {
            if(col_idx[j] == i)  // Check if this is a diagonal element
            {
              values[j] = rng(DT_(0.1), DT_(1));
              diagonal_entry = true;
              break;
            }
          }
        }
      }

      template<typename DT_, typename IT_>
      void _make_diagonal_dominant(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        IT_* col_idx = matrix.col_ind();
        IT_* row_ptr = matrix.row_ptr();
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

    }; // class TestMatrixFactory
  } // namespace LAFEM
} // namespace FEAT
