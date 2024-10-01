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
      return (TestMatrixFlags)(((int)a) & ((int)b));
    }

    /// bit-wise OR operator for TestMatrixFlags
    inline TestMatrixFlags operator|(TestMatrixFlags a, TestMatrixFlags b)
    {
      return (TestMatrixFlags)(((int)a) | ((int)b));
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
      void _fill_generic_values(SparseMatrixCSR<DT_, IT_>& matrix)
      {
        DT_* values = matrix.val();
        Index nnze = matrix.used_elements();
        for(Index i(0); i < nnze; ++i)
          values[i] = rng(-DT_(1), DT_(1));
      }
    }; // class TestMatrixFactory
  } // namespace LAFEM
} // namespace FEAT
