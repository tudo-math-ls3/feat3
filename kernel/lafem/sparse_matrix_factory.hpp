// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_FACTORY_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_FACTORY_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <map>
#include <iterator>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Factory for SparseMatrix construction.
     *
     * The purpose of this class is to make implementing CSR matrices very simple.\n
     * The data is stored in a map and can be transfered to a CSR matrix.
     *
     * \tparam DT_
     * The data type to be used for the sparse matrix.
     *
     * \tparam IT_
     * The index type to be used for the sparse matrix.
     *
     * \author Gesa Pottbrock
     */
    template<typename DT_, typename IT_>
    class SparseMatrixFactory
    {
    private:
      /// Number of rows for the created matrix.
      IT_ _num_row;
      /// Number of columns for the created matrix.
      IT_ _num_col;
      /// Maximal size of matrix.
      static constexpr IT_ max_size= IT_(100000);
      /// Map collecting the input for the matrix.
      std::map<IT_, DT_> memory;

    public:

      /**
       * \brief Constructor
       *
       * \param[in] num_rows The number of rows for the matrix.
       * \param[in] num_cols The number of columns for the matrix.
       */
      explicit SparseMatrixFactory(IT_ num_rows, IT_ num_cols):_num_row(num_rows), _num_col(num_cols)
      {
        XASSERTM(num_rows <= max_size, "User tried to generate a Matrix with more than 100000 rows!");
        XASSERTM(num_cols <= max_size, "User tried to generate a Matrix with more than 100000 colums!");
      }

      /**
       * \brief Adds a new matrix entry to the factory.
       *
       * \param[in] i The row-index of the entry to be added.
       * \param[in] j The column-index of the entry to be added.
       * \param[in] a_ij The value of the entry to be added.
       */
      void add(IT_ i, IT_ j, DT_ a_ij)
      {
        XASSERTM(i < _num_row, "User tried to add input out of bounds of the matrix" );
        XASSERTM(j < _num_col, "User tried to add input out of bounds of the matrix");
        memory.emplace(i * max_size + j, a_ij);
      }

      /**
       * \brief Returns number of columns.
       *
       * \returns Number of matrix columns.
       */
      IT_ columns() const
      {
        return _num_col;
      }
      /**
       * \brief Returns number of rows.
       *
       * \returns Number of matrix rows.
       */
      IT_ rows() const
      {
        return _num_row;
      }
      /**
       * \brief Returns product of number columns and number rows.
       *
       * \returns Number of possible matrix entries.
       */
      IT_ size() const
      {
        return _num_col*_num_row;
      }

      /**
       * \brief Returns number of used elements.
       *
       * \returns  Number of  used elements.
       */
      IT_ used_elements() const
      {
        return IT_(memory.size());
      }

      /**
       * \brief Returns CSR matrix constructed from Map.
       *
       * \returns CSR Matrix.
       */
      SparseMatrixCSR<Mem::Main, DT_, IT_> make_csr() const
      {
        //Creates CSR matrix with dimensions and number of NNZ of map.
        FEAT::LAFEM::SparseMatrixCSR< Mem::Main, DT_, IT_ > matrix (_num_row, _num_col,IT_( memory.size()));
        // pointer on the empty arrays of the CSR matrix structure
        IT_* row_ptr = matrix.row_ptr();
        IT_* col_ind = matrix.col_ind();
        DT_* val = matrix.val();

        row_ptr[0] = IT_(0);
        // numberNotZero elements
        IT_ _nnz = IT_(0);
        // current row count
        IT_ _row_cur = IT_(0);

        //Iterates over the map to extract the data and paste in the corresponding arrays.
        //typename std::map< IT_, DT_>::const_iterator  it = memory.begin();
        for(auto it = memory.begin(); it != memory.end();++it)
        {
          // extract column and row from the map key
          IT_ col = it->first % max_size;
          IT_ row = (it->first - col) / max_size;
          // make sure the row pointer is valid even if there exist preceeding empty rows in the matrix
          while (_row_cur <row )
          {
            row_ptr[++_row_cur] = _nnz;
          }
          val[_nnz] = it->second;
          col_ind[_nnz] = col;
          ++_nnz;
        }
        // make sure that the row pointer is valid even if there are empty rows at the end of the matrix
        while (_row_cur < _num_row)
        {
          row_ptr[++_row_cur] = _nnz;
        }

        /*// For testing Purposes: Printing the three arrays
        std::cout<<"Values: ";
        for (IT_ i = 0; i < _nzv; ++i)
        {
          std::cout << val[i] << "  ";
        }
        std::cout << "  \nColum Indizes:" ;

        for (IT_ i = 0; i < _nzv; ++i)
        {
          std::cout << col_ind[i] << "  ";
        }
        std::cout << "  \nRow Pointers:";
        for (IT_ i = 0; i < _num_row + 1; ++i)
        {
          std::cout << row_ptr[i] << "  ";
        }
        std::cout << " " << std::endl;
        */
        return matrix;
      }
    }; // class SparseMatrixFactory
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_SPARSE_MATRIX_FACTORY_HPP
