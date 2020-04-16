// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_AMAVANKA_HPP
#define KERNEL_SOLVER_AMAVANKA_HPP 1

#include <kernel/lafem/null_matrix.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/util/stop_watch.hpp>

#include <array>

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Vanka Matrix type definition helper class
       *
       * This helper class is used to define the matrix type that the AmaVanka smoother will use
       * based on a given system matrix type. Typically, for any "meta" system matrix type, the
       * corresponding Vanka matrix type will be a LAFEM::TupleMatrix, where each sub-matrix is
       * a LAFEM::SparseMatrixBCSR matrix.
       *
       * \author Peter Zajac
       */
      template<typename SystemMatrix_>
#ifndef DOXYGEN
      struct AmaVankaMatrixHelper;
#else // DOXYGEN
      struct AmaVankaMatrixHelper
      {
        /// The type of the Vanka matrix
        typedef ... VankaMatrix;
      };
#endif // DOXYGEN

      template<typename Mem_, typename DT_, typename IT_, int bh_, int bw_>
      struct AmaVankaMatrixHelper<LAFEM::NullMatrix<Mem_, DT_, IT_, bh_, bw_>>
      {
        typedef LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, bh_, bw_> VankaMatrix;
        static constexpr int num_blocks = 1;
      }; // struct AmaVankaMatrixHelper<LAFEM::NullMatrix<...>>

      template<typename Mem_, typename DT_, typename IT_>
      struct AmaVankaMatrixHelper<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>>
      {
        typedef LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> VankaMatrix;
        static constexpr int num_blocks = 1;
      }; // struct AmaVankaMatrixHelper<LAFEM::SparseMatrixCSR<...>>

      template<typename Mem_, typename DT_, typename IT_, int bh_, int bw_>
      struct AmaVankaMatrixHelper<LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, bh_, bw_>>
      {
        typedef LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, bh_, bw_> VankaMatrix;
        static constexpr int num_blocks = 1;
      }; // struct AmaVankaMatrixHelper<LAFEM::SparseMatrixBCSR<...>>

      template<typename First_, typename... Rest_>
      struct AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<First_, Rest_...>>
      {
        typedef LAFEM::TupleMatrixRow<
          typename AmaVankaMatrixHelper<First_>::VankaMatrix,
          typename AmaVankaMatrixHelper<Rest_>::VankaMatrix...> VankaMatrix;
      }; // struct AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<...>>

      template<typename First_>
      struct AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<First_>>
      {
      public:
        typedef LAFEM::TupleMatrixRow<typename AmaVankaMatrixHelper<First_>::VankaMatrix> VankaMatrix;
      }; // struct AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<First_>>

      template<typename FirstRow_, typename... RestRows_>
      struct AmaVankaMatrixHelper<LAFEM::TupleMatrix<FirstRow_, RestRows_...>>
      {
        typedef LAFEM::TupleMatrix<
          typename AmaVankaMatrixHelper<FirstRow_>::VankaMatrix,
          typename AmaVankaMatrixHelper<RestRows_>::VankaMatrix...> VankaMatrix;
        static constexpr int num_blocks = VankaMatrix::num_col_blocks;
      }; // struct AmaVankaMatrixHelper<LAFEM::TupleMatrix<...>>

      template<typename FirstRow_>
      struct AmaVankaMatrixHelper<LAFEM::TupleMatrix<FirstRow_>>
      {
        typedef LAFEM::TupleMatrix<FirstRow_> MatrixType;
        typedef LAFEM::TupleMatrix<typename AmaVankaMatrixHelper<FirstRow_>::VankaMatrix> VankaMatrix;
        static constexpr int num_blocks = VankaMatrix::num_col_blocks;
      }; // struct AmaVankaMatrixHelper<LAFEM::TupleMatrix<FirstRow_>>

      template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
      struct AmaVankaMatrixHelper<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>>
      {
        typedef typename MatrixA_::MemType MemType;
        typedef typename MatrixA_::DataType DataType;
        typedef typename MatrixA_::IndexType IndexType;
        static constexpr int block_size = MatrixA_::BlockWidth;
        static constexpr int num_blocks = 2;

        typedef LAFEM::TupleMatrix<
          LAFEM::TupleMatrixRow<
            LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, block_size, block_size>,
            LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, block_size, 1>>,
          LAFEM::TupleMatrixRow<
            LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 1, block_size>,
            LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 1, 1>>
        > VankaMatrix;
      }; // struct AmaVankaMatrixHelper<LAFEM::SaddlePointMatrix<...>>

      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */

      /**
       * \brief AmaVanka main operations helper class
       *
       * This class contains the implementations of the main components of the AmaVanka smoother including the
       * necessary template wrapping magic to support a wide variety of meta-matrices.
       *
       * \author Peter Zajac
       */
      class AmaVankaCore
      {
      public:
#ifdef DOXYGEN
        /**
         * \brief Computes the main diagonal block sizes of a Vanka matrix.
         *
         * \note
         * This function is only implemented in specialised overloads,
         * i.e. there exists no generic implementation.
         *
         * \param[in] matrix
         * The Vanka matrix for which the block sizes are to be computed.
         *
         * \param[out] sizes
         * A vector that receives the block sizes.
         */
        template<typename Matrix_>
        static void block_sizes(const Matrix_& matrix, std::vector<Index>& sizes)
        {
        }
#endif // DOXYGEN

        /// specialisation for LAFEM::SparseMatrixCSR
        template<typename Mem_, typename DT_, typename IT_>
        static void block_sizes(const LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>&, std::vector<Index>& sizes)
        {
          sizes.push_back(Index(1));
        }

        /// specialisation for LAFEM::SparseMatrixBCSR
        template<typename Mem_, typename DT_, typename IT_, int bs_>
        static void block_sizes(const LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, bs_, bs_>&, std::vector<Index>& sizes)
        {
          sizes.push_back(Index(bs_));
        }

        /// specialisation for LAFEM::TupleMatrix
        template<typename... Rows_>
        static void block_sizes(const LAFEM::TupleMatrix<Rows_...>& matrix, std::vector<Index>& sizes)
        {
          block_sizes_tpl<0>(matrix, sizes);
        }

        /// auxiliary function for TupleMatrix
        template<int col_, typename FirstRow_, typename SecondRow_, typename... RestRows_>
        static void block_sizes_tpl(const LAFEM::TupleMatrix<FirstRow_, SecondRow_, RestRows_...>& matrix, std::vector<Index>& sizes)
        {
          block_sizes(matrix.first().template at<col_>(), sizes);
          block_sizes_tpl<col_+1>(matrix.rest(), sizes);
        }

        /// auxiliary function for TupleMatrix
        template<int col_, typename FirstRow_>
        static void block_sizes_tpl(const LAFEM::TupleMatrix<FirstRow_>& matrix, std::vector<Index>& sizes)
        {
          block_sizes(matrix.first().template at<col_>(), sizes);
        }

        /* ********************************************************************************************************* */

        /**
         * \brief Automatically deducts the Macro-Dofs graphs for all blocks of a given system matrix.
         *
         * This overload "implements" the default case where no automatic deduction is possible, so this
         * function does nothing except for returning \c false. Currently, the only useful overload for
         * this function is the one for the SaddlePointMatrix.
         *
         * \returns \c false
         */
        template<typename Matrix_>
        static bool deduct_macro_dofs(const Matrix_& DOXY(matrix), std::vector<Adjacency::Graph>& DOXY(macro_dofs))
        {
          return false;
        }

        /**
         * \brief Actually deducts the macro-dofs graph for a SaddlePointMatrix
         *
         * \attention
         * This function works only for discontinuous pressure spaces!
         *
         * \param[in] matrix
         * The SaddlePointMatrix for which the macros-dofs are to be deducted.
         *
         *\param[out] macro_dofs
         * A vector that receives the macro-dofs graphs for the velocity and pressure spaces.
         *
         * \returns \c true
         */
        //template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
        //static bool deduct_macro_dofs(const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix,
        template<typename Mem_, typename DT_, typename IT_, int bs_>
        static bool deduct_macro_dofs(
          const LAFEM::SaddlePointMatrix<
            LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, bs_, bs_>,
            LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, bs_, 1>,
            LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, 1, bs_>>& matrix,
          std::vector<Adjacency::Graph>& macro_dofs)
        {
          typedef IT_ IndexType;

          // download matrices B and D to main memory
          LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, bs_, 1> matrix_b;
          LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, 1, bs_> matrix_d;
          matrix_b.convert(matrix.block_b());
          matrix_d.convert(matrix.block_d());

          // fetch matrix dimensions
          const Index num_dofs_p = Index(matrix_d.rows());

          // fetch matrix array
          const IndexType* row_ptr_b = matrix_b.row_ptr();
          const IndexType* col_idx_b = matrix_b.col_ind();
          const IndexType* row_ptr_d = matrix_d.row_ptr();
          const IndexType* col_idx_d = matrix_d.col_ind();

          // PHASE 1: determine pressure macros
          std::map<IndexType, int> map_s;
          std::vector<int> mask(std::size_t(num_dofs_p), 0);
          std::vector<Index> p_ptr, p_idx;
          p_ptr.push_back(Index(0));

          // loop over all pressure nodes
          for(Index i(0); i < num_dofs_p; ++i)
          {
            // is this node already processed?
            if(mask[i] != 0)
              continue;

            // clear the local map
            map_s.clear();

            // okay, loop over all entries of D
            for(IndexType kd = row_ptr_d[i]; kd < row_ptr_d[i+1]; ++kd)
            {
              // fetch the column index of D
              const IndexType col_d = col_idx_d[kd];

              // loop over the row of B
              for(IndexType kb = row_ptr_b[col_d]; kb < row_ptr_b[col_d+1]; ++kb)
              {
                // fetch the column index of B
                const IndexType col_b = col_idx_b[kb];

                // insert into map
                auto ib = map_s.emplace(col_b, 1);
                if(!ib.second)
                {
                  // node already exists, increment counter then
                  ++(ib.first->second);
                }
              }
            }

            // compute the map degree
            int ideg = 0;
            for(auto it = map_s.begin(); it != map_s.end(); ++it)
              ideg = Math::max(ideg, it->second);

            // loop over all map entries with maximum degree
            for(auto it = map_s.begin(); it != map_s.end(); ++it)
            {
              if(ideg == it->second)
              {
                // insert into map
                p_idx.push_back(Index(it->first));
                mask[it->first] = 1;
              }
            }

            // push row
            p_ptr.push_back(Index(p_idx.size()));
          }

          // PHASE 2: determine velocity macros

          // fetch number of pressure macros
          const Index num_macros = Index(p_ptr.size()-1u);

          // clear block arrays
          std::set<IndexType> set_v;
          std::vector<Index> v_ptr, v_idx;
          v_ptr.reserve(num_macros+1u);
          v_ptr.push_back(Index(0));

          // loop over all pressure blocks
          for(IndexType i(0); i < num_macros; ++i)
          {
            // loop over all pressure dofs in the current block
            for(IndexType j = p_ptr[i]; j < p_ptr[i+1]; ++j)
            {
              // get the pressure dof index
              const Index pix = p_idx[j];

              // now loop over the corresponding row of D
              for(IndexType k = row_ptr_d[pix]; k < row_ptr_d[pix+1]; ++k)
              {
                // insert velocity dof into block set
                set_v.insert(col_idx_d[k]);
              }
            }

            // push velocity block
            for(auto it = set_v.begin(); it != set_v.end(); ++it)
            {
              v_idx.push_back(Index(*it));
            }

            // update pointer
            v_ptr.push_back(Index(v_idx.size()));

            // clear local set
            set_v.clear();
          }

          // create graphs
          macro_dofs.resize(std::size_t(2));
          macro_dofs.front() = Adjacency::Graph(num_macros, matrix_b.rows(), Index(v_idx.size()), v_ptr.data(), v_idx.data());
          macro_dofs.back()  = Adjacency::Graph(num_macros, matrix_d.rows(), Index(p_idx.size()), p_ptr.data(), p_idx.data());
          return true;
        }

        /* ********************************************************************************************************* */

#ifdef DOXYGEN
        /**
         * \brief Gathers the local matrix part of a system matrix w.r.t. a macro
         *
         * \note
         * This function is only implemented in specialised overloads,
         * i.e. there exists no generic implementation.
         *
         * \param[in] matrix
         * The system matrix that the local matrix part is to be gathered from
         *
         * \param[inout] local
         * An array representing the local matrix part stored in row-major order.
         *
         * \param[in] stride
         * The stride of the local matrix part \p local.
         *
         * \param[in] macro
         * The index of the macro whose local matrix part is to be extracted.
         *
         * \param[in] macro_dofs
         * A vector that contains the macro-dof graphs for all matrix blocks.
         *
         * \param[in] row_off
         * The offset of the first row in the local matrix part \p local.
         *
         * \param[in] row_block
         * The row-index of the current meta-matrix block.
         *
         * \param[in] col_off
         * The offset of the first column in the local matrix part \p local.
         *
         * \param[in] col_block
         * The column-index of the current meta-matrix block.
         *
         * \returns
         * An index-pair containing the number of rows and columns of the local matrix part, resp.
         */
        template<typename Matrix_, typename DT_>
        static std::pair<Index,Index> gather(const Matrix_& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
        }
#endif // DOXYGEN

        /// specialisation for LAFEM::NullMatrix
        template<typename DT_, typename IT_, int bh_, int bw_>
        static std::pair<Index,Index> gather(const LAFEM::NullMatrix<Mem::Main, DT_, IT_, bh_, bw_>&,
          DT_*, const Index, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index, const Index row_block, const Index, const Index col_block)
        {
          // get graph domain pointer arrays
          const Index* row_dom_ptr = macro_dofs.at(row_block).get_domain_ptr();
          const Index* col_dom_ptr = macro_dofs.at(col_block).get_domain_ptr();

          // get the number of rows/cols to gather
          const Index num_rows = row_dom_ptr[macro+1] - row_dom_ptr[macro];
          const Index num_cols = col_dom_ptr[macro+1] - col_dom_ptr[macro];

          // return the number of gathered rows
          return std::make_pair(num_rows * Index(bh_), num_cols * Index(bw_));
        }

        /// specialisation for LAFEM::SparseMatrixCSR
        template<typename DT_, typename IT_>
        static std::pair<Index,Index> gather(const LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          // get matrix arrays
          const DT_* vals = matrix.val();
          const IT_* row_ptr = matrix.row_ptr();
          const IT_* col_idx = matrix.col_ind();

          // get the dofs for our row/col blocks
          const Adjacency::Graph& row_dofs = macro_dofs.at(row_block);
          const Adjacency::Graph& col_dofs = macro_dofs.at(col_block);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_dofs.get_domain_ptr();
          const Index* col_dom_ptr = col_dofs.get_domain_ptr();

          // get the number of rows/cols to gather
          const Index num_rows = row_dom_ptr[macro+1] - row_dom_ptr[macro];
          const Index num_cols = col_dom_ptr[macro+1] - col_dom_ptr[macro];

          // get graph image index arrays for our macro
          const Index* row_img_idx = &(row_dofs.get_image_idx()[row_dom_ptr[macro]]);
          const Index* col_img_idx = &(col_dofs.get_image_idx()[col_dom_ptr[macro]]);

          // loop over all rows
          for(Index i(0); i < num_rows; ++i)
          {
            // get the index of this row in the big matrix
            const Index irow = row_img_idx[i];

            // initialise loop variable for local columns
            Index j = Index(0);

            // row over the row
            for(IT_ ai = row_ptr[irow]; ai < row_ptr[irow + Index(1)]; ++ai)
            {
              // get the column index
              const Index icol = Index(col_idx[ai]);

              // now search its position in our local matrix
              while((j < num_cols) && (col_img_idx[j] < icol))
              {
                ++j;
              }

              // did we find our local entry?
              if((j < num_cols) && (col_img_idx[j] == icol))
              {
                // copy macro row
                local[(row_off + i) * stride + col_off + j] = vals[ai];
              }
            }
          }

          // return the number of gathered rows
          return std::make_pair(num_rows, num_cols);
        }

        /// specialisation for LAFEM::SparseMatrixBCSR
        template<typename DT_, typename IT_, int bh_, int bw_>
        static std::pair<Index,Index> gather(const LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, bh_, bw_>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          // get matrix arrays
          const Tiny::Matrix<DT_, bh_, bw_>* vals = matrix.val();
          const IT_* row_ptr = matrix.row_ptr();
          const IT_* col_idx = matrix.col_ind();

          // get the dofs for our row/col blocks
          const Adjacency::Graph& row_dofs = macro_dofs.at(row_block);
          const Adjacency::Graph& col_dofs = macro_dofs.at(col_block);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_dofs.get_domain_ptr();
          const Index* col_dom_ptr = col_dofs.get_domain_ptr();

          // get the number of rows/cols to gather
          const Index num_rows = row_dom_ptr[macro+1] - row_dom_ptr[macro];
          const Index num_cols = col_dom_ptr[macro+1] - col_dom_ptr[macro];

          // get graph image index arrays for our macro
          const Index* row_img_idx = &(row_dofs.get_image_idx()[row_dom_ptr[macro]]);
          const Index* col_img_idx = &(col_dofs.get_image_idx()[col_dom_ptr[macro]]);

          // loop over all rows
          for(Index i(0); i < num_rows; ++i)
          {
            // get the index of this row in the big matrix
            const Index irow = row_img_idx[i];

            // initialise loop variable for local columns
            Index j = Index(0);

            // row over the row
            for(IT_ ai = row_ptr[irow]; ai < row_ptr[irow + Index(1)]; ++ai)
            {
              // get the column index
              const Index icol = Index(col_idx[ai]);

              // now search its position in our local matrix
              while((j < num_cols) && (col_img_idx[j] < icol))
              {
                ++j;
              }

              // did we find our local entry?
              if((j < num_cols) && (col_img_idx[j] == icol))
              {
                // found it, so store it in our local matrix
                const Tiny::Matrix<DT_, bh_, bw_>& mv = vals[ai];

                // copy matrix macro
                for(int ii(0); ii < bh_; ++ii)
                {
                  // compute offset within local matrix
                  const Index lro = (row_off + i * Index(bh_) + Index(ii)) * stride + col_off + j * Index(bw_);

                  // copy macro row
                  for(int jj(0); jj < bw_; ++jj)
                  {
                    local[lro + IT_(jj)] = mv[ii][jj];
                  }
                }
              }
            }
          }

          // return the number of gathered rows
          return std::make_pair(num_rows * Index(bh_), num_cols * Index(bw_));
        }

        /// specialisation for LAFEM::TupleMatrixRow
        template<typename DT_, typename First_, typename Second_, typename... Rest_>
        static std::pair<Index,Index> gather(const LAFEM::TupleMatrixRow<First_, Second_, Rest_...>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          const std::pair<Index,Index> nrc_f = AmaVankaCore::gather(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);

          const std::pair<Index,Index> nrc_r = AmaVankaCore::gather(matrix.rest(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off + nrc_f.second, col_block + Index(1));

          XASSERTM(nrc_f.first == nrc_r.first, "block row count mismatch");

          return std::make_pair(nrc_f.first, nrc_f.second + nrc_r.second);
        }

        /// specialisation for LAFEM::TupleMatrixRow (single column)
        template<typename DT_, typename First_>
        static std::pair<Index,Index> gather(const LAFEM::TupleMatrixRow<First_>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          return AmaVankaCore::gather(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);
        }

        /// specialisation for LAFEM::TupleMatrix
        template<typename DT_, typename FirstRow_, typename SecondRow_, typename... RestRows_>
        static std::pair<Index,Index> gather(const LAFEM::TupleMatrix<FirstRow_, SecondRow_, RestRows_...>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          const std::pair<Index,Index> nrc_f = AmaVankaCore::gather(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);

          const std::pair<Index,Index> nrc_r = AmaVankaCore::gather(matrix.rest(), local, stride, macro,
            macro_dofs, row_off + nrc_f.first, row_block + Index(1), col_off, col_block);

          XASSERTM(nrc_f.second == nrc_r.second, "block column count mismatch");

          return std::make_pair(nrc_f.first + nrc_r.first, nrc_f.second);
        }

        /// specialisation for LAFEM::TupleMatrix (single row)
        template<typename DT_, typename FirstRow_>
        static std::pair<Index,Index> gather(const LAFEM::TupleMatrix<FirstRow_>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          return AmaVankaCore::gather(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);
        }

        /// specialisation for LAFEM::SaddlePointMatrix
        template<typename DT_, typename MatrixA_, typename MatrixB_, typename MatrixD_>
        static std::pair<Index,Index> gather(const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          // gather matrix A
          const std::pair<Index,Index> nrc_a = AmaVankaCore::gather(matrix.block_a(), local, stride,
            macro, macro_dofs, row_off, row_block, col_off, col_block);

          // gather matrix B
          const std::pair<Index,Index> nrc_b = AmaVankaCore::gather(matrix.block_b(), local, stride,
            macro, macro_dofs, row_off, row_block, col_off + nrc_a.second, col_block + Index(1));

          // gather matrix D
          const std::pair<Index,Index> nrc_d = AmaVankaCore::gather(matrix.block_d(), local, stride,
            macro, macro_dofs, row_off + nrc_a.first, row_block + Index(1), col_off, col_block);

          XASSERTM(nrc_a.first  == nrc_b.first, "block row count mismatch");
          XASSERTM(nrc_a.second == nrc_d.second, "block column count mismatch");

          return std::make_pair(nrc_a.first + nrc_d.first, nrc_a.second + nrc_b.second);
        }

        /* ********************************************************************************************************* */

#ifdef DOXYGEN
        /**
         * \brief Allocates the Vanka matrix based on the macro-dof graphs
         *
         * \note
         * This function is only implemented in specialised overloads,
         * i.e. there exists no generic implementation.
         *
         * \param[out] matrix
         * The Vanka matrix that is to be allocated.
         *
         * \param[in] dof_macros
         * The vector of all dof-macro graphs.
         *
         * \param[in] macro_dofs
         * The vector of all macro-dof graphs.
         *
         * \param[in] row_block
         * The row-index of the current meta-matrix block.
         *
         * \param[in] col_block
         * The column-index of the current meta-matrix block.
         */
        template<typename Matrix_>
        static void alloc(Matrix_& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
        }
#endif // DOXYGEN

        /// specialisation for LAFEM::SparseMatrixCSR
        template<typename Mem_, typename DT_, typename IT_>
        static void alloc(LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          XASSERT(row_block < Index(dof_macros.size()));
          XASSERT(col_block < Index(macro_dofs.size()));
          Adjacency::Graph graph(Adjacency::RenderType::injectify_sorted, dof_macros.at(row_block), macro_dofs.at(col_block));
          LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_> matrix_main(graph);
          matrix.convert(matrix_main);
        }

        /// specialisation for LAFEM::SparseMatrixBCSR
        template<typename Mem_, typename DT_, typename IT_, int bh_, int bw_>
        static void alloc(LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, bh_, bw_>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          XASSERT(row_block < Index(dof_macros.size()));
          XASSERT(col_block < Index(macro_dofs.size()));
          Adjacency::Graph graph(Adjacency::RenderType::injectify_sorted, dof_macros.at(row_block), macro_dofs.at(col_block));
          LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, bh_, bw_> matrix_main(graph);
          matrix.convert(matrix_main);
        }

        /// specialisation for LAFEM::TupleMatrixRow
        template<typename First_, typename Second_, typename... Rest_>
        static void alloc(LAFEM::TupleMatrixRow<First_, Second_, Rest_...>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          AmaVankaCore::alloc(matrix.first(), dof_macros, macro_dofs, row_block, col_block);
          AmaVankaCore::alloc(matrix.rest(), dof_macros, macro_dofs, row_block, col_block + Index(1));
        }

        /// specialisation for LAFEM::TupleMatrixRow (single column)
        template<typename First_>
        static void alloc(LAFEM::TupleMatrixRow<First_>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          AmaVankaCore::alloc(matrix.first(), dof_macros, macro_dofs, row_block, col_block);
        }

        /// specialisation for LAFEM::TupleMatrix
        template<typename FirstRow_, typename SecondRow_, typename... RestRows_>
        static void alloc(LAFEM::TupleMatrix<FirstRow_, SecondRow_, RestRows_...>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          AmaVankaCore::alloc(matrix.first(), dof_macros, macro_dofs, row_block, col_block);
          AmaVankaCore::alloc(matrix.rest(), dof_macros, macro_dofs, row_block + Index(1), col_block);
        }

        /// specialisation for LAFEM::TupleMatrix (single row)
        template<typename FirstRow_>
        static void alloc(LAFEM::TupleMatrix<FirstRow_>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          AmaVankaCore::alloc(matrix.first(), dof_macros, macro_dofs, row_block, col_block);
        }

        /* ********************************************************************************************************* */

#ifdef DOXYGEN
        /**
         * \brief Scatter-Adds the local matrix part onto a system matrix w.r.t. a macro
         *
         * \note
         * This function is virtually the counter-part of the gather() function.
         *
         * \note
         * This function is only implemented in specialised overloads,
         * i.e. there exists no generic implementation.
         *
         * \param[inout] matrix
         * The system matrix that the local matrix part is to be scatted-added onto
         *
         * \param[in] local
         * An array representing the local matrix part stored in row-major order.
         *
         * \param[in] stride
         * The stride of the local matrix part \p local.
         *
         * \param[in] macro
         * The index of the macro whose local matrix part is to be scattered.
         *
         * \param[in] macro_dofs
         * A vector that contains the macro-dof graphs for all matrix blocks.
         *
         * \param[in] row_off
         * The offset of the first row in the local matrix part \p local.
         *
         * \param[in] row_block
         * The row-index of the current meta-matrix block.
         *
         * \param[in] col_off
         * The offset of the first column in the local matrix part \p local.
         *
         * \param[in] col_block
         * The column-index of the current meta-matrix block.
         *
         * \returns
         * An index-pair containing the number of rows and columns of the local matrix part, resp.
         */
        template<typename DT_, typename Matrix_>
        static std::pair<Index,Index> scatter_add(Matrix_& matrix,
          const DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
        }
#endif // DOXYGEN

        /// specialisation for LAFEM::SparseMatrixCSR
        template<typename DT_, typename IT_>
        static std::pair<Index,Index> scatter_add(LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>& matrix,
          const DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          // get matrix arrays
          DT_* vals = matrix.val();
          const IT_* row_ptr = matrix.row_ptr();
          const IT_* col_idx = matrix.col_ind();

          // get the dofs for our row/col blocks
          const Adjacency::Graph& row_dofs = macro_dofs.at(row_block);
          const Adjacency::Graph& col_dofs = macro_dofs.at(col_block);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_dofs.get_domain_ptr();
          const Index* col_dom_ptr = col_dofs.get_domain_ptr();

          // get the number of rows/cols to gather
          const Index num_rows = row_dom_ptr[macro+1] - row_dom_ptr[macro];
          const Index num_cols = col_dom_ptr[macro+1] - col_dom_ptr[macro];

          // get graph image index arrays for our macro
          const Index* row_img_idx = &(row_dofs.get_image_idx()[row_dom_ptr[macro]]);
          const Index* col_img_idx = &(col_dofs.get_image_idx()[col_dom_ptr[macro]]);

          // loop over all rows
          for(Index i(0); i < num_rows; ++i)
          {
            // get the index of this row in the big matrix
            const Index irow = row_img_idx[i];

            // initialise loop variable for local columns
            Index j = Index(0);

            // row over the row
            for(IT_ ai = row_ptr[irow]; ai < row_ptr[irow + Index(1)]; ++ai)
            {
              // get the column index
              const Index icol = Index(col_idx[ai]);

              // now search its position in our local matrix
              while((j < num_cols) && (col_img_idx[j] < icol))
              {
                ++j;
              }

              // did we find our local entry?
              if((j < num_cols) && (col_img_idx[j] == icol))
              {
                vals[ai] += local[(row_off + i) * stride + col_off + j];
              }
            }
          }

          // return the number of scattered rows
          return std::make_pair(num_rows, num_cols);
        }

        /// specialisation for LAFEM::SparseMatrixBCSR
        template<typename DT_, typename IT_, int bh_, int bw_>
        static std::pair<Index,Index> scatter_add(LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, bh_, bw_>& matrix,
          const DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          // get matrix arrays
          Tiny::Matrix<DT_, bh_, bw_>* vals = matrix.val();
          const IT_* row_ptr = matrix.row_ptr();
          const IT_* col_idx = matrix.col_ind();

          // get the dofs for our row/col blocks
          const Adjacency::Graph& row_dofs = macro_dofs.at(row_block);
          const Adjacency::Graph& col_dofs = macro_dofs.at(col_block);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_dofs.get_domain_ptr();
          const Index* col_dom_ptr = col_dofs.get_domain_ptr();

          // get the number of rows/cols to gather
          const Index num_rows = row_dom_ptr[macro+1] - row_dom_ptr[macro];
          const Index num_cols = col_dom_ptr[macro+1] - col_dom_ptr[macro];

          // get graph image index arrays for our macro
          const Index* row_img_idx = &(row_dofs.get_image_idx()[row_dom_ptr[macro]]);
          const Index* col_img_idx = &(col_dofs.get_image_idx()[col_dom_ptr[macro]]);

          // loop over all rows
          for(Index i(0); i < num_rows; ++i)
          {
            // get the index of this row in the big matrix
            const Index irow = row_img_idx[i];

            // initialise loop variable for local columns
            Index j = Index(0);

            // row over the row
            for(IT_ ai = row_ptr[irow]; ai < row_ptr[irow + Index(1)]; ++ai)
            {
              // get the column index
              const Index icol = Index(col_idx[ai]);

              // now search its position in our local matrix
              while((j < num_cols) && (col_img_idx[j] < icol))
              {
                ++j;
              }

              // did we find our local entry?
              if((j < num_cols) && (col_img_idx[j] == icol))
              {
                // found it, so store it in our local matrix
                Tiny::Matrix<DT_, bh_, bw_>& mv = vals[ai];

                // copy matrix macro
                for(int ii(0); ii < bh_; ++ii)
                {
                  // compute offset within local matrix
                  const Index lro = (row_off + i * Index(bh_) + Index(ii)) * stride + col_off + j * Index(bw_);

                  // copy macro row
                  for(int jj(0); jj < bw_; ++jj)
                  {
                    mv[ii][jj] += local[lro + IT_(jj)];
                  }
                }
              }
            }
          }

          // return the number of scattered rows
          return std::make_pair(num_rows * Index(bh_), num_cols * Index(bw_));
        }

        /// specialisation for LAFEM::TupleMatrixRow
        template<typename DT_, typename First_, typename Second_, typename... Rest_>
        static std::pair<Index,Index> scatter_add(LAFEM::TupleMatrixRow<First_, Second_, Rest_...>& matrix,
          const DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          const std::pair<Index,Index> nrc_f = AmaVankaCore::scatter_add(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);

          const std::pair<Index,Index> nrc_r = AmaVankaCore::scatter_add(matrix.rest(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off + nrc_f.second, col_block + Index(1));

          XASSERTM(nrc_f.first == nrc_r.first, "block row count mismatch");

          return std::make_pair(nrc_f.first, nrc_f.second + nrc_r.second);
        }

        /// specialisation for LAFEM::TupleMatrixRow (single column)
        template<typename DT_, typename First_>
        static std::pair<Index,Index> scatter_add(LAFEM::TupleMatrixRow<First_>& matrix,
          const DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          return AmaVankaCore::scatter_add(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);
        }

        /// specialisation for LAFEM::TupleMatrix
        template<typename DT_, typename FirstRow_, typename SecondRow_, typename... RestRows_>
        static std::pair<Index,Index> scatter_add(LAFEM::TupleMatrix<FirstRow_, SecondRow_, RestRows_...>& matrix,
          const DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          const std::pair<Index,Index> nrc_f = AmaVankaCore::scatter_add(matrix.first(), local, stride,
            macro, macro_dofs, row_off, row_block, col_off, col_block);

          const std::pair<Index,Index> nrc_r = AmaVankaCore::scatter_add(matrix.rest(), local, stride,
            macro, macro_dofs, row_off + nrc_f.first, row_block + Index(1), col_off, col_block);

          XASSERTM(nrc_f.second == nrc_r.second, "block column count mismatch");

          return std::make_pair(nrc_f.first + nrc_r.first, nrc_f.second);
        }

        /// specialisation for LAFEM::TupleMatrix (single row)
        template<typename DT_, typename FirstRow_>
        static std::pair<Index,Index> scatter_add(LAFEM::TupleMatrix<FirstRow_>& matrix,
          const DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          return AmaVankaCore::scatter_add(matrix.first(), local, stride,
            macro, macro_dofs, row_off, row_block, col_off, col_block);
        }

        /* ********************************************************************************************************* */

#ifdef DOXYGEN
        /**
         * \brief Scales the rows of a Vanka matrix.
         *
         * \note
         * This function is only implemented in specialised overloads,
         * i.e. there exists no generic implementation.
         *
         * \param[inout] matrix
         * The Vanka matrix that is to be scaled.
         *
         * \param[in] omega
         * The damping factor.
         *
         * \param[in] dof_macros
         * The vector of dof-macro graphs.
         *
         * \param[in] row_block
         * The row-index of the current meta-matrix block.
         */
        template<typename Matrix_>
        static void scale_rows(Matrix_& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const Index row_block)
        {
        }
#endif // DOXYGEN

        /// specialisation for LAFEM::SparseMatrixCSR
        template<typename DT_, typename IT_>
        static void scale_rows(LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          // get matrix arrays
          DT_* vals = matrix.val();
          const IT_* row_ptr = matrix.row_ptr();

          // get the dofs for our row blocks
          const Adjacency::Graph& row_blocks = dof_macros.at(row_block);
          const Index num_rows = row_blocks.get_num_nodes_domain();
          XASSERT(matrix.rows() == num_rows);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_blocks.get_domain_ptr();
          const Index* row_img_idx = row_blocks.get_image_idx();

          // loop over all rows/dofs
          for(Index i(0); i < num_rows; ++i)
          {
            XASSERT(row_dom_ptr[i] < row_dom_ptr[i+1]);

            // get number of active macros for this DOF
            Index n(0);
            if(macro_mask.empty())
              n = row_dom_ptr[i+1] - row_dom_ptr[i];
            else
            {
              for(Index j(row_dom_ptr[i]); j < row_dom_ptr[i+1]; ++j)
                n += Index(macro_mask[row_img_idx[j]]);
            }

            const DT_ sc = omega / DT_(Math::max(n,Index(1)));
            for(IT_ j(row_ptr[i]); j < row_ptr[i+1]; ++j)
              vals[j] *= sc;
          }
        }

        /// specialisation for LAFEM::SparseMatrixBCSR
        template<typename DT_, typename IT_, int bh_, int bw_>
        static void scale_rows(LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, bh_, bw_>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          // get matrix arrays
          Tiny::Matrix<DT_, bh_, bw_>* vals = matrix.val();
          const IT_* row_ptr = matrix.row_ptr();

          // get the dofs for our row blocks
          const Adjacency::Graph& row_blocks = dof_macros.at(row_block);
          const Index num_rows = row_blocks.get_num_nodes_domain();
          XASSERT(matrix.rows() == num_rows);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_blocks.get_domain_ptr();
          const Index* row_img_idx = row_blocks.get_image_idx();

          // loop over all rows/dofs
          for(Index i(0); i < num_rows; ++i)
          {
            XASSERT(row_dom_ptr[i] < row_dom_ptr[i+1]);

            // get number of active macros for this DOF
            Index n(0);
            if(macro_mask.empty())
              n = row_dom_ptr[i+1] - row_dom_ptr[i];
            else
            {
              for(Index j(row_dom_ptr[i]); j < row_dom_ptr[i+1]; ++j)
                n += Index(macro_mask[row_img_idx[j]]);
            }

            const DT_ sc = omega / DT_(Math::max(n,Index(1)));
            for(IT_ j(row_ptr[i]); j < row_ptr[i+1]; ++j)
              vals[j] *= sc;
          }
        }

        /// specialisation for LAFEM::TupleMatrixRow
        template<typename DT_, typename First_, typename Second_, typename... Rest_>
        static void scale_rows(LAFEM::TupleMatrixRow<First_, Second_, Rest_...>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          AmaVankaCore::scale_rows(matrix.first(), omega, dof_macros, macro_mask, row_block);
          AmaVankaCore::scale_rows(matrix.rest(), omega, dof_macros, macro_mask, row_block);
        }

        /// specialisation for LAFEM::TupleMatrixRow (single column)
        template<typename DT_, typename First_>
        static void scale_rows(LAFEM::TupleMatrixRow<First_>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          AmaVankaCore::scale_rows(matrix.first(), omega, dof_macros, macro_mask, row_block);
        }

        /// specialisation for LAFEM::TupleMatrix
        template<typename DT_, typename FirstRow_, typename SecondRow_, typename... RestRows_>
        static void scale_rows(LAFEM::TupleMatrix<FirstRow_, SecondRow_, RestRows_...>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          AmaVankaCore::scale_rows(matrix.first(), omega, dof_macros, macro_mask, row_block);
          AmaVankaCore::scale_rows(matrix.rest(), omega, dof_macros, macro_mask, row_block + Index(1));
        }

        /// specialisation for LAFEM::TupleMatrix (single row)
        template<typename DT_, typename FirstRow_>
        static void scale_rows(LAFEM::TupleMatrix<FirstRow_>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          AmaVankaCore::scale_rows(matrix.first(), omega, dof_macros, macro_mask, row_block);
        }

        /* ********************************************************************************************************* */

        /**
         * \brief Calculates the stride for the local matrix
         *
         * \param[in] matrix
         * The Vanka matrix for which the local matrix stride is to be computed.
         *
         * \param[in] macro_dofs
         * The vector of macro-dof graphs for which the stride is to be computed.
         *
         * \returns
         * The stride of the local matrix part.
         */
        template<typename Matrix_>
        static Index calc_stride(const Matrix_& matrix, const std::vector<Adjacency::Graph>& macro_dofs)
        {
          const std::size_t num_graphs = macro_dofs.size();
          const Index num_macros = macro_dofs.front().get_num_nodes_domain();

          // compute macro sizes
          std::vector<Index> block_sizes;
          AmaVankaCore::block_sizes(matrix, block_sizes);

          XASSERT(num_graphs == block_sizes.size());

          // loop over all macros
          Index rtn = Index(0);
          for(Index imacro(0); imacro < num_macros; ++imacro)
          {
            Index count = Index(0);

            // loop over all graphs and summarise degrees
            for(std::size_t igraph(0); igraph < num_graphs; ++igraph)
            {
              count += macro_dofs.at(igraph).degree(imacro) * block_sizes.at(igraph);
            }

            rtn = Math::max(rtn, count);
          }

          return rtn;
        }

      }; // struct AmaVankaCore
    } // namespace Intern
    /// \endcond

    /**
     * \brief Additive Macro-wise Matrix-based Vanka preconditioner/smoother
     *
     * This class implements an additive macro-wise Vanka smoother, which stores its
     * pre-computed operator as a sparse matrix, so that each application of the Vanka
     * smoother consists of only one sparse matrix-vector multiplication.
     *
     * This class supports a whole batch of different matrix types:
     * - LAFEM::SparseMatrixCSR
     * - LAFEM::SparseMatrixBCSR
     * - LAFEM::SaddlePointMatrix with sub-blocks of type
     *   - LAFEM::SparseMatrixBCSR
     * - LAFEM::TupleMatrix with sub-blocks of type
     *   - LAFEM::SparseMatrixCSR
     *   - LAFEM::SparseMatrixBCSR
     *
     * \attention
     * In contrast to the Solver::Vanka class, one must specify the macros for the Vanka smoother for each
     * finite element space which is used to discretise the individual blocks of the solution space by
     * using the #push_macro_dofs() function. The only exception is for the LAFEM::SaddlePointMatrix,
     * as for this matrix type, there exists an automatic macro deduction algorithm, which will compute
     * the element-wise macro graphs based on the sparsity patterns of the B and D matrices, resp.
     *
     * <b>Example:</b> Assume that your solution space is the usual velocity-pressure space pair used in
     * a Stokes system, then you need to add the dof-mapping graphs of both the velocity and the pressure
     * space to the smoother in the following way:
     * \code{.cpp}
     * auto vanka = Solver::new_amavanka(system_matrix, system_filter);
     * vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_velocity));
     * vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_pressure));
     * \endcode
     * If your solution space has other components (such as e.g. temperature, density or stress), then you
     * have to add the rendered dof-mappings of those spaces in the correct order as well.
     *
     * \note
     * In the case of a LAFEM::SaddlePointMatrix, the smoother implemented in this class is mathematically
     * equivalent to a Solver::Vanka smoother of type Solver::VankaType::block_full_add.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class AmaVanka :
      public Solver::SolverBase<typename Matrix_::VectorTypeL>
    {
    public:
      /// our base-class
      typedef Solver::SolverBase<typename Matrix_::VectorTypeL> BaseClass;

      /// our memory type
      typedef typename Matrix_::MemType MemType;
      /// our data type
      typedef typename Matrix_::DataType DataType;
      /// our index type
      typedef typename Matrix_::IndexType IndexType;
      /// our vector type
      typedef typename Matrix_::VectorTypeL VectorType;

    protected:
      /// the type of our Vanka matrix
      typedef typename Intern::AmaVankaMatrixHelper<Matrix_>::VankaMatrix VankaMatrixType;
      /// the type of our Vanka matrix in main memory
      typedef typename VankaMatrixType::template ContainerTypeByMDI<Mem::Main, DataType, IndexType> VankaMatrixMainType;

      /// the system matrix
      const Matrix_& _matrix;
      /// the system filter
      const Filter_& _filter;
      /// the Vanka preconditioner matrix
      VankaMatrixType _vanka;
      /// deduct macro dofs automatically?
      bool _auto_macros;
      /// skip singular macros?
      bool _skip_singular;
      /// the DOF-macro graphs
      std::vector<Adjacency::Graph> _macro_dofs, _dof_macros;
      /// the macro mask
      std::vector<int> _macro_mask;
      /// number of steps
      Index _num_steps;
      /// damping parameter
      DataType _omega;
      /// temporary vectors
      VectorType _vec_c, _vec_d;

      // stop watch for symbolic factorisation
      StopWatch watch_init_symbolic;
      // stop watch for numeric factorisation
      StopWatch watch_init_numeric;
      // stop watch for apply time
      StopWatch watch_apply;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The saddle-point system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] omega
       * The damping parameter.
       *
       * \param[in] num_steps
       * The number of smoothing steps to be performed.
       */
      explicit AmaVanka(const Matrix_& matrix, const Filter_& filter,
        const DataType omega = DataType(1), const Index num_steps = Index(1)) :
        _matrix(matrix),
        _filter(filter),
        _vanka(),
        _auto_macros(true),
        _skip_singular(false),
        _num_steps(num_steps),
        _omega(omega)
      {
      }

      /**
       * \brief Clears the macro dofs graphs.
       */
      void clear_macro_dofs()
      {
        this->_macro_dofs.clear();
        _auto_macros = true;
      }

      /**
       * \brief Pushes the dofs-at-macro graph of the next block.
       *
       * \param[in] dofs
       * The dofs-at-macro graph of the block.
       */
      void push_macro_dofs(Adjacency::Graph&& dofs)
      {
        _auto_macros = false;

        // make sure that we do not push more graphs than we have blocks
        if(int(_macro_dofs.size()) >= Intern::AmaVankaMatrixHelper<Matrix_>::num_blocks)
          XABORTM("all macro-dofs graphs have already been added");

        // make sure that the number of macros matches our previous graph
        if(!_macro_dofs.empty() && (_macro_dofs.back().get_num_nodes_domain() != dofs.get_num_nodes_domain()))
          XABORTM("macro count mismatch");

        // push graph into the list
        _macro_dofs.emplace_back(std::forward<Adjacency::Graph>(dofs));

        // sort the dof indices
        _macro_dofs.back().sort_indices();
      }

      /**
       * \brief Sets the number of smoothing steps.
       *
       * \param[in] num_steps
       * The number of smoothing steps to be performed.
       */
      void set_num_steps(Index num_steps)
      {
        XASSERT(num_steps > Index(0));
        this->_num_steps = num_steps;
      }

      /**
       * \brief Sets the damping parameter omega.
       *
       * \param[in] omega
       * The damping parameter.
       */
      void set_omega(DataType omega)
      {
        XASSERT(omega > DataType(0));
        this->_omega = omega;
      }

      /**
       * \brief Sets whether singular macros are to be skipped.
       *
       * \param[in] skip_sing
       * Specifies whether singular macros are to be skipped.
       */
      void set_skip_singular(bool skip_sing)
      {
        this->_skip_singular = skip_sing;
      }

      /**
       * \brief Returns the total number of bytes currently allocated in this object.
       */
      std::size_t bytes() const
      {
        std::size_t s = _vanka.bytes();
        for(const auto& g : _macro_dofs)
          s += sizeof(Index) * std::size_t(g.get_num_nodes_domain() + g.get_num_indices());
        for(const auto& g : _dof_macros)
          s += sizeof(Index) * std::size_t(g.get_num_nodes_domain() + g.get_num_indices());
        s += _vec_c.bytes();
        s += _vec_d.bytes();
        return s;
      }

      /**
       * \brief Returns the total data size used by the AmaVanka smoother.
       *
       * \returns
       * The total data size, i.e. the total number of floating point values used in the factorisation.
       */
      std::size_t data_size() const
      {
        return std::size_t(_vanka.template used_elements<LAFEM::Perspective::pod>());
      }

      /**
       * \brief Resets the internal stop watches for time measurement.
       */
      void reset_timings()
      {
        watch_init_symbolic.reset();
        watch_init_numeric.reset();
        watch_apply.reset();
      }

      /**
       * \brief Returns the total accumulated time for symbolic initialisation.
       */
      double time_init_symbolic() const
      {
        return watch_init_symbolic.elapsed();
      }

      /**
       * \brief Returns the total accumulated time for numeric initialisation.
       */
      double time_init_numeric() const
      {
        return watch_init_numeric.elapsed();
      }

      /**
       * \brief Returns the total accumulated time for the solver application.
       */
      double time_apply() const
      {
        return watch_apply.elapsed();
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "AmaVanka";
      }

      /// Performs symbolic factorisation
      virtual void init_symbolic() override
      {
        watch_init_symbolic.start();

        BaseClass::init_symbolic();

        // we also need two vectors if we have to perform multiple steps
        if(this->_num_steps > Index(1))
        {
          this->_vec_c = this->_matrix.create_vector_l();
          this->_vec_d = this->_matrix.create_vector_l();
        }

        // automatically deduct macros?
        if(this->_auto_macros)
        {
          // try to deduct macros by pressure matrices in SaddlePointMatrix
          if(!Intern::AmaVankaCore::deduct_macro_dofs(this->_matrix, this->_macro_dofs))
            XABORTM("Cannot auto-deduct macros for this matrix type");
        }

        // make sure we have one macro-dof graph for each matrix block
        XASSERTM(int(_macro_dofs.size()) == Intern::AmaVankaMatrixHelper<Matrix_>::num_blocks,
           "invalid number of macro-dof graphs; did you push all of them?");

        // compute dof-macro graphs by transposing
        this->_dof_macros.resize(this->_macro_dofs.size());
        for(std::size_t i(0); i < this->_macro_dofs.size(); ++i)
        {
          // ensure that we have the same number of macros in all graphs
          XASSERT(this->_macro_dofs.at(i).get_num_nodes_domain() == this->_macro_dofs.front().get_num_nodes_domain());

          // transpose macro-dofs graph
          this->_dof_macros.at(i) = Adjacency::Graph(Adjacency::RenderType::transpose, this->_macro_dofs.at(i));
        }

        // allocate macro skip mask?
        if(this->_skip_singular)
          this->_macro_mask.resize(this->_macro_dofs.front().get_num_nodes_domain(), 0);

        // allocate Vanka matrix
        VankaMatrixMainType vanka_main;
        Solver::Intern::AmaVankaCore::alloc(vanka_main, this->_dof_macros, this->_macro_dofs, Index(0), Index(0));
        this->_vanka.convert(vanka_main);

        watch_init_symbolic.stop();
      }

      /// Releases the symbolic factorisation data
      virtual void done_symbolic() override
      {
        this->_vanka.clear();
        this->_macro_mask.clear();
        this->_dof_macros.clear();
        if(this->_auto_macros)
          this->_macro_dofs.clear();
        if(this->_num_steps > Index(1))
        {
          this->_vec_d.clear();
          this->_vec_c.clear();
        }
        BaseClass::done_symbolic();
      }

      /// Performs numeric factorisation
      virtual void init_numeric() override
      {
        const DataType eps = Math::eps<DataType>();

        watch_init_numeric.start();
        BaseClass::init_numeric();

        // get maximum macro size
        const Index num_macros = Index(this->_macro_dofs.front().get_num_nodes_domain());
        const Index stride = Intern::AmaVankaCore::calc_stride(this->_vanka, this->_macro_dofs);

        // allocate arrays for local matrix
        std::vector<DataType> vec_local(stride*stride, DataType(0)), vec_local_t(stride*stride, DataType(0));
        std::vector<Index> vec_pivot(stride);
        DataType* local = vec_local.data();
        DataType* local_t = vec_local_t.data();
        Index* pivot = vec_pivot.data();

        // convert matrix to main memory
        typename Matrix_::template ContainerTypeByMDI<Mem::Main, DataType, IndexType> matrix_main;
        matrix_main.convert(this->_matrix);

        // get Vanka matrix in main mem and format
        VankaMatrixMainType vanka_main;
        vanka_main.convert(this->_vanka);
        vanka_main.format();

        // loop over all macros
        for(Index imacro(0); imacro < num_macros; ++imacro)
        {
          // gather local matrix
          const std::pair<Index,Index> nrc = Intern::AmaVankaCore::gather(matrix_main,
            local, stride, imacro, this->_macro_dofs, Index(0), Index(0), Index(0), Index(0));

          // make sure we have gathered a square matrix
          XASSERTM(nrc.first == nrc.second, "local matrix is not square");

          // do we check for singular macros?
          if(this->_skip_singular)
          {
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

            // invert local matrix
            Math::invert_matrix(nrc.first, stride, local, pivot);

            // compute (squared) Frobenius norm of (I - A*A^{-1})
            DataType norm = DataType(0);
            for(Index i(0); i < nrc.first; ++i)
            {
              for(Index j(0); j < nrc.first; ++j)
              {
                DataType xij = DataType(i == j ? 1 : 0);
                for(Index k(0); k < nrc.first; ++k)
                  xij -= local_t[i*stride+k] * local[k*stride+j]; // A_ik * (A^{-1})_kj
                norm += xij * xij;
              }
            }

            // is the matrix block singular?
            // Note: we check for !(norm < eps) instead of (norm >= eps),
            // because the latter one evaluates to false if norm is NaN,
            // which would result in a false negative
            const bool singular = !(norm < eps);

            // set macro regularity mask
            this->_macro_mask[imacro] = (singular ? 0 : 1);

            // scatter local matrix
            if(!singular)
            {
              Intern::AmaVankaCore::scatter_add(vanka_main, local, stride, imacro, this->_macro_dofs,
                Index(0), Index(0), Index(0), Index(0));
            }
          }
          else // no singularity check
          {
            // invert local matrix
            Math::invert_matrix(nrc.first, stride, local, pivot);

            // scatter local matrix
            Intern::AmaVankaCore::scatter_add(vanka_main, local, stride, imacro, this->_macro_dofs,
              Index(0), Index(0), Index(0), Index(0));
          }

          // reformat local matrix
          for(Index i(0); i < nrc.first; ++i)
            for(Index j(0); j < nrc.second; ++j)
              local[i*stride+j] = DataType(0);
        }

        // scale rows of Vanka matrix
        Solver::Intern::AmaVankaCore::scale_rows(vanka_main, this->_omega, this->_dof_macros, this->_macro_mask, Index(0));

        // convert back
        /// \todo use copy maybe?
        this->_vanka.convert(vanka_main);

        watch_init_numeric.stop();
      }

      /// applies the preconditioner
      virtual Status apply(VectorType& vec_x, const VectorType& vec_b) override
      {
        watch_apply.start();

        // first step
        this->_vanka.apply(vec_x, vec_b);
        this->_filter.filter_cor(vec_x);

        // steps 2, ..., n   (if any)
        for(Index step(1); step < _num_steps; ++step)
        {
          // compute defect
          this->_matrix.apply(this->_vec_d, vec_x, vec_b, -DataType(1));
          // filter defect
          this->_filter.filter_def(this->_vec_d);
          // apply Vanka matrix
          this->_vanka.apply(this->_vec_c, this->_vec_d);
          // filter correct
          this->_filter.filter_cor(this->_vec_c);
          // update solution
          vec_x.axpy(this->_vec_c, vec_x);
        }

        watch_apply.stop();

        return Status::success;
      }
    }; // class AmaVanka

    /**
     * \brief Creates a new AmaVanka smoother object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter
     *
     * \param[in] omega
     * The damping parameter.
     *
     * \param[in] num_steps
     * The number of Vanka iterations to be performed.
     *
     * \returns
     * A shared pointer to a new AmaVanka object.
     */
    template<typename Matrix_, typename Filter_>
    std::shared_ptr<AmaVanka<Matrix_, Filter_>> new_amavanka(const Matrix_& matrix, const Filter_& filter,
      typename Matrix_::DataType omega = typename Matrix_::DataType(1), Index num_steps = Index(1))
    {
      return std::make_shared<AmaVanka<Matrix_, Filter_>>(matrix, filter, omega, num_steps);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_AMAVANKA_HPP
