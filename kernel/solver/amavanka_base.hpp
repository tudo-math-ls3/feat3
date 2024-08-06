
// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_AMAVANKA_BASE_HPP
#define KERNEL_SOLVER_AMAVANKA_BASE_HPP 1

#include <kernel/lafem/null_matrix.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>

#ifdef __CUDACC__
#include <cuda/std/utility>
#endif

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
       * \author Peter Zajac and Maximilian Esser
       */
      template<typename SystemMatrix_>
#ifndef DOXYGEN
      struct AmaVankaMatrixHelper;
#else // DOXYGEN
      struct AmaVankaMatrixHelper
      {
        /// The type of the Vanka matrix
        typedef ... VankaMatrix;
        /// The number of blocks
        static constexpr int num_blocks = 0;

        /**
         * \brief Recursive function for filling the matrix data part of the wrapper.
         *
         * The function gathers from a square meta matrix all data pointers and required information
         * into n_*n_ arrays in a recursive fashion.
         *
         * \tparam n_ Integer with the meta row_size of the system matrix.
         *
         * \param[out] wrap
         *  Reference to the matrix wrapper to be filled.
         *
         * \param[in] matrix
         *  Matrix to be used to fill the wrapper.
         *
         * \param[in] row
         * Meta row to be filled.
         *
         * \param[in] col
         * Meta col to be filled.
         */
        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<typename SystemMatrix_::DataType, typename SystemMatrix_::IndexType, n_>& wrap, const SystemMatrix_& matrix, int row, int col);

        /**
         * \brief Recursive function for filling the column info part of the wrapper.
         *
         * This function gathers all matrix data, that is redudant over multiple rows,
         * for example blocksizes due to the squareness of the meta matrix.
         *
         * \tparam n_ Integer with the meta row_size of the system matrix.
         *
         * \param[out] wrap
         *  Reference to the matrix wrapper to be filled.
         *
         * \param[in] matrix
         *  Matrix to be used to fill the wrapper.
         *
         * \param[in] row
         * Meta row to be filled.
         *
         * \param[in] col
         * Meta col to be filled.
         */
        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<typename SystemMatrix_::DataType, typename SystemMatrix_::IndexType, n_>& wrap, const SystemMatrix_& matrix, int row, int col);

        /**
         * \brief Compares two matrices by caluclating the normalized Froebnius-Norm of the difference.
         *
         * \param[in] matrix_a
         *  First matrix to compare.
         *
         * \param[in] matrix_b
         *  Second matrix to compare.
         *
         * \returns If matrix are the same.
         */
        static bool compare(const SystemMatrix_& matrix_a, const SystemMatrix_& matrix_b);
      };
#endif // DOXYGEN

      /**
       * \brief TupleMatrixWrapper holding all relavant information
       *
       * This wrapper datastruct is intended to hold all pointer and size information
       * of a square tuple system Consisting of CSR typed matrices or zero matrices.
       * This means that the saved ptr have to be valid through the lifetime of
       * the object, for which in turn the user is responsible.
       */
      template<typename DT_, typename IT_, int n_>
      struct CSRTupleMatrixWrapper
      {
        typedef DT_ DataType;
        typedef IT_ IndexType;
        static constexpr int n = n_;
        //arrays hold pointer in an row major order, i.e. (1,1), (1,2) ...
        DataType* data_arrays[n_*n_];
        IndexType* col_arrays[n_*n_];
        IndexType* row_arrays[n_*n_];
        IndexType used_elements[n_*n_];
        //since the layout is symmetric, we only need to save the row and colum counts of the first
        //meta row, which is entirely determined by the rows of the first entry and the columns
        //of the first meta row
        //required??
        IndexType tensor_counts[n_+1];
        //since the overall layout has to be symmteric and quadratic, we only have to save the
        //the block dimension of the first column
        int blocksizes[n_+1];

        CUDA_HOST FEAT::String print()
        {
          FEAT::String str("--------------------------\n\nMatrixWrapper: \n");
          for(int i = 0; i < n_; ++i)
          {
            for(int j = 0; j < n_; ++j)
            {
              str += FEAT::String("Block (") + stringify(i) + ", " + stringify(j) + ")\n";
              str += FEAT::String("Tensor size (") + stringify(tensor_counts[Index(i)+std::min(Index(i),Index(1))])
                    + ", " + stringify(tensor_counts[j+1]) + "), Blocksize ("
                    + stringify(blocksizes[Index(i)+std::min(Index(i),Index(1))]) + ", " + stringify(blocksizes[j+1])
                    + ")\n";
            }
          }
          str += String("\n--------------------------\n");
          return str;
        }

        // the mapping from the actual quadratic indices to the n_+1 sized row is simply
        // (i,j) -> (i+min(i,1), j+1)
        CUDA_HOST std::pair<IndexType, IndexType> get_tensor_size(Index i, Index j) const
        {
          return std::make_pair(tensor_counts[i+std::min(i,Index(1))], tensor_counts[j+1]);
        }

        CUDA_HOST std::pair<IndexType, IndexType> get_block_sizes(Index i, Index j) const
        {
          return std::make_pair(blocksizes[i+std::min(i,Index(1))], blocksizes[j+1]);
        }

        CUDA_HOST Index get_all_rows_size() const
        {
          Index val = tensor_counts[0];
          for(int i = 1; i < n_; ++i)
            val += tensor_counts[i+1];
          return val;
        }
      };

      template<typename DT_, typename IT_, int bh_, int bw_>
      struct AmaVankaMatrixHelper<LAFEM::NullMatrix<DT_, IT_, bh_, bw_>>
      {
        typedef LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_> VankaMatrix;
        static constexpr int num_blocks = 1;

        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<DT_, IT_, n_>& wrap, const LAFEM::NullMatrix<DT_, IT_, bh_, bw_>&, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          wrap.data_arrays[n_*row + col] = nullptr;
          wrap.col_arrays[n_*row + col] = nullptr;
          wrap.row_arrays[n_*row + col] = nullptr;
          wrap.used_elements[n_*row + col] = IT_(0);
        }

        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<DT_, IT_, n_>& wrap, const LAFEM::NullMatrix<DT_, IT_, bh_, bw_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          ASSERTM(row == 0, "Row counter has to be zero");
          #ifndef DEBUG
          (void)row;
          #endif
          if(col == 0)
          {
            wrap.tensor_counts[0] = matrix.rows();
            wrap.blocksizes[0] = bh_;
          }
          wrap.tensor_counts[col+1] = matrix.columns();
          wrap.blocksizes[col+1] = bw_;
        }

        static void fill_row_helper(const LAFEM::NullMatrix<DT_, IT_, bh_, bw_>& matrix, Index* accum_row_idx, Index* accum_row_ctr, int row_block, int curr_row)
        {
          // fill with list from 0 to num_rows in raw perspective
          std::generate(accum_row_idx+curr_row, accum_row_idx+curr_row + matrix.template rows<LAFEM::Perspective::pod>(), [n=0]() mutable{return n++;});
          // fill with constant value
          std::fill(accum_row_ctr+curr_row, accum_row_ctr+curr_row + matrix.template rows<LAFEM::Perspective::pod>(), Index(row_block));
        }

        static bool compare(const LAFEM::NullMatrix<DT_, IT_, bh_, bw_>& matrix_a, const LAFEM::NullMatrix<DT_, IT_, bh_, bw_>& matrix_b)
        {
          return (matrix_a.rows() == matrix_b.rows()) && (matrix_b.columns() == matrix_b.columns());
        }
      }; // struct AmaVankaMatrixHelper<LAFEM::NullMatrix<...>>

      template<typename DT_, typename IT_>
      struct AmaVankaMatrixHelper<LAFEM::SparseMatrixCSR<DT_, IT_>>
      {
        typedef LAFEM::SparseMatrixCSR<DT_, IT_> VankaMatrix;
        static constexpr int num_blocks = 1;

        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<DT_, IT_, n_>& wrap, const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          wrap.data_arrays[n_*row + col] = const_cast<DT_*>(matrix.val());
          wrap.col_arrays[n_*row + col] = const_cast<IT_*>(matrix.col_ind());
          wrap.row_arrays[n_*row + col] = const_cast<IT_*>(matrix.row_ptr());
          wrap.used_elements[n_*row + col] = IT_(matrix.used_elements());
        }

        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<DT_, IT_, n_>& wrap, const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          ASSERTM(row == 0, "Row counter has to be zero");
          #ifndef DEBUG
          (void)row;
          #endif
          if(col == 0)
          {
            wrap.tensor_counts[0] = IT_(matrix.rows());
            wrap.blocksizes[0] = 1;
          }
          wrap.tensor_counts[col+1] = IT_(matrix.columns());
          wrap.blocksizes[col+1] = 1;
        }

        static void fill_row_helper(const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix, Index* accum_row_idx, Index* accum_row_ctr, int row_block, int curr_row)
        {
          ASSERTM((accum_row_idx+curr_row + matrix.rows()) == std::find(accum_row_idx+curr_row, accum_row_idx+curr_row + matrix.rows(), ~Index(0)), "Row array was already written to");
          ASSERTM((accum_row_ctr+curr_row + matrix.rows()) == std::find(accum_row_ctr+curr_row, accum_row_ctr+curr_row + matrix.rows(), ~Index(0)), "Row array was already written to");
          // fill with list from 0 to num_rows in raw perspective
          std::generate(accum_row_idx+curr_row, accum_row_idx+curr_row + matrix.rows(), [n=0]() mutable{return n++;});
          // fill with constant value
          std::fill(accum_row_ctr+curr_row, accum_row_ctr+curr_row + matrix.rows(), Index(row_block));
        }

        static bool compare(const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix_a, const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix_b)
        {
          auto tmp = matrix_a.clone(LAFEM::CloneMode::Weak);
          tmp.axpy(matrix_a, matrix_b, DT_(-1));
          const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));
          return (tmp.norm_frobenius()/(Math::sqrt(DT_(tmp.size()))) <= tol);
        }
      }; // struct AmaVankaMatrixHelper<LAFEM::SparseMatrixCSR<...>>

      template<typename DT_, typename IT_, int bh_, int bw_>
      struct AmaVankaMatrixHelper<LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>>
      {
        typedef LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_> VankaMatrix;
        static constexpr int num_blocks = 1;

        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<DT_, IT_, n_>& wrap, const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          wrap.data_arrays[n_*row + col] = const_cast<DT_*>(matrix.template val<LAFEM::Perspective::pod>());
          wrap.col_arrays[n_*row + col] = const_cast<IT_*>(matrix.col_ind());
          wrap.row_arrays[n_*row + col] = const_cast<IT_*>(matrix.row_ptr());
          wrap.used_elements[n_*row + col] = IT_(matrix.used_elements());
        }

        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<DT_, IT_, n_>& wrap, const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          ASSERTM(row == 0, "Row counter has to be zero");
          #ifndef DEBUG
          (void)row;
          #endif
          if(col == 0)
          {
            wrap.tensor_counts[0] = IT_(matrix.rows());
            wrap.blocksizes[0] = bh_;
          }
          wrap.tensor_counts[col+1] = IT_(matrix.columns());
          wrap.blocksizes[col+1] = bw_;
        }

        static void fill_row_helper(const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix, Index* accum_row_idx, Index* accum_row_ctr, int row_block, int curr_row)
        {
          ASSERTM((accum_row_idx+curr_row + matrix.template rows<LAFEM::Perspective::pod>()) == std::find(accum_row_idx+curr_row, accum_row_idx+curr_row + matrix.template rows<LAFEM::Perspective::pod>(), ~Index(0)), "Row array was already written to");
          ASSERTM((accum_row_ctr+curr_row + matrix.template rows<LAFEM::Perspective::pod>()) == std::find(accum_row_ctr+curr_row, accum_row_ctr+curr_row + matrix.template rows<LAFEM::Perspective::pod>(), ~Index(0)), "Row array was already written to");
          // fill with list from 0 to num_rows in raw perspective
          std::generate(accum_row_idx+curr_row, accum_row_idx+curr_row + matrix.template rows<LAFEM::Perspective::pod>(), [n=0]() mutable{return n++;});
          // fill with constant value
          std::fill(accum_row_ctr+curr_row, accum_row_ctr+curr_row + matrix.template rows<LAFEM::Perspective::pod>(), Index(row_block));
        }

        static bool compare(const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix_a, const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix_b)
        {
          auto tmp = matrix_a.clone(LAFEM::CloneMode::Weak);
          tmp.axpy(matrix_a, matrix_b, DT_(-1));
          const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));
          return (tmp.norm_frobenius()/(Math::sqrt(DT_(tmp.size()))) <= tol);
        }
      }; // struct AmaVankaMatrixHelper<LAFEM::SparseMatrixBCSR<...>>

      template<typename First_, typename... Rest_>
      struct AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<First_, Rest_...>>
      {
        typedef LAFEM::TupleMatrixRow<
          typename AmaVankaMatrixHelper<First_>::VankaMatrix,
          typename AmaVankaMatrixHelper<Rest_>::VankaMatrix...> VankaMatrix;

        typedef typename First_::DataType DataType;
        typedef typename First_::IndexType IndexType;
        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::TupleMatrixRow<First_, Rest_...>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          AmaVankaMatrixHelper<First_>::fill_ptr_wrapper(wrap, matrix.first(), row, col);
          AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<Rest_...>>::fill_ptr_wrapper(wrap, matrix.rest(), row, col+1);
        }

        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::TupleMatrixRow<First_, Rest_...>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          ASSERTM(row == 0, "Row counter has to be zero");
          AmaVankaMatrixHelper<First_>::fill_column_info(wrap, matrix.first(), row, col);
          AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<Rest_...>>::fill_column_info(wrap, matrix.rest(), row, col+1);
        }

        static void fill_row_helper(const LAFEM::TupleMatrixRow<First_, Rest_...>& matrix, Index* accum_row_idx, Index* accum_row_ctr, int row_block, int curr_row)
        {
          AmaVankaMatrixHelper<First_>::fill_row_helper(matrix.first(), accum_row_idx, accum_row_ctr, row_block, curr_row);
        }

        static bool compare(const LAFEM::TupleMatrixRow<First_, Rest_...>& matrix_a, const LAFEM::TupleMatrixRow<First_, Rest_...>& matrix_b)
        {
          bool compare = AmaVankaMatrixHelper<First_>::compare(matrix_a.first(), matrix_b.first());
          compare &= AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<Rest_...>>::compare(matrix_a.rest(), matrix_b.rest());
          return compare;
        }
      }; // struct AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<...>>

      template<typename First_>
      struct AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<First_>>
      {
      public:
        typedef LAFEM::TupleMatrixRow<typename AmaVankaMatrixHelper<First_>::VankaMatrix> VankaMatrix;

        typedef typename First_::DataType DataType;
        typedef typename First_::IndexType IndexType;
        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::TupleMatrixRow<First_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          AmaVankaMatrixHelper<First_>::fill_ptr_wrapper(wrap, matrix.first(), row, col);
        }

        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::TupleMatrixRow<First_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          ASSERTM(row == 0, "Row counter has to be zero");
          AmaVankaMatrixHelper<First_>::fill_column_info(wrap, matrix.first(), row, col);
        }

        static void fill_row_helper(const LAFEM::TupleMatrixRow<First_>& matrix, Index* accum_row_idx, Index* accum_row_ctr, int row_block, int curr_row)
        {
          AmaVankaMatrixHelper<First_>::fill_row_helper(matrix.first(), accum_row_idx, accum_row_ctr, row_block, curr_row);
        }

        static bool compare(const LAFEM::TupleMatrixRow<First_>& matrix_a, const LAFEM::TupleMatrixRow<First_>& matrix_b)
        {
          bool compare = AmaVankaMatrixHelper<First_>::compare(matrix_a.first(), matrix_b.first());
          return compare;
        }
      }; // struct AmaVankaMatrixHelper<LAFEM::TupleMatrixRow<First_>>

      template<typename FirstRow_, typename... RestRows_>
      struct AmaVankaMatrixHelper<LAFEM::TupleMatrix<FirstRow_, RestRows_...>>
      {
        typedef LAFEM::TupleMatrix<
          typename AmaVankaMatrixHelper<FirstRow_>::VankaMatrix,
          typename AmaVankaMatrixHelper<RestRows_>::VankaMatrix...> VankaMatrix;
        static constexpr int num_blocks = VankaMatrix::num_col_blocks;

        typedef typename FirstRow_::DataType DataType;
        typedef typename FirstRow_::IndexType IndexType;
        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::TupleMatrix<FirstRow_, RestRows_...>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          AmaVankaMatrixHelper<FirstRow_>::fill_ptr_wrapper(wrap, matrix.first(), row, col);
          AmaVankaMatrixHelper<LAFEM::TupleMatrix<RestRows_...>>::fill_ptr_wrapper(wrap, matrix.rest(), row+1, col);
        }

        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::TupleMatrix<FirstRow_, RestRows_...>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          ASSERTM(row == 0, "Row counter has to be zero");
          AmaVankaMatrixHelper<FirstRow_>::fill_column_info(wrap, matrix.first(), row, col);
        }

        static void fill_row_helper(const LAFEM::TupleMatrix<FirstRow_, RestRows_...>& matrix, Index* accum_row_idx, Index* accum_row_ctr, int row_block, int curr_row)
        {
          AmaVankaMatrixHelper<FirstRow_>::fill_row_helper(matrix.first(), accum_row_idx, accum_row_ctr, row_block, curr_row);
          AmaVankaMatrixHelper<LAFEM::TupleMatrix<RestRows_...>>::fill_row_helper(matrix.rest(), accum_row_idx, accum_row_ctr, row_block+1, curr_row+matrix.first().template rows<LAFEM::Perspective::pod>());
        }

        static bool compare(const LAFEM::TupleMatrix<FirstRow_, RestRows_...>& matrix_a, const LAFEM::TupleMatrix<FirstRow_, RestRows_...>& matrix_b)
        {
          bool compare = AmaVankaMatrixHelper<FirstRow_>::compare(matrix_a.first(), matrix_b.first());
          compare &= AmaVankaMatrixHelper<LAFEM::TupleMatrix<RestRows_...>>::compare(matrix_a.rest(), matrix_b.rest());
          return compare;
        }
      }; // struct AmaVankaMatrixHelper<LAFEM::TupleMatrix<...>>

      template<typename FirstRow_>
      struct AmaVankaMatrixHelper<LAFEM::TupleMatrix<FirstRow_>>
      {
        typedef LAFEM::TupleMatrix<FirstRow_> MatrixType;
        typedef LAFEM::TupleMatrix<typename AmaVankaMatrixHelper<FirstRow_>::VankaMatrix> VankaMatrix;
        static constexpr int num_blocks = VankaMatrix::num_col_blocks;

        typedef typename FirstRow_::DataType DataType;
        typedef typename FirstRow_::IndexType IndexType;
        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::TupleMatrix<FirstRow_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          AmaVankaMatrixHelper<FirstRow_>::fill_ptr_wrapper(wrap, matrix.first(), row, col);
        }

        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::TupleMatrix<FirstRow_>& matrix, int row, int col)
        {
          ASSERTM((n_ > row) && (n_ > col), "MatrixWrapper is too small");
          ASSERTM(row == 0, "Row counter has to be zero");
          AmaVankaMatrixHelper<FirstRow_>::fill_column_info(wrap, matrix.first(), row, col);
        }

        static void fill_row_helper(const LAFEM::TupleMatrix<FirstRow_>& matrix, Index* accum_row_idx, Index* accum_row_ctr, int row_block, int curr_row)
        {
          AmaVankaMatrixHelper<FirstRow_>::fill_row_helper(matrix.first(), accum_row_idx, accum_row_ctr, row_block, curr_row);
        }

        static bool compare(const LAFEM::TupleMatrix<FirstRow_>& matrix_a, const LAFEM::TupleMatrix<FirstRow_>& matrix_b)
        {
          bool compare = AmaVankaMatrixHelper<FirstRow_>::compare(matrix_a.first(), matrix_b.first());
          return compare;
        }
      }; // struct AmaVankaMatrixHelper<LAFEM::TupleMatrix<FirstRow_>>

      template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
      struct AmaVankaMatrixHelper<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>>
      {
        typedef typename MatrixA_::DataType DataType;
        typedef typename MatrixA_::IndexType IndexType;
        static constexpr int block_size = MatrixA_::BlockWidth;
        static constexpr int num_blocks = 2;

        typedef LAFEM::TupleMatrix<
          LAFEM::TupleMatrixRow<
            LAFEM::SparseMatrixBCSR<DataType, IndexType, block_size, block_size>,
            LAFEM::SparseMatrixBCSR<DataType, IndexType, block_size, 1>>,
          LAFEM::TupleMatrixRow<
            LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, block_size>,
            LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, 1>>
        > VankaMatrix;

        template<int n_>
        static void fill_ptr_wrapper(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix, int row, int col)
        {
          ASSERTM((n_+1 > row) && (n_+1 > col), "MatrixWrapper is too small");
          AmaVankaMatrixHelper<MatrixA_>::fill_ptr_wrapper(wrap, matrix.template at<0,0>(), row, col);
          AmaVankaMatrixHelper<MatrixB_>::fill_ptr_wrapper(wrap, matrix.template at<0,1>(), row, col+1);
          AmaVankaMatrixHelper<MatrixD_>::fill_ptr_wrapper(wrap, matrix.template at<1,0>(), row+1, col);
          AmaVankaMatrixHelper<LAFEM::NullMatrix<DataType, IndexType>>::fill_ptr_wrapper(wrap, LAFEM::NullMatrix<DataType, IndexType>(matrix.template at<1,0>().rows(), matrix.template at<0,1>().columns()), row+1, col+1);
        }

        template<int n_>
        static void fill_column_info(CSRTupleMatrixWrapper<DataType, IndexType, n_>& wrap, const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix, int row, int col)
        {
          ASSERTM((n_+1 > row) && (n_+1 > col), "MatrixWrapper is too small");
          ASSERTM(row == 0, "Row counter has to be zero");
          AmaVankaMatrixHelper<MatrixA_>::fill_column_info(wrap, matrix.template at<0,0>(), row, col);
          AmaVankaMatrixHelper<MatrixB_>::fill_column_info(wrap, matrix.template at<0,1>(), row, col+1);
        }

        static void fill_row_helper(const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix, Index* accum_row_idx, Index* accum_row_ctr, int row_block, int curr_row)
        {
          AmaVankaMatrixHelper<MatrixA_>::fill_row_helper(matrix.template at<0,0>(), accum_row_idx, accum_row_ctr, row_block, curr_row);
          AmaVankaMatrixHelper<MatrixD_>::fill_row_helper(matrix.template at<1,0>(), accum_row_idx, accum_row_ctr, row_block+1, curr_row + matrix.template at<0,0>().template rows<LAFEM::Perspective::pod>());
        }

        static bool compare(const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix_a, const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix_b)
        {
          bool compare = true;
          compare &= AmaVankaMatrixHelper<MatrixA_>::compare(matrix_a.template at<0,0>(), matrix_b.template at<0,0>());
          compare &= AmaVankaMatrixHelper<MatrixB_>::compare(matrix_a.template at<0,1>(), matrix_b.template at<0,1>());
          compare &= AmaVankaMatrixHelper<MatrixD_>::compare(matrix_a.template at<1,0>(), matrix_b.template at<1,0>());
          return compare;
        }

      }; // struct AmaVankaMatrixHelper<LAFEM::SaddlePointMatrix<...>>

      template<typename Matrix_>
      CSRTupleMatrixWrapper<typename Matrix_::DataType, typename Matrix_::IndexType, Intern::AmaVankaMatrixHelper<Matrix_>::num_blocks>
      get_meta_matrix_wrapper(const Matrix_& matrix)
      {
        CSRTupleMatrixWrapper<typename Matrix_::DataType, typename Matrix_::IndexType, Intern::AmaVankaMatrixHelper<Matrix_>::num_blocks>
          mat_wrapper;
        AmaVankaMatrixHelper<Matrix_>::fill_ptr_wrapper(mat_wrapper, matrix, 0, 0);
        AmaVankaMatrixHelper<Matrix_>::fill_column_info(mat_wrapper, matrix, 0, 0);
        return mat_wrapper;
      }

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
         * This function is only implemented in specialized overloads,
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

        /// specialization for LAFEM::SparseMatrixCSR
        template<typename DT_, typename IT_>
        static void block_sizes(const LAFEM::SparseMatrixCSR<DT_, IT_>&, std::vector<Index>& sizes)
        {
          sizes.push_back(Index(1));
        }

        /// specialization for LAFEM::SparseMatrixBCSR
        template<typename DT_, typename IT_, int bs_>
        static void block_sizes(const LAFEM::SparseMatrixBCSR<DT_, IT_, bs_, bs_>&, std::vector<Index>& sizes)
        {
          sizes.push_back(Index(bs_));
        }

        /// specialization for LAFEM::TupleMatrix
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
        template<typename DT_, typename IT_, int bs_>
        static bool deduct_macro_dofs(
          const LAFEM::SaddlePointMatrix<
            LAFEM::SparseMatrixBCSR<DT_, IT_, bs_, bs_>,
            LAFEM::SparseMatrixBCSR<DT_, IT_, bs_, 1>,
            LAFEM::SparseMatrixBCSR<DT_, IT_, 1, bs_>>& matrix,
          std::vector<Adjacency::Graph>& macro_dofs)
        {
          typedef IT_ IndexType;

          // fetch matrix dimensions
          const Index num_dofs_p = Index(matrix.block_d().rows());

          // fetch matrix array
          const IndexType* row_ptr_b = matrix.block_b().row_ptr();
          const IndexType* col_idx_b = matrix.block_b().col_ind();
          const IndexType* row_ptr_d = matrix.block_d().row_ptr();
          const IndexType* col_idx_d = matrix.block_d().col_ind();

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
            for(Index j = p_ptr[i]; j < p_ptr[i+1]; ++j)
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
          macro_dofs.front() = Adjacency::Graph(num_macros, matrix.block_b().rows(), Index(v_idx.size()), v_ptr.data(), v_idx.data());
          macro_dofs.back()  = Adjacency::Graph(num_macros, matrix.block_d().rows(), Index(p_idx.size()), p_ptr.data(), p_idx.data());
          return true;
        }

        /* ********************************************************************************************************* */

#ifdef DOXYGEN
        /**
         * \brief Gathers the local matrix part of a system matrix w.r.t. a macro
         *
         * \note
         * This function is only implemented in specialized overloads,
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

        /// specialization for LAFEM::NullMatrix
        template<typename DT_, typename IT_, int bh_, int bw_>
        static std::pair<Index,Index> gather(const LAFEM::NullMatrix<DT_, IT_, bh_, bw_>&,
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

        /// specialization for LAFEM::SparseMatrixCSR
        template<typename DT_, typename IT_>
        static std::pair<Index,Index> gather(const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix,
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

            // initialize loop variable for local columns
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

        /// specialization for LAFEM::SparseMatrixBCSR
        template<typename DT_, typename IT_, int bh_, int bw_>
        static std::pair<Index,Index> gather(const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix,
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

            // initialize loop variable for local columns
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

        /// specialization for LAFEM::TupleMatrixRow
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

        /// specialization for LAFEM::TupleMatrixRow (single column)
        template<typename DT_, typename First_>
        static std::pair<Index,Index> gather(const LAFEM::TupleMatrixRow<First_>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          return AmaVankaCore::gather(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);
        }

        /// specialization for LAFEM::TupleMatrix
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

        /// specialization for LAFEM::TupleMatrix (single row)
        template<typename DT_, typename FirstRow_>
        static std::pair<Index,Index> gather(const LAFEM::TupleMatrix<FirstRow_>& matrix,
          DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          return AmaVankaCore::gather(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);
        }

        /// specialization for LAFEM::SaddlePointMatrix
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
         * This function is only implemented in specialized overloads,
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

        /// specialization for LAFEM::SparseMatrixCSR
        template<typename DT_, typename IT_>
        static void alloc(LAFEM::SparseMatrixCSR<DT_, IT_>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          XASSERT(row_block < Index(dof_macros.size()));
          XASSERT(col_block < Index(macro_dofs.size()));
          Adjacency::Graph graph(Adjacency::RenderType::injectify_sorted, dof_macros.at(row_block), macro_dofs.at(col_block));
          LAFEM::SparseMatrixCSR<DT_, IT_> matrix_main(graph);
          matrix.convert(matrix_main);
        }

        /// specialization for LAFEM::SparseMatrixBCSR
        template<typename DT_, typename IT_, int bh_, int bw_>
        static void alloc(LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          XASSERT(row_block < Index(dof_macros.size()));
          XASSERT(col_block < Index(macro_dofs.size()));
          Adjacency::Graph graph(Adjacency::RenderType::injectify_sorted, dof_macros.at(row_block), macro_dofs.at(col_block));
          matrix.convert(graph);
        }

        /// specialization for LAFEM::TupleMatrixRow
        template<typename First_, typename Second_, typename... Rest_>
        static void alloc(LAFEM::TupleMatrixRow<First_, Second_, Rest_...>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          AmaVankaCore::alloc(matrix.first(), dof_macros, macro_dofs, row_block, col_block);
          AmaVankaCore::alloc(matrix.rest(), dof_macros, macro_dofs, row_block, col_block + Index(1));
        }

        /// specialization for LAFEM::TupleMatrixRow (single column)
        template<typename First_>
        static void alloc(LAFEM::TupleMatrixRow<First_>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          AmaVankaCore::alloc(matrix.first(), dof_macros, macro_dofs, row_block, col_block);
        }

        /// specialization for LAFEM::TupleMatrix
        template<typename FirstRow_, typename SecondRow_, typename... RestRows_>
        static void alloc(LAFEM::TupleMatrix<FirstRow_, SecondRow_, RestRows_...>& matrix,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_block, const Index col_block)
        {
          AmaVankaCore::alloc(matrix.first(), dof_macros, macro_dofs, row_block, col_block);
          AmaVankaCore::alloc(matrix.rest(), dof_macros, macro_dofs, row_block + Index(1), col_block);
        }

        /// specialization for LAFEM::TupleMatrix (single row)
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
         * This function is only implemented in specialized overloads,
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

        /// specialization for LAFEM::SparseMatrixCSR
        template<typename DT_, typename IT_>
        static std::pair<Index,Index> scatter_add(LAFEM::SparseMatrixCSR<DT_, IT_>& matrix,
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

            // initialize loop variable for local columns
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

        /// specialization for LAFEM::SparseMatrixBCSR
        template<typename DT_, typename IT_, int bh_, int bw_>
        static std::pair<Index,Index> scatter_add(LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix,
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

            // initialize loop variable for local columns
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

        /// specialization for LAFEM::TupleMatrixRow
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

        /// specialization for LAFEM::TupleMatrixRow (single column)
        template<typename DT_, typename First_>
        static std::pair<Index,Index> scatter_add(LAFEM::TupleMatrixRow<First_>& matrix,
          const DT_* local, const Index stride, const Index macro, const std::vector<Adjacency::Graph>& macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          return AmaVankaCore::scatter_add(matrix.first(), local, stride, macro,
            macro_dofs, row_off, row_block, col_off, col_block);
        }

        /// specialization for LAFEM::TupleMatrix
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

        /// specialization for LAFEM::TupleMatrix (single row)
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
         * This function is only implemented in specialized overloads,
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

        /// specialization for LAFEM::SparseMatrixCSR
        template<typename DT_, typename IT_>
        static void scale_rows(LAFEM::SparseMatrixCSR<DT_, IT_>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          // get matrix arrays
          DT_* vals = matrix.val();
          const IT_* row_ptr = matrix.row_ptr();
          const IT_* col_idx = matrix.col_ind();

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

            if(n > 0u)
            {
              const DT_ sc = omega / DT_(Math::max(n,Index(1)));
              for(IT_ j(row_ptr[i]); j < row_ptr[i+1]; ++j)
                vals[j] *= sc;
            }
            else // null row: replace by unit row
            {
              for(IT_ j(row_ptr[i]); j < row_ptr[i+1]; ++j)
              {
                vals[j] = DT_(col_idx[j] == i ? 1 : 0);
              }
            }
          }
        }

        /// specialization for LAFEM::SparseMatrixBCSR
        template<typename DT_, typename IT_, int bh_, int bw_>
        static void scale_rows(LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          // get matrix arrays
          Tiny::Matrix<DT_, bh_, bw_>* vals = matrix.val();
          const IT_* row_ptr = matrix.row_ptr();
          const IT_* col_idx = matrix.col_ind();

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

            if(n > 0u)
            {
              const DT_ sc = omega / DT_(n);
              for(IT_ j(row_ptr[i]); j < row_ptr[i+1]; ++j)
                vals[j] *= sc;
            }
            else // null row: replace by unit row
            {
              for(IT_ j(row_ptr[i]); j < row_ptr[i+1]; ++j)
              {
                if(col_idx[j] == i)
                  vals[j].set_identity();
              }
            }
          }
        }

        /// specialization for LAFEM::TupleMatrixRow
        template<typename DT_, typename First_, typename Second_, typename... Rest_>
        static void scale_rows(LAFEM::TupleMatrixRow<First_, Second_, Rest_...>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          AmaVankaCore::scale_rows(matrix.first(), omega, dof_macros, macro_mask, row_block);
          AmaVankaCore::scale_rows(matrix.rest(), omega, dof_macros, macro_mask, row_block);
        }

        /// specialization for LAFEM::TupleMatrixRow (single column)
        template<typename DT_, typename First_>
        static void scale_rows(LAFEM::TupleMatrixRow<First_>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          AmaVankaCore::scale_rows(matrix.first(), omega, dof_macros, macro_mask, row_block);
        }

        /// specialization for LAFEM::TupleMatrix
        template<typename DT_, typename FirstRow_, typename SecondRow_, typename... RestRows_>
        static void scale_rows(LAFEM::TupleMatrix<FirstRow_, SecondRow_, RestRows_...>& matrix, const DT_ omega,
          const std::vector<Adjacency::Graph>& dof_macros, const std::vector<int>& macro_mask, const Index row_block)
        {
          AmaVankaCore::scale_rows(matrix.first(), omega, dof_macros, macro_mask, row_block);
          AmaVankaCore::scale_rows(matrix.rest(), omega, dof_macros, macro_mask, row_block + Index(1));
        }

        /// specialization for LAFEM::TupleMatrix (single row)
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

            // loop over all graphs and summarize degrees
            for(std::size_t igraph(0); igraph < num_graphs; ++igraph)
            {
              count += macro_dofs.at(igraph).degree(imacro) * block_sizes.at(igraph);
            }

            rtn = Math::max(rtn, count);
          }

          return rtn;
        }

      }; // struct AmaVankaCore

      /**
       * \brief VoxelAmaVanka main operations helper class
       *
       * This class contains the implementations of the main components of the VoxelAmaVanka smoother including the
       * necessary template wrapping magic to support a wide variety of meta-matrix wrappers on host and device side.
       *
       * \author Maximilian Esser
       */
      struct VoxelAmaVankaCore
      {
        template<typename DT_, typename IT_>
        CUDA_HOST_DEVICE static void gather_loc_cell(DT_* local, const DT_* vals, const IT_* row_ptr, const IT_* col_idx,
        const Index* row_img_idx, const Index* col_img_idx, const int bh_, const int bw_, const Index num_rows, const Index num_cols,
        const Index stride, const Index row_off, const Index col_off)
        {
          if(vals == nullptr)
            return;
          // loop over all rows
          for(Index i(0); i < num_rows; ++i)
          {
            // get the index of this row in the big matrix
            const Index irow = row_img_idx[i];

            // initialize loop variable for local columns
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
                // offset by ai times bw*bh size matrices
                const DT_* mv = vals + ai*IT_(bw_)*IT_(bh_);
                // copy matrix macro
                for(int ii(0); ii < bh_; ++ii)
                {
                  // compute offset within local matrix
                  const Index lro = (row_off + i * Index(bh_) + Index(ii)) * stride + col_off + j * Index(bw_);

                  // copy macro row
                  for(int jj(0); jj < bw_; ++jj)
                  {
                    local[lro + IT_(jj)] = mv[ii*bw_+jj];
                  }
                }
              }
            }
          }
        }

        // ##################### HOST VARIANT ################################
        template<typename DT_, typename IT_, int n_>
        CUDA_HOST static std::pair<Index,Index> gather_loc(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          DT_* local, const Index stride, const Index macro, const Adjacency::Graph** macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          //gather arrays
          const Index cur_ind = Index(n_)*row_block + col_block;
          const DT_* vals = mat_wrap.data_arrays[cur_ind];
          const IT_* row_ptr = mat_wrap.row_arrays[cur_ind];
          const IT_* col_idx = mat_wrap.col_arrays[cur_ind];
          const Adjacency::Graph& row_dofs = *macro_dofs[row_block];
          const Adjacency::Graph& col_dofs = *macro_dofs[col_block];

          // get our internal blocksize
          const int  bh_ = int(mat_wrap.get_block_sizes(row_block, col_block).first);
          const int  bw_ = int(mat_wrap.get_block_sizes(row_block, col_block).second);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_dofs.get_domain_ptr();
          const Index* col_dom_ptr = col_dofs.get_domain_ptr();

          // get the number of rows/cols to gather
          const Index num_rows = row_dom_ptr[macro+1] - row_dom_ptr[macro];
          const Index num_cols = col_dom_ptr[macro+1] - col_dom_ptr[macro];

          // get graph image index arrays for our macro
          const Index* row_img_idx = &(row_dofs.get_image_idx()[row_dom_ptr[macro]]);
          const Index* col_img_idx = &(col_dofs.get_image_idx()[col_dom_ptr[macro]]);

          gather_loc_cell(local, vals, row_ptr, col_idx, row_img_idx, col_img_idx, bh_, bw_, num_rows, num_cols,
                            stride, row_off, col_off);

          // return the number of gathered rows
          return std::make_pair(num_rows * Index(bh_), num_cols * Index(bw_));
        }

        template<typename DT_, typename IT_, int n_>
        CUDA_HOST static std::pair<Index,Index> gather_row(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          DT_* local, const Index stride, const Index macro, const Adjacency::Graph** macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          std::pair<Index, Index> nrc_r{0,0};
          //loop over cols
          for(int l = 0; l < n_; ++l)
          {
            std::pair<Index, Index> nrc_l = gather_loc(mat_wrap, local, stride, macro, macro_dofs, row_off, row_block, col_off+nrc_r.second, col_block+Index(l));
            nrc_r.first = nrc_l.first;
            nrc_r.second += nrc_l.second;
          }

          return nrc_r;
        }

        template<typename DT_, typename IT_, int n_>
        CUDA_HOST static std::pair<Index,Index> gather(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          DT_* local, const Index stride, const Index macro, const Adjacency::Graph** macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          std::pair<Index, Index> nrc_f{0,0};
          //loop over rows
          for(int l = 0; l < n_; ++l)
          {
            std::pair<Index, Index> nrc_r = gather_row(mat_wrap, local, stride, macro, macro_dofs, row_off+nrc_f.first, row_block+Index(l), col_off, col_block);
            nrc_f.first += nrc_r.first;
            nrc_f.second = nrc_r.second;
          }

          return nrc_f;
        }

        #ifdef __CUDACC__
        // ##################### DEVICE VARIANT ################################
        template<typename DT_, typename IT_, int n_, bool uniform_macros_>
        CUDA_HOST_DEVICE static cuda::std::pair<Index,Index> gather_loc(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          DT_* local, const Index stride, const Index macro, Index** macro_dofs, Index* max_degrees_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          //gather arrays
          const Index cur_ind = Index(n_)*row_block + col_block;
          const DT_* vals = mat_wrap.data_arrays[cur_ind];
          const IT_* row_ptr = mat_wrap.row_arrays[cur_ind];
          const IT_* col_idx = mat_wrap.col_arrays[cur_ind];
          const Index* row_dofs = macro_dofs[row_block];
          const Index* col_dofs = macro_dofs[col_block];
          const Index degree_row = max_degrees_dofs[row_block] + Index(!uniform_macros_);
          const Index degree_col = max_degrees_dofs[col_block] + Index(!uniform_macros_);

          // get our internal blocksize
          const int  bh_ = int(mat_wrap.blocksizes[row_block+CudaMath::cuda_min(row_block,Index(1))]);
          const int  bw_ = int(mat_wrap.blocksizes[col_block+Index(1)]);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_dofs + degree_row*macro;
          const Index* col_dom_ptr = col_dofs + degree_col*macro;

          // get the number of rows/cols to gather
          const Index num_rows = uniform_macros_ ? degree_row : row_dom_ptr[0];
          const Index num_cols = uniform_macros_ ? degree_col : col_dom_ptr[0];

          // // get graph image index arrays for our macro
          // const Index* row_img_idx = &(row_dofs.get_image_idx()[row_dom_ptr[macro]]);
          // const Index* col_img_idx = &(col_dofs.get_image_idx()[col_dom_ptr[macro]]);

          gather_loc_cell(local, vals, row_ptr, col_idx, row_dom_ptr + Index(!uniform_macros_), col_dom_ptr + Index(!uniform_macros_), bh_, bw_, num_rows, num_cols,
                            stride, row_off, col_off);

          // return the number of gathered rows
          return cuda::std::make_pair(num_rows * Index(bh_), num_cols * Index(bw_));
        }

        template<typename DT_, typename IT_, int n_, bool uniform_macros_>
        CUDA_HOST_DEVICE static cuda::std::pair<Index,Index> gather_row(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          DT_* local, const Index stride, const Index macro, Index** macro_dofs, Index* max_degrees_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          cuda::std::pair<Index, Index> nrc_r{0,0};
          //loop over cols
          for(int l = 0; l < n_; ++l)
          {
            cuda::std::pair<Index, Index> nrc_l = gather_loc<DT_, IT_, n_, uniform_macros_>(mat_wrap, local, stride, macro, macro_dofs, max_degrees_dofs, row_off, row_block, col_off+nrc_r.second, col_block+Index(l));
            nrc_r.first = nrc_l.first;
            nrc_r.second += nrc_l.second;
          }

          return nrc_r;
        }

        template<typename DT_, typename IT_, int n_, bool uniform_macros_>
        CUDA_HOST_DEVICE static cuda::std::pair<Index,Index> gather(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          DT_* local, const Index stride, const Index macro, Index** macro_dofs, Index* max_degrees_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          cuda::std::pair<Index, Index> nrc_f{0,0};
          //loop over rows
          for(int l = 0; l < n_; ++l)
          {
            cuda::std::pair<Index, Index> nrc_r = gather_row<DT_, IT_, n_, uniform_macros_>(mat_wrap, local, stride, macro, macro_dofs, max_degrees_dofs, row_off+nrc_f.first, row_block+Index(l), col_off, col_block);
            nrc_f.first += nrc_r.first;
            nrc_f.second = nrc_r.second;
          }

          return nrc_f;
        }
      #endif

        template<typename DT_, typename IT_>
        CUDA_HOST_DEVICE static void scatter_add_loc_cell(const DT_* local, DT_* vals, const IT_* row_ptr, const IT_* col_idx,
        const Index* row_img_idx, const Index* col_img_idx, const int bh_, const int bw_, const Index num_rows, const Index num_cols,
        const Index stride, const Index row_off, const Index col_off)
        {
          ASSERT(vals != nullptr);
          // loop over all rows
          for(Index i(0); i < num_rows; ++i)
          {
            // get the index of this row in the big matrix
            const Index irow = row_img_idx[i];

            // initialize loop variable for local columns
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
                // offset by ai times bh*bw site matrices
                DT_* mv = vals + ai*IT_(bh_)*IT_(bw_);
                // copy matrix macro
                for(int ii(0); ii < bh_; ++ii)
                {
                  // compute offset within local matrix
                  const Index lro = (row_off + i * Index(bh_) + Index(ii)) * stride + col_off + j * Index(bw_);
                  // copy macro row
                  for(int jj(0); jj < bw_; ++jj)
                  {
                    mv[ii*bw_+jj] += local[lro + IT_(jj)];
                  }
                }
              }
            }
          }
        }

        //################## Host Version#####################
        template<typename DT_, typename IT_, int n_>
        CUDA_HOST static std::pair<Index,Index> scatter_add_loc(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          const DT_* local, const Index stride, const Index macro, const Adjacency::Graph** macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          //gather arrays
          const Index cur_ind = Index(n_)*row_block + col_block;
          DT_* vals = mat_wrap.data_arrays[cur_ind];
          const IT_* row_ptr = mat_wrap.row_arrays[cur_ind];
          const IT_* col_idx = mat_wrap.col_arrays[cur_ind];
          const Adjacency::Graph& row_dofs = *macro_dofs[row_block];
          const Adjacency::Graph& col_dofs = *macro_dofs[col_block];

          // get our internal blocksize
          const int  bh_ = int(mat_wrap.get_block_sizes(row_block, col_block).first);
          const int  bw_ = int(mat_wrap.get_block_sizes(row_block, col_block).second);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_dofs.get_domain_ptr();
          const Index* col_dom_ptr = col_dofs.get_domain_ptr();

          // get the number of rows/cols to gather
          const Index num_rows = row_dom_ptr[macro+1] - row_dom_ptr[macro];
          const Index num_cols = col_dom_ptr[macro+1] - col_dom_ptr[macro];

          // get graph image index arrays for our macro
          const Index* row_img_idx = &(row_dofs.get_image_idx()[row_dom_ptr[macro]]);
          const Index* col_img_idx = &(col_dofs.get_image_idx()[col_dom_ptr[macro]]);

          scatter_add_loc_cell(local, vals, row_ptr, col_idx, row_img_idx, col_img_idx, bh_, bw_, num_rows, num_cols,
                            stride, row_off, col_off);

          // return the number of gathered rows
          return std::make_pair(num_rows * Index(bh_), num_cols * Index(bw_));
        }

        template<typename DT_, typename IT_, int n_>
        CUDA_HOST static std::pair<Index,Index> scatter_add_row(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          const DT_* local, const Index stride, const Index macro, const Adjacency::Graph** macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          std::pair<Index, Index> nrc_r{0,0};
          //loop over cols
          for(int l = 0; l < n_; ++l)
          {
            std::pair<Index, Index> nrc_l = scatter_add_loc(mat_wrap, local, stride, macro, macro_dofs, row_off, row_block, col_off+nrc_r.second, col_block+Index(l));
            nrc_r.first = nrc_l.first;
            nrc_r.second += nrc_l.second;
          }

          return nrc_r;
        }

        template<typename DT_, typename IT_, int n_>
        CUDA_HOST static std::pair<Index,Index> scatter_add(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          const DT_* local, const Index stride, const Index macro, const Adjacency::Graph** macro_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          std::pair<Index, Index> nrc_f{0,0};
          //loop over rows
          for(int l = 0; l < n_; ++l)
          {
            std::pair<Index, Index> nrc_r = scatter_add_row(mat_wrap, local, stride, macro, macro_dofs, row_off+nrc_f.first, row_block+Index(l), col_off, col_block);
            nrc_f.first += nrc_r.first;
            nrc_f.second = nrc_f.second;
          }

          return nrc_f;
        }

        #ifdef __CUDACC__
        //################# Device version #######################
        template<typename DT_, typename IT_, int n_, bool uniform_macros_>
        CUDA_HOST_DEVICE static cuda::std::pair<Index,Index> scatter_add_loc(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          const DT_* local, const Index stride, const Index macro, Index** macro_dofs, Index* max_degrees_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          //gather arrays
          const Index cur_ind = Index(n_)*row_block + col_block;
          DT_* vals = mat_wrap.data_arrays[cur_ind];
          const IT_* row_ptr = mat_wrap.row_arrays[cur_ind];
          const IT_* col_idx = mat_wrap.col_arrays[cur_ind];
          const Index* row_dofs = macro_dofs[row_block];
          const Index* col_dofs = macro_dofs[col_block];
          const Index degree_row = max_degrees_dofs[row_block] + Index(!uniform_macros_);
          const Index degree_col = max_degrees_dofs[col_block] + Index(!uniform_macros_);

          // get our internal blocksize
          const int  bh_ = int(mat_wrap.blocksizes[row_block+CudaMath::cuda_min(row_block,Index(1))]);
          const int  bw_ = int(mat_wrap.blocksizes[col_block+Index(1)]);

          // get graph domain pointer arrays
          const Index* row_dom_ptr = row_dofs + degree_row*macro;
          const Index* col_dom_ptr = col_dofs + degree_col*macro;

          // get the number of rows/cols to gather
          const Index num_rows = uniform_macros_ ? degree_row : row_dom_ptr[0];
          const Index num_cols = uniform_macros_ ? degree_col : col_dom_ptr[0];

          // // get graph image index arrays for our macro
          // const Index* row_img_idx = &(row_dofs.get_image_idx()[row_dom_ptr[macro]]);
          // const Index* col_img_idx = &(col_dofs.get_image_idx()[col_dom_ptr[macro]]);

          scatter_add_loc_cell(local, vals, row_ptr, col_idx, row_dom_ptr + Index(!uniform_macros_), col_dom_ptr + Index(!uniform_macros_), bh_, bw_, num_rows, num_cols,
                            stride, row_off, col_off);

          // return the number of gathered rows
          return cuda::std::make_pair(num_rows * Index(bh_), num_cols * Index(bw_));
        }

        template<typename DT_, typename IT_, int n_, bool uniform_macros_>
        CUDA_HOST_DEVICE static cuda::std::pair<Index,Index> scatter_add_row(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          const DT_* local, const Index stride, const Index macro, Index** macro_dofs, Index* max_degrees_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          cuda::std::pair<Index, Index> nrc_r{0,0};
          //loop over cols
          for(int l = 0; l < n_; ++l)
          {
            cuda::std::pair<Index, Index> nrc_l = scatter_add_loc<DT_, IT_, n_, uniform_macros_>(mat_wrap, local, stride, macro, macro_dofs, max_degrees_dofs, row_off, row_block, col_off+nrc_r.second, col_block+Index(l));
            nrc_r.first = nrc_l.first;
            nrc_r.second += nrc_l.second;
          }

          return nrc_r;
        }

        template<typename DT_, typename IT_, int n_, bool uniform_macros_>
        CUDA_HOST_DEVICE static cuda::std::pair<Index,Index> scatter_add(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
          const DT_* local, const Index stride, const Index macro, Index** macro_dofs, Index* max_degrees_dofs,
          const Index row_off, const Index row_block, const Index col_off, const Index col_block)
        {
          cuda::std::pair<Index, Index> nrc_f{0,0};
          //loop over rows
          for(int l = 0; l < n_; ++l)
          {
            cuda::std::pair<Index, Index> nrc_r = scatter_add_row<DT_, IT_, n_, uniform_macros_>(mat_wrap, local, stride, macro, macro_dofs, max_degrees_dofs, row_off+nrc_f.first, row_block+Index(l), col_off, col_block);
            nrc_f.first += nrc_r.first;
            nrc_f.second = nrc_f.second;
          }

          return nrc_f;
        }
        #endif

        template<typename DT_, typename IT_, bool skip_singular_>
        CUDA_HOST static void scale_row(DT_* vals, DT_ omega, const IT_* row_ptr, const IT_* col_idx, const Index* row_dom_ptr, const Index* row_img_idx,
                                  int hw, int hb, Index cur_row, const int* m_mask)
        {
          ASSERT(row_dom_ptr[cur_row] < row_dom_ptr[cur_row+1]);
          //get number of active macros for this DOF
          Index n(0);
          if constexpr(skip_singular_)
          {
            for(Index j = row_dom_ptr[cur_row]; j < row_dom_ptr[cur_row+1]; ++j)
            {
              n += Index(m_mask[row_img_idx[j]]);
            }
          }
          else
          {
            n = row_dom_ptr[cur_row+1] - row_dom_ptr[cur_row];
          }

          if(n > 0)
          {
            const DT_ sc = omega/ DT_(Math::max(n, Index(1)));
            for(IT_ j = row_ptr[cur_row]; j < row_ptr[cur_row+1]; ++j)
            {
              DT_* loc_val = vals + j*IT_(hw)*IT_(hb);
              for(int ii = 0; ii < hw; ++ii)
              {
                for(int jj = 0; jj < hb; ++jj)
                {
                  loc_val[ii*hb+jj] *= sc;
                }
              }
            }
          }
          else // null row replace by unit row
          {
            for(IT_ j = row_ptr[cur_row]; j < row_ptr[cur_row+1]; ++j)
            {
              DT_* loc_val = vals + j*IT_(hw)*IT_(hb);
              const DT_ diag_val = DT_(col_idx[j] == cur_row ? 1 : 0);
              for(int ii = 0; ii < hw; ++ii)
              {
                for(int jj = 0; jj < hb; ++jj)
                {
                  loc_val[ii*hb+jj] = DT_(0);
                }
                loc_val[ii*hb +ii] = diag_val;
              }
            }
          }

        }

        #ifdef __CUDACC__
        template<typename DT_, typename IT_, bool skip_singular_>
        CUDA_HOST_DEVICE static void scale_row(DT_* vals, DT_ omega, const IT_* row_ptr, const IT_* col_idx, const Index* row_dom_ptr, const Index row_number,
                                  int hw, int hb, Index cur_row, const int* m_mask)
        {
          //get number of active macros for this DOF
          Index n(0);
          if constexpr(skip_singular_)
          {
            for(Index j = 0; j < row_number; ++j)
            {
              n += Index(m_mask[row_dom_ptr[j]]);
            }
          }
          else
          {
            n = row_number;
          }

          if(n > 0)
          {
            const DT_ sc = omega/ DT_(CudaMath::cuda_max(n, Index(1)));
            for(IT_ j = row_ptr[cur_row]; j < row_ptr[cur_row+1]; ++j)
            {
              DT_* loc_val = vals + j*hw*hb;
              for(int ii = 0; ii < hw; ++ii)
              {
                for(int jj = 0; jj < hb; ++jj)
                {
                  loc_val[ii*hb+jj] *= sc;
                }
              }
            }
          }
          else // null row replace by unit row
          {
            for(IT_ j = row_ptr[cur_row]; j < row_ptr[cur_row+1]; ++j)
            {
              DT_* loc_val = vals + j*hw*hb;
              const DT_ diag_val = DT_(col_idx[j] == cur_row ? 1 : 0);
              for(int ii = 0; ii < hw; ++ii)
              {
                for(int jj = 0; jj < hb; ++jj)
                {
                  loc_val[ii*hb+jj] = DT_(0);
                }
                loc_val[ii*hb +ii] = diag_val;
              }
            }
          }

        }
        #endif
      };

    } // namespace Intern
    /// \endcond
  }
}

#endif // KERNEL_SOLVER_AMAVANKA_HPP