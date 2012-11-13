#pragma once
#ifndef KERNEL_ASSEMBLY_LAFEM_BACKEND_HPP
#define KERNEL_ASSEMBLY_LAFEM_BACKEND_HPP 1

// includes, FEAST
#include <kernel/util/tiny_algebra.hpp>

// includes, FEAST-LAFEM
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /// \todo: all the class templates in this file shall be moved to LAFEM later...
    template<typename Vector_>
    class VectorScatterAxpy DOXY({});

    template<typename Vector_>
    class VectorGatherAxpy DOXY({});

    template<typename Matrix_>
    class MatrixScatterAxpy DOXY({});

    template<typename Matrix_>
    class MatrixGatherAxpy DOXY({});

    template<typename DataType_>
    class VectorScatterAxpy< LAFEM::DenseVector<Mem::Main, DataType_> >
    {
    public:
      typedef LAFEM::DenseVector<Mem::Main, DataType_> VectorType;

    private:
      Index _num_entries;
      DataType_* _data;

    public:
      explicit VectorScatterAxpy(VectorType& vector) :
        _num_entries(vector.size()),
        _data(vector.elements())
      {
      }

      virtual ~VectorScatterAxpy()
      {
      }

      template<typename LocalVectorData_>
      void operator()(
        const LocalVectorData_& loc_vec,
        DataType_ alpha = DataType_(1))
      {
        // loop over all local entries
        for(Index i(0); i < loc_vec.get_num_entries(); ++i)
        {
          // pre-multiply local entry by alpha
          DataType_ dx(alpha * loc_vec(i));

          // loop over all entry contributions
          for(Index ic(0); ic < loc_vec.get_num_contribs(i); ++ic)
          {
            // update vector data
            _data[loc_vec.get_index(i, ic)] += loc_vec.get_weight(i, ic) * dx;
          }
        }
      }
    }; // class VectorScatterAxpy<LAFEM::DenseVector<Mem::Main,...>>

    template<typename DataType_>
    class VectorGatherAxpy< LAFEM::DenseVector<Mem::Main, DataType_> >
    {
    public:
      typedef LAFEM::DenseVector<Mem::Main, DataType_> VectorType;

    private:
      Index _num_entries;
      const DataType_* _data;

    public:
      explicit VectorGatherAxpy(const VectorType& vector) :
        _num_entries(vector.size()),
        _data(vector.elements())
      {
      }

      virtual ~VectorGatherAxpy()
      {
      }

      template<typename LocalVectorData_>
      void operator()(
        LocalVectorData_& loc_vec,
        DataType_ alpha = DataType_(1))
      {
        // loop over all local entries
        for(Index i(0); i < loc_vec.get_num_entries(); ++i)
        {
          // clear accumulation entry
          DataType_ dx(DataType_(0));

          // loop over all entry contributions
          for(Index ic(0); ic < loc_vec.get_num_contribs(i); ++ic)
          {
            // update accumulator
            dx += loc_vec.get_weight(i, ic) * _data[loc_vec.get_index(i, ic)];
          }

          // update local vector data
          loc_vec(i) += alpha * dx;
        }
      }
    }; // class VectorGatherAxpy<LAFEM::DenseVector<Mem::Main,...>>

    template<typename DataType_>
    class MatrixScatterAxpy< LAFEM::SparseMatrixCSR<Mem::Main, DataType_> >
    {
    public:
      typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_> MatrixType;

    private:
#ifdef DEBUG
      Index _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      Index* _row_ptr;
      Index* _row_end;
      Index* _col_idx;
      Index* _col_ptr;
      DataType_ *_data;

    public:
      explicit MatrixScatterAxpy(MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~Index(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _row_ptr(matrix.row_ptr()),
        _row_end(matrix.row_ptr_end()),
        _col_idx(matrix.col_ind()),
        _col_ptr(nullptr),
        _data(matrix.val())
      {
        // allocate column-pointer array
        _col_ptr = new Index[matrix.columns()];
#ifdef DEBUG
        for(Index i(0); i < _num_cols; ++i)
        {
          _col_ptr[i] = _deadcode;
        }
#endif
      }

      virtual ~MatrixScatterAxpy()
      {
        if(_col_ptr != nullptr)
        {
          delete [] _col_ptr;
        }
      }

      template<typename LocalMatrixData_>
      void operator()(
        const LocalMatrixData_& loc_mat,
        DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(Index i(0); i < loc_mat.get_num_rows(); ++i)
        {
          // loop over all row entry contributations
          for(Index ic(0); ic < loc_mat.get_num_row_contribs(i); ++ic)
          {
            // fetch row entry weight and pre-multiply by alpha
            DataType_ iw = alpha * DataType_(loc_mat.get_row_weight(i, ic));

            // fetch row index
            Index ix = loc_mat.get_row_index(i, ic);

            // build column pointer for this row entry contribution
            for(Index k(_row_ptr[ix]); k < _row_end[ix]; ++k)
            {
              _col_ptr[_col_idx[k]] = k;
            }

            // loop over all local column entries
            for(Index j(0); j < loc_mat.get_num_cols(); ++j)
            {
              // loop over all column entry contributions
              for(Index jc(0); jc < loc_mat.get_num_col_contribs(j); ++jc)
              {
                // fetch trial function dof weight
                DataType_ jw = DataType_(loc_mat.get_col_weight(j, jc));

                // fetch column index
                Index jx = loc_mat.get_col_index(j, jc);

#ifdef DEBUG
                // ensure that the column pointer is valid for this index
                ASSERT(_col_ptr[jx] != _deadcode, "invalid column index");
#endif

                // incorporate data into global matrix
                _data[_col_ptr[jx]] += iw * jw * loc_mat(i,j);

                // continue with next column contribution
              }
              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(Index k(_row_ptr[ix]); k < _row_end[ix]; ++k)
            {
              _col_ptr[_col_idx[k]] = _deadcode;
            }
#endif
            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class MatrixScatterAdder<LAFEM::SparseMatrixCSR<Mem::Main,...>>

    template<typename DataType_>
    class MatrixGatherAxpy< LAFEM::SparseMatrixCSR<Mem::Main, DataType_> >
    {
    public:
      typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_> MatrixType;

    private:
#ifdef DEBUG
      Index _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      const Index* _row_ptr;
      const Index* _row_end;
      const Index* _col_idx;
      Index* _col_ptr;
      const DataType_ *_data;

    public:
      explicit MatrixGatherAxpy(const MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~Index(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _row_ptr(matrix.row_ptr()),
        _row_end(matrix.row_ptr_end()),
        _col_idx(matrix.col_ind()),
        _col_ptr(nullptr),
        _data(matrix.val())
      {
        // allocate column-pointer array
        _col_ptr = new Index[matrix.columns()];
#ifdef DEBUG
        for(Index i(0); i < _num_cols; ++i)
        {
          _col_ptr[i] = _deadcode;
        }
#endif
      }

      virtual ~MatrixGatherAxpy()
      {
        if(_col_ptr != nullptr)
        {
          delete [] _col_ptr;
        }
      }

      template<typename LocalMatrixData_>
      void operator()(
        LocalMatrixData_& loc_mat,
        DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(Index i(0); i < loc_mat.get_num_rows(); ++i)
        {
          // loop over all row entry contributations
          for(Index ic(0); ic < loc_mat.get_num_row_contribs(i); ++ic)
          {
            // fetch row index
            Index ix = loc_mat.get_row_index(i, ic);

            // build column pointer for this row entry contribution
            for(Index k(_row_ptr[ix]); k < _row_end[ix]; ++k)
            {
              _col_ptr[_col_idx[k]] = k;
            }

            // loop over all local column entries
            for(Index j(0); j < loc_mat.get_num_cols(); ++j)
            {
              // clear  accumulation entry
              DataType_ dx(DataType_(0));

              // loop over all column entry contributions
              for(Index jc(0); jc < loc_mat.get_num_col_contribs(j); ++jc)
              {
                // fetch column index
                Index jx = loc_mat.get_col_index(j, jc);

#ifdef DEBUG
                // ensure that the column pointer is valid for this index
                ASSERT(_col_ptr[jx] != _deadcode, "invalid column index");
#endif

                // update accumulator
                dx += DataType_(loc_mat.get_col_weight(j, jc)) * _data[_col_ptr[jx]];

                // continue with next column contribution
              }

              // update local matrix data
              loc_mat(i,j) += alpha * DataType_(loc_mat.get_row_weight(i, ic)) * dx;

              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(Index k(_row_ptr[ix]); k < _row_end[ix]; ++k)
            {
              _col_ptr[_col_idx[k]] = _deadcode;
            }
#endif

            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class MatrixGatherAxpy<LAFEM::SparseMatrixCSR<Mem::Main,...>>

    /// \cond internal
    namespace Intern
    {
      template<typename Matrix_>
      class RowScaler;

      template<typename DataType_>
      class RowScaler< LAFEM::SparseMatrixCSR<Mem::Main, DataType_> >
      {
      public:
        typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_> MatrixType;

        static void apply(MatrixType& matrix, const DataType_ x[])
        {
          Index* row_ptr(matrix.row_ptr());
          Index* row_end(matrix.row_ptr_end());
          DataType_* data(matrix.val());

          for(Index i(0); i < matrix.rows(); ++i)
          {
            for(Index j(row_ptr[i]); j < row_end[i]; ++j)
              data[j] *= x[i];
          }
        }
      };
    }
    /// \endcond
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_LAFEM_BACKEND_HPP
