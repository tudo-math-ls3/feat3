#pragma once
#ifndef KERNEL_LAFEM_GATHER_SCATTER_HPP
#define KERNEL_LAFEM_GATHER_SCATTER_HPP 1

// includes, FEAST-LAFEM
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Scatter-Axpy class template
     *
     * \tparam Container_
     * The type of the container that is the target of the scatter-axpy operation.
     *
     * \author Peter Zajac
     */
    template<typename Container_>
#ifndef DOXYGEN
    class ScatterAxpy;
#else
    class ScatterAxpy
    {
    public:
      /**
       * \brief Constructor
       *
       * \param[in] container
       * A reference to the container that serves as the scatter target.
       */
      explicit(Container_& container);

      /**
       * \brief Scatter-Axpy operator
       *
       * \tparam LocalData_
       * The class of the local data container.
       * See the Assembly::LocalVectorData and Assembly::LocalMatrixData class templates for
       * an interface definition and documentation.
       *
       * \param[in] local_data
       * The local data object that serves as the scatter source.
       *
       * \param[in] alpha
       * The scaling factor for the scatter-axpy operation.
       */
      template<typename LocalData_>
      void operator()(
        const LocalData_& local_data,
        DataType alpha = DataType_(1));
    };
#endif // DOXYGEN

    /**
     * \brief Gather-Axpy class template
     *
     * \tparam Container_
     * The type of the container that is the source of the gather-axpy operation.
     *
     * \author Peter Zajac
     */
    template<typename Container_>
#ifndef DOXYGEN
    class GatherAxpy;
#else
    class GatherAxpy
    {
    public:
      /**
       * \brief Constructor
       *
       * \param[in] container
       * A const reference to the container that serves as the gather source.
       */
      explicit(const Container_& container);

      /**
       * \brief Gather-Axpy operator
       *
       * \tparam LocalData_
       * The class of the local data container.
       *
       * \param[in] local_data
       * The local data object that serves as the gather target.
       *
       * \param[in] alpha
       * The scaling factor for the gather-axpy operation.
       */
      template<typename LocalData_>
      void operator()(
        LocalData_& local_data,
        DataType alpha = DataType_(1));
    };
#endif // DOXYGEN

    /**
     * \brief Scatter-Axpy specialisation for DenseVector
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_>
    class ScatterAxpy< LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> VectorType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
      Index _num_entries;
      DataType_* _data;

    public:
      explicit ScatterAxpy(VectorType& vector) :
        _num_entries(vector.size()),
        _data(vector.elements())
      {
      }

      virtual ~ScatterAxpy()
      {
      }

      template<typename LocalData_>
      void operator()(
        const LocalData_& loc_vec,
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
    }; // class ScatterAxpy<LAFEM::DenseVector<Mem::Main,...>>

    /**
     * \brief Gather-Axpy specialisation for DenseVector
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_>
    class GatherAxpy< LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> VectorType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
      Index _num_entries;
      const DataType_* _data;

    public:
      explicit GatherAxpy(const VectorType& vector) :
        _num_entries(vector.size()),
        _data(vector.elements())
      {
      }

      virtual ~GatherAxpy()
      {
      }

      template<typename LocalData_>
      void operator()(
        LocalData_& loc_vec,
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
    }; // class GatherAxpy<LAFEM::DenseVector<Mem::Main,...>>

    /**
     * \brief Scatter-Axpy specialisation for SparseMatrixCSR
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_>
    class ScatterAxpy< LAFEM::SparseMatrixCSR<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      IndexType_* _row_ptr;
      IndexType_* _col_idx;
      IndexType_* _col_ptr;
      DataType_ *_data;

    public:
      explicit ScatterAxpy(MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~IndexType_(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _row_ptr(matrix.row_ptr()),
        _col_idx(matrix.col_ind()),
        _col_ptr(nullptr),
        _data(matrix.val())
      {
        // allocate column-pointer array
        _col_ptr = new IndexType_[matrix.columns()];
#ifdef DEBUG
        for(Index i(0); i < _num_cols; ++i)
        {
          _col_ptr[i] = _deadcode;
        }
#endif
      }

      virtual ~ScatterAxpy()
      {
        if(_col_ptr != nullptr)
        {
          delete [] _col_ptr;
        }
      }

      template<typename LocalData_>
      void operator()(
        const LocalData_& loc_mat,
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
            for(IndexType_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
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
            for(IndexType_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
            {
              _col_ptr[_col_idx[k]] = _deadcode;
            }
#endif
            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class ScatterAxpy<LAFEM::SparseMatrixCSR<Mem::Main,...>>

    /**
     * \brief Gather-Axpy specialisation for SparseMatrixCSR
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_>
    class GatherAxpy< LAFEM::SparseMatrixCSR<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      const IndexType_* _row_ptr;
      const IndexType_* _col_idx;
      IndexType_* _col_ptr;
      const DataType_ *_data;

    public:
      explicit GatherAxpy(const MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~IndexType_(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _row_ptr(matrix.row_ptr()),
        _col_idx(matrix.col_ind()),
        _col_ptr(nullptr),
        _data(matrix.val())
      {
        // allocate column-pointer array
        _col_ptr = new IndexType_[matrix.columns()];
#ifdef DEBUG
        for(Index i(0); i < _num_cols; ++i)
        {
          _col_ptr[i] = _deadcode;
        }
#endif
      }

      virtual ~GatherAxpy()
      {
        if(_col_ptr != nullptr)
        {
          delete [] _col_ptr;
        }
      }

      template<typename LocalData_>
      void operator()(
        LocalData_& loc_mat,
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
            for(IndexType_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
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
            for(IndexType_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
            {
              _col_ptr[_col_idx[k]] = _deadcode;
            }
#endif

            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class GatherAxpy<LAFEM::SparseMatrixCSR<Mem::Main,...>>

    /**
     * \brief Scatter-Axpy specialisation for SparseMatrixCOO
     *
     * \author Christoph Lohmann
     */
    template<typename DataType_, typename IndexType_>
    class ScatterAxpy< LAFEM::SparseMatrixCOO<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::SparseMatrixCOO<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      Index _used_elements;
      IndexType_* _row_idx;
      IndexType_* _col_idx;
      IndexType_* _col_ptr;
      DataType_ *_data;

    public:
      explicit ScatterAxpy(MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~IndexType_(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _used_elements(matrix.used_elements()),
        _row_idx(matrix.row_indices()),
        _col_idx(matrix.column_indices()),
        _col_ptr(nullptr),
        _data(matrix.val())
      {
        // allocate column-pointer array
        _col_ptr = new IndexType_[matrix.columns()];
#ifdef DEBUG
        for(Index i(0); i < _num_cols; ++i)
        {
          _col_ptr[i] = _deadcode;
        }
#endif
      }

      virtual ~ScatterAxpy()
      {
        if(_col_ptr != nullptr)
        {
          delete [] _col_ptr;
        }
      }

      template<typename LocalData_>
      void operator()(const LocalData_& loc_mat,
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
            IndexType_ k(0);
            while (_row_idx[k] < ix)
            {
              ++k;
            }
            while (k < _used_elements && _row_idx[k] <= ix)
            {
              _col_ptr[_col_idx[k]] = k;
              ++k;
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
            k = IndexType_(0);
            while (_row_idx[k] < ix)
            {
              ++k;
            }
            while (k < _used_elements && _row_idx[k] <= ix)
            {
              _col_ptr[_col_idx[k]] = _deadcode;
              ++k;
            }
#endif
            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class ScatterAxpy<LAFEM::SparseMatrixCOO<Mem::Main,...>>

    /**
     * \brief Gather-Axpy specialisation for SparseMatrixCOO
     *
     * \author Christoph Lohmann
     */
    template<typename DataType_, typename IndexType_>
    class GatherAxpy< LAFEM::SparseMatrixCOO<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::SparseMatrixCOO<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      Index _used_elements;
      const IndexType_* _row_idx;
      const IndexType_* _col_idx;
      IndexType_* _col_ptr;
      const DataType_ *_data;

    public:
      explicit GatherAxpy(const MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~IndexType_(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _used_elements(matrix.used_elements()),
        _row_idx(matrix.row_indices()),
        _col_idx(matrix.column_indices()),
        _col_ptr(nullptr),
        _data(matrix.val())
      {
        // allocate column-pointer array
        _col_ptr = new IndexType_[matrix.columns()];
#ifdef DEBUG
        for(Index i(0); i < _num_cols; ++i)
        {
          _col_ptr[i] = _deadcode;
        }
#endif
      }

      virtual ~GatherAxpy()
      {
        if(_col_ptr != nullptr)
        {
          delete [] _col_ptr;
        }
      }

      template<typename LocalData_>
      void operator()(
                      LocalData_& loc_mat,
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
            IndexType_ k(0);
            while (_row_idx[k] < ix)
            {
              ++k;
            }
            while (k < _used_elements && _row_idx[k] <= ix)
            {
              _col_ptr[_col_idx[k]] = k;
              ++k;
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
            k = IndexType_(0);
            while (_row_idx[k] < ix)
            {
              ++k;
            }
            while (k < _used_elements && _row_idx[k] <= ix)
            {
              _col_ptr[_col_idx[k]] = _deadcode;
              ++k;
            }
#endif

            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class GatherAxpy<LAFEM::SparseMatrixCOO<Mem::Main,...>>

    /**
     * \brief Scatter-Axpy specialisation for SparseMatrixELL
     *
     * \author Christoph Lohmann
     */
    template<typename DataType_, typename IndexType_>
    class ScatterAxpy< LAFEM::SparseMatrixELL<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::SparseMatrixELL<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      IndexType_ _C;
      const IndexType_ * _cs;
      const IndexType_ * _col_ind;
      IndexType_ * _col_ptr;
      DataType_ * _val;

    public:
      explicit ScatterAxpy(MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~IndexType_(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _C(IndexType_(matrix.C())),
        _cs(matrix.cs()),
        _col_ind(matrix.col_ind()),
        _col_ptr(nullptr),
        _val(matrix.val())
      {
        // allocate column-pointer array
        _col_ptr = new IndexType_[matrix.columns()];
#ifdef DEBUG
        for(Index i(0); i < _num_cols; ++i)
        {
          _col_ptr[i] = _deadcode;
        }
#endif
      }

      virtual ~ScatterAxpy()
      {
        if(_col_ptr != nullptr)
        {
          delete [] _col_ptr;
        }
      }

      template<typename LocalData_>
      void operator()(const LocalData_& loc_mat,
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
            IndexType_ ix = IndexType_(loc_mat.get_row_index(i, ic));

            // build column pointer for this row entry contribution
            for(IndexType_ k(_cs[ix/_C] + ix%_C); k < _cs[ix/_C + 1]; k += _C)
            {
              _col_ptr[_col_ind[k]] = k;
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
                _val[_col_ptr[jx]] += iw * jw * loc_mat(i,j);

                // continue with next column contribution
              }
              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(IndexType_ k(_cs[ix/_C] + ix%_C); k < _cs[ix/_C + 1]; k += _C)
            {
              _col_ptr[_col_ind[k]] = _deadcode;
            }
#endif
            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class ScatterAxpy<LAFEM::SparseMatrixELL<Mem::Main,...>>

    /**
     * \brief Gather-Axpy specialisation for SparseMatrixELL
     *
     * \author Christoph Lohmann
     */
    template<typename DataType_, typename IndexType_>
    class GatherAxpy< LAFEM::SparseMatrixELL<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::SparseMatrixELL<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      IndexType_ _C;
      const IndexType_ * _cs;
      const IndexType_ * _col_ind;
      IndexType_ * _col_ptr;
      const DataType_ * _val;

    public:
      explicit GatherAxpy(const MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~IndexType_(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _C(IndexType_(matrix.C())),
        _cs(matrix.cs()),
        _col_ind(matrix.col_ind()),
        _col_ptr(nullptr),
        _val(matrix.val())
      {
        // allocate column-pointer array
        _col_ptr = new IndexType_[matrix.columns()];
#ifdef DEBUG
        for(Index i(0); i < _num_cols; ++i)
        {
          _col_ptr[i] = _deadcode;
        }
#endif
      }

      virtual ~GatherAxpy()
      {
        if(_col_ptr != nullptr)
        {
          delete [] _col_ptr;
        }
      }

      template<typename LocalData_>
      void operator()(LocalData_& loc_mat,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(Index i(0); i < loc_mat.get_num_rows(); ++i)
        {
          // loop over all row entry contributations
          for(Index ic(0); ic < loc_mat.get_num_row_contribs(i); ++ic)
          {
            // fetch row index
            IndexType_ ix = IndexType_(loc_mat.get_row_index(i, ic));

            // build column pointer for this row entry contribution
            for(IndexType_ k(_cs[ix/_C] + ix%_C); k < _cs[ix/_C + 1]; k += _C)
            {
              _col_ptr[_col_ind[k]] = k;
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
                dx += DataType_(loc_mat.get_col_weight(j, jc)) * _val[_col_ptr[jx]];

                // continue with next column contribution
              }

              // update local matrix data
              loc_mat(i,j) += alpha * DataType_(loc_mat.get_row_weight(i, ic)) * dx;

              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(IndexType_ k(_cs[ix/_C] + ix%_C); k < _cs[ix/_C + 1]; k += _C)
            {
              _col_ptr[_col_ind[k]] = _deadcode;
            }
#endif

            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class GatherAxpy<LAFEM::SparseMatrixELL<Mem::Main,...>>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_GATHER_SCATTER_HPP
