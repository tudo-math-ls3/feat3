#pragma once
#ifndef KERNEL_LAFEM_GATHER_SCATTER_HPP
#define KERNEL_LAFEM_GATHER_SCATTER_HPP 1

// includes, FEAST-LAFEM
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>

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
       * \brief Vector Scatter-Axpy operator
       *
       * \tparam LocalVector_
       * The class of the local vector container. Usually a Tiny::Vector instance.
       *
       * \tparam Mapping_
       * The class of the mapping for the scatter-axpy operation.
       *
       * \param[in] local_vector
       * The local vector object that serves as the scatter source.
       *
       * \param[in] mapping
       * The mapping to be used for scattering.
       *
       * \param[in] alpha
       * The scaling factor for the scatter-axpy operation.
       */
      template<typename LocalVector_, typename Mapping_>
      void operator()(
                      const LocalVector_& local_vector,
                      const Mapping_& mapping,
                      DataType alpha = DataType_(1));

      /**
       * \brief Matrix Scatter-Axpy operator
       *
       * \tparam LocalMatrix_
       * The class of the local matrix container. Usually a Tiny::Matrix instance.
       *
       * \tparam RowMapping_
       * The class of the row-mapping for the scatter-axpy operation.
       *
       * \tparam ColMapping_
       * The class of the column-mapping for the scatter-axpy operation.
       *
       * \param[in] local_matrix
       * The local matrix object that serves as the scatter source.
       *
       * \param[in] row_map
       * The row-mapping to be used for scattering.
       *
       * \param[in] col_map
       * The col-mapping to be used for scattering.
       *
       * \param[in] alpha
       * The scaling factor for the scatter-axpy operation.
       */
      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      const LocalMatrix_& local_matrix,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
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
       * \brief Vector Gather-Axpy operator
       *
       * \tparam LocalVector_
       * The class of the local vector container. Usually a Tiny::Vector
       *
       * \tparam Mapping_
       * The class of the mapping for the gather-axpy operation.
       *
       * \param[out] local_vector
       * The local vector object that serves as the gather target.
       *
       * \param[in] mapping
       * The mapping to be used for scattering.
       *
       * \param[in] alpha
       * The scaling factor for the gather-axpy operation.
       */
      template<typename LocalData_, typename Mapping_>
      void operator()(
                      LocalData_& local_data,
                      const Mapping_& mapping,
                      DataType alpha = DataType_(1));

      /**
       * \brief Matrix Gather-Axpy operator
       *
       * \tparam LocalMatrix_
       * The class of the local matrix container. Usually a Tiny::Matrix
       *
       * \tparam RowMapping_
       * The class of the row-mapping for the gather-axpy operation.
       *
       * \tparam ColMapping_
       * The class of the column-mapping for the gather-axpy operation.
       *
       * \param[out] local_matrix
       * The local matrix object that serves as the gather target.
       *
       * \param[in] row_map
       * The row-mapping to be used for scattering.
       *
       * \param[in] col_map
       * The col-mapping to be used for scattering.
       *
       * \param[in] alpha
       * The scaling factor for the gather-axpy operation.
       */
      template<typename LocalData_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      LocalMatrix_& local_matrix,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
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

      template<typename LocalVector_, typename Mapping_>
      void operator()(
                      const LocalVector_& loc_vec,
                      const Mapping_& mapping,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local entries
        for(int i(0); i < mapping.get_num_local_dofs(); ++i)
        {
          // pre-multiply local entry by alpha
          DataType_ dx(alpha * loc_vec(i));

          // loop over all entry contributions
          for(int ic(0); ic < mapping.get_num_contribs(i); ++ic)
          {
            // update vector data
            _data[mapping.get_index(i, ic)] += DataType_(mapping.get_weight(i, ic)) * dx;
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

      template<typename LocalVector_, typename Mapping_>
      void operator()(
                      LocalVector_& loc_vec,
                      const Mapping_& mapping,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local entries
        for(int i(0); i < mapping.get_num_local_dofs(); ++i)
        {
          // clear accumulation entry
          DataType_ dx(DataType_(0));

          // loop over all entry contributions
          for(int ic(0); ic < mapping.get_num_contribs(i); ++ic)
          {
            // update accumulator
            dx += DataType_(mapping.get_weight(i, ic)) * _data[mapping.get_index(i, ic)];
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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      const LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row entry weight and pre-multiply by alpha
            DataType_ iw = alpha * DataType_(row_map.get_weight(i, ic));

            // fetch row index
            Index ix = row_map.get_index(i, ic);

            // build column pointer for this row entry contribution
            for(IndexType_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
            {
              _col_ptr[_col_idx[k]] = k;
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch trial function dof weight
                DataType_ jw = DataType_(col_map.get_weight(j, jc));

                // fetch column index
                Index jx = col_map.get_index(j, jc);

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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row index
            Index ix = row_map.get_index(i, ic);

            // build column pointer for this row entry contribution
            for(IndexType_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
            {
              _col_ptr[_col_idx[k]] = k;
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // clear  accumulation entry
              DataType_ dx(DataType_(0));

              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch column index
                Index jx = col_map.get_index(j, jc);

#ifdef DEBUG
                // ensure that the column pointer is valid for this index
                ASSERT(_col_ptr[jx] != _deadcode, "invalid column index");
#endif

                // update accumulator
                dx += DataType_(col_map.get_weight(j, jc)) * _data[_col_ptr[jx]];

                // continue with next column contribution
              }

              // update local matrix data
              loc_mat(i,j) += alpha * DataType_(row_map.get_weight(i, ic)) * dx;

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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      const LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row entry weight and pre-multiply by alpha
            DataType_ iw = alpha * DataType_(row_map.get_weight(i, ic));

            // fetch row index
            Index ix = row_map.get_index(i, ic);

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
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch trial function dof weight
                DataType_ jw = DataType_(col_map.get_weight(j, jc));

                // fetch column index
                Index jx = col_map.get_index(j, jc);

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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row index
            Index ix = row_map.get_index(i, ic);

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
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // clear  accumulation entry
              DataType_ dx(DataType_(0));

              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch column index
                Index jx = col_map.get_index(j, jc);

#ifdef DEBUG
                // ensure that the column pointer is valid for this index
                ASSERT(_col_ptr[jx] != _deadcode, "invalid column index");
#endif

                // update accumulator
                dx += DataType_(col_map.get_weight(j, jc)) * _data[_col_ptr[jx]];

                // continue with next column contribution
              }

              // update local matrix data
              loc_mat(i,j) += alpha * DataType_(row_map.get_weight(i, ic)) * dx;

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
      const IndexType_ * _rl;
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
        _rl(matrix.rl()),
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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      const LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row entry weight and pre-multiply by alpha
            DataType_ iw = alpha * DataType_(row_map.get_weight(i, ic));

            // fetch row index
            IndexType_ ix = IndexType_(row_map.get_index(i, ic));

            // build column pointer for this row entry contribution
            for(IndexType_ k(_cs[ix/_C] + ix%_C); k < _cs[ix/_C] + ix%_C + _rl[ix]*_C; k += _C)
            {
              _col_ptr[_col_ind[k]] = k;
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch trial function dof weight
                DataType_ jw = DataType_(col_map.get_weight(j, jc));

                // fetch column index
                Index jx = col_map.get_index(j, jc);

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
            for(IndexType_ k(_cs[ix/_C] + ix%_C); k < _cs[ix/_C] + ix%_C + _rl[ix]*_C; k += _C)
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
      const IndexType_ * _rl;
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
        _rl(matrix.rl()),
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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row index
            IndexType_ ix = IndexType_(row_map.get_index(i, ic));

            // build column pointer for this row entry contribution
            for(IndexType_ k(_cs[ix/_C] + ix%_C); k < _cs[ix/_C] + ix%_C + _rl[ix]*_C; k += _C)
            {
              _col_ptr[_col_ind[k]] = k;
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // clear  accumulation entry
              DataType_ dx(DataType_(0));

              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch column index
                Index jx = col_map.get_index(j, jc);

#ifdef DEBUG
                // ensure that the column pointer is valid for this index
                ASSERT(_col_ptr[jx] != _deadcode, "invalid column index");
#endif

                // update accumulator
                dx += DataType_(col_map.get_weight(j, jc)) * _val[_col_ptr[jx]];

                // continue with next column contribution
              }

              // update local matrix data
              loc_mat(i,j) += alpha * DataType_(row_map.get_weight(i, ic)) * dx;

              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(IndexType_ k(_cs[ix/_C] + ix%_C); k < _cs[ix/_C] + ix%_C + _rl[ix]*_C; k += _C)
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

    /**
     * \brief Scatter-Axpy specialisation for SparseMatrixBanded
     *
     * \author Christoph Lohmann
     */
    template<typename DataType_, typename IndexType_>
    class ScatterAxpy< LAFEM::SparseMatrixBanded<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::SparseMatrixBanded<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      Index _num_of_offsets;
      IndexType_* _offsets;
      IndexType_* _col_ptr;
      DataType_ *_data;

    public:
      explicit ScatterAxpy(MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~IndexType_(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _num_of_offsets(matrix.num_of_offsets()),
        _offsets(matrix.offsets()),
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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      const LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row entry weight and pre-multiply by alpha
            DataType_ iw = alpha * DataType_(row_map.get_weight(i, ic));

            // fetch row index
            Index ix = row_map.get_index(i, ic);

            // build column pointer for this row entry contribution
            for(IndexType_ k(0); k < _num_of_offsets; ++k)
            {
              if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
              {
                _col_ptr[_offsets[k] + ix + 1 - _num_rows] = IndexType_(k * _num_rows + ix);
              }
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch trial function dof weight
                DataType_ jw = DataType_(col_map.get_weight(j, jc));

                // fetch column index
                Index jx = col_map.get_index(j, jc);

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
            for(IndexType_ k(0); k < _num_of_offsets; ++k)
            {
              if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
              {
                _col_ptr[_offsets[k] + ix + 1 - _num_rows] = _deadcode;
              }
            }
#endif
            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class ScatterAxpy<LAFEM::SparseMatrixBanded<Mem::Main,...>>

    /**
     * \brief Gather-Axpy specialisation for SparseMatrixBanded
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_>
    class GatherAxpy< LAFEM::SparseMatrixBanded<Mem::Main, DataType_, IndexType_> >
    {
    public:
      typedef LAFEM::SparseMatrixBanded<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      Index _num_of_offsets;
      const IndexType_* _offsets;
      IndexType_* _col_ptr;
      const DataType_ *_data;

    public:
      explicit GatherAxpy(const MatrixType& matrix) :
#ifdef DEBUG
        _deadcode(~IndexType_(0)),
#endif
        _num_rows(matrix.rows()),
        _num_cols(matrix.columns()),
        _num_of_offsets(matrix.num_of_offsets()),
        _offsets(matrix.offsets()),
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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row index
            Index ix = row_map.get_index(i, ic);

            // build column pointer for this row entry contribution
            for(IndexType_ k(0); k < _num_of_offsets; ++k)
            {
              if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
              {
                _col_ptr[_offsets[k] + ix + 1 - _num_rows] = IndexType_(k * _num_rows + ix);
              }
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // clear  accumulation entry
              DataType_ dx(DataType_(0));

              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch column index
                Index jx = col_map.get_index(j, jc);

#ifdef DEBUG
                // ensure that the column pointer is valid for this index
                ASSERT(_col_ptr[jx] != _deadcode, "invalid column index");
#endif

                // update accumulator
                dx += DataType_(col_map.get_weight(j, jc)) * _data[_col_ptr[jx]];

                // continue with next column contribution
              }

              // update local matrix data
              loc_mat(i,j) += alpha * DataType_(row_map.get_weight(i, ic)) * dx;

              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(IndexType_ k(0); k < _num_of_offsets; ++k)
            {
              if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
              {
                _col_ptr[_offsets[k] + ix + 1 - _num_rows] = _deadcode;
              }
            }
#endif

            // continue with next row contribution
          }
          // continue with next row entry
        }
      }
    }; // class GatherAxpy<LAFEM::SparseMatrixBanded<Mem::Main,...>>

    /**
     * \brief Scatter-Axpy specialisation for SparseMatrixCSRBlocked
     *
     * Apart from the MatrixType and usage of DataTypeBlocked c&p from SparseMatrixCSR.
     *
     * \author Jordi Paul
     */
    template<typename DataType_, typename IndexType_, int BlockHeight_, int BlockWidth_>
    class ScatterAxpy< LAFEM::SparseMatrixCSRBlocked<Mem::Main, DataType_, IndexType_, BlockHeight_, BlockWidth_> >
    {
    public:
      typedef LAFEM::SparseMatrixCSRBlocked<Mem::Main, DataType_, IndexType_, BlockHeight_, BlockWidth_> MatrixType;
      typedef Mem::Main MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;
      typedef Tiny::Matrix<DataType_, BlockHeight_, BlockWidth_> DataTypeBlocked;

    private:
#ifdef DEBUG
      IndexType_ _deadcode;
#endif
      Index _num_rows;
      Index _num_cols;
      IndexType_* _row_ptr;
      IndexType_* _col_idx;
      IndexType_* _col_ptr;
      DataTypeBlocked *_data;

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

      template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
      void operator()(
                      const LocalMatrix_& loc_mat,
                      const RowMapping_& row_map,
                      const ColMapping_& col_map,
                      DataType_ alpha = DataType_(1))
      {
        // loop over all local row entries
        for(int i(0); i < row_map.get_num_local_dofs(); ++i)
        {
          // loop over all row entry contributations
          for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
          {
            // fetch row entry weight and pre-multiply by alpha
            DataType_ iw = alpha * DataType_(row_map.get_weight(i, ic));

            // fetch row index
            Index ix = row_map.get_index(i, ic);

            // build column pointer for this row entry contribution
            for(IndexType_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
            {
              _col_ptr[_col_idx[k]] = k;
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // loop over all column entry contributions
              for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
              {
                // fetch trial function dof weight
                DataType_ jw = DataType_(col_map.get_weight(j, jc));

                // fetch column index
                Index jx = col_map.get_index(j, jc);

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
    }; // class ScatterAxpy<LAFEM::SparseMatrixCSRBlocked<Mem::Main,...>>

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_GATHER_SCATTER_HPP
