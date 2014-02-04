#pragma once
#ifndef KERNEL_LAFEM_UNIT_FILTER_HPP
#define KERNEL_LAFEM_UNIT_FILTER_HPP 1

// includes, FEAST
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Unit Filter class template.
     *
     * \todo Replace the internal data storage by SparseVector once it is implemented.
     *
     * \author Peter Zajac
     */
    template<
      typename MemType_,
      typename DataType_>
    class UnitFilter
    {
    public:
      /// mem-type typedef
      typedef MemType_ MemType;
      /// data-type typedef
      typedef DataType_ DataType;

    protected:
      /// total number of entries to be filtered
      Index _num_entries;
      /// index array
      Index* _indices;
      /// value array
      DataType* _values;

    public:
      /// default constructor
      UnitFilter() :
        _num_entries(0),
        _indices(nullptr),
        _values(nullptr)
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] num_entries
       * The total number of entries for the unit filter.
       */
      explicit UnitFilter(Index num_entries) :
        _num_entries(num_entries),
        _indices(new Index[num_entries]),
        _values(new DataType[num_entries])
      {
      }

      /// move-ctor
      UnitFilter(UnitFilter&& other) :
        _num_entries(other._num_entries),
        _indices(other._indices),
        _values(other._values)
      {
        other._num_entries = Index(0);
        other._indices = nullptr;
        other._values = nullptr;
      }

      /// move-assignment operator
      UnitFilter& operator=(UnitFilter&& other)
      {
        if(this == &other)
          return *this;

        if(_indices != nullptr)
          delete [] _indices;
        if(_values != nullptr)
          delete [] _values;

        _num_entries = other._num_entries;
        _indices = other._indices;
        _values = other._values;
        other._num_entries = Index(0);
        other._indices = nullptr;
        other._values = nullptr;

        return *this;
      }

      /// virtual destructor
      virtual ~UnitFilter()
      {
        if(_values != nullptr)
        {
          delete [] _values;
        }
        if(_indices != nullptr)
        {
          delete [] _indices;
        }
      }

      UnitFilter clone() const
      {
        // create copy
        UnitFilter other(size());

        // copy arrays
        Index* indices(other.get_indices());
        DataType* values(other.get_values());
        for(Index i(0); i < _num_entries; ++i)
        {
          indices[i] = _indices[i];
          values[i] = _values[i];
        }

        // return
        return std::move(other);
      }

      /// \returns The number of entries in the filter.
      Index size() const
      {
        return _num_entries;
      }

      /// \returns The index array.
      Index* get_indices()
      {
        return _indices;
      }

      /// \returns The index array.
      const Index* get_indices() const
      {
        return _indices;
      }

      /// \returns The value array.
      DataType* get_values()
      {
        return _values;
      }

      /// \returns The value array.
      const DataType* get_values() const
      {
        return _values;
      }

      /**
       * \brief Applies the filter onto a system matrix.
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       */
      template<typename Algo_>
      void filter_mat(SparseMatrixCSR<MemType, DataType>& matrix) const
      {
        const Index* row_ptr(matrix.row_ptr());
        const Index* col_idx(matrix.col_ind());
        DataType* v(matrix.val());

        for(Index i(0); i < _num_entries; ++i)
        {
          Index ix(_indices[i]);
          // replace by unit row
          for(Index j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            v[j] = (col_idx[j] == ix) ? DataType(1) : DataType(0);
          }
        }
      }

      template<typename Algo_>
      void filter_offdiag_row_mat(SparseMatrixCSR<MemType, DataType>& matrix) const
      {
        const Index* row_ptr(matrix.row_ptr());
        DataType* v(matrix.val());

        for(Index i(0); i < _num_entries; ++i)
        {
          Index ix(_indices[i]);
          // replace by null row
          for(Index j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            v[j] = DataType(0);
          }
        }
      }

      template<typename Algo_>
      void filter_offdiag_col_mat(SparseMatrixCSR<MemType, DataType>&) const
      {
        // nothing to do here
      }

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      template<typename Algo_>
      void filter_rhs(DenseVector<MemType, DataType>& vector) const
      {
        DataType* v(vector.elements());
        for(Index i(0); i < _num_entries; ++i)
        {
          v[_indices[i]] = _values[i];
        }
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      template<typename Algo_>
      void filter_sol(DenseVector<MemType, DataType>& vector) const
      {
        // same as rhs
        filter_rhs<Algo_>(vector);
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      template<typename Algo_>
      void filter_def(DenseVector<MemType, DataType>& vector) const
      {
        DataType* v(vector.elements());
        for(Index i(0); i < _num_entries; ++i)
        {
          v[_indices[i]] = DataType(0);
        }
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      template<typename Algo_>
      void filter_cor(DenseVector<MemType, DataType>& vector) const
      {
        // same as def
        filter_def<Algo_>(vector);
      }
    }; // class UnitFilter<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_UNIT_FILTER_HPP
