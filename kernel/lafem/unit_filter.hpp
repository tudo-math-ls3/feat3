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
     * \author Peter Zajac
     */
    template<
      typename MemType_,
      typename DataType_>
    class UnitFilter
    {
    public:
      /// arch typedef
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

      /// copy-ctor
      UnitFilter(const UnitFilter& other) :
        _num_entries(other._num_entries),
        _indices(new Index[_num_entries]),
        _values(new DataType[_num_entries])
      {
        const Index* indices(other.get_indices());
        const DataType* values(other.get_values());
        for(Index i(0); i < _num_entries; ++i)
        {
          _indices[i] = indices[i];
          _values[i] = values[i];
        }
      }

      /// assignment operator
      UnitFilter& operator=(const UnitFilter& other)
      {
        if(this == &other)
          return *this;

        if(_indices != nullptr)
          delete [] _indices;
        if(_values != nullptr)
          delete [] _values;
        _num_entries = other._num_entries;
        _indices = new Index[_num_entries];
        _values = new DataType[_num_entries];
        const Index* indices(other.get_indices());
        const DataType* values(other.get_values());
        for(Index i(0); i < _num_entries; ++i)
        {
          _indices[i] = indices[i];
          _values[i] = values[i];
        }
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
      template<typename DataType2_>
      void filter_mat(SparseMatrixCSR<MemType, DataType2_>& matrix) const
      {
        const Index* row_ptr(matrix.row_ptr());
        const Index* col_idx(matrix.col_ind());
        DataType2_* v(matrix.val());

        for(Index i(0); i < _num_entries; ++i)
        {
          Index ix(_indices[i]);
          // replace by unit row
          for(Index j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            v[j] = (col_idx[j] == ix) ? DataType2_(1) : DataType2_(0);
          }
        }
      }

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      template<typename DataType2_>
      void filter_rhs(DenseVector<MemType, DataType2_>& vector) const
      {
        DataType2_* v(vector.elements());
        for(Index i(0); i < _num_entries; ++i)
        {
          v[_indices[i]] = DataType2_(_values[i]);
        }
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      template<typename DataType2_>
      void filter_sol(DenseVector<MemType, DataType2_>& vector) const
      {
        DataType2_* v(vector.elements());
        for(Index i(0); i < _num_entries; ++i)
        {
          v[_indices[i]] = DataType2_(_values[i]);
        }
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      template<typename DataType2_>
      void filter_def(DenseVector<MemType, DataType2_>& vector) const
      {
        DataType2_* v(vector.elements());
        for(Index i(0); i < _num_entries; ++i)
        {
          v[_indices[i]] = DataType2_(0);
        }
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      template<typename DataType2_>
      void filter_cor(DenseVector<MemType, DataType2_>& vector) const
      {
        DataType2_* v(vector.elements());
        for(Index i(0); i < _num_entries; ++i)
        {
          v[_indices[i]] = DataType2_(0);
        }
      }
    }; // class UnitFilter<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_UNIT_FILTER_HPP
