#pragma once
#ifndef KERNEL_LAFEM_UNIT_FILTER_HPP
#define KERNEL_LAFEM_UNIT_FILTER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_vector.hpp>

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
      typename Mem_,
      typename DT_>
    class UnitFilter
    {
    public:
      /// mem-type typedef
      typedef Mem_ MemType;
      /// data-type typedef
      typedef DT_ DataType;

    private:
      /// SparseVector, containing all entries of the unit filter
      SparseVector<Mem_, DT_> _sv;

    public:
      /// default constructor
      UnitFilter() :
        _sv()
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] num_entries
       * The total number of entries for the unit filter.
       */
      explicit UnitFilter(Index num_entries) :
        _sv(num_entries)
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] values DenseVector containing element values
       * \param[in] indices DenseVector containing element indices
       */
      explicit UnitFilter(DenseVector<Mem_, DT_> & values, DenseVector<Mem_, Index> & indices) :
        _sv(values.size(), values, indices)
      {
        if (values.size() != indices.size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size missmatch!");
      }

      /// move-ctor
      UnitFilter(UnitFilter && other) :
        _sv(std::move(other._sv))
      {
      }

      /// move-assignment operator
      UnitFilter & operator=(UnitFilter && other)
      {
        if(this != &other)
        {
          _sv = std::forward<UnitFilter>(other._sv);
        }
        return *this;
      }

      /// virtual destructor
      virtual ~UnitFilter()
      {
      }

      void clone(const UnitFilter & other)
      {
        _sv.clone(other._sv);
      }

      /// \returns The number of entries in the filter.
      Index size() const
      {
        return _sv.used_elements();
      }

      /// \returns The index array.
      Index* get_indices()
      {
        return _sv.indices();
      }

      /// \returns The index array.
      const Index* get_indices() const
      {
        return _sv.indices;
      }

      /// \returns The value array.
      DT_* get_values()
      {
        return _sv.elements();
      }

      /// \returns The value array.
      const DT_* get_values() const
      {
        return _sv.elements();
      }

      /**
       * \brief Applies the filter onto a system matrix.
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       */
      template<typename Algo_>
      void filter_mat(SparseMatrixCSR<Mem_, DT_> & matrix) const
      {
        const Index* row_ptr(matrix.row_ptr());
        const Index* col_idx(matrix.col_ind());
        DT_* v(matrix.val());

        for(Index i(0); i < _sv.used_elements(); ++i)
        {
          Index ix(_sv.indices()[i]);
          // replace by unit row
          for(Index j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            v[j] = (col_idx[j] == ix) ? DT_(1) : DT_(0);
          }
        }
      }

      template<typename Algo_>
      void filter_offdiag_row_mat(SparseMatrixCSR<Mem_, DT_> & matrix) const
      {
        const Index* row_ptr(matrix.row_ptr());
        DT_* v(matrix.val());

        for(Index i(0); i < _sv.used_elements(); ++i)
        {
          Index ix(_sv.indices()[i]);
          // replace by null row
          for(Index j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            v[j] = DT_(0);
          }
        }
      }

      template<typename Algo_>
      void filter_offdiag_col_mat(SparseMatrixCSR<Mem_, DT_> &) const
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
      void filter_rhs(DenseVector<Mem_, DT_> & vector) const
      {
        DT_* v(vector.elements());
        for(Index i(0); i < _sv.used_elements(); ++i)
        {
          v[_sv.indices()[i]] = _sv.elements()[i];
        }
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      template<typename Algo_>
      void filter_sol(DenseVector<Mem_, DT_> & vector) const
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
      void filter_def(DenseVector<Mem_, DT_> & vector) const
      {
        DT_* v(vector.elements());
        for(Index i(0); i < _sv.used_elements(); ++i)
        {
          v[_sv.indices()[i]] = DT_(0);
        }
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      template<typename Algo_>
      void filter_cor(DenseVector<Mem_, DT_> & vector) const
      {
        // same as def
        filter_def<Algo_>(vector);
      }
    }; // class UnitFilter<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_UNIT_FILTER_HPP
