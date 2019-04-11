// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_UNIT_FILTER_BLOCKED_HPP
#define KERNEL_LAFEM_UNIT_FILTER_BLOCKED_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/lafem/arch/unit_filter_blocked.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Unit Filter Blocked class template.
     *
     * \tparam Mem_
     * Memory architecture
     *
     * \tparam DT_
     * Data type, i.e. double
     *
     * \tparam IT_
     * Index type, i.e. unsigned int
     *
     * \tparam BlockSize_
     * Size of the blocks, i.e. 2 for a filter acting on a velocity field in 2d
     *
     * Mostly c&p from UnitFilter
     *
     * \author Jordi Paul
     */
    template<
      typename Mem_,
      typename DT_,
      typename IT_,
      int BlockSize_>
    class UnitFilterBlocked
    {
    public:
      /// mem-type typedef
      typedef Mem_ MemType;
      /// data-type typedef
      typedef DT_ DataType;
      /// index-type typedef
      typedef IT_ IndexType;
      /// The block size
      static constexpr int BlockSize = BlockSize_;
      /// Value type
      typedef Tiny::Vector<DataType, BlockSize> ValueType;
      /// Our supported vector type
      typedef DenseVectorBlocked<MemType, DataType, IndexType, BlockSize> VectorType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using FilterType = UnitFilterBlocked<Mem2_, DT2_, IT2_, BlockSize_>;

      /// this typedef lets you create a filter with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using FilterTypeByMDI = FilterType<Mem2_, DataType2_, IndexType2_>;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

      static_assert(BlockSize > 1, "BlockSize has to be >= 2 in UnitFilterBlocked!");

    private:
      /// SparseVector, containing all entries of the unit filter
      SparseVectorBlocked<MemType, DataType, IndexType, BlockSize> _sv;

    public:
      /// default constructor
      UnitFilterBlocked() :
        _sv()
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] num_entries
       * The total number of entries for the unit filter.
       */
      explicit UnitFilterBlocked(Index num_entries) :
        _sv(num_entries)
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] size_in The size of the created filter.
       * \param[in] values DenseVector containing element values
       * \param[in] indices DenseVector containing element indices
       */
      explicit UnitFilterBlocked(Index size_in,
                                 DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_> & values,
                                 DenseVector<Mem_, IT_, IT_> & indices) :
        _sv(size_in, values, indices)
      {
        XASSERTM(values.size() == indices.size(), "Vector size mismatch!");
      }

      /// move-ctor
      UnitFilterBlocked(UnitFilterBlocked && other) :
        _sv(std::move(other._sv))
      {
      }

      /// move-assignment operator
      UnitFilterBlocked & operator=(UnitFilterBlocked && other)
      {
        if(this != &other)
        {
          _sv = std::forward<decltype(other._sv)>(other._sv);
        }
        return *this;
      }

      /// virtual destructor
      virtual ~UnitFilterBlocked()
      {
      }

      /// \brief Creates a clone of itself
      UnitFilterBlocked clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        UnitFilterBlocked other;
        other.clone(*this, clone_mode);
        return std::move(other);
      }

      /// \brief Clones data from another UnitFilterBlocked
      void clone(const UnitFilterBlocked & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _sv.clone(other.get_filter_vector(), clone_mode);
      }

      /// \brief Converts data from another UnitFilter
      template<typename Mem2_, typename DT2_, typename IT2_, int BS_>
      void convert(const UnitFilterBlocked<Mem2_, DT2_, IT2_, BS_>& other)
      {
        _sv.convert(other.get_filter_vector());
      }

      /// \brief Clears the underlying data (namely the SparseVector)
      void clear()
      {
        _sv.clear();
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _sv.bytes();
      }

      /// \cond internal
      SparseVectorBlocked<Mem_, DT_, IT_, BlockSize>& get_filter_vector()
      {
        return _sv;
      }
      const SparseVectorBlocked<Mem_, DT_, IT_, BlockSize>& get_filter_vector() const
      {
        return _sv;
      }
      /// \endcond

      /**
       * \brief Adds one element to the filter
       *
       * \param[in] idx Index where to add
       * \param[in] val Value to add
       *
       **/
      void add(IndexType idx, ValueType val)
      {
        _sv(idx, val);
      }

      /// \returns The total size of the filter.
      Index size() const
      {
        return _sv.size();
      }

      /// \returns The number of entries in the filter.
      Index used_elements() const
      {
        return _sv.used_elements();
      }

      /// \returns The index array.
      IT_* get_indices()
      {
        return _sv.indices();
      }

      /// \returns The index array.
      const IT_* get_indices() const
      {
        return _sv.indices();
      }

      /// \returns The value array.
      ValueType* get_values()
      {
        return _sv.elements();
      }

      /// \returns The value array.
      const ValueType* get_values() const
      {
        return _sv.elements();
      }

      /// Permutate internal vector according to the given Permutation
      void permute(Adjacency::Permutation & perm)
      {
        _sv.permute(perm);
      }

#ifdef DOXYGEN
      // The following documentation block is visible to Doxygen only. The actual implementation is matrix type
      // specific and provided below.

      /**
       * \brief Applies the filter onto a system matrix.
       *
       * \tparam BlockWidth_
       * The input matrix' block width
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       * The input matrix has to have a block(ed) structure and its BlockHeight has to agree with the filter's
       * blocksize.
       *
       */
      void filter_mat(MatrixType & matrix) const
      {
      }

      /**
       * \brief Filter the non-diagonal row entries
       *
       * \tparam BlockWidth_
       * The input matrix' block width
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       * The input matrix has to have a block(ed) structure and its BlockHeight has to agree with the filter's
       * blocksize.
       *
       */
      void filter_offdiag_row_mat(MatrixType & matrix) const
      {
      }

      /**
       * \brief Filter the non-diagonal column entries
       *
       * \tparam BlockWidth_
       * The input matrix' block width
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       * The input matrix has to have a block(ed) structure and its BlockHeight has to agree with the filter's
       * blocksize.
       *
       */
      void filter_offdiag_col_mat(MatrixType & matrix) const
      {
      }

#endif
      /// \cond internal
      template<int BlockWidth_>
      void filter_mat(SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockSize_, BlockWidth_> & matrix) const
      {
        XASSERTM(_sv.size() == matrix.rows(), "Matrix size does not match!");

        const IT_* row_ptr(matrix.row_ptr());
        const IT_* col_idx(matrix.col_ind());
        typename SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockSize_, BlockWidth_>::ValueType* v(matrix.val());

        for(Index i(0); i < _sv.used_elements(); ++i)
        {
          IT_ ix(_sv.indices()[i]);
          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            v[j].format(DT_(0));
            if(col_idx[j] == ix)
            {
              for(int l(0); l < Math::min(BlockSize_, BlockWidth_); ++l)
                v[j](l,l) = DT_(1);
            }
          }
        }
      }

      template<int BlockWidth_>
      void filter_offdiag_row_mat(SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockSize_, BlockWidth_> & matrix) const
      {
        XASSERTM(_sv.size() == matrix.rows(), "Matrix size does not match!");

        const IT_* row_ptr(matrix.row_ptr());
        auto* v(matrix.val());

        for(Index i(0); i < _sv.used_elements(); ++i)
        {
          IT_ ix(_sv.indices()[i]);
          // replace by null row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            v[j].format(DT_(0));
          }
        }
      }

      template<int BlockWidth_>
      void filter_offdiag_col_mat(SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockSize_, BlockWidth_> &) const
      {
        // nothing to do here
      }
      /// \endcond

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType& vector) const
      {
        XASSERTM(_sv.size() == vector.size(), "Vector size does not match!");
        if(_sv.used_elements() > Index(0))
          Arch::UnitFilterBlocked<Mem_>::template filter_rhs<DT_, IT_, BlockSize_>
            (vector.template elements<Perspective::pod>(), _sv.template elements<Perspective::pod>(), _sv.indices(), _sv.used_elements());
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      void filter_sol(VectorType& vector) const
      {
        // same as rhs
        filter_rhs(vector);
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      void filter_def(VectorType& vector) const
      {
        XASSERTM(_sv.size() == vector.size(), "Vector size does not match!");
        if(_sv.used_elements() > Index(0))
          Arch::UnitFilterBlocked<Mem_>::template filter_def<DT_, IT_, BlockSize_>
            (vector.template elements<Perspective::pod>(), _sv.indices(), _sv.used_elements() );
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      void filter_cor(VectorType& vector) const
      {
        // same as def
        filter_def(vector);
      }
    }; // class UnitFilterBlocked<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_UNIT_FILTER_BLOCKED_HPP
