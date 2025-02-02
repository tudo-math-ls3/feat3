// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/math.hpp>
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
     * \note
     * This class allows to emulate a "slip-filter-like" behavior by disabling individual filter
     * components on a per-DOF basis by turning on the "ignore NaNs" functionality and setting all
     * component values, which should be ignored by the filter, to NaN. Example: if you want to
     * emulate a 3D slip-filter for the Y-component, then set the filter value (NaN, 0, NaN) for
     * all DOFs that should be affected by the filter.
     *
     * \author Jordi Paul, Peter Zajac
     */
    template<
      typename DT_,
      typename IT_,
      int BlockSize_>
    class UnitFilterBlocked
    {
    public:
      /// data-type typedef
      typedef DT_ DataType;
      /// index-type typedef
      typedef IT_ IndexType;
      /// The block size
      static constexpr int BlockSize = BlockSize_;
      /// Value type
      typedef Tiny::Vector<DataType, BlockSize> ValueType;
      /// Our supported vector type
      typedef DenseVectorBlocked<DataType, IndexType, BlockSize> VectorType;

      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using FilterType = UnitFilterBlocked<DT2_, IT2_, BlockSize_>;

      /// this typedef lets you create a filter with new Datatype and Index types
      template <typename DataType2_, typename IndexType2_>
      using FilterTypeByDI = FilterType<DataType2_, IndexType2_>;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

      static_assert(BlockSize > 1, "BlockSize has to be >= 2 in UnitFilterBlocked!");

    private:
      /// SparseVector, containing all entries of the unit filter
      SparseVectorBlocked<DataType, IndexType, BlockSize> _sv;

      /// ignore NaNs in filter values
      bool _ignore_nans;

    public:
      /// default constructor
      UnitFilterBlocked() :
        _sv(),
        _ignore_nans(false)
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] num_entries
       * The total number of entries for the unit filter.
       *
       * \param[in] ingore_nans
       * Specifies whether the filter should ignore NaNs filter values
       */
      explicit UnitFilterBlocked(Index num_entries, bool ignore_nans = false) :
        _sv(num_entries),
        _ignore_nans(ignore_nans)
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
                                 DenseVectorBlocked<DT_, IT_, BlockSize_> & values,
                                 DenseVector<IT_, IT_> & indices) :
        _sv(size_in, values, indices)
      {
        XASSERTM(values.size() == indices.size(), "Vector size mismatch!");
      }

      /// move-ctor
      UnitFilterBlocked(UnitFilterBlocked && other) :
        _sv(std::move(other._sv)),
        _ignore_nans(other._ignore_nans)
      {
      }

      /// move-assignment operator
      UnitFilterBlocked & operator=(UnitFilterBlocked && other)
      {
        if(this != &other)
        {
          _sv = std::forward<decltype(other._sv)>(other._sv);
          _ignore_nans = other._ignore_nans;
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
        return other;
      }

      /// \brief Clones data from another UnitFilterBlocked
      void clone(const UnitFilterBlocked & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _sv.clone(other.get_filter_vector(), clone_mode);
        _ignore_nans = other.get_ignore_nans();
      }

      /// \brief Converts data from another UnitFilter
      template<typename DT2_, typename IT2_, int BS_>
      void convert(const UnitFilterBlocked<DT2_, IT2_, BS_>& other)
      {
        _sv.convert(other.get_filter_vector());
        _ignore_nans = other.get_ignore_nans();
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
      SparseVectorBlocked<DT_, IT_, BlockSize>& get_filter_vector()
      {
        return _sv;
      }
      const SparseVectorBlocked<DT_, IT_, BlockSize>& get_filter_vector() const
      {
        return _sv;
      }
      /// \endcond

      /**
       * \brief Specifies whether the filter should ignore NaN filter values
       *
       * \param[in] ingore_nans
       * Specifies whether the filter should ignore NaNs filter values
       */
      void set_ignore_nans(bool ignore_nans)
      {
        _ignore_nans = ignore_nans;
      }

      /// Specifies whether the filter should ignore NaNs filter values
      bool get_ignore_nans() const
      {
        return _ignore_nans;
      }

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
      void filter_mat(SparseMatrixBCSR<DT_, IT_, BlockSize_, BlockWidth_> & matrix) const
      {
        if(_sv.used_elements() == Index(0))
          return;

        XASSERTM(_sv.size() == matrix.rows(), "Matrix size does not match!");

        Arch::UnitFilterBlocked::template filter_unit_mat
          (matrix.template val<LAFEM::Perspective::pod>(), matrix.row_ptr(), matrix.col_ind(), BlockSize_, BlockWidth_,
          _sv.template elements<LAFEM::Perspective::pod>(), _sv.indices(), _sv.used_elements(), _ignore_nans);
      }

      template<int BlockWidth_>
      void filter_offdiag_row_mat(SparseMatrixBCSR<DT_, IT_, BlockSize_, BlockWidth_> & matrix) const
      {
        if(_sv.used_elements() == Index(0))
          return;

        XASSERTM(_sv.size() == matrix.rows(), "Matrix size does not match!");

        Arch::UnitFilterBlocked::template filter_offdiag_row_mat
          (matrix.template val<LAFEM::Perspective::pod>(), matrix.row_ptr(), BlockSize_, BlockWidth_,
          _sv.template elements<LAFEM::Perspective::pod>(), _sv.indices(), _sv.used_elements(), _ignore_nans);
      }

      template<int BlockWidth_>
      void filter_offdiag_col_mat(SparseMatrixBCSR<DT_, IT_, BlockSize_, BlockWidth_> &) const
      {
        // nothing to do here
      }
      /// \endcond


      /**
       * \brief Replaces the rows of the system matrix by scaled rows of another matrix
       *
       * This function replaces all rows of the system matrix A, whose row index is included in this filter's
       * indices, by the corresponding rows of a given donor matrix M, which is typically a mass matrix, scaled
       * by the values which are stored in this filter.
       * This functionality can be used to implement a weak form of Dirichlet boundary conditions, which is also
       * used to employ mass-conserving fictitious boundary conditions in a Stokes system.
       *
       * \param[inout] matrix_a
       * The system matrix whose rows are to be replaced
       *
       * \param[in] matrix_m
       * The donor matrix whose rows are to be copied into the system matrix
       */
      template<int BlockWidth_>
      void filter_weak_matrix_rows(SparseMatrixBCSR<DT_, IT_, BlockSize_, BlockWidth_> & matrix_a,
        const SparseMatrixBCSR<DT_, IT_, BlockSize_, BlockWidth_>& matrix_m) const
      {
        if(_sv.used_elements() == Index(0))
          return;

        if(matrix_a.val() == matrix_m.val())
        {
          XABORTM("Matrices are not allowed to hold the same data");
        }

        XASSERTM(_sv.size() == matrix_a.rows(), "Matrix size does not match!");
        XASSERTM(_sv.size() == matrix_m.rows(), "Matrix size does not match!");

        const IndexType* row_ptr(matrix_a.row_ptr());
        const IndexType* col_idx(matrix_m.col_ind());
        XASSERTM(row_ptr == matrix_m.row_ptr(), "matrix A and M must share their layout");
        XASSERTM(col_idx == matrix_m.col_ind(), "matrix A and M must share their layout");

        Arch::UnitFilterBlocked::template filter_weak_matrix_rows
          (matrix_a.template val<LAFEM::Perspective::pod>(), matrix_m.template val<LAFEM::Perspective::pod>(), row_ptr, BlockSize_, BlockWidth_,
          _sv.template elements<LAFEM::Perspective::pod>(), _sv.indices(), _sv.used_elements());
      }


      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType& vector) const
      {
        if(_sv.used_elements() == Index(0))
          return;

        XASSERTM(_sv.size() == vector.size(), "Vector size does not match!");
        Arch::UnitFilterBlocked::template filter_rhs
          (vector.template elements<Perspective::pod>(), BlockSize_, _sv.template elements<Perspective::pod>(),
          _sv.indices(), _sv.used_elements(), _ignore_nans);
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
        if(_sv.used_elements() == Index(0))
          return;

        XASSERTM(_sv.size() == vector.size(), "Vector size does not match!");
        Arch::UnitFilterBlocked::template filter_def
          (vector.template elements<Perspective::pod>(), BlockSize_, _sv.template elements<Perspective::pod>(),
          _sv.indices(), _sv.used_elements(), _ignore_nans);
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
