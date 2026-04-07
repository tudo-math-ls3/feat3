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
#include <kernel/lafem/arch/unit_filter_block_bcsr.hpp>
#include <kernel/lafem/arch/unit_filter_block_vec.hpp>
#include <kernel/lafem/arch/unit_filter_block_weak_bcsr.hpp>

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
     * \tparam block_size_
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
      int block_size_>
    class UnitFilterBlocked
    {
    public:
      /// data-type typedef
      typedef DT_ DataType;
      /// index-type typedef
      typedef IT_ IndexType;
      /// The block size
      static constexpr int block_size = block_size_;
      /// Value type
      typedef Tiny::Vector<DataType, block_size> ValueType;
      /// Our supported vector type
      typedef DenseVectorBlocked<DataType, IndexType, block_size> VectorType;

      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using FilterType = UnitFilterBlocked<DT2_, IT2_, block_size_>;

      /// this typedef lets you create a filter with new Datatype and Index types
      template <typename DataType2_, typename IndexType2_>
      using FilterTypeByDI = FilterType<DataType2_, IndexType2_>;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

      static_assert(block_size > 1, "block_size has to be >= 2 in UnitFilterBlocked!");

    private:
      /// SparseVector, containing all entries of the unit filter
      SparseVectorBlocked<DataType, IndexType, block_size> _sv;

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
      explicit UnitFilterBlocked(Index num_entries, Index num_nonzeros, bool ignore_nans = false) :
        _sv(num_entries, num_nonzeros),
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
                                 DenseVectorBlocked<DT_, IT_, block_size_> & values,
                                 DenseVector<IT_, IT_> & indices, bool ignore_nans = false) :
        _sv(size_in, values, indices),
        _ignore_nans(ignore_nans)
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
      SparseVectorBlocked<DT_, IT_, block_size>& get_filter_vector()
      {
        return _sv;
      }
      const SparseVectorBlocked<DT_, IT_, block_size>& get_filter_vector() const
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

      /// Returns the total native size of the filter
      Index size() const
      {
        return _sv.size();
      }

      /// Returns the total raw size of the filter
      Index size_raw() const
      {
        return _sv.size_raw();
      }

      /// Returns the number of non-zero entries in the filter
      Index num_nzes() const
      {
        return _sv.num_nzes();
      }

      /// Returns the raw number of non-zero entries in the filter
      Index num_nzes_raw() const
      {
        return _sv.num_nzes_raw();
      }

      /// Checks whether the size of the filter is zero
      bool empty() const
      {
        return _sv.empty();
      }

      /// Checks whether the filter does not contain any non-zero elements
      bool hollow() const
      {
        return _sv.hollow();
      }

      /// Returns a reference to the elements memory arbiter
      Memory::Arbiter& elements_arbiter()
      {
        return _sv.elements_arbiter();
      }

      /// Returns a const reference to the elements memory arbiter
      const Memory::Arbiter& elements_arbiter() const
      {
        return _sv.elements_arbiter();
      }

      /// Creates a read-only memory view for the elements array for a given memory location
      Memory::TypedView<ValueType> elements_view_r(Memory::Location loc = Memory::Location::main) const
      {
        return _sv.elements_view_r(loc);
      }

      /// Creates a write-only memory view for the elements array for a given memory location
      Memory::TypedView<ValueType> elements_view_w(Memory::Location loc = Memory::Location::main)
      {
        return _sv.elements_view_w(loc);
      }

      /// Creates a read-write memory view for the elements array for a given memory location
      Memory::TypedView<ValueType> elements_view_rw(Memory::Location loc = Memory::Location::main)
      {
        return _sv.elements_view_rw(loc);
      }

      /**
       * \brief Creates a memory view for the elements array
       *
       * \param[in] loc
       * The memory location for which the view is to be created
       *
       * \param[in] acc
       * A combination of access rights to grant for the view
       */
      Memory::TypedView<ValueType> elements_view(Memory::Location loc, Memory::Access acc)
      {
        return _sv.elements_view(loc, acc);
      }

      /// Returns a reference to the indices memory arbiter
      Memory::Arbiter& indices_arbiter()
      {
        return _sv.indices_arbiter();
      }

      /// Returns a const reference to the indices memory arbiter
      const Memory::Arbiter& indices_arbiter() const
      {
        return _sv.indices_arbiter();
      }

      /// Creates a read-only memory view for the indices array for a given memory location
      Memory::TypedView<IT_> indices_view_r(Memory::Location loc = Memory::Location::main) const
      {
        return _sv.indices_view_r(loc);
      }

      /// Creates a write-only memory view for the indices array for a given memory location
      Memory::TypedView<IT_> indices_view_w(Memory::Location loc = Memory::Location::main)
      {
        return _sv.indices_view_w(loc);
      }

      /// Creates a read-write memory view for the indices array for a given memory location
      Memory::TypedView<IT_> indices_view_rw(Memory::Location loc = Memory::Location::main)
      {
        return _sv.indices_view_rw(loc);
      }

      /**
       * \brief Creates a memory view for the indices array
       *
       * \param[in] loc
       * The memory location for which the view is to be created
       *
       * \param[in] acc
       * A combination of access rights to grant for the view
       */
      Memory::TypedView<IT_> indices_view(Memory::Location loc, Memory::Access acc)
      {
        return _sv.indices_view(loc, acc);
      }

#ifdef DOXYGEN
      // The following documentation block is visible to Doxygen only. The actual implementation is matrix type
      // specific and provided below.

      /**
       * \brief Applies the filter onto a system matrix.
       *
       * \tparam block_width_
       * The input matrix' block width
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       * The input matrix has to have a block(ed) structure and its block_height has to agree with the filter's
       * blocksize.
       *
       */
      void filter_mat(MatrixType & matrix) const
      {
      }

      /**
       * \brief Filter the non-diagonal row entries
       *
       * \tparam block_width_
       * The input matrix' block width
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       * The input matrix has to have a block(ed) structure and its block_height has to agree with the filter's
       * blocksize.
       *
       */
      void filter_offdiag_row_mat(MatrixType & matrix) const
      {
      }

      /**
       * \brief Filter the non-diagonal column entries
       *
       * \tparam block_width_
       * The input matrix' block width
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       * The input matrix has to have a block(ed) structure and its block_height has to agree with the filter's
       * blocksize.
       *
       */
      void filter_offdiag_col_mat(MatrixType & matrix) const
      {
      }

#endif
      /// \cond internal
      template<int block_width_>
      void filter_mat(SparseMatrixBCSR<DT_, IT_, block_size_, block_width_> & matrix) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == matrix.num_rows(), "Matrix size does not match!");

        Arch::UnitFilterBlockBCSR::template exec<DT_, IT_, block_size_, block_width_>(matrix.val_arbiter(), matrix.row_ptr_arbiter(),
          matrix.col_idx_arbiter(), _sv.elements_arbiter(), _sv.indices_arbiter(), _sv.num_nzes(), true, _ignore_nans);
      }

      template<int block_width_>
      void filter_offdiag_row_mat(SparseMatrixBCSR<DT_, IT_, block_size_, block_width_> & matrix) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == matrix.num_rows(), "Matrix size does not match!");

        Arch::UnitFilterBlockBCSR::template exec<DT_, IT_, block_size_, block_width_>(
          matrix.val_arbiter(), matrix.row_ptr_arbiter(), matrix.col_idx_arbiter(),
          _sv.elements_arbiter(), _sv.indices_arbiter(), _sv.num_nzes(), false, _ignore_nans);
      }

      template<int block_width_>
      void filter_offdiag_col_mat(SparseMatrixBCSR<DT_, IT_, block_size_, block_width_> &) const
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
      template<int block_width_>
      void filter_weak_matrix_rows(SparseMatrixBCSR<DT_, IT_, block_size_, block_width_> & matrix_a,
        const SparseMatrixBCSR<DT_, IT_, block_size_, block_width_>& matrix_m) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(matrix_a.val_arbiter() != matrix_m.val_arbiter(), "Matrices are not allowed to hold the same data");

        XASSERTM(_sv.size() == matrix_a.num_rows(), "Matrix size does not match!");
        XASSERTM(_sv.size() == matrix_m.num_rows(), "Matrix size does not match!");
        XASSERTM(matrix_a.num_cols() == matrix_m.num_cols(), "matrix A and M must share their layout");
        XASSERTM(matrix_a.num_nzes() == matrix_m.num_nzes(), "matrix A and M must share their layout");

        Arch::UnitFilterBlockWeakBCSR::template exec<DT_, IT_, block_size_, block_width_>(
          matrix_a.val_arbiter(), matrix_m.val_arbiter(), matrix_a.row_ptr_arbiter(),
          _sv.elements_arbiter(), _sv.indices_arbiter(), _sv.num_nzes(), _ignore_nans);
      }


      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType& vector) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == vector.size(), "Vector size does not match!");

        Arch::UnitFilterBlockVec::template exec<DT_, IT_, block_size_>(vector.elements_arbiter(),
          _sv.elements_arbiter(), _sv.indices_arbiter(), _sv.num_nzes(), false, _ignore_nans);
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
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == vector.size(), "Vector size does not match!");

        Arch::UnitFilterBlockVec::template exec<DT_, IT_, block_size_>(vector.elements_arbiter(),
          _sv.elements_arbiter(), _sv.indices_arbiter(), _sv.num_nzes(), true, _ignore_nans);
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
