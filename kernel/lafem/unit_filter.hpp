// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_vector.hpp>
#include <kernel/lafem/arch/unit_filter_dense_vec.hpp>
#include <kernel/lafem/arch/unit_filter_dense_csr.hpp>
#include <kernel/lafem/arch/unit_filter_dense_weak_csr.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Unit Filter class template.
     *
     * \author Peter Zajac
     */
    template<
      typename DT_,
      typename IT_ = Index>
    class UnitFilter
    {
    public:
      /// data-type typedef
      typedef DT_ DataType;
      /// index-type typedef
      typedef IT_ IndexType;

      /// our supported vector type
      typedef DenseVector<DataType, IndexType> VectorType;

      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using FilterType = UnitFilter<DT2_, IT2_>;

      /// this typedef lets you create a filter with new Datatype and Index types
      template <typename DataType2_, typename IndexType2_>
      using FilterTypeByDI = FilterType<DataType2_, IndexType2_>;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

    private:
      /// SparseVector, containing all entries of the unit filter
      SparseVector<DT_, IT_> _sv;

    public:
      /// default constructor
      UnitFilter() :
        _sv()
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] size_in
       * The total size of the unit filter.
       */
      explicit UnitFilter(Index size_in, Index num_nzes_in) :
        _sv(size_in, num_nzes_in)
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] size_in The size of the created filter.
       * \param[in] values DenseVector containing element values
       * \param[in] indices DenseVector containing element indices
       */
      explicit UnitFilter(Index size_in, DenseVector<DT_, IT_> & values, DenseVector<IT_, IT_> & indices) :
        _sv(size_in, values, indices)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] sv
       * The sparse vector containing the filter entries
       */
      explicit UnitFilter(SparseVector<DT_, IT_> && sv) :
        _sv(std::forward<SparseVector<DT_, IT_>>(sv))
      {
      }

      /// move-ctor
      UnitFilter(UnitFilter && other) :
        _sv(std::forward<SparseVector<DT_, IT_>>(other._sv))
      {
      }

      /// move-assignment operator
      UnitFilter & operator=(UnitFilter && other)
      {
        if(this != &other)
        {
          _sv = std::forward<SparseVector<DT_, IT_>>(other._sv);
        }
        return *this;
      }

      /// virtual destructor
      virtual ~UnitFilter()
      {
      }

      /// \brief Creates a clone of itself
      UnitFilter clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        UnitFilter other;
        other.clone(*this, clone_mode);
        return other;
      }

      /// \brief Clones data from another UnitFilter
      void clone(const UnitFilter & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _sv.clone(other.get_filter_vector(), clone_mode);
      }

      /// \brief Converts data from another UnitFilter
      template<typename DT2_, typename IT2_>
      void convert(const UnitFilter<DT2_, IT2_>& other)
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
      SparseVector<DT_, IT_>& get_filter_vector()
      {
        return _sv;
      }

      const SparseVector<DT_, IT_>& get_filter_vector() const
      {
        return _sv;
      }
      /// \endcond

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
      Memory::TypedView<DT_> elements_view_r(Memory::Location loc = Memory::Location::main) const
      {
        return _sv.elements_view_r(loc);
      }

      /// Creates a write-only memory view for the elements array for a given memory location
      Memory::TypedView<DT_> elements_view_w(Memory::Location loc = Memory::Location::main)
      {
        return _sv.elements_view_w(loc);
      }

      /// Creates a read-write memory view for the elements array for a given memory location
      Memory::TypedView<DT_> elements_view_rw(Memory::Location loc = Memory::Location::main)
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
      Memory::TypedView<DT_> elements_view(Memory::Location loc, Memory::Access acc)
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

      /// Permute internal vector according to the given Permutation
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
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       */
      void filter_mat(MatrixType & matrix) const
      {
      }

      /**
       * \brief Filter the non-diagonal entries, row wise
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       */
      void filter_offdiag_row_mat(MatrixType & matrix) const
      {
      }

      /**
       * \brief Filter the non-diagonal entries, column wise
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       *
       */
      void filter_offdiag_col_mat(MatrixType & matrix) const
      {
      }

#endif
      ///\cond internal
      void filter_mat(SparseMatrixCSR<DT_, IT_> & matrix) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == matrix.num_rows(), "Matrix size does not match!");

        Arch::UnitFilterDenseCSR::template exec<DT_, IT_>(matrix.val_arbiter(), matrix.row_ptr_arbiter(), matrix.col_idx_arbiter(),
          _sv.indices_arbiter(), _sv.num_nzes(), true, 1);
      }

      void filter_offdiag_row_mat(SparseMatrixCSR<DT_, IT_> & matrix) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == matrix.num_rows(), "Matrix size does not match!");

        Arch::UnitFilterDenseCSR::template exec<DT_, IT_>(matrix.val_arbiter(), matrix.row_ptr_arbiter(), matrix.col_idx_arbiter(),
          _sv.indices_arbiter(), _sv.num_nzes(), false, 1);
      }

      void filter_offdiag_col_mat(SparseMatrixCSR<DT_, IT_> &) const
      {
        // nothing to do here
      }

      template<int block_width_>
      void filter_offdiag_row_mat(SparseMatrixBCSR<DT_, IT_, 1, block_width_> & matrix) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == matrix.num_rows(), "Matrix size does not match!");

        Arch::UnitFilterDenseCSR::template exec<DT_, IT_>(matrix.val_arbiter(), matrix.row_ptr_arbiter(), matrix.col_idx_arbiter(),
          _sv.indices_arbiter(), _sv.num_nzes(), false, block_width_);
      }

      template<int block_height_>
      void filter_offdiag_row_mat(SparseMatrixBCSR<DT_, IT_, block_height_, 1> &) const
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
      void filter_weak_matrix_rows(SparseMatrixCSR<DT_, IT_> & matrix_a, const SparseMatrixCSR<DT_, IT_>& matrix_m) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == matrix_a.num_rows(), "Matrix size does not match!");
        XASSERTM(_sv.size() == matrix_m.num_rows(), "Matrix size does not match!");
        XASSERTM(matrix_a.num_cols() == matrix_m.num_cols(), "matrix A and M must share their layout");
        XASSERTM(matrix_a.num_nzes() == matrix_m.num_nzes(), "matrix A and M must share their layout");

        Arch::UnitFilterDenseWeakCSR::template exec<DT_, IT_>(
          matrix_a.val_arbiter(), matrix_m.val_arbiter(), matrix_a.row_ptr_arbiter(),
          _sv.elements_arbiter(), _sv.indices_arbiter(), _sv.num_nzes());
      }


      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(DenseVector<DT_, IT_> & vector) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == vector.size(), "Vector size does not match!");

        Arch::UnitFilterDenseVec::template exec<DT_, IT_>(vector.elements_arbiter(),
          _sv.elements_arbiter(), _sv.indices_arbiter(), _sv.num_nzes(), false);
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      void filter_sol(DenseVector<DT_, IT_> & vector) const
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
      void filter_def(DenseVector<DT_, IT_> & vector) const
      {
        if(_sv.hollow())
          return;

        XASSERTM(_sv.size() == vector.size(), "Vector size does not match!");

        Arch::UnitFilterDenseVec::template exec<DT_, IT_>(vector.elements_arbiter(),
          _sv.elements_arbiter(), _sv.indices_arbiter(), _sv.num_nzes(), true);
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      void filter_cor(DenseVector<DT_, IT_> & vector) const
      {
        // same as def
        filter_def(vector);
      }
    }; // class UnitFilter<...>
  } // namespace LAFEM
} // namespace FEAT
