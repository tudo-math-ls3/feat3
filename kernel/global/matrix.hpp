// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GLOBAL_MATRIX_HPP
#define KERNEL_GLOBAL_MATRIX_HPP 1

#include <kernel/lafem/container.hpp> // required for LAFEM::CloneMode
#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/synch_mat.hpp>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Global Matrix wrapper class template
     *
     * This class implements a wrapper that contains a LAFEM matrix as its core data object and
     * provides the necessary synchronization functions required in an MPI-parallel simulation
     * based on the overlapping domain decomposition approach. Effectively, this class only
     * couples a local LAFEM matrix with its corresponding pair of row/column Global::Gate objects,
     * which are required to define the compatible L/R vector types, which then take care of the
     * actual synchronization dirty work.
     *
     * \tparam LocalVector_
     * The type of the local matrix container; may be any valid combination of LAFEM (meta-)matrix types
     *
     * \tparam RowMirror_
     * The type of the row mirror; must be compatible to the L-vector type of the local matrix
     *
     * \tparam ColMirror_
     * The type of the column mirror; must be compatible to the R-vector type of the local matrix
     *
     * \author Peter Zajac
     */
    template<typename LocalMatrix_, typename RowMirror_, typename ColMirror_>
    class Matrix
    {
    public:
      typedef LocalMatrix_ LocalMatrixType;
      typedef RowMirror_ RowMirrorType;
      typedef ColMirror_ ColMirrorType;

      typedef typename LocalMatrix_::MemType MemType;
      typedef typename LocalMatrix_::DataType DataType;
      typedef typename LocalMatrix_::IndexType IndexType;

      typedef typename LocalMatrix_::VectorTypeL LocalVectorTypeL;
      typedef typename LocalMatrix_::VectorTypeR LocalVectorTypeR;

      typedef Vector<LocalVectorTypeL, RowMirror_> VectorTypeL;
      typedef Vector<LocalVectorTypeR, ColMirror_> VectorTypeR;

      typedef Gate<LocalVectorTypeL, RowMirror_> GateRowType;
      typedef Gate<LocalVectorTypeR, ColMirror_> GateColType;

      /// Our 'base' class type
      template <typename LocalMatrix2_, typename RowMirror2_ = RowMirror_, typename ColMirror2_ = ColMirror_>
      using ContainerType = Matrix<LocalMatrix2_, RowMirror2_, ColMirror2_>;

      /// this typedef lets you create a matrix container with new Memory, Datatype and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = Matrix<
        typename LocalMatrix_::template ContainerType<Mem2_, DataType2_, IndexType2_>,
        typename RowMirror_::template MirrorType<Mem2_, DataType2_, IndexType2_>,
        typename ColMirror_::template MirrorType<Mem2_, DataType2_, IndexType2_> >;

      /// this is a global matrix class
      static constexpr bool is_global = true;
      /// this is not a local matrix class
      static constexpr bool is_local = false;

    protected:
      /// a pointer to the row gate responsible for synchronization
      GateRowType* _row_gate;
      /// a pointer to the column gate responsible for synchronization
      GateColType* _col_gate;
      /// the internal local matrix object
      LocalMatrix_ _matrix;

    public:
      /// standard constructor
      Matrix() :
        _row_gate(nullptr),
        _col_gate(nullptr),
        _matrix()
      {
      }

      /**
       * \brief Forwarding constructor
       *
       * \param[in] row_gate
       * A \resident pointer to the row gate to be used for synchronization
       *
       * \param[in] col_gate
       * A \resident pointer to the column gate to be used for synchronization
       *
       * \param[in] args
       * The arguments that are to be passed to the local matrix object constructor
       *
       * \attention
       * The two gates may be different objects, which is e.g. required for matrices which are
       * defined by different test- and trial-spaces in the finite element context, but they
       * always must use the same internal communicator!
       */
      template<typename... Args_>
      explicit Matrix(GateRowType* row_gate, GateColType* col_gate, Args_&&... args) :
        _row_gate(row_gate),
        _col_gate(col_gate),
        _matrix(std::forward<Args_>(args)...)
      {
        if((_row_gate != nullptr) && (_col_gate != nullptr))
        {
          // the gates must be defined on the same communicator
          XASSERT(_row_gate->get_comm() == _col_gate->get_comm());
        }
      }

      /**
       * \brief Returns a reference to the internal local LAFEM matrix object.
       *
       * \returns A reference to the internal local LAFEM matrix object.
       */
      LocalMatrix_& local()
      {
        return _matrix;
      }

      /**
       * \brief Returns a const reference to the internal local LAFEM matrix object.
       *
       * \returns A const reference to the internal local LAFEM matrix object.
       */
      const LocalMatrix_& local() const
      {
        return _matrix;
      }

      template<typename OtherGlobalMatrix_>
      void convert(GateRowType* row_gate, GateColType* col_gate, const OtherGlobalMatrix_ & other)
      {
        this->_row_gate = row_gate;
        this->_col_gate = col_gate;
        this->_matrix.convert(other.local());
      }

      /**
       * \brief Returns a const pointer to the internal row gate of the matrix.
       *
       * \returns A const pointer to the internal row gate of the matrix.
       */
      const GateRowType* get_row_gate() const
      {
        return _row_gate;
      }

      /**
       * \brief Returns a const pointer to the internal column gate of the matrix.
       *
       * \returns A const pointer to the internal column gate of the matrix.
       */
      const GateColType* get_col_gate() const
      {
        return _col_gate;
      }

      /**
       * \brief Returns a const pointer to the internal communicator of the gates of the matrix.
       *
       * \returns a const pointer to the internal communicator of the gates of the matrix.
       */
      const Dist::Comm* get_comm() const
      {
        return _row_gate != nullptr ? _row_gate->get_comm() : nullptr;
      }

      /**
       * \brief Creates and returns a clone of this global matrix
       *
       * \param[in] mode
       * Specifies the clone mode for the internal matrix object.
       *
       * \returns The created clone object
       */
      Matrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return Matrix(_row_gate, _col_gate, _matrix.clone(mode));
      }

      /**
       * \brief Creates and returns a new L-compatible global vector object
       *
       * \returns A new L-compatible global vector that is linked to the row gate of this matrix
       */
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(_row_gate, _matrix.create_vector_l());
      }

      /**
       * \brief Creates and returns a new R-compatible global vector object
       *
       * \returns A new R-compatible global vector that is linked to the column gate of this matrix
       */
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(_col_gate, _matrix.create_vector_r());
      }

      /**
       * \brief Gets the total number of columns in this distributed matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       *
       * \returns The number of columns of this matrix
       */
      Index columns() const
      {
        // Compute total number of columns
        auto vec_r = create_vector_r();
        vec_r.format(DataType(1));

        return Index(vec_r.norm2sqr());
      }

      /**
       * \brief Returns the total number of rows in this distributed matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       *
       * \returns The number of rows of this matrix
       */
      Index rows() const
      {
        // Compute total number of rows
        auto vec_l = create_vector_l();
        vec_l.format(DataType(1));

        return Index(vec_l.norm2sqr());
      }

      /**
       * \brief Returns the total number of non-zeros in this distributed matrix
       *
       * \note This always returns the raw (or POD - Plain Old Data) count, as everything else is ambiguous.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       *
       * \returns The total number of non-zeros in this matrix
       */
      Index used_elements() const
      {
        Index my_used_elements(_matrix.template used_elements<LAFEM::Perspective::pod>());
        const Dist::Comm* comm = (_row_gate != nullptr ? _row_gate->get_comm() : (_col_gate != nullptr ? _col_gate->get_comm() : nullptr));
        if((comm != nullptr) && (comm->size() > 1))
          comm->allreduce(&my_used_elements, &my_used_elements, std::size_t(1), Dist::op_sum);
        return my_used_elements;
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        size_t my_bytes(_matrix.bytes());
        const Dist::Comm* comm = (_row_gate != nullptr ? _row_gate->get_comm() : (_col_gate != nullptr ? _col_gate->get_comm() : nullptr));
        if((comm != nullptr) && (comm->size() > 1))
          comm->allreduce(&my_bytes, &my_bytes, std::size_t(1), Dist::op_sum);
        return my_bytes;
      }

      /**
       * \brief Extracts the main diagonal of the matrix as a vector
       *
       * \param[in] diag
       * A \transient reference to the vector that shall receive the main diagonal elements
       *
       * \param[in] sync
       * Specifies whether the main diagonal vector is to be synchronized to obtain a type-1 vector.
       * If set to \c false, the resulting vector will be a type-0 vector.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application might deadlock.
       */
      void extract_diag(VectorTypeL& diag, bool sync = true) const
      {
        _matrix.extract_diag(diag.local());
        if(sync)
        {
          diag.sync_0();
        }
      }

      /**
       * \brief Performs a matrix-vector multiplication: r <- A*x
       *
       * \param[inout] r
       * A \transient reference to the vector that should receive the result of the matrix-vector product.
       * Must be allocated to the correct sizes and must have a valid gate assigned, but its numerical
       * contents upon entry are ignored. Must not be the same object as \p x.
       *
       * \param[in] x
       * A \transient reference to the vector that is to be multiplied by this matrix.
       * Must not be the same object as \p r.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        _matrix.apply(r.local(), x.local());
        r.sync_0();
      }

      /**
       * \brief Performs a matrix-vector multiplication: r <- A*x
       *
       * \param[inout] r
       * A \transient reference to the vector that should receive the result of the matrix-vector product.
       * Must be allocated to the correct sizes and must have a valid gate assigned, but its numerical
       * contents upon entry are ignored. Must not be the same object as \p x.
       *
       * \param[in] x
       * A \transient reference to the vector that is to be multiplied by this matrix.
       * Must not be the same object as \p r.
       *
       * \returns A SynchVectorTicket object that waits for the operation to complete.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      auto apply_async(VectorTypeL& r, const VectorTypeR& x) const -> decltype(r.sync_0_async())
      {
        _matrix.apply(r.local(), x.local());
        return r.sync_0_async();
      }

      /**
       * \brief Performs a matrix-vector multiplication: r <- y + alpha*A*x
       *
       * \param[inout] r
       * A \transient reference to the vector that should receive the result of the matrix-vector product.
       * Must be allocated to the correct sizes and must have a valid gate assigned, but its numerical
       * contents upon entry are ignored. Must not be the same object as \p x, but it may be
       * the same object as \p y.
       *
       * \param[in] x
       * A \transient reference to the vector that is to be multiplied by this matrix.
       * Must not be the same object as \p r.
       *
       * \param[in] y
       * A \transient reference to the vector that is to be added onto the result of the product.
       * May be the same object as \p r.
       *
       * \param[in] alpha
       * A scaling factor for the matrix-vector product A*x.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const
      {
        // copy y to r
        r.copy(y);

        // convert from type-1 to type-0
        r.from_1_to_0();

        // r <- r + alpha*A*x
        _matrix.apply(r.local(), x.local(), r.local(), alpha);

        // synchronize r
        r.sync_0();
      }

      /**
       * \brief Performs a matrix-vector multiplication: r <- y + alpha*A*x
       *
       * \param[inout] r
       * A \transient reference to the vector that should receive the result of the matrix-vector product.
       * Must be allocated to the correct sizes and must have a valid gate assigned, but its numerical
       * contents upon entry are ignored. Must not be the same object as \p x, but it may be
       * the same object as \p y.
       *
       * \param[in] x
       * A \transient reference to the vector that is to be multiplied by this matrix.
       * Must not be the same object as \p r.
       *
       * \param[in] y
       * A \transient reference to the vector that is to be added onto the result of the product.
       * May be the same object as \p r.
       *
       * \param[in] alpha
       * A scaling factor for the matrix-vector product A*x.
       *
       * \returns A SynchVectorTicket object that waits for the operation to complete.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      auto apply_async(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const -> decltype(r.sync_0_async())
      {
        // copy y to r
        r.copy(y);

        // convert from type-1 to type-0
        r.from_1_to_0();

        // r <- r + alpha*A*x
        _matrix.apply(r.local(), x.local(), r.local(), alpha);

        // synchronize r
        return r.sync_0_async();
      }

      /**
       * \brief Computes the lumped rows of the matrix as a vector
       *
       * Each entry of the lumped rows vector contains the sum of all matrix elements in the
       * corresponding row of this matrix.
       *
       * \param[in] lump
       * A \transient reference to the vector that shall receive the lumped row elements
       *
       * \param[in] sync
       * Specifies whether the lumped rows vector is to be synchronized to obtain a type-1 vector.
       * If set to \c false, the resulting vector will be a type-0 vector.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application might deadlock.
       */
      void lump_rows(VectorTypeL& lump, bool sync = true) const
      {
        XASSERTM(lump.local().size() == _matrix.rows(), "lump vector size does not match matrix row count!");

        _matrix.lump_rows(lump.local());

        if(sync)
        {
          lump.sync_0();
        }
      }

      /**
       * \brief Computes and returns the lumped rows of the matrix as a vector
       *
       * Each entry of the lumped rows vector contains the sum of all matrix elements in the
       * corresponding row of this matrix.
       *
       * \param[in] sync
       * Specifies whether the lumped rows vector is to be synchronized to obtain a type-1 vector.
       * If set to \c false, the resulting vector will be a type-0 vector.
       *
       * \returns
       * A new vector containing the lumped row elements of this matrix.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application might deadlock.
       */
      VectorTypeL lump_rows(bool sync = true) const
      {
        VectorTypeL lump = create_vector_l();

        lump_rows(lump, sync);

        return lump;
      }

      /**
       * \brief Computes and returns the type-1 conversion of this matrix as a local matrix
       *
       * This function performs a type-0 to type-1 conversion of this matrix, which is effectively
       * the matrix counterpart of the type-0 synchronization of a vector, i.e. all matrix entries,
       * which are shared by multiple processes, are replaced by the sum of all contributions from
       * all processes which share the corresponding row/column DOFs. This converted matrix is
       * required for Schwarz-like solver approaches.
       *
       * \attention
       * The resulting matrix may (and usually will) have a larger stencil than the original underlying
       * patch-local LAFEM matrix, which is stored in this object's _matrix variable. The reason is that
       * a neighbor process may have additional couplings between DOFs that this process does not have,
       * and these additional coupling are also included in this process's type-1 matrix. However, the
       * stencil of the returned matrix is never smaller than the original local matrix stencil, i.e.
       * the stencil of the patch-local type-0 matrix is always a sub-stencil of the returned type-1
       * matrix.
       *
       * \returns
       * A new local matrix object containing the type-1 matrix for this process patch.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application might deadlock.
       */
      LocalMatrix_ convert_to_1() const
      {
        LocalMatrix_ locmat = _matrix.clone(LAFEM::CloneMode::Weak);
        if((_row_gate != nullptr) && (_col_gate != nullptr))
          synch_matrix(locmat, *_row_gate->_comm, _row_gate->_ranks, _row_gate->_mirrors, _col_gate->_mirrors);
        return locmat;
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(LAFEM::SerialConfig& config)
      {
        return _matrix.get_checkpoint_size(config);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char>& data)
      {
        _matrix.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, LAFEM::SerialConfig& config)
      {
        return _matrix.set_checkpoint_data(data, config);
      }
    }; // class Matrix<...>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_MATRIX_HPP
