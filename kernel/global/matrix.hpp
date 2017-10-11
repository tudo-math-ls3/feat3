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
     * \author Peter Zajac
     */
    template<typename LocalMatrix_, typename RowMirror_, typename ColMirror_>
    class Matrix
    {
    public:
      typedef LocalMatrix_ LocalMatrix;
      typedef RowMirror_ RowMirror;
      typedef ColMirror_ ColMirror;

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
      using ContainerType = class Matrix<LocalMatrix2_, RowMirror2_, ColMirror2_>;

      /// this typedef lets you create a matrix container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = class Matrix<typename LocalMatrix_::template ContainerType<Mem2_, DataType2_, IndexType2_>, typename RowMirror_::template MirrorType<Mem2_, DataType2_, IndexType2_>,
            typename ColMirror_::template MirrorType<Mem2_, DataType2_, IndexType2_> >;

    protected:
      GateRowType* _row_gate;
      GateColType* _col_gate;
      LocalMatrix_ _matrix;

    public:
      Matrix() :
        _row_gate(nullptr),
        _col_gate(nullptr),
        _matrix()
      {
      }

      template<typename... Args_>
      explicit Matrix(GateRowType* row_gate, GateColType* col_gate, Args_&&... args) :
        _row_gate(row_gate),
        _col_gate(col_gate),
        _matrix(std::forward<Args_>(args)...)
      {
      }

      LocalMatrix_& operator*()
      {
        return _matrix;
      }

      const LocalMatrix_& operator*() const
      {
        return _matrix;
      }

      LocalMatrix_& local()
      {
        return _matrix;
      }

      const LocalMatrix_& local() const
      {
        return _matrix;
      }

      template<typename OtherGlobalMatrix_>
      void convert(GateRowType* row_gate, GateColType* col_gate, const OtherGlobalMatrix_ & other)
      {
        this->_row_gate = row_gate;
        this->_col_gate = col_gate;
        this->_matrix.convert(*other);
      }

      const GateRowType* get_row_gate() const
      {
        return _row_gate;
      }

      const GateColType* get_col_gate() const
      {
        return _col_gate;
      }

      Matrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return Matrix(_row_gate, _col_gate, _matrix.clone(mode));
      }

      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(_row_gate, _matrix.create_vector_l());
      }

      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(_col_gate, _matrix.create_vector_r());
      }

      /**
       * \brief Gets the total number of columns in this matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of colums
       */
      Index columns() const
      {
        // Compute total number of rows and columns
        auto vec_r = create_vector_r();
        vec_r.format(DataType(1));

        return Index(vec_r.norm2sqr());
      }

      /**
       * \brief Gets the total number of rows in this matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of colums
       */
      Index rows() const
      {
        // Compute total number of rows and rows
        auto vec_l = create_vector_l();
        vec_l.format(DataType(1));

        return Index(vec_l.norm2sqr());
      }

      /**
       * \brief Returns the total number of non-zeros in this matrix.
       *
       * \note This always returns the raw (or POD - Plain Old Data) count, as everything else is ambiguous.
       *
       * \returns The total number of nonzeros in this matrix
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

      void extract_diag(VectorTypeL& diag, bool sync = true) const
      {
        _matrix.extract_diag(*diag);
        if(sync)
        {
          diag.sync_0();
        }
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        _matrix.apply(*r, *x);
        r.sync_0();
      }

      auto apply_async(VectorTypeL& r, const VectorTypeR& x) const -> decltype(r.sync_0_async())
      {
        _matrix.apply(*r, *x);
        return r.sync_0_async();
      }

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const
      {
        // copy y to r
        r.copy(y);

        // convert from type-1 to type-0
        r.from_1_to_0();

        // r <- r + alpha*A*x
        _matrix.apply(*r, *x, *r, alpha);

        // synchronise r
        r.sync_0();
      }

      auto apply_async(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const -> decltype(r.sync_0_async())
      {
        // copy y to r
        r.copy(y);

        // convert from type-1 to type-0
        r.from_1_to_0();

        // r <- r + alpha*A*x
        _matrix.apply(*r, *x, *r, alpha);

        // synchronise r
        r.sync_0_async();
      }

      /**
       * \brief Computes the lumped rows vector
       *
       * \param[out] lump
       * The vector receiving the lumped rows
       *
       * \param[in] sync
       * Synchronise the lumped rows vector?
       *
       * Each entry in the returned lumped rows vector contains the the sum of all matrix elements in the
       * corresponding row.
       *
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
       * \brief Returns the lumped rows vector
       *
       * Each entry in the returned lumped rows vector contains the the sum of all matrix elements in the
       * corresponding row.
       *
       * \returns
       * The lumped vector.
       */
      VectorTypeL lump_rows(bool sync = true) const
      {
        VectorTypeL lump = create_vector_l();

        lump_rows(lump, sync);

        return lump;
      }

      LocalMatrix convert_to_1() const
      {
        LocalMatrix locmat = _matrix.clone(LAFEM::CloneMode::Weak);
        if((_row_gate != nullptr) && (_col_gate != nullptr))
          synch_matrix(locmat, *_row_gate->_comm, _row_gate->_ranks, _row_gate->_mirrors, _col_gate->_mirrors);
        return locmat;
      }
    };
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_MATRIX_HPP
