#pragma once
#ifndef KERNEL_GLOBAL_MATRIX_HPP
#define KERNEL_GLOBAL_MATRIX_HPP 1

#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>

namespace FEAST
{
  namespace Global
  {
    /**
     * \brief Global Matrix wrapper class template
     *
     * \author Peter Zajac
     */
    template<typename LocalMatrix_>
    class Matrix
    {
    public:
      typedef typename LocalMatrix_::MemType MemType;
      typedef typename LocalMatrix_::DataType DataType;
      typedef typename LocalMatrix_::IndexType IndexType;

      typedef typename LocalMatrix_::VectorTypeL LocalVectorTypeL;
      typedef typename LocalMatrix_::VectorTypeR LocalVectorTypeR;

      typedef Vector<LocalVectorTypeL> VectorTypeL;
      typedef Vector<LocalVectorTypeR> VectorTypeR;

      typedef Gate<LocalVectorTypeL> GateRowType;
      typedef Gate<LocalVectorTypeR> GateColType;;

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

      template<typename OtherLocalMatrix_>
      void convert(GateRowType* row_gate, GateColType* col_gate, const Global::Matrix<OtherLocalMatrix_>& other)
      {
        this->_row_gate = row_gate;
        this->_col_gate = col_gate;
        this->_matrix.convert(*other);
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
    };
  } // namespace Global
} // namespace FEAST

#endif // KERNEL_GLOBAL_MATRIX_HPP
