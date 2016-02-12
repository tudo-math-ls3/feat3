#pragma once
#ifndef KERNEL_GLOBAL_NONLINEAR_FUNCTIONAL_HPP
#define KERNEL_GLOBAL_NONLINEAR_FUNCTIONAL_HPP 1

#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>

namespace FEAST
{
  namespace Global
  {
    /**
     * \brief Global NonlinearFunctional wrapper class template
     *
     * \author Jordi Paul
     */
    template<typename LocalNonlinearFunctional_>
    class NonlinearFunctional
    {
    public:
      typedef typename LocalNonlinearFunctional_::MemType MemType;
      typedef typename LocalNonlinearFunctional_::DataType DataType;
      typedef typename LocalNonlinearFunctional_::IndexType IndexType;

      typedef typename LocalNonlinearFunctional_::VectorTypeL LocalVectorTypeL;
      typedef typename LocalNonlinearFunctional_::VectorTypeR LocalVectorTypeR;

      typedef Vector<LocalVectorTypeL> VectorTypeL;
      typedef Vector<LocalVectorTypeR> VectorTypeR;

      typedef Gate<LocalVectorTypeL> GateRowType;
      typedef Gate<LocalVectorTypeR> GateColType;

      typedef VectorTypeL GradientType;

    private:
      Index _columns;
      Index _rows;


    protected:
      GateRowType* _row_gate;
      GateColType* _col_gate;
      LocalNonlinearFunctional_ _nonlinear_functional;

    public:
      NonlinearFunctional() :
        _row_gate(nullptr),
        _col_gate(nullptr),
        _nonlinear_functional()
      {
      }

      template<typename... Args_>
      explicit NonlinearFunctional(GateRowType* row_gate, GateColType* col_gate, Args_&&... args) :
        _row_gate(row_gate),
        _col_gate(col_gate),
        _nonlinear_functional(std::forward<Args_>(args)...)
      {
      }

      LocalNonlinearFunctional_& operator*()
      {
        return _nonlinear_functional;
      }

      const LocalNonlinearFunctional_& operator*() const
      {
        return _nonlinear_functional;
      }

      //template<typename OtherLocalNonlinearFunctional_>
      //void convert(GateRowType* row_gate, GateColType* col_gate, const Global::NonlinearFunctional<OtherLocalNonlinearFunctional_>& other)
      //{
      //  this->_row_gate = row_gate;
      //  this->_col_gate = col_gate;
      //  this->_nonlinear_functional.convert(*other);
      //}

      //NonlinearFunctional clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      //{
      //  return NonlinearFunctional(_row_gate, _col_gate, _nonlinear_functional.clone(mode));
      //}

      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(_row_gate, _nonlinear_functional.create_vector_l());
      }

      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(_col_gate, _nonlinear_functional.create_vector_r());
      }

      //void extract_diag(VectorTypeL& diag, bool sync = true) const
      //{
      //  _nonlinear_functional.extract_diag(*diag);
      //  if(sync)
      //  {
      //    diag.sync_0();
      //  }
      //}

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        _nonlinear_functional.apply(*r, *x);
        r.sync_0();
      }

        /**
         * \brief Prepares the operator for evaluation by setting the current state
         *
         * \param[in] vec_state
         * The current state
         *
         */
        void prepare(VectorTypeR& vec_state)
        {
          _nonlinear_functional.prepare(*vec_state);
        }

        const Index& columns()
        {
          _columns = _nonlinear_functional.columns();
          Foundation::Comm::allreduce(&_columns, 1, &_columns);
          return _columns;
        }

        DataType compute_func()
        {
          DataType my_fval(_nonlinear_functional.compute_func());
          Foundation::Comm::allreduce(&my_fval, 1, &my_fval);
          return my_fval;
        }

        void compute_grad(VectorTypeL& grad)
        {
          _nonlinear_functional.compute_grad(*grad);
          grad.sync_0();
        }



      //void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const
      //{
      //  // copy y to r
      //  r.copy(y);

      //  // convert from type-1 to type-0
      //  r.from_1_to_0();

      //  // r <- r + alpha*A*x
      //  _nonlinear_functional.apply(*r, *x, *r, alpha);

      //  // synchronise r
      //  r.sync_0();
      //}
    };
  } // namespace Global
} // namespace FEAST

#endif // KERNEL_GLOBAL_NONLINEAR_FUNCTIONAL
