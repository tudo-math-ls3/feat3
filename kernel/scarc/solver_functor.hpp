#pragma once
#ifndef SCARC_GUARD_SOLVER_FUNCTOR_HH
#define SCARC_GUARD_SOLVER_FUNCTOR_HH 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/foundation/functor.hpp>
#include<kernel/lafem/defect.hpp>
#include<kernel/lafem/sum.hpp>

using namespace FEAST::Foundation;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {
    template<typename Algo_, typename VT_, typename MT_>
    class DefectFunctor : public FunctorBase
    {
      public:
        DefectFunctor(VT_& d, const VT_& b, const MT_& A, const VT_& x) :
          _d(d),
          _b(b),
          _A(A),
          _x(x)
        {
        }

        virtual const std::string type_name()
        {
          return "DefectFunctor";
        }

        virtual void execute()
        {
          LAFEM::Defect<Algo_>::value(_d, _b, _A, _x);
        }

        virtual void undo()
        {
          throw ScaRCError("Error: Numerical functors can not be undone!");
        }

        DefectFunctor& operator=(const DefectFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_d = rhs._d;
          this->_b = rhs._b;
          this->_A = rhs._A;
          this->_x = rhs._x;
          return *this;
        }

        DefectFunctor(const DefectFunctor& other) :
          _d(other._d),
          _b(other._b),
          _A(other._A),
          _x(other._x)
        {
        }

      private:
        VT_& _d;
        const VT_& _b;
        const MT_& _A;
        const VT_& _x;
    };


    template<typename Algo_, typename VT_>
    class SumFunctor : public FunctorBase
    {
      public:
        SumFunctor(VT_& y, const VT_& r, const VT_& l) :
          _y(y),
          _r(r),
          _l(l)
        {
        }

        virtual const std::string type_name()
        {
          return "SumFunctor";
        }

        virtual void execute()
        {
          LAFEM::Sum<Algo_>::value(_y, _r, _l);
        }

        virtual void undo()
        {
          throw ScaRCError("Error: Numerical functors can not be undone!");
        }

        SumFunctor& operator=(const SumFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_r = rhs._r;
          this->_l = rhs._l;
          return *this;
        }

        SumFunctor(const SumFunctor& other) :
          _y(other._y),
          _r(other._r),
          _l(other._l)
        {
        }

      private:
        VT_& _y;
        const VT_& _r;
        const VT_& _l;
    };
  }
}

#endif
