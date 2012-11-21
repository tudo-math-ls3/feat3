#pragma once
#ifndef SCARC_GUARD_SOLVER_FUNCTOR_HH
#define SCARC_GUARD_SOLVER_FUNCTOR_HH 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/foundation/functor.hpp>
#include<kernel/lafem/defect.hpp>

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
  }
}

#endif
