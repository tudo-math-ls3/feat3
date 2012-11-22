#pragma once
#ifndef SCARC_GUARD_SOLVER_FUNCTOR_HH
#define SCARC_GUARD_SOLVER_FUNCTOR_HH 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/foundation/functor.hpp>
#include<kernel/lafem/defect.hpp>
#include<kernel/lafem/sum.hpp>
#include<kernel/lafem/product.hpp>

using namespace FEAST::Foundation;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {
    ///ApplicableAsPrecon interface
    template<typename MyFunctorType_, typename VT_>
    class ApplicableAsPrecon
    {
      public:
        ///needed in substitution of CSF
        typedef VT_ vector_type_;

        ///named copy-CTOR
        virtual MyFunctorType_ substitute(VT_& arg) = 0;
    };


    template<typename VT_>
    class ProxyPreconApplyFunctor : public FunctorBase
    {
      public:
        ProxyPreconApplyFunctor(VT_& x) :
          _x(x)
        {
        }

        virtual const std::string type_name()
        {
          return "ProxyPreconApplyFunctor";
        }

        virtual void execute()
        {
          throw ScaRCError("Error: Proxy functors can not be executed - substitute first!");
        }

        virtual void undo()
        {
          throw ScaRCError("Error: Numerical functors can not be undone!");
        }

        ProxyPreconApplyFunctor& operator=(const ProxyPreconApplyFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_x = rhs._x;
          return *this;
        }

        ProxyPreconApplyFunctor(const ProxyPreconApplyFunctor& other) :
          _x(other._x)
        {
        }

        VT_& get_argument()
        {
          return _x;
        }

      private:
        ///apply preconditioner to what?
        VT_& _x;
    };

    template<template<typename, typename> class StorageType_ = std::vector>
    class CompoundSolverFunctor : public CompoundFunctor<StorageType_>
    {
      public:
        ///overwrite execute and undo functions
        virtual void execute()
        {
          for(Index i(0) ; i < (this->_functors).size() ; ++i)
          {
            (this->_functors).at(i)->execute();
          }
        }

        virtual const std::string type_name()
        {
          return "CompoundSolverFunctor";
        }

        virtual void undo()
        {
          throw ScaRCError("Error: Numerical functor lists can not be undone!");
        }

        ///substitute first ProxyPreconApplyFunctor by new_precon
        template<typename T_>
        void substitute_first(T_& new_precon)
        {
          if(new_precon.type_name() == "CompoundSolverFunctor")
          {
            for(Index i(0) ; i < (this->_functors).size() ; ++i)
            {
              if( (this->_functors).at(i)->type_name() == "ProxyPreconApplyFunctor")
              {
                ///TODO
              }
            }
          }
          else
          {
            for(Index i(0) ; i < (this->_functors).size() ; ++i)
            {
              if( (this->_functors).at(i)->type_name() == "ProxyPreconApplyFunctor")
              {
                (this->_functors).at(i) = std::shared_ptr<FunctorBase>(new T_(new_precon.substitute(((ProxyPreconApplyFunctor<typename T_::vector_type_>*)(this->_functors.at(i).get()))->get_argument())));
                return;
              }
            }
          }
          throw ScaRCError("Error: No ProxyPreconApplyFunctor in functor list!");
        }
    };

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


    template<typename Algo_, typename VT_, typename MT_>
    class ProductFunctor : public FunctorBase, public ApplicableAsPrecon<ProductFunctor<Algo_, VT_, MT_>, VT_>
    {
      public:
        ProductFunctor(VT_& y, const MT_& A, const VT_& x) :
          _y(y),
          _A(A),
          _x(x)
        {
        }

        virtual const std::string type_name()
        {
          return "ProductFunctor";
        }

        virtual void execute()
        {
          LAFEM::Product<Algo_>::value(_y, _A, _x);
        }

        virtual void undo()
        {
          throw ScaRCError("Error: Numerical functors can not be undone!");
        }

        ProductFunctor& operator=(const ProductFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_A = rhs._A;
          this->_x = rhs._x;
          return *this;
        }

        ProductFunctor(const ProductFunctor& other) :
          _y(other._y),
          _A(other._A),
          _x(other._x)
        {
        }

        ///implementation of ApplicableAsPrecon interface
        virtual ProductFunctor<Algo_, VT_, MT_> substitute(VT_& arg)
        {
          return ProductFunctor<Algo_, VT_, MT_>(_y, _A, arg);
        }

      private:
        VT_& _y;
        const MT_& _A;
        const VT_& _x;
    };
  }
}

#endif
