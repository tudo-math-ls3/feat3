#pragma once
#ifndef SCARC_GUARD_SOLVER_FUNCTOR_HH
#define SCARC_GUARD_SOLVER_FUNCTOR_HH 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/foundation/functor.hpp>
#include<kernel/lafem/defect.hpp>
#include<kernel/lafem/sum.hpp>
#include<kernel/lafem/product_matvec.hpp>
#include<kernel/lafem/norm.hpp>

using namespace FEAST::Foundation;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {
    template<typename VT_>
    class SolverFunctorBase
    {
      public:
        ///needed in substitution of CSF
        typedef VT_ solver_vector_type_;

        virtual void substitute(VT_& arg) = 0;

        virtual void execute() = 0;
        virtual const std::string type_name() = 0;

        virtual ~SolverFunctorBase()
        {
        }

        virtual bool is_preconditioner()
        {
          return false;
        }

        virtual bool is_complete()
        {
          return _complete;
        }

      protected:
        bool _complete;
    };


    ///always substitute left argument
    template<typename Algo_, typename VT_>
    class SumFunctorProxyLeft : public SolverFunctorBase<VT_>
    {
      public:
        SumFunctorProxyLeft(VT_& y, VT_& l, const VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "SumFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete SumFunctor can not be executed!");

          LAFEM::Sum<Algo_>::value(_y, _l, _r);
        }

        SumFunctorProxyLeft& operator=(const SumFunctorProxyLeft& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        SumFunctorProxyLeft(const SumFunctorProxyLeft& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _l = arg;
          this->_complete = true;
        }

      private:
        VT_& _y;
        VT_& _l;
        const VT_& _r;
    };


    ///always substitute result and left arguments
    template<typename Algo_, typename VT_>
    class SumFunctorProxyResultLeft : public SolverFunctorBase<VT_>
    {
      public:
        SumFunctorProxyResultLeft(VT_& y, VT_& l, const VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "SumFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete SumFunctor can not be executed!");

          LAFEM::Sum<Algo_>::value(_y, _l, _r);
        }

        SumFunctorProxyResultLeft& operator=(const SumFunctorProxyResultLeft& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        SumFunctorProxyResultLeft(const SumFunctorProxyResultLeft& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _y = arg;
          _l = arg;
          this->_complete = true;
        }

      private:
        VT_& _y;
        VT_& _l;
        const VT_& _r;
    };

    ///always substitute all arguments
    template<typename Algo_, typename VT_>
    class SumFunctorProxyAll : public SolverFunctorBase<VT_>
    {
      public:
        SumFunctorProxyAll(VT_& y, VT_& l, VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "SumFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete SumFunctor can not be executed!");

          LAFEM::Sum<Algo_>::value(_y, _l, _r);
        }

        SumFunctorProxyAll& operator=(const SumFunctorProxyAll& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        SumFunctorProxyAll(const SumFunctorProxyAll& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _y = arg;
          _l = arg;
          _r = arg;
          this->_complete = true;
        }

      private:
        VT_& _y;
        VT_& _l;
        VT_& _r;
    };


    template<typename Algo_, typename VT_>
    class SumFunctor : public SolverFunctorBase<VT_>
    {
      public:
        SumFunctor(VT_& y, const VT_& l, const VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "SumFunctor";
        }

        virtual void execute()
        {
          LAFEM::Sum<Algo_>::value(_y, _l, _r);
        }

        SumFunctor& operator=(const SumFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        SumFunctor(const SumFunctor& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
        }

      private:
        VT_& _y;
        const VT_& _l;
        const VT_& _r;
    };

    ///always substitute right argument
    template<typename Algo_, typename VT_, typename MT_>
    class ProductFunctorProxyRight : public SolverFunctorBase<VT_>
    {
      public:
        ProductFunctorProxyRight(VT_& y, const MT_& l, VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "ProductFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete ProductFunctor can not be executed!");

          LAFEM::ProductMatVec<Algo_>::value(_y, _l, _r);
        }

        ProductFunctorProxyRight& operator=(const ProductFunctorProxyRight& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        ProductFunctorProxyRight(const ProductFunctorProxyRight& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _r = arg;
          this->_complete = true;
        }

      private:
        VT_& _y;
        const MT_& _l;
        VT_& _r;
    };

    ///always substitute result and right arguments
    template<typename Algo_, typename VT_, typename MT_>
    class ProductFunctorProxyResultRight : public SolverFunctorBase<VT_>
    {
      public:
        ProductFunctorProxyResultRight(VT_& y, const MT_& l, VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "ProductFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete ProductFunctor can not be executed!");

          LAFEM::ProductMatVec<Algo_>::value(_y, _l, _r);
        }

        ProductFunctorProxyResultRight& operator=(const ProductFunctorProxyResultRight& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        ProductFunctorProxyResultRight(const ProductFunctorProxyResultRight& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _y = arg;
          _r = arg;
          this->_complete = true;
        }

      private:
        VT_& _y;
        const MT_& _l;
        VT_& _r;
    };

    template<typename Algo_, typename VT_, typename MT_>
    class ProductFunctor : public SolverFunctorBase<VT_>
    {
      public:
        ProductFunctor(VT_& y, const MT_& l, const VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "ProductFunctor";
        }

        virtual void execute()
        {
          LAFEM::ProductMatVec<Algo_>::value(_y, _l, _r);
        }

        ProductFunctor& operator=(const ProductFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        ProductFunctor(const ProductFunctor& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
        }

      private:
        VT_& _y;
        const MT_& _l;
        const VT_& _r;
    };

    template<typename Algo_, typename VT_, typename MT_>
    class DefectFunctorProxyRight : public SolverFunctorBase<VT_>
    {
      public:
        DefectFunctorProxyRight(VT_& y, const VT_& l, const MT_& m, VT_& r) :
          _y(y),
          _l(l),
          _m(m),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "DefectFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete DefectFunctor can not be executed!");

          LAFEM::Defect<Algo_>::value(_y, _l, _m, _r);
        }

        DefectFunctorProxyRight& operator=(const DefectFunctorProxyRight& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_m = rhs._r;
          this->_r = rhs._r;
          return *this;
        }

        DefectFunctorProxyRight(const DefectFunctorProxyRight& other) :
          _y(other._y),
          _l(other._l),
          _m(other._m),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _r = arg;
          this->_complete = true;
        }

      private:
        VT_& _y;
        const VT_& _l;
        const MT_& _m;
        VT_& _r;
    };

    template<typename Algo_, typename VT_, typename MT_>
    class DefectFunctorProxyResultRight : public SolverFunctorBase<VT_>
    {
      public:
        DefectFunctorProxyResultRight(VT_& y, const VT_& l, const MT_& m, VT_& r) :
          _y(y),
          _l(l),
          _m(m),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "DefectFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete DefectFunctor can not be executed!");

          LAFEM::Defect<Algo_>::value(_y, _l, _m, _r);
        }

        DefectFunctorProxyResultRight& operator=(const DefectFunctorProxyResultRight& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_m = rhs._r;
          this->_r = rhs._r;
          return *this;
        }

        DefectFunctorProxyResultRight(const DefectFunctorProxyResultRight& other) :
          _y(other._y),
          _l(other._l),
          _m(other._m),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _y = arg;
          _r = arg;
          this->_complete = true;
        }

      private:
        VT_& _y;
        const VT_& _l;
        const MT_& _m;
        VT_& _r;
    };

    template<typename Algo_, typename VT_, typename MT_>
    class DefectFunctor : public SolverFunctorBase<VT_>
    {
      public:
        DefectFunctor(VT_& y, const VT_& l, const MT_& m, const VT_& r) :
          _y(y),
          _l(l),
          _m(m),
          _r(r)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "DefectFunctor";
        }

        virtual void execute()
        {
          LAFEM::Defect<Algo_>::value(_y, _l, _m, _r);
        }

        DefectFunctor& operator=(const DefectFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_m = rhs._r;
          this->_r = rhs._r;
          return *this;
        }

        DefectFunctor(const DefectFunctor& other) :
          _y(other._y),
          _l(other._l),
          _m(other._m),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
        }

      private:
        VT_& _y;
        const VT_& _l;
        const MT_& _m;
        const VT_& _r;
    };


    template<typename Algo_, typename VT_>
    class PreconFunctor : public SolverFunctorBase<VT_>
    {
      public:
        PreconFunctor(VT_& y) :
          _y(y),
          _functor()
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          if(!this->_complete)
            return "PreconFunctor[]";
          else
          {
            std::string result("PreconFunctor[");

            result.append((const std::string)((this->_functor)->type_name()));
            result.append("]");
            return result;
          }
        }

        virtual bool is_preconditioner()
        {
          return true;
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete PreconFunctor can not be executed!");

          _functor->execute();
        }

        PreconFunctor& operator=(const PreconFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_functor = rhs._functor;
          return *this;
        }

        PreconFunctor(const PreconFunctor& other) :
          _y(other._y),
          _functor(other._functor)
        {
        }

        virtual void substitute(VT_& arg)
        {
        }

        void set_precon_functor(std::shared_ptr<SolverFunctorBase<VT_> >& functor)
        {
          functor->substitute(_y);
          _functor = functor;
          this->_complete = true;
        }

      private:
        VT_& _y;
        std::shared_ptr<SolverFunctorBase<VT_> > _functor;
    };

    template<typename Algo_, typename VT_>
    class PreconFunctorProxy : public SolverFunctorBase<VT_>
    {
      public:
        PreconFunctorProxy(VT_& y) :
          _y(y),
          _functor(),
          _precon_complete(false),
          _substitution_complete(false)
        {
          this->_complete = false;
        }

        virtual bool is_preconditioner()
        {
          return true;
        }

        virtual const std::string type_name()
        {
          if(!_precon_complete)
            return "PreconFunctor[]";
          else
          {
            std::string result("PreconFunctor[");
            result.append((const std::string)((this->_functor)->type_name()));
            result.append("]");
            return result;
          }
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete PreconFunctor can not be executed!");

          _functor->execute();
        }

        PreconFunctorProxy& operator=(const PreconFunctorProxy& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_functor = rhs._functor;
          this->_precon_complete = rhs._precon_complete;
          this->_substitution_complete = rhs._substitution_complete;
          return *this;
        }

        PreconFunctorProxy(const PreconFunctorProxy& other) :
          _y(other._y),
          _functor(other._functor),
          _precon_complete(other._precon_complete),
          _substitution_complete(other._substitution_complete)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _y = arg;
          _substitution_complete = true;
          this->_complete = this->_substitution_complete && _precon_complete;
        }

        void set_precon_functor(std::shared_ptr<SolverFunctorBase<VT_> >& functor)
        {
          functor->substitute(_y);
          _functor = functor;
          _precon_complete = true;
          this->_complete = this->_substitution_complete && _precon_complete;
        }

      private:
        VT_& _y;
        std::shared_ptr<SolverFunctorBase<VT_> > _functor;
        bool _precon_complete;
        bool _substitution_complete;
    };

    enum ComparisonOpCode
    {
      coc_eq = 0,
      coc_neq,
      coc_less,
      coc_greater,
      coc_leq,
      coc_geq,
      coc_countloop
    };

    template<typename Algo_, typename VT_, typename ST_>
    class IterateFunctor : public SolverFunctorBase<VT_>
    {
      public:
        IterateFunctor(std::shared_ptr<SolverFunctorBase<VT_> >& functor,
                       ST_& arg1,
                       ST_& arg2,
                       Index& count,
                       Index max_iter = 100,
                       ComparisonOpCode opcode = coc_countloop) :
          _functor(functor),
          _arg1(arg1),
          _arg2(arg2),
          _count(count),
          _maxiter(max_iter),
          _opcode(opcode)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          std::string result("IterateFunctor[");

          result.append((const std::string)((this->_functor)->type_name()));
          result.append("]");
          return result;
        }

        virtual void execute()
        {
          Index count(0);
          while(true)
          {
            _functor->execute();
            ++count;

            bool break_condition(
                                 _opcode == coc_eq ? (_arg1 == _arg2) :
                                 (_opcode == coc_neq ? (_arg1 != _arg2) :
                                 (_opcode == coc_less ? (_arg1 < _arg2) :
                                 (_opcode == coc_greater ? (_arg1 > _arg2) :
                                 (_opcode == coc_leq ? (_arg1 <= _arg2) :
                                 (_opcode == coc_geq ? (_arg1 >= _arg2) :
                                 (count == _maxiter))))))
                                );
            if(break_condition || count == _maxiter)
            {
              _count = count;
              break;
            }
          }
        }

        IterateFunctor& operator=(const IterateFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_functor = rhs._functor;
          this->_arg1 = rhs._arg1;
          this->_arg2 = rhs._arg2;
          this->_count = rhs._count;
          this->_maxiter = rhs._maxiter;
          this->_opcode = rhs._opcode;
          return *this;
        }

        IterateFunctor(const IterateFunctor& other) :
          _functor(other._functor),
          _arg1(other._arg1),
          _arg2(other._arg2),
          _count(other._count),
          _maxiter(other._maxiter),
          _opcode(other._opcode)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _functor->substitute(arg);
        }

      private:
        std::shared_ptr<SolverFunctorBase<VT_> > _functor;
        ST_& _arg1;
        ST_& _arg2;
        Index& _count;
        Index _maxiter;
        ComparisonOpCode _opcode;
    };

    template<typename Algo_, typename VT_>
    class CopyFunctor : public SolverFunctorBase<VT_>
    {
      public:
        CopyFunctor(VT_& y, const VT_& x) :
          _y(y),
          _x(x)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "CopyFunctor";
        }

        virtual void execute()
        {
          copy(_y, _x);
        }

        CopyFunctor& operator=(const CopyFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          return *this;
        }

        CopyFunctor(const CopyFunctor& other) :
          _y(other._y),
          _x(other._x)
        {
        }

        virtual void substitute(VT_& arg)
        {
        }

      private:
        VT_& _y;
        const VT_& _x;
    };

    template<typename Algo_, typename VT_, typename ST_>
    class NormFunctor : public SolverFunctorBase<VT_>
    {
      public:
        NormFunctor(ST_& y, const VT_& x) :
          _y(y),
          _x(x)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "NormFunctor";
        }

        virtual void execute()
        {
          _y = LAFEM::Norm2<Algo_>::value(_x);
        }

        NormFunctor& operator=(const NormFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          return *this;
        }

        NormFunctor(const NormFunctor& other) :
          _y(other._y),
          _x(other._x)
        {
        }

        virtual void substitute(VT_& arg)
        {
        }

      private:
        ST_& _y;
        const VT_& _x;
    };

    template<typename Algo_, typename VT_, typename ST_>
    class NormFunctorProxyRight : public SolverFunctorBase<VT_>
    {
      public:
        NormFunctorProxyRight(ST_& y, VT_& x) :
          _y(y),
          _x(x)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "NormFunctor";
        }

        virtual void execute()
        {
          _y = LAFEM::Norm2<Algo_>::value(_x);
        }

        NormFunctorProxyRight& operator=(const NormFunctorProxyRight& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          return *this;
        }

        NormFunctorProxyRight(const NormFunctorProxyRight& other) :
          _y(other._y),
          _x(other._x)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _x = arg;
          this->_complete = true;
        }

      private:
        ST_& _y;
        VT_& _x;
    };

    template<typename Algo_, typename VT_, typename HaloStorageType_>
    class SynchFunctor : public SolverFunctorBase<VT_>
    {
      public:
        SynchFunctor(VT_& y, const VT_& x, const HaloStorageType_& halos) :
          _y(y),
          _x(x),
          _halos(halos)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "SynchFunctor";
        }

        virtual void execute()
        {
          ///TODO
          ///Foundation::Synchronisation(_y, _x, _halos);
        }

        SynchFunctor& operator=(const SynchFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          this->_halos = rhs._halos;
          return *this;
        }

        SynchFunctor(const SynchFunctor& other) :
          _y(other._y),
          _x(other._x),
          _halos(other._halos)
        {
        }

        virtual void substitute(VT_& arg)
        {
        }

      private:
        VT_& _y;
        const VT_& _x;
        const HaloStorageType_& _halos;
    };

    ///in scalar functors, we still need a vector type in order to make the functor compatible with SolverFunctorBase<VT_>
    template<typename VT_, typename DT_>
    class DivFunctor : public SolverFunctorBase<VT_>
    {
      public:
        DivFunctor(DT_& y, const DT_& l, const DT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "DivFunctor";
        }

        virtual void execute()
        {
          _y = _l / _r;
        }

        DivFunctor& operator=(const DivFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        DivFunctor(const DivFunctor& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& arg)
        {
        }

      private:
        DT_& _y;
        const DT_& _l;
        const DT_& _r;
    };

    template<typename Algo_, typename VT_, template<typename, typename> class StorageType_ = std::vector>
    class CompoundSolverFunctor : public SolverFunctorBase<VT_>
    {
      public:
        typedef StorageType_<std::shared_ptr<SolverFunctorBase<VT_> >, std::allocator<std::shared_ptr<SolverFunctorBase<VT_> > > > storage_type_;

        CompoundSolverFunctor() :
          _functors()
        {
          this->_complete = false;
        }

        virtual void substitute(VT_& arg)
        {
          for(Index i(0) ; i < (this->_functors).size() ; ++i)
          {
            (this->_functors).at(i)->substitute(arg);
          }
        }

        virtual void execute()
        {
          for(Index i(0) ; i < (this->_functors).size() ; ++i)
          {
            (this->_functors).at(i)->execute();
          }
        }

        virtual const std::string type_name()
        {
          std::string result("[");

          for(Index i(0) ; i < (this->_functors).size() ; ++i)
          {
            if(i != 0)
              result.append(", ");

            result.append((const std::string)((this->_functors).at(i)->type_name()));
          }
          result.append("]");
          return result;
        }

        void add_functor(SolverFunctorBase<VT_>* functor)
        {
          _functors.push_back(std::shared_ptr<SolverFunctorBase<VT_> >(functor));
        }

        void add_functor(const std::shared_ptr<SolverFunctorBase<VT_> >& functor)
        {
          _functors.push_back(functor);
        }

        storage_type_& get_functors()
        {
          return _functors;
        }

        const Index size()
        {
          return Index(_functors.size());
        }

        const storage_type_& get_functors() const
        {
          return _functors;
        }

        Index size() const
        {
          return Index(_functors.size());
        }

        CompoundSolverFunctor& operator=(const CompoundSolverFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_functors = rhs._functors;
          return *this;
        }

        CompoundSolverFunctor(const CompoundSolverFunctor& other) :
          _functors(other._functors)
        {
        }

        template<typename OVT_>
        void set_preconditioner(std::shared_ptr<SolverFunctorBase<OVT_> >& functor)
        {
          for(Index i(0) ; i < (this->_functors).size() ; ++i)
          {
            if( (this->_functors).at(i)->is_preconditioner())
            {
              ((PreconFunctor<Algo_, VT_>*)((this->_functors).at(i).get()))->set_precon_functor(functor);
              return;
            }
          }
          throw ScaRCError("Error: No PreconFunctor in functor list!");
        }

      protected:
        storage_type_ _functors;

    };
  }
}

#endif
