/**
 * \file
 * \brief FEAST milestone 1 ScaRC functor implementations
 * \author Markus Geveler
 * \date 2012 - 2013
 *
 * See class documentation.
 *
 */

#pragma once
#ifndef SCARC_GUARD_SOLVER_FUNCTOR_HH
#define SCARC_GUARD_SOLVER_FUNCTOR_HH 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/foundation/synch.hpp>

using namespace FEAST::Foundation;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {
    /**
     * \brief Base class for all solver functors
     *
     * Solver functors store references to variable data containers and values for constants.
     * The specific SolverFunctor subclass specifies the operation by overwriting its execute() member function.
     * Hence, no function pointers and bindings are required for the functor library implementation.
     *
     * \tparam VT_
     * Vector type of the solution space, necesssary for substitution function
     *
     * \author Markus Geveler
     */
    template<typename VT_>
    class SolverFunctorBase
    {
      public:
        /**
         * \brief substitutes specific references by renference to arg
         *
         * ProxyFunctors can overwrite this function in order to 'replace' references to dummy containers.
         *
         * \param[in] arg
         * the substitute
         *
         */
        virtual void substitute(VT_& arg) = 0;

        /**
         * \brief executes the functor on the stored arguments
         *
         */
        virtual void execute() = 0;

        ///returns a string that describes the functor
        virtual const std::string type_name() = 0;

        ///virtual DTOR
        virtual ~SolverFunctorBase()
        {
        }

        /**
         * \brief substitutes a preconditioner by a functor
         *
         * Preconditioner functors can overwrite this function. With this mechanism, solver layers can be created
         * independently and plugged into each other dynamically.
         *
         * \param[in] precon
         * the substitute
         *
         */
        virtual void set_preconditioner(std::shared_ptr<SolverFunctorBase<VT_> >& /*precon*/)
        {
        }

        ///better visualisation of type_names
        static String pretty_printer(String input)
        {
          String output;
          String d_indent(">-<");
          String indent;
          Index last_token(0);

          for (Index i(0) ; i < input.length() ; ++i)
          {
            if(input[i] == '[')
            {
              output += indent + input.substr(last_token, i - last_token) + "\n";
              last_token = i+1;
              d_indent.replace(1, 1, stringify(indent.size() / d_indent.size() + 1));
              indent += d_indent;
            }
            else if(input[i] == ']')
            {
              if (input[i-1] != ']')
                output += indent + input.substr(last_token, i - last_token) + "\n";
              last_token = i+1;
              indent.erase(indent.size() - d_indent.size(), d_indent.size());
            }
            else if(input[i] == ',')
            {
              if (input[i-1] != ']')
                output += indent + input.substr(last_token, i - last_token) + "\n";
              last_token = i+2;
            }

          }
          return output;
        }

      protected:
        ///flag for signalling that all substitutions are done - ProxyFunctor CTORs must init it with false
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

          _y.template sum<Algo_>(_l, _r);
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

          _y.template sum<Algo_>(_l, _r);
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

          _y.template sum<Algo_>(_l, _r);
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
          _y.template sum<Algo_>(_l, _r);
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

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        VT_& _y;
        const VT_& _l;
        const VT_& _r;
    };

    ///always substitute left argument
    template<typename Algo_, typename VT_>
    class DifferenceFunctorProxyLeft : public SolverFunctorBase<VT_>
    {
      public:
        DifferenceFunctorProxyLeft(VT_& y, VT_& l, const VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "DifferenceFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete DifferenceFunctor can not be executed!");

          _y.template difference<Algo_>( _l, _r);
        }

        DifferenceFunctorProxyLeft& operator=(const DifferenceFunctorProxyLeft& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        DifferenceFunctorProxyLeft(const DifferenceFunctorProxyLeft& other) :
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
    class DifferenceFunctorProxyResultLeft : public SolverFunctorBase<VT_>
    {
      public:
        DifferenceFunctorProxyResultLeft(VT_& y, VT_& l, const VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "DifferenceFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete DifferenceFunctor can not be executed!");

          _y.template difference<Algo_>( _l, _r);
        }

        DifferenceFunctorProxyResultLeft& operator=(const DifferenceFunctorProxyResultLeft& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        DifferenceFunctorProxyResultLeft(const DifferenceFunctorProxyResultLeft& other) :
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
    class DifferenceFunctorProxyAll : public SolverFunctorBase<VT_>
    {
      public:
        DifferenceFunctorProxyAll(VT_& y, VT_& l, VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "DifferenceFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete DifferenceFunctor can not be executed!");

          _y.template difference<Algo_>( _l, _r);
        }

        DifferenceFunctorProxyAll& operator=(const DifferenceFunctorProxyAll& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        DifferenceFunctorProxyAll(const DifferenceFunctorProxyAll& other) :
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
    class DifferenceFunctor : public SolverFunctorBase<VT_>
    {
      public:
        DifferenceFunctor(VT_& y, const VT_& l, const VT_& r) :
          _y(y),
          _l(l),
          _r(r)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "DifferenceFunctor";
        }

        virtual void execute()
        {
          _y.template difference<Algo_>( _l, _r);
        }

        DifferenceFunctor& operator=(const DifferenceFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        DifferenceFunctor(const DifferenceFunctor& other) :
          _y(other._y),
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& /*arg*/)
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

          _y.template product_matvec<Algo_>(_l, _r);
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

          _y.template product_matvec<Algo_>(_l, _r);
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
          _y.template product_matvec<Algo_>(_l, _r);
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

        virtual void substitute(VT_& /*arg*/)
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

          _y.template defect<Algo_>(_l, _m, _r);
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

          _y.template defect<Algo_>(_l, _m, _r);
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
          _y.template defect<Algo_>(_l, _m, _r);
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

        virtual void substitute(VT_& /*arg*/)
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

        virtual void substitute(VT_& /*arg*/)
        {
        }

        void set_preconditioner(std::shared_ptr<SolverFunctorBase<VT_> >& functor)
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

        void set_preconditioner(std::shared_ptr<SolverFunctorBase<VT_> >& functor)
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
                       const ST_& arg1,
                       ST_ arg2,
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

        virtual bool is_solver_loop()
        {
          return true;
        }

        virtual void execute()
        {
          if(_maxiter == 0)
            return;

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
            if(break_condition || (count == _maxiter))
            {
              //if(_opcode != coc_countloop)
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

        virtual void set_preconditioner(std::shared_ptr<SolverFunctorBase<VT_> >& precon)
        {
          _functor->set_preconditioner(precon);
        }

      private:
        std::shared_ptr<SolverFunctorBase<VT_> > _functor;
        const ST_& _arg1;
        ST_ _arg2;
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

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        VT_& _y;
        const VT_& _x;
    };

    template<typename Algo_, typename VT_>
    class CopyFunctorProxyResult : public SolverFunctorBase<VT_>
    {
      public:
        CopyFunctorProxyResult(VT_& y, const VT_& x) :
          _y(y),
          _x(x)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "CopyFunctorProxyResult";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete SumFunctor can not be executed!");

          copy(_y, _x);
        }

        CopyFunctorProxyResult& operator=(const CopyFunctorProxyResult& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          return *this;
        }

        CopyFunctorProxyResult(const CopyFunctorProxyResult& other) :
          _y(other._y),
          _x(other._x)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _y = arg;
          this->_complete = true;
        }

      private:
        VT_& _y;
        const VT_& _x;
    };

    template<typename Algo_, typename VT_, typename ST_>
    class NormFunctor2wosqrt : public SolverFunctorBase<VT_>
    {
      public:
        NormFunctor2wosqrt(ST_& y, const VT_& x) :
          _y(y),
          _x(x)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "NormFunctor2wosqrt";
        }

        virtual void execute()
        {
          _y = _x.template norm2wosqrt<Algo_>();
        }

        NormFunctor2wosqrt& operator=(const NormFunctor2wosqrt& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          return *this;
        }

        NormFunctor2wosqrt(const NormFunctor2wosqrt& other) :
          _y(other._y),
          _x(other._x)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        ST_& _y;
        const VT_& _x;
    };

    template<typename Algo_, typename VT_, typename ST_>
    class NormFunctor2wosqrtProxyRight : public SolverFunctorBase<VT_>
    {
      public:
        NormFunctor2wosqrtProxyRight(ST_& y, VT_& x) :
          _y(y),
          _x(x)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "NormFunctor2wosqrt";
        }

        virtual void execute()
        {
          _y = _x.template norm2wosqrt<Algo_>();
        }

        NormFunctor2wosqrtProxyRight& operator=(const NormFunctor2wosqrtProxyRight& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          return *this;
        }

        NormFunctor2wosqrtProxyRight(const NormFunctor2wosqrtProxyRight& other) :
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

    template<typename Algo_, typename VT_, typename ST_>
    class NormFunctor2 : public SolverFunctorBase<VT_>
    {
      public:
        NormFunctor2(ST_& y, const VT_& x) :
          _y(y),
          _x(x)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "NormFunctor2";
        }

        virtual void execute()
        {
          _y = _x.template norm2<Algo_>();
        }

        NormFunctor2& operator=(const NormFunctor2& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          return *this;
        }

        NormFunctor2(const NormFunctor2& other) :
          _y(other._y),
          _x(other._x)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        ST_& _y;
        const VT_& _x;
    };

    template<typename Algo_, typename VT_, typename ST_>
    class NormFunctor2ProxyRight : public SolverFunctorBase<VT_>
    {
      public:
        NormFunctor2ProxyRight(ST_& y, VT_& x) :
          _y(y),
          _x(x)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "NormFunctor2";
        }

        virtual void execute()
        {
          _y = _x.template norm2<Algo_>();
        }

        NormFunctor2ProxyRight& operator=(const NormFunctor2ProxyRight& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_x = rhs._x;
          return *this;
        }

        NormFunctor2ProxyRight(const NormFunctor2ProxyRight& other) :
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

        virtual void substitute(VT_& /*arg*/)
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
          std::string result("CompoundSolverFunctor[");

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

        Index size()
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

        void set_preconditioner(std::shared_ptr<SolverFunctorBase<VT_> >& functor)
        {
          for(Index i(0) ; i < (this->_functors).size() ; ++i)
          {
              ((this->_functors).at(i).get())->set_preconditioner(functor);
          }
        }

      protected:
        storage_type_ _functors;

    };

    template<typename Algo_,
             typename VT_,
             typename VMT_,
             Tier2CommModes cm_,
             template<typename, typename> class StoreT_ = std::vector>
    class SynchVecFunctor : public SolverFunctorBase<VT_>
    {
      public:
        SynchVecFunctor(VT_& l,
                        StoreT_<VMT_, std::allocator<VMT_> >& mirrors,
                        StoreT_<VT_, std::allocator<VT_> >& sendbufs,
                        StoreT_<VT_, std::allocator<VT_> >& recvbufs,
                        StoreT_<Index, std::allocator<Index> >& dest_ranks,
                        StoreT_<Index, std::allocator<Index> >& source_ranks) :
          _l(l),
          _mirrors(mirrors),
          _sendbufs(sendbufs),
          _recvbufs(recvbufs),
          _dest_ranks(dest_ranks),
          _source_ranks(source_ranks)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "SynchVecFunctor";
        }

        virtual void execute()
        {
#ifndef SERIAL
          SynchVec<Algo_, Parallel, cm_>::execute(_l, _mirrors, _sendbufs, _recvbufs, _dest_ranks, _source_ranks);
#else
          SynchVec<Algo_, Serial, cm_>::execute(_l, _mirrors, _sendbufs, _recvbufs, _dest_ranks, _source_ranks);
#endif
        }

        SynchVecFunctor& operator=(const SynchVecFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_l = rhs._l;
          this->_mirrors = rhs._mirrors;
          this->_sendbufs = rhs._sendbufs;
          this->_recvbufs = rhs._recvbufs;
          this->_dest_ranks(_dest_ranks);
          this->_source_ranks(_source_ranks);

          return *this;
        }

        SynchVecFunctor(const SynchVecFunctor& other) :
          _l(other._l),
          _mirrors(other._mirrors),
          _sendbufs(other._sendbufs),
          _recvbufs(other._recvbufs),
          _dest_ranks(other._dest_ranks),
          _source_ranks(other._source_ranks)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        VT_& _l;
        StoreT_<VMT_, std::allocator<VMT_> >& _mirrors;
        StoreT_<VT_, std::allocator<VT_> >& _sendbufs;
        StoreT_<VT_, std::allocator<VT_> >& _recvbufs;
        StoreT_<Index, std::allocator<Index> >& _dest_ranks;
        StoreT_<Index, std::allocator<Index> >& _source_ranks;
    };

    template<typename Algo_,
             typename VT_,
             typename VMT_,
             Tier2CommModes cm_,
             template<typename, typename> class StoreT_ = std::vector>
    class SynchVecFunctorProxy : public SolverFunctorBase<VT_>
    {
      public:
        SynchVecFunctorProxy(VT_& l,
                             StoreT_<VMT_, std::allocator<VMT_> >& mirrors,
                             StoreT_<VT_, std::allocator<VT_> >& sendbufs,
                             StoreT_<VT_, std::allocator<VT_> >& recvbufs,
                             StoreT_<Index, std::allocator<Index> >& dest_ranks,
                             StoreT_<Index, std::allocator<Index> >& source_ranks) :
          _l(l),
          _mirrors(mirrors),
          _sendbufs(sendbufs),
          _recvbufs(recvbufs),
          _dest_ranks(dest_ranks),
          _source_ranks(source_ranks)
        {
          this->_complete = false;
        }

        virtual const std::string type_name()
        {
          return "SynchVecFunctor";
        }

        virtual void execute()
        {
          if(!this->_complete)
            throw ScaRCError("Error: Incomplete SynchVecFunctor can not be executed!");

#ifndef SERIAL
          SynchVec<Algo_, Parallel, cm_>::execute(_l, _mirrors, _sendbufs, _recvbufs, _dest_ranks, _source_ranks);
#else
          SynchVec<Algo_, Serial, cm_>::execute(_l, _mirrors, _sendbufs, _recvbufs, _dest_ranks, _source_ranks);
#endif
        }

        SynchVecFunctorProxy& operator=(const SynchVecFunctorProxy& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_l = rhs._l;
          this->_mirrors = rhs._mirrors;
          this->_sendbufs = rhs._sendbufs;
          this->_recvbufs = rhs._recvbufs;
          this->_dest_ranks(_dest_ranks);
          this->_source_ranks(_source_ranks);

          return *this;
        }

        SynchVecFunctorProxy(const SynchVecFunctorProxy& other) :
          _l(other._l),
          _mirrors(other._mirrors),
          _sendbufs(other._sendbufs),
          _recvbufs(other._recvbufs),
          _dest_ranks(other._dest_ranks),
          _source_ranks(other._source_ranks)
        {
        }

        virtual void substitute(VT_& arg)
        {
          _l = arg;
          this->_complete = true;
        }

      private:
        VT_& _l;
        StoreT_<VMT_, std::allocator<VMT_> >& _mirrors;
        StoreT_<VT_, std::allocator<VT_> >& _sendbufs;
        StoreT_<VT_, std::allocator<VT_> >& _recvbufs;
        StoreT_<Index, std::allocator<Index> >& _dest_ranks;
        StoreT_<Index, std::allocator<Index> >& _source_ranks;
    };

    template<typename Algo_,
             typename VT_,
             typename DT_,
             Tier2CommModes cm_>
    class SynchScalFunctor : public SolverFunctorBase<VT_>
    {
      public:
        SynchScalFunctor(DT_& l,
                         DT_& sendbuf,
                         DT_& recvbuf) :
          _l(l),
          _sendbuf(sendbuf),
          _recvbuf(recvbuf)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "SynchScalFunctor";
        }

        virtual void execute()
        {
#ifndef SERIAL
          SynchScal<Parallel, cm_>::execute(_l, _sendbuf, _recvbuf);
#else
          SynchScal<Serial, cm_>::execute(_l, _sendbuf, _recvbuf);
#endif
        }

        SynchScalFunctor& operator=(const SynchScalFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_l = rhs._l;
          this->_sendbuf = rhs._sendbuf;
          this->_recvbuf = rhs._recvbuf;
          return *this;
        }

        SynchScalFunctor(const SynchScalFunctor& other) :
          _l(other._l),
          _sendbuf(other._sendbuf),
          _recvbuf(other._recvbuf)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        DT_& _l;
        DT_& _sendbuf;
        DT_& _recvbuf;
    };

    ///filter functors
    template<typename Algo_, typename VT_, typename FilterType_>
    class FilterDefectFunctor : public SolverFunctorBase<VT_>
    {
      public:
        FilterDefectFunctor(VT_& l, const FilterType_& r) :
          _l(l),
          _r(r)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "FilterDefectFunctor";
        }

        virtual void execute()
        {
          _r.template filter_def<Algo_>(_l);
        }

        FilterDefectFunctor& operator=(const FilterDefectFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        FilterDefectFunctor(const FilterDefectFunctor& other) :
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        VT_& _l;
        const FilterType_& _r;
    };

    template<typename Algo_, typename VT_, typename FilterType_>
    class FilterCorrectionFunctor : public SolverFunctorBase<VT_>
    {
      public:
        FilterCorrectionFunctor(VT_& l, const FilterType_& r) :
          _l(l),
          _r(r)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "FilterCorrectionFunctor";
        }

        virtual void execute()
        {
          _r.template filter_cor<Algo_>(_l);
        }

        FilterCorrectionFunctor& operator=(const FilterCorrectionFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_l = rhs._l;
          this->_r = rhs._r;
          return *this;
        }

        FilterCorrectionFunctor(const FilterCorrectionFunctor& other) :
          _l(other._l),
          _r(other._r)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        VT_& _l;
        const FilterType_& _r;
    };

    template<typename Algo_, typename VT_, typename T_>
    class InspectionFunctor : public SolverFunctorBase<VT_>
    {
      public:
        InspectionFunctor(const T_& y, const std::string tracetag) :
          _y(y),
          _tracetag(tracetag)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "InspectionFunctor";
        }

        virtual void execute()
        {
          std::cout << _tracetag << "| element of size " << sizeof(T_) << " at address " << &_y << " with value " << T_(_y) << " !" << std::endl;
        }

        InspectionFunctor& operator=(const InspectionFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_tracetag = rhs._tracetag;
          return *this;
        }

        InspectionFunctor(const InspectionFunctor& other) :
          _y(other._y),
          _tracetag(other._tracetag)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        const T_& _y;
        const std::string _tracetag;
    };

    template<typename Algo_, typename VT_, typename T_, template<typename, typename> class StoT_ = std::vector>
    class InspectionFunctorSTL : public SolverFunctorBase<VT_>
    {
      public:
        InspectionFunctorSTL(const StoT_<T_, std::allocator<T_> >& y, Index index, const std::string tracetag) :
          _y(y),
          _index(index),
          _tracetag(tracetag)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "InspectionFunctorSTL";
        }

        virtual void execute()
        {
          std::cout << _tracetag << "| element of size " << sizeof(T_) << " at index " << _index << " with value " << _y.at(_index) << " !" << std::endl;
        }

        InspectionFunctorSTL& operator=(const InspectionFunctorSTL& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_index = rhs._index;
          this->_tracetag = rhs._tracetag;
          return *this;
        }

        InspectionFunctorSTL(const InspectionFunctorSTL& other) :
          _y(other._y),
          _index(other._index),
          _tracetag(other._tracetag)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        const StoT_<T_, std::allocator<T_> >& _y;
        Index _index;
        const std::string _tracetag;
    };

    template<typename Algo_, typename VT_, typename T_>
    class InspectionConstFunctor : public SolverFunctorBase<VT_>
    {
      public:
        InspectionConstFunctor(const T_& y, const std::string tracetag) :
          _y(y),
          _tracetag(tracetag)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "InspectionConstFunctor";
        }

        virtual void execute()
        {
          std::cout << _tracetag << "| element of size " << sizeof(T_) << " at address " << &_y << " with value " << T_(_y) << " !" << std::endl;
        }

        InspectionConstFunctor& operator=(const InspectionConstFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          this->_tracetag = rhs._tracetag;
          return *this;
        }

        InspectionConstFunctor(const InspectionConstFunctor& other) :
          _y(other._y),
          _tracetag(other._tracetag)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        const T_ _y;
        const std::string _tracetag;
    };

    template<typename Algo_, typename VT_>
    class SetFunctor : public SolverFunctorBase<VT_>
    {
      public:
        SetFunctor(VT_& y, typename VT_::DataType val) :
          _y(y),
          _val(val)
        {
          this->_complete = true;
        }

        virtual const std::string type_name()
        {
          return "SetFunctor";
        }

        virtual void execute()
        {
          for(Index i(0) ; i < _y.size() ; ++i)
            _y(i, _val);
        }

        SetFunctor& operator=(const SetFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_y = rhs._y;
          return *this;
        }

        SetFunctor(const SetFunctor& other) :
          _y(other._y),
          _val(other._val)
        {
        }

        virtual void substitute(VT_& /*arg*/)
        {
        }

      private:
        VT_& _y;
        const typename VT_::DataType _val;
    };
  }
}

#endif
