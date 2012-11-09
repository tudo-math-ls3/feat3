#pragma once
#ifndef SCARC_GUARD_SOLVER_FUNCTOR_HH
#define SCARC_GUARD_SOLVER_FUNCTOR_HH 1

#include<kernel/base_header.hpp>
#include<kernel/foundation/functor.hpp>
#include<kernel/scarc/matrix.hpp>
#include<kernel/scarc/vector.hpp>

using namespace FEAST::Foundation;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {
    template<typename T_>
    class ProxyMatrix : public FunctorBase
    {
      public:
        ProxyMatrix()
        {
        }

        T_& cast()
        {
          return static_cast<T_&>(*this);
        }

        const T_& cast() const
        {
          return static_cast<const T_&>(*this);
        }

        virtual const std::string type_name()
        {
          return "ProxyMatrix";
        }

        virtual void execute()
        {
        }

        virtual void undo()
        {
        }

        ProxyMatrix& operator=(const ProxyMatrix& rhs)
        {
          if(this == &rhs)
            return *this;

          return *this;
        }

        ProxyMatrix(const ProxyMatrix& other)
        {
        }
    };

    class MatrixData : public ProxyMatrix<MatrixData>
    {
    };

    template<typename T_>
    class ProxyVector : public FunctorBase
    {
      public:
        ProxyVector()
        {
        }

        T_& cast()
        {
          return static_cast<T_&>(*this);
        }

        const T_& cast() const
        {
          return static_cast<const T_&>(*this);
        }

        virtual const std::string type_name()
        {
          return "ProxyVector";
        }

        virtual void execute()
        {
        }

        virtual void undo()
        {
        }

        ProxyVector& operator=(const ProxyVector& rhs)
        {
          if(this == &rhs)
            return *this;

          return *this;
        }

        ProxyVector(const ProxyVector& other)
        {
        }
    };

    class VectorData : public ProxyVector<VectorData>
    {
    };


    template<typename T1_, typename T2_>
    class ProxyVectorSum : public ProxyVector<ProxyVectorSum<T1_, T2_> >
    {
      public:
        ProxyVectorSum(T1_& left, T2_& right) :
          _left(left),
          _right(right)
        {
        }

        ProxyVectorSum& operator=(const ProxyVectorSum& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_left = rhs._left;
          this->_right = rhs._right;

          return *this;
        }

        ProxyVectorSum(const ProxyVectorSum& other) :
          _left(other._left),
          _right(other._right)
        {
        }

        const std::string type_name()
        {
          return _left.cast().type_name() + " + " + _right.cast().type_name();
        }

      private:
        T1_& _left;
        T2_& _right;
    };

    class UninitializedProxyPreconApply : public ProxyVector<UninitializedProxyPreconApply >
    {
      public:
        UninitializedProxyPreconApply()
        {
        }

        UninitializedProxyPreconApply& operator=(const UninitializedProxyPreconApply& rhs)
        {
          if(this == &rhs)
            return *this;

          return *this;
        }

        UninitializedProxyPreconApply(const UninitializedProxyPreconApply& other)
        {
        }

        virtual const std::string type_name()
        {
          return "__UNINITIALIZED_PRECONDITIONER__()";
        }
    };

    class ProxyPreconApply : public ProxyVector<ProxyPreconApply>
    {
      public:
        ProxyPreconApply() :
          _right(std::shared_ptr<FunctorBase>(new UninitializedProxyPreconApply()))
        {
        }

        ProxyPreconApply(std::shared_ptr<FunctorBase>& right) :
          _right(right)
        {
        }

        ProxyPreconApply& operator=(const ProxyPreconApply& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_right = rhs._right;

          return *this;
        }

        ProxyPreconApply(const ProxyPreconApply& other) :
          _right(other._right)
        {
        }

        virtual const std::string type_name()
        {
          return "__precon__(" + _right.get()->type_name() + ")";
        }

        std::shared_ptr<FunctorBase>& get()
        {
          return _right;
        }

        const std::shared_ptr<FunctorBase>& get() const
        {
          return _right;
        }

      private:
        std::shared_ptr<FunctorBase> _right;
    };


    template<typename T1_, typename T2_, typename T3_>
    class ProxyDefect : public ProxyVector<ProxyDefect<T1_, T2_, T3_> >
    {
      public:
        ProxyDefect(T1_& left, T2_& middle, T3_& right) :
          _left(left),
          _middle(middle),
          _right(right)
        {
        }

        ProxyDefect& operator=(const ProxyDefect& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_left = rhs._left;
          this->_middle = rhs._middle;
          this->_right = rhs._right;

          return *this;
        }

        ProxyDefect(const ProxyDefect& other) :
          _left(other._left),
          _middle(other._middle),
          _right(other._right)
        {
        }

        virtual const std::string type_name()
        {
          return "__defect__(" + _left.cast().type_name() + "," + _middle.cast().type_name() + "," + _right.cast().type_name() + ")";
        }

      private:
        T1_& _left;
        T2_& _middle;
        T3_& _right;
    };

    class Richardson
    {
    };

    class RichardsonLocal
    {
    };

    template<typename Pattern_>
    class SolverPatternGeneration
    {
    };

    template<>
    class SolverPatternGeneration<Richardson>
    {
      public:

        //typedef std::shared_ptr<ProxyVectorSum<VectorData, ProxyPreconApply> > result_type_;

        static std::shared_ptr<FunctorBase> execute(std::shared_ptr<VectorData>& x,
            std::shared_ptr<FunctorBase>& p)
        {
          p = std::shared_ptr<FunctorBase>(new ProxyPreconApply);
          return std::shared_ptr<FunctorBase>(new ProxyVectorSum<VectorData, ProxyPreconApply>(*(x.get()), *((ProxyPreconApply*)(p.get()))));
        }
    };

    template<>
      class SolverPatternGeneration<RichardsonLocal>
      {
        public:

          //typedef std::shared_ptr<ProxyVectorSum<VectorData, ProxyDefect<VectorData, MatrixData, VectorData> > > result_type_;

          static std::shared_ptr<FunctorBase> execute(std::shared_ptr<MatrixData>& A,
              std::shared_ptr<VectorData>& x,
              std::shared_ptr<VectorData>& b,
              std::shared_ptr<FunctorBase>& p)
          {
            std::shared_ptr<FunctorBase> defect(new ProxyDefect<VectorData, MatrixData, VectorData>(*(b.get()), *(A.get()), *(x.get())));
            p = std::shared_ptr<FunctorBase>(new ProxyPreconApply(defect));

            return std::shared_ptr<FunctorBase>(new ProxyVectorSum<VectorData, ProxyPreconApply>(*(x.get()), *((ProxyPreconApply*)(p.get()))));
          }
      };

    /*template<
      template<typename, typename> class StorageType_ = std::vector
    >
    class SolverProxyFunctor : public FunctorBase
    {
      public:
        SolverProxyFunctor(CompoundFunctor<StorageType_>& pattern) :
          _pattern(pattern),
        {
        }

        virtual const std::string type_name()
        {
          return "SolverProxyFunctor";
        }

        virtual void execute()
        {
        }

        virtual void undo()
        {
        }

        SolverProxyFunctor& operator=(const SolverProxyFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_pattern = rhs._pattern;

          return *this;
        }

        SolverProxyFunctor(const SolverProxyFunctor& other) :
          _pattern(other._pattern),
        {
        }

      private:
        CompoundFunctor<StorageType_> _pattern;
    };*/

/*
    template<
      typename DT_ = double,
      template<typename, typename> class StorageType_ = std::vector,
      typename MT_ = DynamicAOSMatrix<StorageType_>,
      typename DT_ = DynamicVector<StorageType_>,
      typename IndexType_ = Index
    >
    class OpFunctor
    {
      private:
        std::shared_ptr<FunctorBase> _preconditioner_functor;
        std::shared_ptr<MT_> _A;
        std::shared_ptr<VT_> _x;
        std::shared_ptr<VT_> _b;
    } */
  }
}

#endif
