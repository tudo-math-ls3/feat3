#pragma once
#ifndef SCARC_GUARD_SOLVER_EXPR_HH
#define SCARC_GUARD_SOLVER_EXPR_HH 1

#include <kernel/util/cpp11_smart_pointer.hpp>

using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {

    ///basic expressions
    template<typename T_>
    class VecExpr
    {
      public:
        T_& cast()
        {
          return static_cast<T_&>(*this);
        }

        const T_& cast() const
        {
          return static_cast<const T_&>(*this);
        }

        std::string get_id()
        {
          return this->cast().get_id();
        }

        const std::string get_id() const
        {
          return this->cast().get_id();
        }
    };

    template<typename T_>
    class MatExpr
    {
      public:
        T_& cast()
        {
          return static_cast<T_&>(*this);
        }

        const T_& cast() const
        {
          return static_cast<const T_&>(*this);
        }

        std::string get_id()
        {
          return this->cast().get_id();
        }

        const std::string get_id() const
        {
          return this->cast().get_id();
        }
    };

    template<typename T_>
    class ScalExpr
    {
      public:
        T_& cast()
        {
          return static_cast<T_&>(*this);
        }

        const T_& cast() const
        {
          return static_cast<const T_&>(*this);
        }

        std::string get_id()
        {
          return this->cast().get_id();
        }

        const std::string get_id() const
        {
          return this->cast().get_id();
        }
    };

    template<typename T_>
    class BoolExpr
    {
      public:
        T_& cast()
        {
          return static_cast<T_&>(*this);
        }

        const T_& cast() const
        {
          return static_cast<const T_&>(*this);
        }

        std::string get_id()
        {
          return this->cast().get_id();
        }

        const std::string get_id() const
        {
          return this->cast().get_id();
        }
    };

    ///VirtualGlobalVector
    class VGVector : public VecExpr<VGVector>
    {
      public:
        ///construct from any VecExpr
        template<typename T_>
        VGVector(VecExpr<T_> const& other)
        {
          _id = other.cast().get_id();
        }

        ///construct from id
        VGVector(const std::string id) :
          _id(id)
        {
        }

        std::string get_id()
        {
          return _id;
        }

        const std::string get_id() const
        {
          return _id;
        }

      private:
        std::string _id;
    };

    class Vector : public VecExpr<Vector>
    {
      public:
        ///construct from any VecExpr
        template<typename T_>
        Vector(VecExpr<T_> const& other)
        {
          _id = "[" + other.cast().get_id() +"]_chunk";
        }

        std::string get_id()
        {
          return _id;
        }

        const std::string get_id() const
        {
          return _id;
        }

      private:
        std::string _id;
    };

    ///VirtualGlobalMatrix
    class VGMatrix : public MatExpr<VGMatrix>
    {
      public:
        ///construct from any MatExpr
        template<typename T_>
        VGMatrix(MatExpr<T_> const& other)
        {
          _id = other.cast().get_id();
        }

        ///construct from id
        VGMatrix(const std::string id) :
          _id(id)
        {
        }

        std::string get_id()
        {
          return _id;
        }

        const std::string get_id() const
        {
          return _id;
        }

      private:
        std::string _id;
    };

    class Matrix : public MatExpr<Matrix>
    {
      public:
        ///construct from any MatExpr
        template<typename T_>
        Matrix(MatExpr<T_> const& other)
        {
          _id = "[" + other.cast().get_id() +"]_chunk";
        }

        std::string get_id()
        {
          return _id;
        }

        const std::string get_id() const
        {
          return _id;
        }

      private:
        std::string _id;
    };

    ///VirtualGlobalScalar
    class VGScalar : public ScalExpr<VGScalar>
    {
      public:
        ///construct from any ScalExpr
        template<typename T_>
        VGScalar(ScalExpr<T_> const& other)
        {
          _id = other.cast().get_id();
        }

        ///construct from id
        VGScalar(const std::string id) :
          _id(id)
        {
        }

        std::string get_id()
        {
          return _id;
        }

        const std::string get_id() const
        {
          return _id;
        }

      private:
        std::string _id;
    };

    class Scalar : public ScalExpr<Scalar>
    {
      public:
        ///construct from any ScalExpr
        template<typename T_>
        Scalar(ScalExpr<T_> const& other)
        {
          _id = "[" + other.cast().get_id() +"]_chunk";
        }

        std::string get_id()
        {
          return _id;
        }

        const std::string get_id() const
        {
          return _id;
        }

      private:
        std::string _id;
    };

    class Bool : public BoolExpr<Bool>
    {
      public:
        ///construct from any BoolExpr
        template<typename T_>
        Bool(BoolExpr<T_> const& other)
        {
          _id = "[" + other.cast().get_id() +"]_chunk";
        }

        std::string get_id()
        {
          return _id;
        }

        const std::string get_id() const
        {
          return _id;
        }

      private:
        std::string _id;
    };

    ///VirtualGlobalBool
    class VGBool : public BoolExpr<VGBool>
    {
      public:
        ///construct from any ScalExpr
        template<typename T_>
        VGBool(BoolExpr<T_> const& other)
        {
          _id = other.cast().get_id();
        }

        ///construct from id
        VGBool(const std::string id) :
          _id(id)
        {
        }

        std::string get_id()
        {
          return _id;
        }

        const std::string get_id() const
        {
          return _id;
        }

      private:
        std::string _id;
    };

    ///operation expressions
    template<typename T1_, typename T2_>
    class VecSum : public VecExpr<VecSum<T1_, T2_> >
    {
      public:
        VecSum(const VecExpr<T1_>& l, const VecExpr<T2_>& r) :
          _left(l.cast()),
          _right(r.cast())
        {
        }

        std::string get_id()
        {
          return _left.cast().get_id() + " + " + _right.cast().get_id();
        }

        const std::string get_id() const
        {
          return _left.cast().get_id() + " + " + _right.cast().get_id();
        }

      private:
          const T1_& _left;
          const T2_& _right;
    };

    template<typename T1_>
    class NormVec : public ScalExpr<NormVec<T1_> >
    {
      public:
        NormVec(const VecExpr<T1_>& l) :
          _left(l.cast())
        {
        }

        std::string get_id()
        {
          return "NORMVEC(" + _left.cast().get_id() + ")";
        }

        const std::string get_id() const
        {
          return "NORMVEC(" + _left.cast().get_id() + ")";
        }

      private:
          const T1_& _left;
    };

    template<typename T1_>
    class VecSynch : public VecExpr<VecSynch<T1_> >
    {
      public:
        VecSynch(const VecExpr<T1_>& l) :
          _left(l.cast())
        {
        }

        std::string get_id()
        {
          return "SYNCHVEC(" + _left.cast().get_id() + ")";
        }

        const std::string get_id() const
        {
          return "SYNCHVEC(" + _left.cast().get_id() + ")";
        }

      private:
          const T1_& _left;
    };

    template<typename T_>
    class VecPreconApply : public VecExpr<VecPreconApply<T_> >
    {
      public:
        VecPreconApply(const VecExpr<T_>& l) :
          _left(l.cast())
        {
        }

        std::string get_id()
        {
          return "PRECONAPPLY(" + _left.cast().get_id() + ")";
        }

        const std::string get_id() const
        {
          return "PRECONAPPLY(" + _left.cast().get_id() + ")";
        }

      private:
          const T_& _left;
    };

    template<typename T1_, typename T2_>
    class ScalSum : public ScalExpr<ScalSum<T1_, T2_> >
    {
      public:
        ScalSum(const ScalExpr<T1_>& l, const ScalExpr<T2_>& r) :
          _left(l.cast()),
          _right(r.cast())
        {
        }

        std::string get_id()
        {
          return _left.cast().get_id() + " + " + _right.cast().get_id();
        }

        const std::string get_id() const
        {
          return _left.cast().get_id() + " + " + _right.cast().get_id();
        }

      private:
          const T1_& _left;
          const T2_& _right;
    };

    template<typename T1_, typename T2_>
    class BoolLess : public BoolExpr<BoolLess<T1_, T2_> >
    {
      public:
        BoolLess(const BoolExpr<T1_>& l, const BoolExpr<T2_>& r) :
          _left(l.cast()),
          _right(r.cast())
        {
        }

        std::string get_id()
        {
          return _left.cast().get_id() + " < " + _right.cast().get_id();
        }

        const std::string get_id() const
        {
          return _left.cast().get_id() + " < " + _right.cast().get_id();
        }

      private:
          const T1_& _left;
          const T2_& _right;
    };

    template<typename BT_, typename AT_, typename XT_>
    class Defect : public VecExpr<Defect<BT_, AT_, XT_> >
    {
      public:
        Defect(const VecExpr<BT_>& b, const MatExpr<AT_>& A, const VecExpr<XT_>& x) :
          _b(b.cast()),
          _A(A.cast()),
          _x(x.cast())
        {
        }

        std::string get_id()
        {
          return "DEFECT(" + _b.cast().get_id() + ", " + _A.cast().get_id() + ", " + _x.cast().get_id() + ")";
        }

        const std::string get_id() const
        {
          return "DEFECT(" + _b.cast().get_id() + ", " +_A.cast().get_id() + ", " + _x.cast().get_id() + ")";
        }

      private:
          const BT_& _b;
          const AT_& _A;
          const XT_& _x;
    };

    template<typename T1_, typename T2_>
    class VecIterate : public VecExpr<VecIterate<T1_, T2_> >
    {
      public:
        VecIterate(const VecExpr<T1_>& l, const BoolExpr<T2_>& u) :
          _left(l.cast()),
          _u(u)
        {
        }

        std::string get_id()
        {
          return "ITERATE(" + _left.cast().get_id() + " UNTIL " + _u.cast().get_id() + ")";
        }

        const std::string get_id() const
        {
          return "ITERATE(" + _left.cast().get_id() + " UNTIL " + _u.cast().get_id() + ")";
        }

      private:
          const T1_& _left;
          const T2_& _u;
    };

    template<typename T1_, typename T2_>
    class MatVecProd : public VecExpr<MatVecProd<T1_, T2_> >
    {
      public:
        MatVecProd(const MatExpr<T1_>& l, const VecExpr<T2_>& r) :
          _left(l.cast()),
          _right(r.cast())
        {
        }

        std::string get_id()
        {
          return _left.cast().get_id() + " * " + _right.cast().get_id();
        }

        const std::string get_id() const
        {
          return _left.cast().get_id() + " * " + _right.cast().get_id();
        }

      private:
          const T1_& _left;
          const T2_& _right;
    };

    ///operator overloads and operation shortcuts
    template<typename T1_, typename T2_>
    VecSum<T1_, T2_> operator+(const VecExpr<T1_>& l, const VecExpr<T2_>& r)
    {
      return VecSum<T1_, T2_>(l, r);
    }

    template<typename T1_, typename T2_>
    ScalSum<T1_, T2_> operator+(const ScalExpr<T1_>& l, const ScalExpr<T2_>& r)
    {
      return ScalSum<T1_, T2_>(l, r);
    }

    template<typename T1_>
    NormVec<T1_> normvec_expr(const VecExpr<T1_>& l)
    {
      return NormVec<T1_>(l);
    }

    template<typename T1_>
    VecSynch<T1_> vec_synch_expr(const VecExpr<T1_>& l)
    {
      return VecSynch<T1_>(l);
    }

    template<typename T1_, typename T2_>
    BoolLess<T1_, T2_> operator<(const BoolExpr<T1_>& l, const BoolExpr<T2_>& r)
    {
      return BoolLess<T1_, T2_>(l, r);
    }

    template<typename T2_, typename T3_, typename T4_>
    Defect<T2_, T3_, T4_> defect_expr(const VecExpr<T2_>& b, const MatExpr<T3_>& A, const VecExpr<T4_>& x)
    {
      return Defect<T2_, T3_, T4_>(b, A, x);
    }

    template<typename T1_, typename T2_>
    VecIterate<T1_, T2_> vec_iterate_expr(const VecExpr<T1_>& l, const BoolExpr<T2_>& u)
    {
      return VecIterate<T1_, T2_>(l, u);
    }

    template<typename T_>
    VecPreconApply<T_> vec_preconapply_expr(const VecExpr<T_>& l)
    {
      return VecPreconApply<T_>(l);
    }

    template<typename T1_, typename T2_>
    MatVecProd<T1_, T2_> operator*(const MatExpr<T1_>& l, const VecExpr<T2_>& r)
    {
      return MatVecProd<T1_, T2_>(l, r);
    }
  }
}

#endif
