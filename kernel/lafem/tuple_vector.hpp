#pragma once
#ifndef KERNEL_LAFEM_TUPLE_VECTOR_HPP
#define KERNEL_LAFEM_TUPLE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#include <iostream>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Variadic base class for TupleVector
     *
     * \attention Do not use this class template directly; use the TupleVector class template instead!
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class TupleVectorBase
    {
      template<typename,typename...>
      friend class TupleVectorBase;

      typedef TupleVectorBase<Rest_...> RestClass;

    public:
      /// dummy enum
      enum
      {
        /// number of vector blocks
        num_blocks = TupleVectorBase<Rest_...>::num_blocks + 1
      };

      typedef typename First_::MemType MemType;
      typedef typename First_::DataType DataType;

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(SAME_TYPE(MemType, typename RestClass::MemType), "sub-vectors have different mem-types");
      static_assert(SAME_TYPE(DataType, typename RestClass::DataType), "sub-vectors have different data-types");

    protected:
      First_ _first;
      RestClass _rest;

      explicit TupleVectorBase(First_&& first, TupleVectorBase<Rest_...>&& rest) :
        _first(std::move(first)),
        _rest(std::move(rest))
      {
      }

      explicit TupleVectorBase(First_&& first, Rest_&&... rest) :
        _first(std::move(first)),
        _rest(std::move(rest...))
      {
      }

      TupleVectorBase()
      {
      }

      TupleVectorBase(TupleVectorBase&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      TupleVectorBase& operator=(TupleVectorBase&& other)
      {
        _first = std::move(other._first);
        _rest = std::move(other._rest);
      }

      /// deleted copy-ctor
      TupleVectorBase(const TupleVectorBase&) = delete;
      /// deleted copy-assign operator
      TupleVectorBase& operator=(const TupleVectorBase&) = delete;

      TupleVectorBase clone() const
      {
        return TupleVectorBase(_first.clone(), _rest.clone());
      }

      static String sub_name_list()
      {
        return First_::type_name() + "," + TupleVectorBase<Rest_...>::sub_name_list();
      }

    public:
      /// \cond internal
      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      TupleVectorBase<Rest_...>& rest()
      {
        return _rest;
      }

      const TupleVectorBase<Rest_...>& rest() const
      {
        return _rest;
      }
      /// \endcond

      Index size() const
      {
        return _first.size() + _rest.size();
      }

      void clear(DataType value = DataType(0))
      {
        first().clear(value);
        rest().clear(value);
      }

      //template<typename First2_, typename... Rest2_>
      void copy(const TupleVectorBase/*<First2_, Rest2_...>*/& x)
      {
        first().copy(x.first());
        rest().copy(x.rest());
      }

      template<typename Algo_>
      void axpy(const TupleVectorBase& x, const TupleVectorBase& y, DataType alpha = DataType(1))
      {
        first().template axpy<Algo_>(x.first(), y.first(), alpha);
        rest().template axpy<Algo_>(x.rest(), y.rest(), alpha);
      }

      template <typename Algo_>
      void component_product(const TupleVectorBase & x, const TupleVectorBase & y)
      {
        first().template component_product<Algo_>(x.first(), y.first());
        rest().template component_product<Algo_>(x.rest(), y.rest());
      }

      template <typename Algo_>
      void component_product(const TupleVectorBase & x, const TupleVectorBase & y, const TupleVectorBase& z)
      {
        first().template component_product<Algo_>(x.first(), y.first(), z.first());
        rest().template component_product<Algo_>(x.rest(), y.rest(), z.rest());
      }

      template<typename Algo_>
      void scale(const TupleVectorBase& x, DataType alpha)
      {
        first().template scale<Algo_>(x.first(), alpha);
        rest().template scale<Algo_>(x.rest(), alpha);
      }

      template<typename Algo_>
      DataType dot(const TupleVectorBase& x) const
      {
        return first().template dot<Algo_>(x.first()) + rest().template dot<Algo_>(x.rest());
      }

      template<typename Algo_>
      DataType norm2sqr() const
      {
        return first().template norm2sqr<Algo_>() + rest().template norm2sqr<Algo_>();
      }

      template<typename Algo_>
      DataType norm2() const
      {
        return Math::sqrt(norm2sqr<Algo_>());
      }
    };

    /// \cond internal
    template<typename First_>
    class TupleVectorBase<First_>
    {
      template<typename,typename...>
      friend class TupleVectorBase;

    public:
      enum
      {
        num_blocks = 1
      };
      typedef typename First_::MemType MemType;
      typedef typename First_::DataType DataType;

    protected:
      First_ _first;

      explicit TupleVectorBase(First_&& first) :
        _first(std::move(first))
      {
      }

      TupleVectorBase()
      {
      }

      TupleVectorBase(TupleVectorBase&& other) :
        _first(std::move(other._first))
      {
      }

      TupleVectorBase& operator=(TupleVectorBase&& other)
      {
        _first = std::move(other._first);
      }

      /// deleted copy-ctor
      TupleVectorBase(const TupleVectorBase&) = delete;
      /// deleted copy-assign operator
      TupleVectorBase& operator=(const TupleVectorBase&) = delete;

      TupleVectorBase clone() const
      {
        return TupleVectorBase(_first.clone());
      }

      static String sub_name_list()
      {
        return First_::type_name();
      }

    public:
      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      Index size() const
      {
        return _first.size();
      }

      void clear(DataType value = DataType(0))
      {
        _first.clear(value);
      }

      //template<typename First2_>
      void copy(const TupleVectorBase/*<First2_>*/& x)
      {
        first().copy(x.first());
      }

      template<typename Algo_>
      void axpy(const TupleVectorBase& x, const TupleVectorBase& y, DataType alpha = DataType(1))
      {
        first().template axpy<Algo_>(x.first(), y.first(), alpha);
      }

      template <typename Algo_>
      void component_product(const TupleVectorBase & x, const TupleVectorBase & y)
      {
        first().template component_product<Algo_>(x.first(), y.first());
      }

      template <typename Algo_>
      void component_product(const TupleVectorBase & x, const TupleVectorBase & y, const TupleVectorBase& z)
      {
        first().template component_product<Algo_>(x.first(), y.first(), z.first());
      }

      template<typename Algo_>
      void scale(const TupleVectorBase& x, DataType alpha)
      {
        first().template scale<Algo_>(x.first(), alpha);
      }

      template<typename Algo_>
      DataType dot(const TupleVectorBase& x) const
      {
        return first().template dot<Algo_>(x.first());
      }

      template<typename Algo_>
      DataType norm2sqr() const
      {
        return first().template norm2sqr<Algo_>();
      }

      template<typename Algo_>
      DataType norm2() const
      {
        return first().template norm2<Algo_>();
      }
    };
    /// \endcond

    /**
     * \brief TupleVector element helper class template
     *
     * This class template is a helper to realise the "at" member function of the TupleVector class template.
     *
     * \author Peter Zajac
     */
    template<
      Index i_,
      typename First_,
      typename... Rest_>
    struct TupleVectorElement
    {
      typedef typename TupleVectorElement<i_-1, Rest_...>::Type Type;

      static Type& get(TupleVectorBase<First_, Rest_...>& meta)
      {
        return TupleVectorElement<i_-1, Rest_...>::get(meta.rest());
      }

      static const Type& get(const TupleVectorBase<First_, Rest_...>& meta)
      {
        return TupleVectorElement<i_-1, Rest_...>::get(meta.rest());
      }
    };

    /// \cond internal
    template<typename First_, typename... Rest_>
    struct TupleVectorElement<0, First_, Rest_...>
    {
      typedef First_ Type;

      static Type& get(TupleVectorBase<First_, Rest_...>& meta)
      {
        return meta.first();
      }
      static const Type& get(const TupleVectorBase<First_, Rest_...>& meta)
      {
        return meta.first();
      }
    };
    /// \endcond

    /**
     * \brief Variadic TupleVector class template
     *
     * This class template implements a composition of sub-vectors of arbitrary classes.
     *
     * \tparam First_, ...Rest_
     * A sequence of (meta) vector classes which are to be composed.
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class TupleVector :
      public TupleVectorBase<First_, Rest_...>
    {
      /// base-class typedef
      typedef TupleVectorBase<First_, Rest_...> BaseClass;

    public:
      /// dummy enum
      enum
      {
        // number of blocks in this vector
        num_blocks = BaseClass::num_blocks
      };

    protected:
      /// data-move CTOR; this one is protected for a reason
      explicit TupleVector(TupleVectorBase<First_, Rest_...>&& other) :
        BaseClass(std::move(other))
      {
      }

    public:
      /// default CTOR
      TupleVector()
      {
      }

      /// Sub-Vector emplacement constructor
      explicit TupleVector(First_&& first, Rest_&&... rest) :
        BaseClass(std::move(first), std::move(rest...))
      {
      }

      /// move CTOR
      TupleVector(TupleVector&& other) :
        BaseClass(static_cast<BaseClass&&>(other))
      {
      }

      /// move-assign operator
      TupleVector& operator=(TupleVector&& other)
      {
        static_cast<BaseClass&>(*this).operator=(other);
        return *this;
      }

      /// deleted copy-ctor
      TupleVector(const TupleVector&) = delete;
      /// deleted copy-assign operator
      TupleVector& operator=(const TupleVector&) = delete;

      /// virtual destructor
      virtual ~TupleVector()
      {
      }

      /**
       * \brief Creates and returns a deep copy of this vector.
       */
      TupleVector clone() const
      {
        return TupleVector(BaseClass::clone());
      }

      /**
       * \brief Returns a sub-vector block.
       *
       * \tparam i_
       * The index of the sub-vector block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-vector at position \p i_.
       */
      template<Index i_>
      typename TupleVectorElement<i_, First_, Rest_...>::Type& at()
      {
        static_assert(i_ < Index(num_blocks), "invalid sub-vector index");
        return TupleVectorElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<Index i_>
      typename TupleVectorElement<i_, First_, Rest_...>::Type const& at() const
      {
        static_assert(i_ < Index(num_blocks), "invalid sub-vector index");
        return TupleVectorElement<i_, First_, Rest_...>::get(*this);
      }

      /// Returns the number of blocks in this tuple-vector.
      Index blocks() const
      {
        return Index(num_blocks);
      }

      /// Returns a descriptive string for this container.
      static String type_name()
      {
        return String("TupleVector<") + BaseClass::sub_name_list() + ">";
      }
    }; // class BlockVector<...>

    /// \cond internal
    template <typename First_>
    inline void dump_tuple_vector(std::ostream & os, const TupleVectorBase<First_>& x)
    {
      os << x.first();
    }

    template <typename First_, typename... Rest_>
    inline void dump_tuple_vector(std::ostream & os, const TupleVectorBase<First_, Rest_...>& x)
    {
      os << x.first() << ",";
      dump_tuple_vector<Rest_...>(os, x.rest());
    }

    template <typename First_, typename... Rest_>
    inline std::ostream& operator<< (std::ostream & os, const TupleVector<First_, Rest_...>& x)
    {
      os << "[";
      dump_tuple_vector(os, x);
      os << "]";
      return os;
    }
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_VECTOR_HPP
