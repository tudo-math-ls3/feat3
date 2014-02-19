#pragma once
#ifndef KERNEL_LAFEM_TUPLE_FILTER_HPP
#define KERNEL_LAFEM_TUPLE_FILTER_HPP 1

// includes, FEAST
#include <kernel/lafem/tuple_vector.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief TupleVector meta-filter class template
     *
     * This class template implements a composition of sub-filters of arbitrary classes,
     * which can be applied onto TupleVector meta-container objects.
     *
     * \tparam First_, Rest_...
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class TupleFilter
    {
      template<typename,typename...>
      friend class TupleFilter;

      typedef TupleFilter<Rest_...> RestClass;

    public:
      /// sub-filter mem-type
      typedef typename First_::MemType MemType;
      /// sub-filter data-type
      typedef typename First_::DataType DataType;

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value,
        "sub-filters have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value,
        "sub-filters have different data-types");

    protected:
      /// the first sub-filter
      First_ _first;
      /// all remaining sub-filters
      RestClass _rest;

      /// data-emplacement ctor; this one is protected for a reason
      explicit TupleFilter(First_&& first, RestClass&& rest) :
        _first(std::move(first)),
        _rest(std::move(rest))
      {
      }

    public:
      /// default ctor
      TupleFilter()
      {
      }

      /// sub-filter emplacement ctor
      explicit TupleFilter(First_&& first, Rest_&&... rest) :
        _first(std::move(first)),
        _rest(std::move(rest...))
      {
      }

      /// move-ctor
      TupleFilter(TupleFilter&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      TupleFilter& operator=(TupleFilter&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      TupleFilter clone() const
      {
        return TupleFilter(_first.clone(), _rest.clone());
      }

      /// \cond internal
      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      RestClass& rest()
      {
        return _rest;
      }

      const RestClass& rest() const
      {
        return _rest;
      }
      /// \endcond

      /** \copydoc UnitFilter::filter_rhs() */
      template<typename Algo_, typename Ty_, typename... Tv_>
      void filter_rhs(TupleVector<Ty_, Tv_...>& vector) const
      {
        first().template filter_rhs<Algo_>(vector.first());
        rest().template filter_rhs<Algo_>(vector.rest());
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Algo_, typename Ty_, typename... Tv_>
      void filter_sol(TupleVector<Ty_, Tv_...>& vector) const
      {
        first().template filter_sol<Algo_>(vector.first());
        rest().template filter_sol<Algo_>(vector.rest());
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Algo_, typename Ty_, typename... Tv_>
      void filter_def(TupleVector<Ty_, Tv_...>& vector) const
      {
        first().template filter_def<Algo_>(vector.first());
        rest().template filter_def<Algo_>(vector.rest());
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Algo_, typename Ty_, typename... Tv_>
      void filter_cor(TupleVector<Ty_, Tv_...>& vector) const
      {
        first().template filter_cor<Algo_>(vector.first());
        rest().template filter_cor<Algo_>(vector.rest());
      }
    }; // class TupleFilter

    /// \cond internal
    template<typename First_>
    class TupleFilter<First_>
    {
      template<typename,typename...>
      friend class TupleFilter;

    public:
      /// sub-filter mem-type
      typedef typename First_::MemType MemType;
      /// sub-filter data-type
      typedef typename First_::DataType DataType;

    protected:
      /// the first sub-filter
      First_ _first;

    public:
      /// default ctor
      TupleFilter()
      {
      }

      /// sub-filter emplacement ctor
      explicit TupleFilter(First_&& first) :
        _first(std::move(first))
      {
      }

      /// move-ctor
      TupleFilter(TupleFilter&& other) :
        _first(std::move(other._first))
      {
      }

      /// move-assign operator
      TupleFilter& operator=(TupleFilter&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      TupleFilter clone() const
      {
        return TupleFilter(_first.clone());
      }

      /// \cond internal
      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }
      /// \endcond

#ifdef FEAST_COMPILER_MICROSOFT
      /** \copydoc UnitFilter::filter_rhs() */
      template<typename Algo_, typename... Tv_>
      void filter_rhs(TupleVector<Tv_...>& vector) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        first().template filter_rhs<Algo_>(vector.first());
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Algo_, typename... Tv_>
      void filter_sol(TupleVector<Tv_...>& vector) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        first().template filter_sol<Algo_>(vector.first());
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Algo_, typename... Tv_>
      void filter_def(TupleVector<Tv_...>& vector) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        first().template filter_def<Algo_>(vector.first());
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Algo_, typename... Tv_>
      void filter_cor(TupleVector<Tv_...>& vector) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        first().template filter_cor<Algo_>(vector.first());
      }
#else // any other compiler
      /** \copydoc UnitFilter::filter_rhs() */
      template<typename Algo_, typename Tv_>
      void filter_rhs(TupleVector<Tv_>& vector) const
      {
        first().template filter_rhs<Algo_>(vector.first());
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Algo_, typename Tv_>
      void filter_sol(TupleVector<Tv_>& vector) const
      {
        first().template filter_sol<Algo_>(vector.first());
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Algo_, typename Tv_>
      void filter_def(TupleVector<Tv_>& vector) const
      {
        first().template filter_def<Algo_>(vector.first());
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Algo_, typename Tv_>
      void filter_cor(TupleVector<Tv_>& vector) const
      {
        first().template filter_cor<Algo_>(vector.first());
      }
#endif // FEAST_COMPILER_MICROSOFT
    }; // class TupleFilter
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_FILTER_HPP
