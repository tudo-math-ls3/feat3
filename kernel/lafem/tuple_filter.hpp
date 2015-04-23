#pragma once
#ifndef KERNEL_LAFEM_TUPLE_FILTER_HPP
#define KERNEL_LAFEM_TUPLE_FILTER_HPP 1

// includes, FEAST
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/meta_element.hpp>

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
      /// number of vector blocks
      static constexpr Index num_blocks = TupleFilter<Rest_...>::num_blocks + 1;

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
      explicit TupleFilter(First_&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

    public:
      /// default ctor
      TupleFilter()
      {
      }

      /// sub-filter emplacement ctor
      explicit TupleFilter(First_&& the_first, Rest_&&... the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest...))
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

      template<typename... SubFilter2_>
      void convert(const TupleFilter<SubFilter2_...>& other)
      {
        _first.convert(other._first);
        _rest.convert(other._rest);
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

      template<Index i_>
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        static_assert(i_ < Index(num_blocks), "invalid sub-filter index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<Index i_>
      typename TupleElement<i_, First_, Rest_...>::Type const& at() const
      {
        static_assert(i_ < Index(num_blocks), "invalid sub-filter index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc UnitFilter::filter_rhs() */
      template<typename Ty_, typename... Tv_>
      void filter_rhs(TupleVector<Ty_, Tv_...>& vector) const
      {
        first().filter_rhs(vector.first());
        rest().filter_rhs(vector.rest());
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Ty_, typename... Tv_>
      void filter_sol(TupleVector<Ty_, Tv_...>& vector) const
      {
        first().filter_sol(vector.first());
        rest().filter_sol(vector.rest());
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Ty_, typename... Tv_>
      void filter_def(TupleVector<Ty_, Tv_...>& vector) const
      {
        first().filter_def(vector.first());
        rest().filter_def(vector.rest());
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Ty_, typename... Tv_>
      void filter_cor(TupleVector<Ty_, Tv_...>& vector) const
      {
        first().filter_cor(vector.first());
        rest().filter_cor(vector.rest());
      }
    }; // class TupleFilter

    /// \cond internal
    template<typename First_>
    class TupleFilter<First_>
    {
      template<typename,typename...>
      friend class TupleFilter;

    public:
      static constexpr Index num_blocks = 1;

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
      explicit TupleFilter(First_&& the_first) :
        _first(std::move(the_first))
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

#ifdef FEAST_COMPILER_MICROSOFT
      template<typename... SubFilter2_>
      void convert(const TupleFilter<SubFilter2_...>& other)
      {
        static_assert(sizeof...(SubFilter2_) == std::size_t(1), "invalid TupleFilter size");
        _first.convert(other._first);
      }
#else
      template<typename SubFilter2_>
      void convert(const TupleFilter<SubFilter2_>& other)
      {
        _first.convert(other._first);
      }
#endif

      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      template<Index i_>
      typename TupleElement<i_, First_>::Type& at()
      {
        static_assert(i_ == 0, "invalid sub-filter index");
        return first();
      }

      template<Index i_>
      typename TupleElement<i_, First_>::Type const& at() const
      {
        static_assert(i_ == 0, "invalid sub-filter index");
        return first();
      }

#ifdef FEAST_COMPILER_MICROSOFT
      /** \copydoc UnitFilter::filter_rhs() */
      template<typename... Tv_>
      void filter_rhs(TupleVector<Tv_...>& vector) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        first().filter_rhs(vector.first());
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename... Tv_>
      void filter_sol(TupleVector<Tv_...>& vector) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        first().filter_sol(vector.first());
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename... Tv_>
      void filter_def(TupleVector<Tv_...>& vector) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        first().filter_def(vector.first());
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename... Tv_>
      void filter_cor(TupleVector<Tv_...>& vector) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        first().filter_cor(vector.first());
      }
#else // any other compiler
      /** \copydoc UnitFilter::filter_rhs() */
      template<typename Tv_>
      void filter_rhs(TupleVector<Tv_>& vector) const
      {
        first().filter_rhs(vector.first());
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Tv_>
      void filter_sol(TupleVector<Tv_>& vector) const
      {
        first().filter_sol(vector.first());
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Tv_>
      void filter_def(TupleVector<Tv_>& vector) const
      {
        first().filter_def(vector.first());
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Tv_>
      void filter_cor(TupleVector<Tv_>& vector) const
      {
        first().filter_cor(vector.first());
      }
#endif // FEAST_COMPILER_MICROSOFT
    }; // class TupleFilter
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_FILTER_HPP
