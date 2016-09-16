#pragma once
#ifndef KERNEL_LAFEM_TUPLE_FILTER_HPP
#define KERNEL_LAFEM_TUPLE_FILTER_HPP 1

// includes, FEAT
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/meta_element.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief TupleVector meta-filter class template
     *
     * This class template implements a composition of sub-filters of arbitrary classes,
     * which can be applied onto TupleVector meta-container objects.
     *
     * \note If you are looking for a way to define a filter which is a composition of
     * filters successively applied onto the same vector, see the FilterChain class template.
     *
     * \tparam First_, Rest_...
     * The filters to be composed.
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
      static constexpr int num_blocks = TupleFilter<Rest_...>::num_blocks + 1;

      /// sub-filter mem-type
      typedef typename First_::MemType MemType;
      /// sub-filter data-type
      typedef typename First_::DataType DataType;
      /// sub-filter index-type
      typedef typename First_::IndexType IndexType;

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value,
                    "sub-filters have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value,
                    "sub-filters have different data-types");
      static_assert(std::is_same<IndexType, typename RestClass::IndexType>::value,
                    "sub-filters have different index-types");

      /// corresponding vector
      typedef TupleVector<typename First_::VectorType, typename Rest_::VectorType...> VectorType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using FilterType = TupleFilter<
        typename First_::template FilterType<Mem2_, DT2_, IT2_>,
        typename Rest_::template FilterType<Mem2_, DT2_, IT2_>...>;

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

      /// \brief Creates a clone of itself
      TupleFilter clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return TupleFilter(_first.clone(clone_mode), _rest.clone(clone_mode));
      }

      /// \brief Clones data from another TupleFilter
      void clone(const TupleFilter & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _first.clone(other._first.get_filter_vector(), clone_mode);
        _rest.clone(other._rest, clone_mode);
      }

      template<typename... SubFilter2_>
      void convert(const TupleFilter<SubFilter2_...>& other)
      {
        _first.convert(other._first);
        _rest.convert(other._rest);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
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

      template<int i_>
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        static_assert((0 <= i_) && (i_ < num_blocks), "invalid sub-filter index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<int i_>
      typename TupleElement<i_, First_, Rest_...>::Type const& at() const
      {
        static_assert((0 <= i_) && (i_ < num_blocks), "invalid sub-filter index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc UnitFilter::filter_rhs() */
      void filter_rhs(VectorType& vector) const
      {
        first().filter_rhs(vector.first());
        rest().filter_rhs(vector.rest());
      }

      /** \copydoc UnitFilter::filter_sol() */
      void filter_sol(VectorType& vector) const
      {
        first().filter_sol(vector.first());
        rest().filter_sol(vector.rest());
      }

      /** \copydoc UnitFilter::filter_def() */
      void filter_def(VectorType& vector) const
      {
        first().filter_def(vector.first());
        rest().filter_def(vector.rest());
      }

      /** \copydoc UnitFilter::filter_cor() */
      void filter_cor(VectorType& vector) const
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
      static constexpr int num_blocks = 1;

      /// sub-filter mem-type
      typedef typename First_::MemType MemType;
      /// sub-filter data-type
      typedef typename First_::DataType DataType;
      /// sub-filter index-type
      typedef typename First_::IndexType IndexType;

      typedef TupleVector<typename First_::VectorType> VectorType;

      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using FilterType = class TupleFilter<typename First_::template FilterType<Mem2_, DT2_, IT2_> >;

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

      /// \brief Creates a clone of itself
      TupleFilter clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return TupleFilter(_first.clone(clone_mode));
      }

      /// \brief Clones data from another TupleFilter
      void clone(const TupleFilter & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _first.clone(other._first.get_filter_vector(), clone_mode);
      }

      template<typename SubFilter2_>
      void convert(const TupleFilter<SubFilter2_>& other)
      {
        _first.convert(other._first);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes();
      }

      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      template<int i_>
      typename TupleElement<i_, First_>::Type& at()
      {
        static_assert(i_ == 0, "invalid sub-filter index");
        return first();
      }

      template<int i_>
      typename TupleElement<i_, First_>::Type const& at() const
      {
        static_assert(i_ == 0, "invalid sub-filter index");
        return first();
      }

      /** \copydoc UnitFilter::filter_rhs() */
      void filter_rhs(VectorType& vector) const
      {
        first().filter_rhs(vector.first());
      }

      /** \copydoc UnitFilter::filter_sol() */
      void filter_sol(VectorType& vector) const
      {
        first().filter_sol(vector.first());
      }

      /** \copydoc UnitFilter::filter_def() */
      void filter_def(VectorType& vector) const
      {
        first().filter_def(vector.first());
      }

      /** \copydoc UnitFilter::filter_cor() */
      void filter_cor(VectorType& vector) const
      {
        first().filter_cor(vector.first());
      }
    }; // class TupleFilter
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_TUPLE_FILTER_HPP
