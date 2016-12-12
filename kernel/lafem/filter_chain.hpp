#pragma once
#ifndef KERNEL_LAFEM_FILTER_CHAIN_HPP
#define KERNEL_LAFEM_FILTER_CHAIN_HPP 1

// includes, FEAT
#include <kernel/lafem/meta_element.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Filter Chainclass template
     *
     * This class template implements a chain of filters which are successively applied
     * onto a vector.
     *
     * \tparam First_, Rest_...
     * The filters to be chained.
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class FilterChain
    {
      template<typename,typename...>
      friend class FilterChain;

      typedef FilterChain<Rest_...> RestClass;

    public:
      /// number of vector blocks
      static constexpr int num_blocks = FilterChain<Rest_...>::num_blocks + 1;

      /// sub-filter mem-type
      typedef typename First_::MemType MemType;
      /// sub-filter data-type
      typedef typename First_::DataType DataType;
      /// sub-filter index-type
      typedef typename First_::IndexType IndexType;

      typedef typename First_::VectorType VectorType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using FilterType = FilterChain<
        typename First_::template FilterType<Mem2_, DT2_, IT2_>,
        typename Rest_::template FilterType<Mem2_, DT2_, IT2_>...>;

      /// this typedef lets you create a matrix container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using FilterTypeByMDI = FilterType<Mem2_, DT2_, IT2_>;

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value,
                    "sub-filters have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value,
                    "sub-filters have different data-types");
      static_assert(std::is_same<IndexType, typename RestClass::IndexType>::value,
                    "sub-filters have different index-types");

    protected:
      /// the first sub-filter
      First_ _first;
      /// all remaining sub-filters
      RestClass _rest;

      /// data-emplacement ctor; this one is protected for a reason
      explicit FilterChain(First_&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

    public:
      /// default ctor
      FilterChain()
      {
      }

      /// sub-filter emplacement ctor
      explicit FilterChain(First_&& the_first, Rest_&&... the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest...))
      {
      }

      /// move-ctor
      FilterChain(FilterChain&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      FilterChain& operator=(FilterChain&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      /// \brief Creates a clone of itself
      FilterChain clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return FilterChain(_first.clone(clone_mode), _rest.clone(clone_mode));
      }

      /// \brief Clones data from another FilterChain
      void clone(const FilterChain & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _first.clone(other._first.get_filter_vector(), clone_mode);
        _rest.clone(other._rest, clone_mode);
      }

      template<typename... SubFilter2_>
      void convert(const FilterChain<SubFilter2_...>& other)
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
      template<typename Vector_>
      void filter_rhs(Vector_& vector) const
      {
        first().filter_rhs(vector);
        rest().filter_rhs(vector);
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Vector_>
      void filter_sol(Vector_& vector) const
      {
        first().filter_sol(vector);
        rest().filter_sol(vector);
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Vector_>
      void filter_def(Vector_& vector) const
      {
        first().filter_def(vector);
        rest().filter_def(vector);
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Vector_>
      void filter_cor(Vector_& vector) const
      {
        first().filter_cor(vector);
        rest().filter_cor(vector);
      }
    }; // class FilterChain

    /// \cond internal
    template<typename First_>
    class FilterChain<First_>
    {
      template<typename,typename...>
      friend class FilterChain;

    public:
      static constexpr int num_blocks = 1;

      /// sub-filter mem-type
      typedef typename First_::MemType MemType;
      /// sub-filter data-type
      typedef typename First_::DataType DataType;
      /// sub-filter index-type
      typedef typename First_::IndexType IndexType;

      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using FilterType = class FilterChain<typename First_::template FilterType<Mem2_, DT2_, IT2_> >;

    protected:
      /// the first sub-filter
      First_ _first;

    public:
      /// default ctor
      FilterChain()
      {
      }

      /// sub-filter emplacement ctor
      explicit FilterChain(First_&& the_first) :
        _first(std::move(the_first))
      {
      }

      /// move-ctor
      FilterChain(FilterChain&& other) :
        _first(std::move(other._first))
      {
      }

      /// move-assign operator
      FilterChain& operator=(FilterChain&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      /// \brief Creates a clone of itself
      FilterChain clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return FilterChain(_first.clone(clone_mode));
      }

      /// \brief Clones data from another FilterChain
      void clone(const FilterChain & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _first.clone(other._first.get_filter_vector(), clone_mode);
      }

      template<typename SubFilter2_>
      void convert(const FilterChain<SubFilter2_>& other)
      {
        _first.convert(other._first);
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
      template<typename Vector_>
      void filter_rhs(Vector_& vector) const
      {
        first().filter_rhs(vector);
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Vector_>
      void filter_sol(Vector_& vector) const
      {
        first().filter_sol(vector);
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Vector_>
      void filter_def(Vector_& vector) const
      {
        first().filter_def(vector);
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Vector_>
      void filter_cor(Vector_& vector) const
      {
        first().filter_cor(vector);
      }
    }; // class FilterChain
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_FILTER_CHAIN_HPP
