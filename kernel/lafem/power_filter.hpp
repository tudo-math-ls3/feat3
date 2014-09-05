#pragma once
#ifndef KERNEL_LAFEM_POWER_FILTER_HPP
#define KERNEL_LAFEM_POWER_FILTER_HPP 1

#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/meta_element.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief PowerVector meta-filter class template
     *
     * This class template implements a composition of \e n sub-filters of the same type,
     * which can be applied onto PowerVector meta-container objects.
     *
     * \note For a composed meta-filter of different sub-filters, see the TupleFilter class template.
     *
     * \tparam SubFilter_
     * The type of the sub-filter.
     *
     * \tparam count_
     * The number of sub-filter blocks.
     *
     * \author Peter Zajac
     */
    template<
      typename SubFilter_,
      Index count_>
    class PowerFilter
    {
      // declare this class template as a friend for recursive inheritance
      template<typename,Index>
      friend class PowerFilter;

      /// base-class typedef
      typedef PowerFilter<SubFilter_, count_-1> RestClass;

    public:
      /// sub-filter type
      typedef SubFilter_ SubFilterType;
      /// sub-filter mem-type
      typedef typename SubFilter_::MemType MemType;
      /// sub-filter data-type
      typedef typename SubFilter_::DataType DataType;

      /// dummy enum
      enum
      {
        /// number of filter blocks
        num_blocks = count_
      };

    protected:
      /// the first sub-filter
      SubFilterType _first;
      /// the remaining part
      RestClass _rest;

      /// base-class ctor; for internal use only
      explicit PowerFilter(SubFilterType&& the_first, RestClass&& the_rest) :
        _first(the_first),
        _rest(the_rest)
      {
      }

    public:
      /// default ctor
      PowerFilter()
      {
      }

      /// move CTOR
      PowerFilter(PowerFilter&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      PowerFilter& operator=(PowerFilter&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      /// Creates and returns a deep copy of this filter.
      PowerFilter clone() const
      {
        return PowerFilter(first().clone(), rest().clone());
      }

      /// Conversion method
      template<typename SubFilter2_>
      void convert(const PowerFilter<SubFilter2_, count_>& other)
      {
        _first.convert(other._first);
        _rest.convert(other._rest);
      }

      /// \cond internal
      SubFilterType& first()
      {
        return _first;
      }

      const SubFilterType& first() const
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
      SubFilterType& at()
      {
        static_assert(i_ < count_, "invalid sub-filter index");
        return PowerElement<i_, SubFilterType>::get(*this);
      }

      template<Index i_>
      const SubFilterType& at() const
      {
        static_assert(i_ < count_, "invalid sub-filter index");
        return PowerElement<i_, SubFilterType>::get(*this);
      }

      /** \copydoc UnitFilter::filter_rhs() */
      template<typename Algo_, typename SubVector_>
      void filter_rhs(PowerVector<SubVector_, count_>& vector) const
      {
        first().template filter_rhs<Algo_>(vector.first());
        rest().template filter_rhs<Algo_>(vector.rest());
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Algo_, typename SubVector_>
      void filter_sol(PowerVector<SubVector_, count_>& vector) const
      {
        first().template filter_sol<Algo_>(vector.first());
        rest().template filter_sol<Algo_>(vector.rest());
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Algo_, typename SubVector_>
      void filter_def(PowerVector<SubVector_, count_>& vector) const
      {
        first().template filter_def<Algo_>(vector.first());
        rest().template filter_def<Algo_>(vector.rest());
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Algo_, typename SubVector_>
      void filter_cor(PowerVector<SubVector_, count_>& vector) const
      {
        first().template filter_cor<Algo_>(vector.first());
        rest().template filter_cor<Algo_>(vector.rest());
      }
    }; // class PowerFilter<...>

    /// \cond internal
    template<typename SubFilter_>
    class PowerFilter<SubFilter_, 1>
    {
      template<typename,Index>
      friend class PowerFilter;

    public:
      typedef SubFilter_ SubFilterType;

      enum
      {
        num_blocks = 1
      };

    protected:
      SubFilterType _first;

      explicit PowerFilter(SubFilterType&& the_first) :
        _first(the_first)
      {
      }

    public:
      PowerFilter()
      {
      }

      /// move CTOR
      PowerFilter(PowerFilter&& other) :
        _first(std::move(other._first))
      {
      }

      PowerFilter& operator=(PowerFilter&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      PowerFilter clone() const
      {
        return PowerFilter(first().clone());
      }

      template<typename SubFilter2_>
      void convert(const PowerFilter<SubFilter2_, 1>& other)
      {
        _first.convert(other._first);
      }

      SubFilterType& first()
      {
        return _first;
      }

      const SubFilterType& first() const
      {
        return _first;
      }

      template<Index i_>
      SubFilterType& at()
      {
        static_assert(i_ != 0, "invalid sub-filter index");
        return _first;
      }

      template<Index i_>
      const SubFilterType& at() const
      {
        static_assert(i_ != 0, "invalid sub-filter index");
        return _first;
      }

      template<typename Algo_, typename SubVector_>
      void filter_rhs(PowerVector<SubVector_, 1>& vector) const
      {
        first().template filter_rhs<Algo_>(vector.first());
      }

      template<typename Algo_, typename SubVector_>
      void filter_sol(PowerVector<SubVector_, 1>& vector) const
      {
        first().template filter_sol<Algo_>(vector.first());
      }

      template<typename Algo_, typename SubVector_>
      void filter_def(PowerVector<SubVector_, 1>& vector) const
      {
        first().template filter_def<Algo_>(vector.first());
      }

      template<typename Algo_, typename SubVector_>
      void filter_cor(PowerVector<SubVector_, 1>& vector) const
      {
        first().template filter_cor<Algo_>(vector.first());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_FILTER_HPP
