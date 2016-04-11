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
      int count_>
    class PowerFilter
    {
      // Note: the case = 1 is specialised below
      static_assert(count_ > 1, "invalid block size");

      // declare this class template as a friend for recursive inheritance
      template<typename,int>
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
      /// sub-filter index-type
      typedef typename SubFilter_::IndexType IndexType;

      /// number of filter blocks
      static constexpr int num_blocks = count_;

      /// vector type
      typedef PowerVector<typename SubFilter_::VectorType, count_> VectorType;

    protected:
      /// the first sub-filter
      SubFilterType _first;
      /// the remaining part
      RestClass _rest;

      /// base-class ctor; for internal use only
      explicit PowerFilter(SubFilterType&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
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

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
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

      template<int i_>
      SubFilterType& at()
      {
        static_assert((0 <= i_) && (i_ < count_), "invalid sub-filter index");
        return PowerElement<i_, SubFilterType>::get(*this);
      }

      template<int i_>
      const SubFilterType& at() const
      {
        static_assert((0 <= i_) && (i_ < count_), "invalid sub-filter index");
        return PowerElement<i_, SubFilterType>::get(*this);
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
    }; // class PowerFilter<...>

    /// \cond internal
    template<typename SubFilter_>
    class PowerFilter<SubFilter_, 1>
    {
      template<typename,int>
      friend class PowerFilter;

    public:
      typedef SubFilter_ SubFilterType;
      typedef typename SubFilter_::MemType MemType;
      typedef typename SubFilter_::DataType DataType;
      typedef typename SubFilter_::IndexType IndexType;

      static constexpr int num_blocks = 1;

      typedef PowerVector<typename SubFilter_::VectorType, 1> VectorType;

    protected:
      SubFilterType _first;

      explicit PowerFilter(SubFilterType&& the_first) :
        _first(std::move(the_first))
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

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes();
      }

      SubFilterType& first()
      {
        return _first;
      }

      const SubFilterType& first() const
      {
        return _first;
      }

      template<int i_>
      SubFilterType& at()
      {
        static_assert(i_ != 0, "invalid sub-filter index");
        return _first;
      }

      template<int i_>
      const SubFilterType& at() const
      {
        static_assert(i_ != 0, "invalid sub-filter index");
        return _first;
      }

      void filter_rhs(VectorType& vector) const
      {
        first().filter_rhs(vector.first());
      }

      void filter_sol(VectorType& vector) const
      {
        first().filter_sol(vector.first());
      }

      void filter_def(VectorType& vector) const
      {
        first().filter_def(vector.first());
      }

      void filter_cor(VectorType& vector) const
      {
        first().filter_cor(vector.first());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_FILTER_HPP
