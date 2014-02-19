#pragma once
#ifndef KERNEL_LAFEM_POWER_FILTER_HPP
#define KERNEL_LAFEM_POWER_FILTER_HPP 1

#include <kernel/lafem/power_vector.hpp>

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
    class PowerFilter :
      protected PowerFilter<SubFilter_, count_-1>
    {
      // declare this class template as a friend for recursive inheritance
      template<typename,Index>
      friend class PowerFilter;

      /// base-class typedef
      typedef PowerFilter<SubFilter_, count_-1> BaseClass;

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
      /// the last sub-filter
      SubFilterType _last;

      /// base-class ctor; for internal use only
      explicit PowerFilter(BaseClass&& other_base, SubFilterType&& other_last) :
        BaseClass(std::move(other_base)),
        _last(std::move(other_last))
      {
      }

    public:
      /// default ctor
      PowerFilter()
      {
      }

      /// move CTOR
      PowerFilter(PowerFilter&& other) :
        BaseClass(static_cast<BaseClass&&>(other)),
        _last(std::move(other._last))
      {
      }

      /// move-assign operator
      PowerFilter& operator=(PowerFilter&& other)
      {
        if(this != &other)
        {
          base().operator=(static_cast<BaseClass&&>(other));
          _last = std::move(other._last);
        }
        return *this;
      }

      /// Creates and returns a deep copy of this filter.
      PowerFilter clone() const
      {
        return PowerFilter(base().clone(), last().clone());
      }

      /// \cond internal
      BaseClass& base()
      {
        return static_cast<BaseClass&>(*this);
      }

      const BaseClass& base() const
      {
        return static_cast<const BaseClass&>(*this);
      }

      SubFilterType& last()
      {
        return _last;
      }

      const SubFilterType& last() const
      {
        return _last;
      }
      /// \endcond

      template<Index i_>
      SubFilterType& at()
      {
        static_assert(i_ < count_, "invalid sub-filter index");
        return static_cast<PowerFilter<SubFilter_, i_+1>&>(*this)._last;
      }

      template<Index i_>
      const SubFilterType& at() const
      {
        static_assert(i_ < count_, "invalid sub-filter index");
        return static_cast<const PowerFilter<SubFilter_, i_+1>&>(*this)._last;
      }

      /** \copydoc UnitFilter::filter_rhs() */
      template<typename Algo_, typename SubVector_>
      void filter_rhs(PowerVector<SubVector_, count_>& vector) const
      {
        base().template filter_rhs<Algo_>(vector.base());
        last().template filter_rhs<Algo_>(vector.last());
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Algo_, typename SubVector_>
      void filter_sol(PowerVector<SubVector_, count_>& vector) const
      {
        base().template filter_sol<Algo_>(vector.base());
        last().template filter_sol<Algo_>(vector.last());
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Algo_, typename SubVector_>
      void filter_def(PowerVector<SubVector_, count_>& vector) const
      {
        base().template filter_def<Algo_>(vector.base());
        last().template filter_def<Algo_>(vector.last());
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Algo_, typename SubVector_>
      void filter_cor(PowerVector<SubVector_, count_>& vector) const
      {
        base().template filter_cor<Algo_>(vector.base());
        last().template filter_cor<Algo_>(vector.last());
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
      SubFilterType _last;

      explicit PowerFilter(SubFilterType&& other_last) :
        _last(std::move(other_last))
      {
      }

    public:
      PowerFilter()
      {
      }

      /// move CTOR
      PowerFilter(PowerFilter&& other) :
        _last(std::move(other._last))
      {
      }

      PowerFilter& operator=(PowerFilter&& other)
      {
        _last = std::move(other._last);
        return *this;
      }

      PowerFilter clone() const
      {
        return PowerFilter(last().clone());
      }

      SubFilterType& last()
      {
        return _last;
      }

      const SubFilterType& last() const
      {
        return _last;
      }

      template<Index i_>
      SubFilterType& at()
      {
        static_assert(i_ != 0, "invalid sub-filter index");
        return _last;
      }

      template<Index i_>
      const SubFilterType& at() const
      {
        static_assert(i_ != 0, "invalid sub-filter index");
        return _last;
      }

      template<typename Algo_, typename SubVector_>
      void filter_rhs(PowerVector<SubVector_, 1>& vector) const
      {
        last().template filter_rhs<Algo_>(vector.last());
      }

      template<typename Algo_, typename SubVector_>
      void filter_sol(PowerVector<SubVector_, 1>& vector) const
      {
        last().template filter_sol<Algo_>(vector.last());
      }

      template<typename Algo_, typename SubVector_>
      void filter_def(PowerVector<SubVector_, 1>& vector) const
      {
        last().template filter_def<Algo_>(vector.last());
      }

      template<typename Algo_, typename SubVector_>
      void filter_cor(PowerVector<SubVector_, 1>& vector) const
      {
        last().template filter_cor<Algo_>(vector.last());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_FILTER_HPP
