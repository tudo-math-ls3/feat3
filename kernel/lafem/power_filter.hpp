// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_POWER_FILTER_HPP
#define KERNEL_LAFEM_POWER_FILTER_HPP 1

#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/meta_element.hpp>

namespace FEAT
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

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using FilterType = PowerFilter<typename SubFilterType::template FilterType<Mem2_, DT2_, IT2_>, count_>;

      /// this typedef lets you create a filter with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using FilterTypeByMDI = FilterType<Mem2_, DataType2_, IndexType2_>;

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

      /// \brief Creates a clone of itself
      PowerFilter clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return PowerFilter(first().clone(clone_mode), rest().clone(clone_mode));
      }

      /// \brief Clones data from another PowerFilter
      void clone(const PowerFilter & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _first.clone(other._first.get_filter_vector(), clone_mode);
        _rest.clone(other._rest, clone_mode);
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

      /**
       * \brief Returns a sub-filter block.
       *
       * \param[in] i
       * The index of the sub-filter block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-filter at position \p i.
       */
      SubFilterType& get(int i)
      {
        XASSERTM((0 <= i) && (i < count_), "invalid sub-filter index");
        return (i == 0) ? _first : _rest.get(i-1);
      }

      /** \copydoc get() */
      const SubFilterType& get(int i) const
      {
        XASSERTM((0 <= i) && (i < count_), "invalid sub-filter index");
        return (i == 0) ? _first : _rest.get(i-1);
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

      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using FilterType = PowerFilter<typename SubFilterType::template FilterType<Mem2_, DT2_, IT2_>, Index(1)>;

      /// this typedef lets you create a filter with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using FilterTypeByMDI = FilterType<Mem2_, DataType2_, IndexType2_>;

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

      /// \brief Creates a clone of itself
      PowerFilter clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return PowerFilter(first().clone(clone_mode));
      }

      /// \brief Clones data from another PowerFilter
      void clone(const PowerFilter & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _first.clone(other._first.get_filter_vector(), clone_mode);
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

      SubFilterType& get(int i)
      {
        XASSERTM(i == 0, "invalid sub-filter index");
        return _first;
      }

      const SubFilterType& get(int i) const
      {
        XASSERTM(i == 0, "invalid sub-filter index");
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
} // namespace FEAT

#endif // KERNEL_LAFEM_POWER_FILTER_HPP
