// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/global/vector.hpp>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Global Filter wrapper class template
     *
     * \author Peter Zajac
     */
    template<typename LocalFilter_, typename Mirror_>
    class Filter
    {
    public:
      typedef Global::Vector<typename LocalFilter_::VectorType, Mirror_> VectorType;
      typedef Mirror_ MirrorType;
      typedef LocalFilter_ LocalFilterType;

      /// Our 'base' class type
      template <typename LocalFilter2_, typename Mirror2_ = Mirror_>
      using FilterType = Filter<LocalFilter2_, Mirror2_>;

      /// this typedef lets you create a filter with different Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using FilterTypeByDI = Filter<typename LocalFilter_::template FilterType<DataType2_, IndexType2_>, typename Mirror_::template MirrorType<DataType2_, IndexType2_> >;

      static constexpr bool is_global = true;
      static constexpr bool is_local = false;

    protected:
      LocalFilter_ _filter;

    public:
      template<typename... Args_>
      explicit Filter(Args_&&... args) :
        _filter(std::forward<Args_>(args)...)
      {
      }

      LocalFilter_& local()
      {
        return _filter;
      }

      const LocalFilter_& local() const
      {
        return _filter;
      }

      template<typename OtherGlobalFilter_>
      void convert(const OtherGlobalFilter_ & other)
      {
        this->_filter.convert(other.local());
      }

      Filter clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return Filter(_filter.clone(mode));
      }

      void clone(const Filter& other, LAFEM::CloneMode mode = LAFEM::CloneMode::Weak)
      {
        XASSERTM(&(other.local()) != &(this->local()), "Trying to self-clone a Global::Filter!");

        this->local() = other.local().clone(mode);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _filter.bytes();
      }

      void filter_rhs(VectorType& vector) const
      {
        _filter.filter_rhs(vector.local());
      }

      void filter_sol(VectorType& vector) const
      {
        _filter.filter_sol(vector.local());
      }

      void filter_def(VectorType& vector) const
      {
        _filter.filter_def(vector.local());
      }

      void filter_cor(VectorType& vector) const
      {
        _filter.filter_cor(vector.local());
      }
    }; // class Filter<...>
  } // namespace Global
} // namespace FEAT
