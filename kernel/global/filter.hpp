#pragma once
#ifndef KERNEL_GLOBAL_FILTER_HPP
#define KERNEL_GLOBAL_FILTER_HPP 1

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
      typedef LocalFilter_ LocalFilter;

      /// Our 'base' class type
      template <typename LocalFilter2_, typename Mirror2_ = Mirror_>
      using FilterType = class Filter<LocalFilter2_, Mirror2_>;

      /// this typedef lets you create a filter with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using FilterTypeByMDI = class Filter<typename LocalFilter_::template FilterType<Mem2_, DataType2_, IndexType2_>, typename Mirror_::template MirrorType<Mem2_, DataType2_, IndexType2_> >;

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

        this->local() = other->local().clone(mode);
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


#endif // KERNEL_GLOBAL_FILTER_HPP
