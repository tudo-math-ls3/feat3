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

    protected:
      LocalFilter_ _filter;

    public:
      template<typename... Args_>
      explicit Filter(Args_&&... args) :
        _filter(std::forward<Args_>(args)...)
      {
      }

      LocalFilter_& operator*()
      {
        return _filter;
      }

      const LocalFilter_& operator*() const
      {
        return _filter;
      }

      template<typename OtherGlobalFilter_>
      void convert(const OtherGlobalFilter_ & other)
      {
        this->_filter.convert(*other);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _filter.bytes();
      }

      void filter_rhs(VectorType& vector) const
      {
        _filter.filter_rhs(*vector);
      }

      void filter_sol(VectorType& vector) const
      {
        _filter.filter_sol(*vector);
      }

      void filter_def(VectorType& vector) const
      {
        _filter.filter_def(*vector);
      }

      void filter_cor(VectorType& vector) const
      {
        _filter.filter_cor(*vector);
      }
    }; // class Filter<...>
  } // namespace Global
} // namespace FEAT


#endif // KERNEL_GLOBAL_FILTER_HPP
