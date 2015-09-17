#pragma once
#ifndef KERNEL_GLOBAL_FILTER_HPP
#define KERNEL_GLOBAL_FILTER_HPP 1

#include <kernel/global/vector.hpp>

namespace FEAST
{
  namespace Global
  {
    /**
     * \brief Global Filter wrapper class template
     *
     * \author Peter Zajac
     */
    template<typename LocalFilter_>
    class Filter
    {
    public:
      typedef Global::Vector<typename LocalFilter_::VectorType> VectorType;

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

      template<typename OtherLocalFilter_>
      void convert(const Global::Filter<OtherLocalFilter_>& other)
      {
        this->_filter.convert(*other);
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
} // namespace FEAST


#endif // KERNEL_GLOBAL_FILTER_HPP