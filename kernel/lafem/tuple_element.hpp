#pragma once
#ifndef KERNEL_LAFEM_TUPLE_ELEMENT_HPP
#define KERNEL_LAFEM_TUPLE_ELEMENT_HPP 1

#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    template<
      Index i_,
      typename First_,
      typename... Rest_>
    struct TupleElement
    {
      typedef typename TupleElement<i_-1, Rest_...>::Type Type;

      template<typename Meta_>
      static Type& get(Meta_& meta)
      {
        return TupleElement<i_-1, Rest_...>::get(meta.rest());
      }

      template<typename Meta_>
      static const Type& get(const Meta_& meta)
      {
        return TupleElement<i_-1, Rest_...>::get(meta.rest());
      }
    };

    template<typename First_, typename... Rest_>
    struct TupleElement<Index(0), First_, Rest_...>
    {
      typedef First_ Type;

      template<typename Meta_>
      static Type& get(Meta_& meta)
      {
        return meta.first();
      }

      template<typename Meta_>
      static const Type& get(const Meta_& meta)
      {
        return meta.first();
      }
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_ELEMENT_HPP
