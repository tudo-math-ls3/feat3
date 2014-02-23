#pragma once
#ifndef KERNEL_LAFEM_VECTOR_BASE_HPP
#define KERNEL_LAFEM_VECTOR_BASE_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    template <typename IT_>
    class VectorBase
    {
      public:
        //virtual const Index & size() const = 0;

        virtual ~VectorBase()
        {
        }

    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
