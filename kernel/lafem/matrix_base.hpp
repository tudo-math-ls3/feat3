#pragma once
#ifndef KERNEL_LAFEM_MATRIX_BASE_HPP
#define KERNEL_LAFEM_MATRIX_BASE_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    class MatrixBase
    {
    public:
      virtual const Index & rows() const = 0;
      virtual const Index & columns() const = 0;

      virtual ~MatrixBase()
      {
      }

    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
