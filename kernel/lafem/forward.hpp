#pragma once
#ifndef KERNEL_LAFEM_FORWARD_HPP
#define KERNEL_LAFEM_FORWARD_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    //forward declarations
    template <typename Mem_, typename DT_>
    class DenseVector;

    template <typename Mem_, typename DT_>
    class SparseMatrixCSR;

    template <typename Mem_, typename DT_>
    class SparseMatrixCOO;

    template <typename Mem_, typename DT_>
    class SparseMatrixELL;

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_FORWARD_HPP
