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
    template <typename Mem_, typename DT_, typename IT_>
    class DenseVector;

    template <typename Mem_, typename DT_, typename IT_>
    class SparseVector;

    template <typename Mem_, typename DT_, typename IT_, int BlockSize_>
    class DenseVectorBlocked;

    template <typename Mem_, typename DT_, typename IT_, int BlockSize_>
    class SparseVectorBlocked;

    template <typename Mem_, typename DT_, typename IT_>
    class DenseMatrix;

    template <typename Mem_, typename DT_, typename IT_>
    class SparseMatrixCSR;

    template <typename Mem_, typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
    class SparseMatrixBCSR;

    template <typename Mem_, typename DT_, typename IT_>
    class SparseMatrixCOO;

    template <typename Mem_, typename DT_, typename IT_>
    class SparseMatrixELL;

    template <typename Mem_, typename DT_, typename IT_>
    class SparseMatrixBanded;

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_FORWARD_HPP
