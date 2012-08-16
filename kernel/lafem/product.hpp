#pragma once
#ifndef KERNEL_LAFEM_PRODUCT_HPP
#define KERNEL_LAFEM_PRODUCT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>



namespace FEAST
{
  template <typename Arch_, typename BType_>
  struct Product
  {
  };

  template <>
  struct Product<Archs::CPU, Archs::Generic>
  {
    template <typename DT_>
    static void value(DenseVector<Archs::CPU, DT_> & r, const SparseMatrixCSR<Archs::CPU, DT_> & a, const DenseVector<Archs::CPU, DT_> & b)
    {
      if (b.size() != a.columns())
        throw InternalError("Vector size does not match!");
      if (a.rows() != r.size())
        throw InternalError("Vector size does not match!");

      const DT_ * bp(b.elements());
      const Index * Aj(a.Aj());
      const DT_ * Ax(a.Ax());
      const Index * Ar(a.Ar());
      DT_ * rp(r.elements());
      const Index rows(a.rows());

      for (Index row(0) ; row < rows ; ++row)
      {
        DT_ sum(0);
        const Index end(Ar[row + 1]);
        for (Index i(Ar[row]) ; i < end ; ++i)
        {
          sum += Ax[i] * bp[Aj[i]];
        }
        rp[row] = sum;
      }
    }
  };

} // namespace FEAST

#endif // KERNEL_LAFEM_PRODUCT_HPP
