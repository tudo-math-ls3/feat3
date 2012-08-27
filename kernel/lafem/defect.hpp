#pragma once
#ifndef KERNEL_LAFEM_DEFECT_HPP
#define KERNEL_LAFEM_DEFECT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Arch_, typename BType_>
    struct Defect
    {
    };

    template <>
    struct Defect<Archs::CPU, Archs::Generic>
    {
      template <typename DT_>
      static void value(DenseVector<Archs::CPU, DT_> & r, const DenseVector<Archs::CPU, DT_> & rhs, const SparseMatrixCSR<Archs::CPU, DT_> & a, const DenseVector<Archs::CPU, DT_> & b)
      {
        if (b.size() != a.columns())
          throw InternalError("Vector size does not match!");
        if (a.rows() != r.size())
          throw InternalError("Vector size does not match!");
        if (a.rows() != rhs.size())
          throw InternalError("Vector size does not match!");

        const DT_ * bp(b.elements());
        const DT_ * rhsp(rhs.elements());
        const Index * Aj(a.Aj());
        const DT_ * Ax(a.Ax());
        const Index * Ar(a.Ar());
        DT_ * rp(r.elements());
        const Index rows(a.rows());

        for (Index row(0) ; row < rows ; ++row)
        {
          DT_ sum(0);
          const Index end(Ar[row * 2 + 1]);
          for (Index i(Ar[row * 2]) ; i < end ; ++i)
          {
            sum += Ax[i] * bp[Aj[i]];
          }
          rp[row] = rhsp[row] - sum;
        }
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DEFECT_HPP
