#pragma once
#ifndef KERNEL_LAFEM_RICHARDSON_HPP
#define KERNEL_LAFEM_RICHARDSON_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sum.hpp>
#include <kernel/lafem/product.hpp>
#include <kernel/lafem/defect.hpp>
#include <kernel/lafem/norm.hpp>
#include <kernel/lafem/element_product.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename BType_>
    struct Richardson
    {
      template <typename Arch_, typename DT_>
      static void value(DenseVector<Arch_, DT_> & x, const SparseMatrixCSR<Arch_, DT_> & A, const DenseVector<Arch_, DT_> & b, Index max_iters, DT_ eps_relative)
      {
        DenseVector<Arch_, DT_> temp_0(x.size());
        DenseVector<Arch_, DT_> temp_1(x.size());

        Defect<Arch_, BType_>::value(temp_0, b, A, x);
        DT_ initial_defect = Norm2<Arch_, BType_>::value(temp_0);

        Index used_iters = 0;
        while(true)
        {
          ++used_iters;

          Defect<Arch_, BType_>::value(temp_0, b, A, x);
          DT_ current_defect = Norm2<Arch_, BType_>::value(temp_0);

          DenseVector<Arch_, DT_> P(x.size());
          for (Index i(0) ; i < x.size() ; ++i)
            P(i, 0.7 / A(i, i));

          ElementProduct<Arch_, BType_>::value(temp_1, P, temp_0);
          Sum<Arch_, BType_>::value(x, x, temp_1);


          if(current_defect < eps_relative * initial_defect || current_defect < eps_relative || used_iters >= max_iters)
          {
            std::cout<<"Initial defect norm: " << initial_defect << ", Final defect norm: " << current_defect << ", used iters: " << used_iters << std::endl;
            break;
          }
        }
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_RICHARDSON_HPP
