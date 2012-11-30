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
#include <kernel/lafem/product_matvec.hpp>
#include <kernel/lafem/defect.hpp>
#include <kernel/lafem/norm.hpp>
#include <kernel/lafem/component_product.hpp>
#include <kernel/lafem/preconditioner.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct Richardson
    {
      template <typename MT_, typename VT_>
      static void value(VT_ & x, const MT_ & A, const VT_ & b, Preconditioner<Algo_, MT_, VT_> & precon, Index max_iters, typename VT_::DataType eps_relative)
      {
        typedef typename VT_::DataType DT_;
        typedef typename VT_::MemType Arch_;

        DenseVector<Arch_, DT_> temp_0(x.size());
        DenseVector<Arch_, DT_> temp_1(x.size());

        Defect<Algo_>::value(temp_0, b, A, x);
        DT_ initial_defect = Norm2<Algo_>::value(temp_0);

        Index used_iters = 0;
        while(true)
        {
          ++used_iters;

          Defect<Algo_>::value(temp_0, b, A, x);
          DT_ current_defect = Norm2<Algo_>::value(temp_0);
          precon.apply(temp_1, temp_0);
          Sum<Algo_>::value(x, x, temp_1);

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
