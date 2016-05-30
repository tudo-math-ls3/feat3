#pragma once
#ifndef KERNEL_LAFEM_RICHARDSON_HPP
#define KERNEL_LAFEM_RICHARDSON_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/preconditioner.hpp>



namespace FEAT
{
  namespace LAFEM
  {
    struct Richardson
    {
      template <typename MT_, typename VT_>
      static void value(VT_ & x, const MT_ & A, const VT_ & b, Preconditioner<MT_, VT_> & precon, Index max_iters, typename VT_::DataType eps_relative)
      {
        typedef typename VT_::DataType DT_;
        typedef typename VT_::MemType Mem_;

        DenseVector<Mem_, DT_> temp_0(x.size());
        DenseVector<Mem_, DT_> temp_1(x.size());

        //temp_0.defect(b, A, x);
        A.apply(temp_0, x, b, -DT_(1));
        DT_ initial_defect = temp_0.norm2();

        Index used_iters = 0;
        while(true)
        {
          ++used_iters;

          //temp_0.defect(b, A, x);
          A.apply(temp_0, x, b, -DT_(1));
          DT_ current_defect = temp_0.norm2();
          precon.apply(temp_1, temp_0);
          x.axpy(x, temp_1);

          if(current_defect < eps_relative * initial_defect || current_defect < eps_relative || used_iters >= max_iters)
          {
            std::cout<<"Initial defect norm: " << initial_defect << ", Final defect norm: " << current_defect << ", used iters: " << used_iters << std::endl;
            break;
          }
        }
      }
    };

  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_RICHARDSON_HPP
