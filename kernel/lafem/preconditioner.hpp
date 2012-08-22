#pragma once
#ifndef KERNEL_LAFEM_PRECONDITIONER_HPP
#define KERNEL_LAFEM_PRECONDITIONER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/element_product.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename BType_, typename MT_, typename VT_>
    class Preconditioner
    {
      public:
        virtual void apply(VT_ & out, const VT_ & in) = 0;
    };

    template <typename BType_, typename MT_, typename VT_>
    class JacobiPreconditioner : public virtual Preconditioner<BType_, MT_, VT_>
    {
      private:
        VT_ _jac;

      public:
        JacobiPreconditioner(const MT_ & A, typename VT_::data_type damping) :
          _jac(A.rows())
        {
          for (Index i(0) ; i < _jac.size() ; ++i)
            _jac(i, damping / A(i, i));
        }

        virtual void apply(VT_ & out, const VT_ & in)
        {
          ElementProduct<typename VT_::arch_type, BType_>::value(out, _jac, in);
        };
    };


  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_PRECONDITIONER_HPP
