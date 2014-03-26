#pragma once
#ifndef KERNEL_LAFEM_ARCH_SYNCH_HPP
#define KERNEL_LAFEM_ARCH_SYNCH_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_, typename Algo_, typename VectorT_>
      class SynchVec0GatewayBase
      {
        public:
          virtual VectorT_& value(VectorT_& x) const = 0;

          virtual ~SynchVec0GatewayBase()
          {
          }
      };

      template <typename Mem_, typename Algo_, typename VectorT_>
      class SynchVec1GatewayBase
      {
        public:
          virtual VectorT_& value(VectorT_& x) const = 0;

          virtual ~SynchVec1GatewayBase()
          {
          }
      };
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
