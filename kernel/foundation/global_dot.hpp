#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_DOT_HPP
#define FOUNDATION_GUARD_GLOBAL_DOT_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/lafem/arch/scale.hpp>
#include<kernel/lafem/arch/component_product.hpp>
#include<kernel/lafem/arch/dot_product.hpp>

namespace FEAST
{
  namespace Foundation
  {
      /// \todo add communicators
      template <typename Mem_>
      struct GlobalDot
      {
      };

      template <>
      struct GlobalDot<Mem::Main>
      {
        public:
#ifndef SERIAL
          template<typename VectorT_>
          static typename VectorT_::DataType value(typename VectorT_::DataType& r,
                                                   const VectorT_& a,
                                                   const VectorT_& b,
                                                   const VectorT_& frequencies)
          {
            if(frequencies.size() == 0)
              return a.dot(b);

            // assumes that a and b are type-1 vectors (full entries at inner boundaries)
            // Compute a^T diag(frequencies) b
            typename VectorT_::DataType sum = frequencies.triple_dot(a, b);

            typename VectorT_::DataType sendbuf(sum), recvbuf;

            Status stat;

            Request req;

            Comm::iallreduce(&sendbuf, Index(1), &recvbuf, req);

            Comm::wait(req, stat);

            r = recvbuf;
            return r;
          }
#else
          template<typename VectorT_>
          static typename VectorT_::DataType value(typename VectorT_::DataType& r,
                                                   const VectorT_& a,
                                                   const VectorT_& b,
                                                   const VectorT_&)
          {
            r = a.dot(b);
            return r;
          }
#endif
      };
  }
}


#endif
