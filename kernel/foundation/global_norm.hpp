#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_NORM_HPP
#define FOUNDATION_GUARD_GLOBAL_NORM_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/global_dot.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

///TODO add communicators
namespace FEAST
{
  namespace Foundation
  {
      template <typename Mem_, typename Algo_>
      struct GlobalNorm2
      {
        public:

          template<typename VectorT_>
          static typename VectorT_::DataType value(typename VectorT_::DataType& r,
                                                   const VectorT_& a,
                                                   const VectorT_& frequencies)
          {
            ///assumes type-1 vector (full entries at inner boundaries)
            GlobalDot<Mem_, Algo_>::value(r, a, a , frequencies);
            r = typename VectorT_::DataType(std::sqrt(r));
            return r;
          }
      };

      template <typename Mem_, typename Algo_>
      struct GlobalNorm2Squared
      {
        public:

          template<typename VectorT_>
          static typename VectorT_::DataType value(typename VectorT_::DataType& r,
                                                   const VectorT_& a,
                                                   const VectorT_& frequencies)
          {
            ///assumes type-1 vector (full entries at inner boundaries)
            GlobalDot<Mem_, Algo_>::value(r, a, a , frequencies);
            return r;
          }
      };
  }
}


#endif
