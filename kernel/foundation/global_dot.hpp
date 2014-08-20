#pragma once
#ifndef FOUNDATION_GUARD_GLOBAL_DOT_HPP
#define FOUNDATION_GUARD_GLOBAL_DOT_HPP 1

#include<kernel/foundation/comm_base.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/lafem/arch/scale.hpp>
#include<kernel/lafem/arch/component_product.hpp>
#include<kernel/lafem/arch/dot_product.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

///TODO add communicators
namespace FEAST
{
  namespace Foundation
  {
      template <typename Mem_, typename Algo_>
      struct GlobalDot
      {
      };

      template <>
      struct GlobalDot<Mem::Main, Algo::Generic>
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
              return a.template dot<Algo::Generic>(b);

            ///assumes type-1 vector (full entries at inner boundaries)

            const typename VectorT_::DataType* a_data(a.elements());
            const typename VectorT_::DataType* b_data(b.elements());
            const typename VectorT_::DataType* f_data(frequencies.elements());
            typename VectorT_::DataType sum(0);
            for(Index i(0) ; i < a.size() ; ++i)
            {
              sum += a_data[i] * b_data[i] * typename VectorT_::DataType(1.)/f_data[i];
            }

            typename VectorT_::DataType sendbuf(sum), recvbuf;

            Status stat;

// The current MS-MPI implementation does not offer the MPI_Iallreduce function...
/// \todo Remove this workaround once MS does its homework.
#ifdef MSMPI_VER
            Comm::allreduce(&sendbuf, Index(1), &recvbuf);
#else
            Request req;

            Comm::iallreduce(&sendbuf, Index(1), &recvbuf, req);

            Comm::wait(req, stat);
#endif // MSMPI_VER
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
            r = a.template dot<Algo::Generic>(b);
            return r;
          }
#endif
      };
  }
}


#endif
