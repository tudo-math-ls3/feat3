#pragma once
#ifndef KERNEL_FOUNDATION_PSYNCH_HPP
#define KERNEL_FOUNDATION_PSYNCH_HPP 1

#include<cstring>
#include<kernel/base_header.hpp>
#include<kernel/foundation/base.hpp>
#include<kernel/util/comm_base.hpp>
#include<kernel/foundation/pgraph.hpp>
#include<kernel/foundation/pexecutor.hpp>

#ifdef FEAT_HAVE_PARMETIS
FEAT_DISABLE_WARNINGS
#include<parmetis.h>
FEAT_RESTORE_WARNINGS
#endif

namespace FEAT
{
  namespace Foundation
  {

    template<typename PExecutorT_>
    class PSynch
    {
      public:
        static typename PExecutorT_::PResult exec(const typename PExecutorT_::PResult& local_partitioning_result,
                                                  typename PExecutorT_::PResult::IndexType num_global_elems
                                                  )
        {
#ifdef FEAT_HAVE_MPI

          //typename PExecutorT_::PResult result(num_global_elems);
          typename PExecutorT_::PResult result(local_partitioning_result.clone());
          result.reset(num_global_elems);

          typename PExecutorT_::PResult::IndexType* sendbuf = local_partitioning_result.get();

          const Index commsize(Util::Comm::size(local_partitioning_result.get_comm()));

          typename PExecutorT_::PResult::IndexType* recvbuf(result.get());

          int* recvcounts = new int[commsize];
          for(Index i(0) ; i < commsize ; ++i)
          {
            recvcounts[i] = int(local_partitioning_result.get_vtxdist()[i + 1] - local_partitioning_result.get_vtxdist()[i]);
          }
          //recvcounts[Util::Comm::rank(local_partitioning_result.get_comm())] = int(0);

          int* rdispls = new int[commsize];
          for(Index i(0) ; i < commsize ; ++i)
          {
            rdispls[i] = int(local_partitioning_result.get_vtxdist()[i]);
          }

          int err(Util::Comm::allgatherv(sendbuf, int(local_partitioning_result.size()), recvbuf, recvcounts, rdispls, local_partitioning_result.get_comm()));

          if(err != MPI_SUCCESS)
            throw(InternalError("Synch during partitioning failed!"));

          delete[] recvcounts;
          delete[] rdispls;

          return result;
#else
          typename PExecutorT_::PResult result(num_global_elems);
          result = local_partitioning_result.clone();
          return result;
#endif

        }

    };

  }
}
#endif
