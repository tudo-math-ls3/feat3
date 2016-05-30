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
#ifndef SERIAL

          //typename PExecutorT_::PResult result(num_global_elems);
          typename PExecutorT_::PResult result(local_partitioning_result.clone());
          result.reset(num_global_elems);

          typename PExecutorT_::PResult::IndexType* sendbuf = local_partitioning_result.get();

          const Index commsize(Util::Comm::size(local_partitioning_result.get_comm()));

          int* sendcounts = new int[commsize];
          for(Index i(0) ; i < commsize ; ++i)
          {
            sendcounts[i] = int(local_partitioning_result.size());
          }
          //sendcounts[Util::Comm::rank(local_partitioning_result.get_comm())] = int(0);

          int* sdispls = new int[commsize];
          for(Index i(0) ; i < commsize ; ++i)
          {
            sdispls[i] = int(0);
          }

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

          int err(Util::Comm::alltoallv(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, local_partitioning_result.get_comm()));

          if(err != MPI_SUCCESS)
            throw(InternalError("Synch during partitioning failed!"));

          delete[] sendcounts;
          delete[] sdispls;
          delete[] recvcounts;
          delete[] rdispls;

          return result;
#else
          typename PExecutorT_::PResult result(num_global_elems);
          result = local_partitioning_result.clone();
          return result;
#endif

        }

#ifndef SERIAL
        static void exec(std::stringstream& iss, Util::Communicator comm = Util::Communicator(MPI_COMM_WORLD))
        {
          Index me(Util::Comm::rank(comm));
          Index size;
          std::string str;

          //bcast size
          if(me == 0)
          {
            str = (iss.str());
            size = Index(str.length());
          }
          // synchronize length
          Util::Comm::bcast(&size, 1, 0, comm);

          //allocate
          char* buf = new char[size + 1];

          //fill
          if(me == 0) //master
          {
            std::strcpy(buf, str.c_str());
          }

          //bcast data
          Util::Comm::bcast(buf, size, 0, comm);

          //convert
          if(me != 0)
          {
            std::string res_str(buf, size);
            iss << res_str;
          }

          delete[] buf;
        }
#endif
    };

  }
}
#endif
