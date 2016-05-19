#pragma once
#ifndef GLOBAL_GUARD_MATRIX_CONVERSION_HPP
#define GLOBAL_GUARD_MATRIX_CONVERSION_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/util/comm_base.hpp>
#include<kernel/lafem/matrix_mirror.hpp>
#include<kernel/assembly/mirror_assembler.hpp>

namespace FEAT
{
  namespace Global
  {
#ifdef FEAT_HAVE_MPI
    ///type-0 to type-1 matrix conversion
    struct SynchMat1
    {
      template<typename MT_, typename VMT_, template<typename, typename> class ST_>
      static MT_ exec(const MT_& origin,
                      const ST_<VMT_, std::allocator<VMT_> >& vec_mirrors_row,
                      const ST_<VMT_, std::allocator<VMT_> >& vec_mirrors_column,
                      const ST_<typename MT_::IndexType, std::allocator<typename MT_::IndexType> >& other_ranks,
                      const ST_<typename MT_::IndexType, std::allocator<typename MT_::IndexType> >& tags,
                      Util::Communicator communicator = Util::Communicator(MPI_COMM_WORLD) )
      {
        MT_ result;
        result.clone(origin);

        ST_<std::vector<char>, std::allocator<std::vector<char> > > recv_buf;
        ST_<std::vector<char>, std::allocator<std::vector<char> > > send_buf;

        ST_<Util::CommRequest, std::allocator<Util::CommRequest> > recvrequests(vec_mirrors_row.size());
        ST_<Util::CommRequest, std::allocator<Util::CommRequest> > sendrequests(vec_mirrors_row.size());
        ST_<Util::CommRequest, std::allocator<Util::CommRequest> > presendrequests(vec_mirrors_row.size());

        ST_<Util::CommStatus, std::allocator<Util::CommStatus> > recvstatus;
        ST_<Util::CommStatus, std::allocator<Util::CommStatus> > sendstatus;
        ST_<Util::CommStatus, std::allocator<Util::CommStatus> > prerecvstatus;
        ST_<Util::CommStatus, std::allocator<Util::CommStatus> > presendstatus;

        for(Index i(0) ; i < vec_mirrors_row.size() ; ++i)
        {
          LAFEM::MatrixMirror<VMT_> mat_mirror(vec_mirrors_row.at(i), vec_mirrors_column.at(i));
          MT_ sendbuf_mat;
          Assembly::MirrorAssembler::assemble_buffer_matrix(sendbuf_mat, mat_mirror, result);
          mat_mirror.gather(sendbuf_mat, result);

          send_buf.push_back(sendbuf_mat.serialise());
        }

        std::vector<Index> recv_msg_size(vec_mirrors_row.size());
        for(Index i(0); i < vec_mirrors_row.size();i++)
        {
          Util::CommStatus ss;

          presendstatus.push_back(ss);
          Index send_size(send_buf.at(i).size());
          Util::Comm::isend(&send_size,
                      Index(1),
                      other_ranks.at(i),
                      presendrequests.at(i),
                      tags.at(i),
                      communicator
              );
        }
        for(Index i(0) ; i < vec_mirrors_row.size() ; ++i)
        {

          Util::CommStatus rs;

          prerecvstatus.push_back(rs);

          Util::Comm::recv(&recv_msg_size.at(i),
                      Index(1),
                      other_ranks.at(i),
                      prerecvstatus.at(i),
                      tags.at(i),
                      communicator
              );
        }

        for(Index i(0) ; i < vec_mirrors_row.size() ; ++i)
        {

          Util::CommStatus rs;

          recvstatus.push_back(rs);

          recv_buf.push_back(std::vector<char>(recv_msg_size.at(i)));

          Util::Comm::irecv(recv_buf.at(i).data(),
                      Index(recv_msg_size.at(i)),
                      other_ranks.at(i),
                      recvrequests.at(i),
                      tags.at(i),
                      communicator
              );
        }

        for(Index i(0) ; i < vec_mirrors_row.size() ; ++i)
        {
          Util::CommStatus ss;

          sendstatus.push_back(ss);

          Util::Comm::isend(send_buf.at(i).data(),
                      Index(send_buf.at(i).size()),
                      other_ranks.at(i),
                      sendrequests.at(i),
                      tags.at(i),
                      communicator
              );
        }

        int* recvflags = new int[recvrequests.size()];
        int* taskflags = new int[recvrequests.size()];
        for(Index i(0) ; i < recvrequests.size() ; ++i)
        {
          recvflags[i] = 0;
          taskflags[i] = 0;
        }

        Index count(0);
        while(count != recvrequests.size())
        {
          //go through all requests round robin
          for(Index i(0) ; i < recvrequests.size() ; ++i)
          {
            if(taskflags[i] == 0)
            {
              Util::Comm::test(recvrequests.at(i), recvflags[i], recvstatus.at(i));
              if(recvflags[i] != 0)
              {
                MT_ other_mat(recv_buf.at(i));
                LAFEM::MatrixMirror<VMT_> mat_mirror(vec_mirrors_row.at(i), vec_mirrors_row.at(i));
                mat_mirror.scatter_axpy(result, other_mat);
                ++count;
                taskflags[i] = 1;
              }
            }
          }
        }

        for(Index i(0) ; i < sendrequests.size() ; ++i)
        {
          Util::CommStatus ws,ws2;
          Util::Comm::wait(sendrequests.at(i), ws);
          Util::Comm::wait(presendrequests.at(i), ws2);
        }

        delete[] recvflags;
        delete[] taskflags;

        return result;
      }

#else
    struct SynchMat1
    {
      template<typename MT_, typename VMT_, template<typename, typename> class ST_>
      static MT_ exec(const MT_& origin,
                     const ST_<VMT_, std::allocator<VMT_> >&,
                     const ST_<VMT_, std::allocator<VMT_> >&,
                     const ST_<typename MT_::IndexType, std::allocator<typename MT_::IndexType> >&,
                     const ST_<typename MT_::IndexType, std::allocator<typename MT_::IndexType> >&,
                     Util::Communicator = Util::Communicator(0))
      {
        MT_ result;
        result.clone(origin);
        return result;
      }

#endif
    };
  }
}

#endif
