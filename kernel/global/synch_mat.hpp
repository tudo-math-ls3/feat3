#pragma once
#ifndef SCARC_GUARD_MATRIX_CONVERSION_HPP
#define SCARC_GUARD_MATRIX_CONVERSION_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/scarc/scarc_data.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/matrix_mirror.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/archs.hpp>
#include <kernel/assembly/mirror_assembler.hpp>

namespace FEAT
{
  namespace Global
  {
#ifndef SERIAL
    ///type-0 to type-1 matrix conversion
    template<typename Mem_>
    struct SynchMat0
    {
      template<typename DT_,
               typename IT_,
               template<typename, typename, typename> class MT_,
               template<typename, typename> class ST_, typename VMT_>
      static MT_<Mem_, DT_, IT_> value(const MT_<Mem_, DT_, IT_>& origin,
                                       const ST_<VMT_, std::allocator<VMT_> >& vec_mirrors,
                                       const ST_<IT_, std::allocator<IT_> >& other_ranks,
                                       const ST_<IT_, std::allocator<IT_> >& tags,
                                       Foundation::Communicator communicator = Foundation::Communicator(MPI_COMM_WORLD) )
      {
        MT_<Mem_, DT_, IT_> result;
        result.clone(origin);

        ST_<std::vector<char>, std::allocator<std::vector<char> > > recv_buf;
        ST_<std::vector<char>, std::allocator<std::vector<char> > > send_buf;

        ST_<Foundation::Request, std::allocator<Foundation::Request> > recvrequests(vec_mirrors.size());
        ST_<Foundation::Request, std::allocator<Foundation::Request> > sendrequests(vec_mirrors.size());
        ST_<Foundation::Request, std::allocator<Foundation::Request> > presendrequests(vec_mirrors.size());

        ST_<Foundation::Status, std::allocator<Foundation::Status> > recvstatus;
        ST_<Foundation::Status, std::allocator<Foundation::Status> > sendstatus;
        ST_<Foundation::Status, std::allocator<Foundation::Status> > prerecvstatus;
        ST_<Foundation::Status, std::allocator<Foundation::Status> > presendstatus;

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {
          LAFEM::MatrixMirror<VMT_> mat_mirror(vec_mirrors.at(i), vec_mirrors.at(i));
          LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> sendbuf_mat;
          Assembly::MirrorAssembler::assemble_buffer_matrix(sendbuf_mat, mat_mirror, result);
          mat_mirror.gather(sendbuf_mat, result);

          send_buf.push_back(sendbuf_mat.serialise());
        }

        std::vector<Index> recv_msg_size(vec_mirrors.size());
        for(Index i(0); i < vec_mirrors.size();i++)
        {
          Foundation::Status ss;

          presendstatus.push_back(ss);
          Index send_size(send_buf.at(i).size());
          Foundation::Comm::isend(&send_size,
                      Index(1),
                      other_ranks.at(i),
                      presendrequests.at(i),
                      tags.at(i),
                      communicator
              );
        }
        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {

          Foundation::Status rs;

          prerecvstatus.push_back(rs);

          Foundation::Comm::recv(&recv_msg_size.at(i),
                      Index(1),
                      other_ranks.at(i),
                      prerecvstatus.at(i),
                      tags.at(i),
                      communicator
              );
        }

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {

          Foundation::Status rs;

          recvstatus.push_back(rs);

          recv_buf.push_back(std::vector<char>(recv_msg_size.at(i)));

          Foundation::Comm::irecv(recv_buf.at(i).data(),
                      Index(recv_msg_size.at(i)),
                      other_ranks.at(i),
                      recvrequests.at(i),
                      tags.at(i),
                      communicator
              );
        }

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {
          Foundation::Status ss;

          sendstatus.push_back(ss);

          Foundation::Comm::isend(send_buf.at(i).data(),
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
              Foundation::Comm::test(recvrequests.at(i), recvflags[i], recvstatus.at(i));
              if(recvflags[i] != 0)
              {
                LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> other_mat(recv_buf.at(i));
                LAFEM::MatrixMirror<VMT_> mat_mirror(vec_mirrors.at(i), vec_mirrors.at(i));
                mat_mirror.scatter_axpy(result, other_mat);
                ++count;
                taskflags[i] = 1;
              }
            }
          }
        }

        for(Index i(0) ; i < sendrequests.size() ; ++i)
        {
          Foundation::Status ws,ws2;
          Foundation::Comm::wait(sendrequests.at(i), ws);
          Foundation::Comm::wait(presendrequests.at(i), ws2);
        }

        delete[] recvflags;
        delete[] taskflags;

        return result;
      }

#else
    template<typename Mem_>
    struct SynchMat0
    {
      template<typename DT_,
               typename IT_,
               template<typename, typename, typename> class MT_,
               template<typename, typename> class ST_, typename VMT_>
      static MT_<Mem_, DT_, IT_> value(const MT_<Mem_, DT_, IT_>& origin,
                                       const ST_<VMT_, std::allocator<VMT_> >&,
                                       const ST_<IT_, std::allocator<IT_> >&,
                                       const ST_<IT_, std::allocator<IT_> >&,
                                       Foundation::Communicator = Foundation::Communicator(0))
      {
        MT_<Mem_, DT_, IT_> result;
        result.clone(origin);
        return result;
      }

#endif
    };
  }
}

#endif
