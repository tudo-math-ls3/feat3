#pragma once
#ifndef SCARC_GUARD_MATRIX_CONVERSION_HPP
#define SCARC_GUARD_MATRIX_CONVERSION_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/scarc/scarc_data.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/matrix_mirror.hpp>
#include <kernel/archs.hpp>
#include <kernel/assembly/mirror_assembler.hpp>

namespace FEAST
{
  namespace ScaRC
  {
#ifndef SERIAL
    ///type-0 to type-1 matrix conversion
    template<typename Mem_,
             typename Algo_>
    struct MatrixConversion
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

        ST_<Foundation::Request, std::allocator<Foundation::Request> > recvrequests;
        ST_<Foundation::Request, std::allocator<Foundation::Request> > sendrequests;

        ST_<Foundation::Status, std::allocator<Foundation::Status> > recvstatus;
        ST_<Foundation::Status, std::allocator<Foundation::Status> > sendstatus;

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {
          LAFEM::MatrixMirror<VMT_> mat_mirror(vec_mirrors.at(i), vec_mirrors.at(i));
          MT_<Mem_, DT_, IT_> sendbuf_mat;
          Assembly::MirrorAssembler::assemble_buffer_matrix(sendbuf_mat, mat_mirror, result);
          mat_mirror.template gather<Algo_>(sendbuf_mat, result);

          send_buf.push_back(sendbuf_mat.serialise());
        }

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {

          Foundation::Request rr;
          Foundation::Status rs;

          recvrequests.push_back(rr);
          recvstatus.push_back(rs);

          recv_buf.push_back(std::vector<char>(send_buf.at(i).size()));

          Foundation::Comm::irecv(recv_buf.at(i).data(),
                      Index(recv_buf.at(i).size()),
                      other_ranks.at(i),
                      recvrequests.at(i),
                      tags.at(i),
                      communicator
              );

        }

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {
          Foundation::Request sr;
          Foundation::Status ss;

          sendrequests.push_back(sr);
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
                MT_<Mem_, DT_, IT_> other_mat(recv_buf.at(i));
                LAFEM::MatrixMirror<VMT_> mat_mirror(vec_mirrors.at(i), vec_mirrors.at(i));
                LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> buf_mat;
                Assembly::MirrorAssembler::assemble_buffer_matrix(buf_mat, mat_mirror, result);

                mat_mirror.template gather<Algo_>(buf_mat, result);
                buf_mat.template axpy<Algo_>(buf_mat, other_mat);
                mat_mirror.template scatter<Algo_>(result, buf_mat);
                ++count;
                taskflags[i] = 1;
              }
            }
          }
        }

        for(Index i(0) ; i < sendrequests.size() ; ++i)
        {
          Foundation::Status ws;
          Foundation::Comm::wait(sendrequests.at(i), ws);
        }

        delete[] recvflags;
        delete[] taskflags;

        return result;
      }
#else
    template<typename Mem_,
             typename Algo_>
    struct MatrixConversion
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
