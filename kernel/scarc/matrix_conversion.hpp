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

using namespace FEAST;
using namespace FEAST::Foundation;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;

namespace FEAST
{
  namespace ScaRC
  {
#ifndef SERIAL
    ///type-0 to type-1 matrix conversion
    template<
             typename Mem_,
             typename DT_,
             typename IT_,
             template<typename, typename, typename> class MT_>
    struct MatrixConversion
    {
      template<template<typename, typename> class ST_, typename VMT_>
      static MT_<Mem_, DT_, IT_> value(const MT_<Mem_, DT_, IT_>& origin,
                                       const ST_<VMT_, std::allocator<VMT_> >& vec_mirrors,
                                       const ST_<IT_, std::allocator<IT_> >& other_ranks,
                                       const ST_<IT_, std::allocator<IT_> >& tags,
                                       Communicator communicator = Communicator(MPI_COMM_WORLD) )
      {
        MT_<Mem_, DT_, IT_> result;
        result.clone(origin);

        ST_<std::pair<Index, char*>, std::allocator<std::pair<Index, char*> > > recv_buf;
        ST_<std::pair<Index, char*>, std::allocator<std::pair<Index, char*> > > send_buf;

        ST_<Request, std::allocator<Request> > recvrequests;
        ST_<Request, std::allocator<Request> > sendrequests;

        ST_<Status, std::allocator<Status> > recvstatus;
        ST_<Status, std::allocator<Status> > sendstatus;

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {
          MatrixMirror<Mem_, DT_> mat_mirror(vec_mirrors.at(i), vec_mirrors.at(i));
          MT_<Mem_, DT_, IT_> sendbuf_mat;
          Assembly::MirrorAssembler::assemble_buffer_matrix(sendbuf_mat, mat_mirror, result);
          mat_mirror.gather_op(sendbuf_mat, result);

          send_buf.push_back(sendbuf_mat.serialize());
        }

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {

          Request rr;
          Status rs;

          recvrequests.push_back(rr);
          recvstatus.push_back(rs);

          char* recv_buf_p(new char[send_buf.at(i).first]);
          recv_buf.push_back(std::pair<Index, char*>(send_buf.at(i).first, recv_buf_p));

          Comm::irecv(recv_buf.at(i).second,
                      recv_buf.at(i).first,
                      other_ranks.at(i),
                      recvrequests.at(i),
                      tags.at(i),
                      communicator
              );

        }

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {
          Request sr;
          Status ss;

          sendrequests.push_back(sr);
          sendstatus.push_back(ss);

          Comm::isend(send_buf.at(i).second,
                      send_buf.at(i).first,
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
              Comm::test(recvrequests.at(i), recvflags[i], recvstatus.at(i));
              if(recvflags[i] != 0)
              {
                MT_<Mem_, DT_, IT_> other_mat(recv_buf.at(i));
                MatrixMirror<Mem::Main, double> mat_mirror(vec_mirrors.at(i), vec_mirrors.at(i));
                SparseMatrixCSR<Mem::Main, double> buf_mat;
                Assembly::MirrorAssembler::assemble_buffer_matrix(buf_mat, mat_mirror, result);

                mat_mirror.gather_op(buf_mat, result);

                buf_mat.axpy<Algo::Generic>(buf_mat, other_mat);
                mat_mirror.scatter_op(result, buf_mat);
                ++count;
                taskflags[i] = 1;
              }
            }
          }
        }

        for(Index i(0) ; i < sendrequests.size() ; ++i)
        {
          Status ws;
          Comm::wait(sendrequests.at(i), ws);
        }

        delete[] recvflags;
        delete[] taskflags;
        for(auto recv_buf_p_i : recv_buf)
          delete[] recv_buf_p_i.second;
        for(auto send_buf_p_i : send_buf)
          delete[] send_buf_p_i.second;

        return result;
      }
#else
    template<
             typename Mem_,
             typename DT_,
             typename IT_,
             template<typename, typename, typename> class MT_>
    struct MatrixConversion
    {
      template<template<typename, typename> class ST_, typename VMT_>
      static MT_<Mem_, DT_, IT_> value(const MT_<Mem_, DT_, IT_>& origin,
                                       const ST_<VMT_, std::allocator<VMT_> >&,
                                       const ST_<IT_, std::allocator<IT_> >&,
                                       const ST_<IT_, std::allocator<IT_> >&,
                                       Communicator = Communicator(0)
                                       )
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
