#pragma once
#ifndef SCARC_GUARD_SCARC_LOG_HPP
#define SCARC_GUARD_SCARC_LOG_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/foundation/communication.hpp>

using namespace FEAST;
using namespace FEAST::Foundation;

namespace FEAST
{
  namespace ScaRC
  {
    template<template<typename, typename> class ST_ = std::vector>
    struct ScaRCLog
    {
      Index basetag;
      std::ostringstream stream;
      std::string msg;

      ScaRCLog(Index bt = 100000) :
        basetag(bt),
        stream(),
        msg()
      {
      }

#ifndef SERIAL
      std::string synch(Communicator c = Communicator(MPI_COMM_WORLD))
      {
        msg = std::string(stream.str());
        Index len_sb;
        len_sb = msg.size() + 1;
        Index len_rb[Comm::size(c)];

        Comm::allgather(&len_sb,
                        1,
                        &len_rb[0],
                        1,
                        c);

        ST_<Request, std::allocator<Request> >recv_req;
        ST_<Request, std::allocator<Request> >send_req;
        ST_<Status, std::allocator<Status> >recv_stat;
        ST_<Status, std::allocator<Status> >send_stat;
        char** recvbufs(new char*[Comm::size(c)]);
        char* sendbuf(new char[len_sb]);

        for(Index i(0) ; i < Comm::size(c) ; ++i)
        {
          Status rs;
          recv_stat.push_back(std::move(rs));
          Request rr;
          recv_req.push_back(std::move(rr));
          recvbufs[i] = new char[len_rb[i]];

          Comm::irecv(recvbufs[i],
                      len_rb[i],
                      i,
                      recv_req.at(i),
                      basetag + i,
                      c);
        }

        for(Index i(0) ; i < Comm::size(c) ; ++i)
        {
          Status ss;
          send_stat.push_back(std::move(ss));
          Request sr;
          send_req.push_back(std::move(sr));
          std::copy(msg.c_str(), msg.c_str() + len_sb, sendbuf);

          Comm::isend(sendbuf,
                      len_sb,
                      i,
                      send_req.at(i),
                      basetag + Comm::rank(c),
                      c);
        }

        int recvflags[recv_req.size()];
        int taskflags[recv_req.size()];
        for(Index i(0) ; i < recv_req.size() ; ++i)
        {
          recvflags[i] = 0;
          taskflags[i] = 0;
        }

        Index count(0);
        std::ostringstream res;
        while(count != recv_req.size())
        {
          for(Index i(0) ; i < recv_req.size() ; ++i)
          {
            if(taskflags[i] == 0)
            {
              Comm::test(recv_req.at(i), recvflags[i], recv_stat.at(i));
              if(recvflags[i] != 0)
              {
                res << recvbufs[i] << std::endl;
                ++count;
                taskflags[i] = 1;
              }
            }
          }
        }

        for(Index i(0) ; i < Comm::size(c) ; ++i)
        {
          delete[] recvbufs[i];
        }
        delete[] recvbufs;
        delete[] sendbuf;

        msg = res.str();
        return msg;
      }
#else
      std::string synch(Communicator c = Communicator(0))
      {
        msg = std::string(stream.str());
        return msg;
      }
#endif

      template<typename...DT_>
      int checkin_line(DT_&...data)
      {
        stream << Comm::rank() << "| " ;
        int r0[sizeof...(data)] = {(stream << data, 0)...};

        stream << std::endl;

        if(sizeof...(data) >= 0)
          return r0[0];
        else
          return 0;
      }
    };

  }
}

#endif
