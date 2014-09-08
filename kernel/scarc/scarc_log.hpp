#pragma once
#ifndef SCARC_GUARD_SCARC_LOG_HPP
#define SCARC_GUARD_SCARC_LOG_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/environment.hpp>

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

      ScaRCLog() :
        basetag(Foundation::Environment::reserve_tag()),
        stream(),
        msg()
      {
      }

#ifndef SERIAL
      std::string synch(Foundation::Communicator c = Foundation::Communicator(MPI_COMM_WORLD))
      {
        msg = std::string(stream.str());
        Index len_sb;
        len_sb = msg.size() + 1;
        Index* len_rb = new Index[Foundation::Comm::size(c)];

        Foundation::Comm::allgather(&len_sb,
                        1,
                        &len_rb[0],
                        1,
                        c);

        ST_<Foundation::Request, std::allocator<Foundation::Request> >recv_req;
        ST_<Foundation::Request, std::allocator<Foundation::Request> >send_req;
        ST_<Foundation::Status, std::allocator<Foundation::Status> >recv_stat;
        ST_<Foundation::Status, std::allocator<Foundation::Status> >send_stat;
        char** recvbufs = new char*[Foundation::Comm::size(c)];
        char* sendbuf(new char[len_sb]);

        for(Index i(0) ; i < Foundation::Comm::size(c) ; ++i)
        {
          Foundation::Status rs;
          recv_stat.push_back(std::move(rs));
          Foundation::Request rr;
          recv_req.push_back(std::move(rr));
          recvbufs[i] = new char[len_rb[i]];

          Foundation::Comm::irecv(recvbufs[i],
                      len_rb[i],
                      i,
                      recv_req.at(i),
                      basetag + i,
                      c);
        }

        for(Index i(0) ; i < Foundation::Comm::size(c) ; ++i)
        {
          Foundation::Status ss;
          send_stat.push_back(std::move(ss));
          Foundation::Request sr;
          send_req.push_back(std::move(sr));
          std::copy(msg.c_str(), msg.c_str() + len_sb, sendbuf);

          Foundation::Comm::isend(sendbuf,
                      len_sb,
                      i,
                      send_req.at(i),
                      basetag + Foundation::Comm::rank(c),
                      c);
        }

        int* recvflags = new int[recv_req.size()];
        int* taskflags = new int[recv_req.size()];
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
              Foundation::Comm::test(recv_req.at(i), recvflags[i], recv_stat.at(i));
              if(recvflags[i] != 0)
              {
                res << recvbufs[i] << std::endl;
                ++count;
                taskflags[i] = 1;
              }
            }
          }
        }

        for(Index i(0) ; i < Foundation::Comm::size(c) ; ++i)
        {
          delete[] recvbufs[i];
        }
        delete[] recvbufs;
        delete[] sendbuf;
        delete[] recvflags;
        delete[] taskflags;
        delete[] len_rb;

        msg = res.str();
        return msg;
      }
#else
      std::string synch(Foundation::Communicator c = Foundation::Communicator(0))
      {
        msg = std::string(stream.str());
        return msg;
      }
#endif

      template<typename...DT_>
      int checkin_line(DT_&...data)
      {
        stream << Foundation::Comm::rank() << "| " ;
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
