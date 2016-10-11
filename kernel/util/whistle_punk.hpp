#pragma once
#ifndef UTIL_GUARD_WHISTLE_PUNK_HPP
#define UTIL_GUARD_WHISTLE_PUNK_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/util/comm_base.hpp>
#include<kernel/util/environment.hpp>

namespace FEAT
{
  namespace Util
  {

    template<template<typename, typename> class ST_ = std::vector>
      struct WhistlePunk
      {
        Index basetag;
        std::ostringstream stream;
        std::string msg;

        WhistlePunk() :
          basetag(Util::Environment::reserve_tag()),
          stream(),
          msg()
        {
        }

#ifdef FEAT_HAVE_MPI
        std::string synch(Util::Communicator c = Util::Communicator(MPI_COMM_WORLD))
        {
          msg = std::string(stream.str());
          Index len_sb;
          len_sb = msg.size() + 1;
          Index* len_rb = new Index[Util::Comm::size(c)];

          Util::Comm::allgather(&len_sb,
                                1,
                                &len_rb[0],
                                1,
                                c);

          ST_<Util::CommRequest, std::allocator<Util::CommRequest> >recv_req(Util::Comm::size(c));
          ST_<Util::CommRequest, std::allocator<Util::CommRequest> >send_req(Util::Comm::size(c));
          ST_<Util::CommStatus, std::allocator<Util::CommStatus> >recv_stat;
          ST_<Util::CommStatus, std::allocator<Util::CommStatus> >send_stat;
          char** recvbufs = new char*[Util::Comm::size(c)];
          char* sendbuf(new char[len_sb]);

          for(Index i(0) ; i < Util::Comm::size(c) ; ++i)
          {
            Util::CommStatus rs;
            recv_stat.push_back(std::move(rs));
            recvbufs[i] = new char[len_rb[i]];

            Util::Comm::irecv(recvbufs[i],
                len_rb[i],
                i,
                recv_req.at(i),
                basetag + i,
                c);
          }

          for(Index i(0) ; i < Util::Comm::size(c) ; ++i)
          {
            Util::CommStatus ss;
            send_stat.push_back(std::move(ss));
            std::copy(msg.c_str(), msg.c_str() + len_sb, sendbuf);

            Util::Comm::isend(sendbuf,
                len_sb,
                i,
                send_req.at(i),
                basetag + Util::Comm::rank(c),
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
                Util::Comm::test(recv_req.at(i), recvflags[i], recv_stat.at(i));
                if(recvflags[i] != 0)
                {
                  res << recvbufs[i] << std::endl;
                  ++count;
                  taskflags[i] = 1;
                }
              }
            }
          }

          for(Index i(0) ; i < Util::Comm::size(c) ; ++i)
          {
            delete[] recvbufs[i];
          }

          for(Index i(0) ; i < send_req.size() ; ++i)
          {
            Util::CommStatus ws;
            Util::Comm::wait(send_req.at(i), ws);
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
        std::string synch()
        {
          msg = std::string(stream.str());
          return msg;
        }
#endif

        template<typename...DT_>
          int log(DT_&...data)
          {
            stream << Util::Comm::rank() << "| " ;
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
