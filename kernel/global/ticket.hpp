#pragma once
#ifndef KERNEL_GLOBAL_TICKET_HPP
#define KERNEL_GLOBAL_TICKET_HPP 1

#include<kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/comm_base.hpp>
#include<kernel/util/time_stamp.hpp>
#include<kernel/util/statistics.hpp>
#include <vector>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Ticket class for asynchronous global operations on scalars
     */
    template <typename DT_>
    class ScalTicket
    {
      public:
        /// Our request for the corresponding iallreduce mpi call
        std::shared_ptr<Util::CommRequest> req;
        /// buffer containing the send data
        DT_ x;
        /// buffer containing the received data
        DT_ r;
        /// should we compute the sqrt of the result
        bool sqrt;
        /// signals, whether wait was already called
        bool finished;

        /**
         * \brief Constructor
         *
         * Ticket constructor
         *
         * \param[in] in The value to be accumulated by this rank
         *
         */
        ScalTicket(DT_ & in) :
          req(std::make_shared<Util::CommRequest>()),
          x(in),
          sqrt(false),
          finished(false)
        {
        }

        /// Unwanted copy constructor: Do not implement!
        ScalTicket(const ScalTicket &) = delete;
        /// Unwanted copy assignment operator: Do not implement!
        ScalTicket & operator=(const ScalTicket &) = delete;

        /**
         * \brief wait method
         *
         * wait for completion barrier
         *
         * \returns the accumulated data
         *
         */
        DT_ wait()
        {
          TimeStamp ts_start;

          Util::CommStatusIgnore stat;
          Util::Comm::wait(*req, stat);

          finished = true;

          TimeStamp ts_stop;
          Statistics::add_time_mpi_wait(ts_stop.elapsed(ts_start));

          if (sqrt)
            return Math::sqrt(r);
          else
            return r;
        }

        /// Destructor
        ~ScalTicket()
        {
          XASSERT(finished);
        }
    };

#ifdef FEAT_HAVE_MPI
    /**
     * \brief Ticket class for asynchronous global operations on vectors
     */
    template <typename VT_, typename MTs_>
    class VecTicket
    {
      public:
        using BufferVectorType = LAFEM::DenseVector<Mem::Main, typename VT_::DataType, typename VT_::IndexType>;

        VT_& target;
        const VT_& frequencies;
        const MTs_& mirrors;
        bool scale;
        bool finished;
        std::vector<BufferVectorType> sendbufs, recvbufs;
        std::vector<std::shared_ptr<Util::CommRequest>> recvrequests;
        std::vector<std::shared_ptr<Util::CommRequest>> sendrequests;
        std::vector<int> recvflags;

        /**
         * \brief Constructor
         *
         * VecTicket constructor
         *
         * \param[out] target_in the target vector to be updated
         * \param[in] frequencies_in the frequencies vector
         * \param[in] mirrors_in a list of corresponding mirrors for the neighbour2neighbour synchronisation
         */
        VecTicket(VT_ & target_in, const VT_ & frequencies_in, const MTs_ & mirrors_in) :
          target(target_in),
          frequencies(frequencies_in),
          mirrors(mirrors_in),
          scale(false),
          finished(false),
          sendbufs(mirrors.size()),
          recvbufs(mirrors.size()),
          recvflags(mirrors.size())
        {
          for (Index i(0) ; i < mirrors.size() ; ++i)
          {
            sendrequests.push_back(std::make_shared<Util::CommRequest>());
            recvrequests.push_back(std::make_shared<Util::CommRequest>());
            sendbufs.at(i) = mirrors.at(i).create_buffer_vector();
            recvbufs.at(i) = mirrors.at(i).create_buffer_vector();
          }

          for(Index i(0) ; i < recvrequests.size() ; ++i)
          {
            recvflags.at(i) = 0;
          }
        }

        /// Unwanted copy constructor: Do not implement!
        VecTicket(const VecTicket &) = delete;
        /// Unwanted copy assignment operator: Do not implement!
        VecTicket & operator=(const VecTicket &) = delete;

        /**
         * \brief wait method
         *
         * wait for completion barrier
         * on return, the target vector has been updated
         *
         */
        void wait()
        {
          TimeStamp ts_start;

          Util::CommStatusIgnore s;
          Index count(0);
          //handle receives
          while(count != recvrequests.size())
          {
            //go through all requests round robin
            for(Index i(0) ; i < recvrequests.size() ; ++i)
            {
              if(recvflags.at(i) == 0)
              {
                Util::Comm::test(*recvrequests.at(i), recvflags.data()[i], s);
                if(recvflags.at(i) != 0)
                {
                  mirrors.at(i).scatter_axpy_dual(target, recvbufs.at(i));
                  ++count;
                }
              }
            }
          }

          for(Index i(0) ; i < sendrequests.size() ; ++i)
          {
            Util::Comm::wait(*sendrequests.at(i), s);
          }

          finished = true;

          if (scale)
            target.component_product(target, frequencies);

          TimeStamp ts_stop;
          Statistics::add_time_mpi_wait(ts_stop.elapsed(ts_start));
        }

        /// Destructor
        ~VecTicket()
        {
          XASSERT(finished);
        }
    };
#else
    template <typename VT_, typename MTs_>
    class VecTicket
    {
      public:
        void wait()
        {
        }
    };
#endif //FEAT_HAVE_MPI
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_TICKET_HPP
