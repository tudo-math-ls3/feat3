// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/string.hpp>

#include <numeric>
#include <type_traits>
#include <vector>

namespace FEAT::Global
{
  namespace Intern
  {
    /**
     * \brief Determine the MPI pack size of a scalar type \c T_
     *
     * \tparam T_ Scalar type to determine pack size of
     *
     * \param[in] comm MPI communicator
     *
     * \returns The pack size of \c T_
     *
     * \note This is just a convenience wrapper around \c comm.pack_size() for use in fold-expressions
     */
    template<typename T_>
    int pack_size(const Dist::Comm& comm)
    {
      int size(0);
      comm.pack_size<T_>(1, &size);
      return size;
    }

    /**
     * \brief Determine an upper bound on the buffer size required to pack a tuple using Dist::Comm::pack
     *
     * \tparam Ts_ Tuple types
     *
     * \param[in] comm MPI communicator
     *
     * \returns An upper bound on the buffer size required to pack a tuple of \c Ts_
     */
    template<typename... Ts_>
    int tuple_pack_size(const Dist::Comm& comm)
    {
      // NOTE: std::decay_t and std::remove_reference_t don't seem to work in fold expressions
      return (pack_size<typename std::decay<typename std::remove_reference<Ts_>::type>::type>(comm) + ...);
    }

    /**
     * \brief Pack a tuple into a buffer
     *
     * \tparam i Index into parameter pack
     * \tparam Ts_ Tuple types
     *
     * \param[in] comm MPI communicator
     * \param[in] tuple Tuple to pack
     * \param[out] outbuf Output buffer
     * \param[in] outsize Size of output buffer
     * \param[inout] position Offset into buffer for packing tuple value at position \c i
     */
    template<std::size_t i = 0, typename... Ts_>
    void pack_tuple(
      const Dist::Comm& comm,
      const std::tuple<Ts_...>& tuple,
      void* outbuf,
      std::size_t outsize,
      int* position)
    {
      using T = std::decay_t<std::remove_reference_t<decltype(std::get<i>(tuple))>>;

      const T& t = std::get<i>(tuple);
      comm.pack<T>(&t, 1, outbuf, outsize, position);

      if constexpr(i + 1 != sizeof...(Ts_))
      {
        pack_tuple<i + 1>(comm, tuple, outbuf, outsize, position);
      }
    }

    /**
     * \brief Unpack a tuple from a buffer
     *
     * \tparam i Index into parameter pack
     * \tparam Ts_ Tuple types
     *
     * \param[in] comm MPI communicator
     * \param[in] tuple Tuple to unpack into
     * \param[out] inbuf Input buffer
     * \param[in] insize Size of input buffer
     * \param[inout] position Offset into buffer for unpacking tuple value at position \c i
     */
    template<std::size_t i = 0, typename... Ts_>
    void unpack_tuple(
      const Dist::Comm& comm,
      std::tuple<Ts_...>& tuple,
      const void* inbuf,
      std::size_t insize,
      int* position)
    {
      using T = std::decay_t<std::remove_reference_t<decltype(std::get<i>(tuple))>>;

      T& t = std::get<i>(tuple);
      comm.unpack(inbuf, insize, position, &t, 1);

      if constexpr(i + 1 != sizeof...(Ts_))
      {
        unpack_tuple<i + 1>(comm, tuple, inbuf, insize, position);
      }
    }

    /**
     * \brief Retrieves the nth element of a range
     *
     * \tparam Range_ Range type
     *
     * \param[in] r Range to retrieve element of
     * \param[in] n Index of element to retrieve
     *
     * \returns The nth element of the range
     */
    template<typename Range_>
    auto nth_of_range(const Range_& r, Index n) -> decltype(*r.begin())&
    {
      auto it = r.begin();
      std::advance(it, n);
      return *it;
    }

    /**
     * \brief Normalizes a parameter-pack by removing references and constness
     *
     * \tparam T_ Parameter pack
     *
     * For a given parameter pack \c T_ the normalized tuple is the type of the tuples
     * that actually get sent via MPI in the RemoteLambda class.
     */
    template<typename... T_>
    using NormalizedTuple =
      std::tuple<typename std::decay<typename std::remove_reference<T_>::type>::type...>; // NOLINT
  } // namespace Intern

#ifdef FEAT_HAVE_MPI
  /**
   * \brief MPI-abstraction for running scalar lambdas on other ranks
   *
   * This is the base template declaration and must not be instantiated.
   * See the function-type specialization for the actual implementation.
   */
  template<typename T_>
  class ScalarRemoteLambda;

  /**
   * \brief MPI-abstraction for running scalar lambdas on other ranks
   *
   * \tparam R_ Return type of the lambda function
   * \tparam Args_ Parameter-pack of argument types of the lambda function
   *
   * This class-template allows calling scalar lambda functions on other ranks.
   * We call a lambda function scalar if it returns a single MPI-transmittable type.
   *
   * \see RemoteLambda<Range_<R_>(Args_...)> for the specialization for vector-valued lambdas
   */
  template<typename R_, typename... Args_>
  class ScalarRemoteLambda<R_(Args_...)>
  {
  public:
    /// Return type of the lambda function
    using ReturnType = R_;

    /// Tuple of arguments for the lambda function
    using ArgumentTuple = Intern::NormalizedTuple<Args_...>;

  private:
    /// Tag for messages carrying an argument tuple
    static constexpr int tag_argument = 1;
    /// Tag for messages carrying a return value
    static constexpr int tag_return = 2;
    /// The communicator used by the remote lambda
    Dist::Comm _comm;

  public:
    /// Constructor
    explicit ScalarRemoteLambda(const Dist::Comm& comm) : _comm(comm.comm_dup())
    {
    }

    /**
     * \brief Remote lambda call
     *
     * \tparam RankRange_
     * Range of `int`s, must have `begin` and `end` methods.
     * \tparam IndexRange_
     * Range of `Index`, must have `begin` and `end` methods.
     * \tparam ArgsRange_ Range of argument tuples
     * Range of `ArgumentTuples`s, must have `begin` and `end` methods.
     * \tparam Lambda_
     * Lambda function type. Must be compatible with \c ReturnType and \c ArgumentTuple
     *
     * \param[in] ranks
     * Range of ranks, for each entry a remote call will be made to that rank
     * \param[in] indices
     * Range of indices, must have same size as \c ranks and contain valid indices into the \c args range
     * \param[in] args
     * Range of argument tuples
     * \param[in] lambda
     * Lambda function to call on this rank
     *
     * \returns A vector containing the return values of all calls described by `ranks`, `indices`, and
     * `args`, in the same order.
     *
     * Assume we have ranges `ranks`, `indices`, and `args`.
     * Assume `result = remote_lambda.call(ranks, indices, args, lambda)` for some `remote_lambda` and some `lambda`.
     * Let `n` be the size of ranks. Then for 0 <= i < n:
     * `result[i] = lambda(args[indices[i]])` executed on rank `rank[i]`
     *
     * \note
     * This method works by exchanging argument tuples between ranks, calling the local lambdas, and then exchanging
     * return values back. Note that only arguments and return values are exchanged, not the lambdas themselves. This is
     * both easier/possible to implement and allows the lambda to capture from its environment.
     *
     * \note
     * This is a collective operation.
     */
    template<typename RankRange_, typename IndexRange_, typename ArgsRange_, typename Lambda_>
    std::vector<ReturnType>
    call(const RankRange_& ranks, const IndexRange_& indices, const ArgsRange_& args, Lambda_ lambda)
    {
      XASSERTM(
        std::distance(ranks.begin(), ranks.end()) == std::distance(indices.begin(), indices.end()),
        "size mismatch between ranks and indices");

      const auto num_ranks = (std::size_t)_comm.size();
      const auto num_calls = (std::size_t)std::distance(ranks.begin(), ranks.end());
      const auto pack_size = (std::size_t)Intern::tuple_pack_size<Args_...>(_comm);

      ////////////////////////////////////////////
      // Phase 1 - Exchange number of messages
      ////////////////////////////////////////////

      std::vector<int> receivecounts(num_ranks, 0);
      {
        std::vector<int> sendcounts(num_ranks, 0);
        for(auto i : ranks)
        {
          sendcounts[(std::size_t)i] += 1;
        }
        _comm.alltoall(sendcounts.data(), 1, receivecounts.data(), 1);
      }

      ////////////////////////////////////////////
      // Phase 2 - Send args and receive results
      ////////////////////////////////////////////

      Dist::RequestVector arg_sends;
      Dist::RequestVector result_receives;
      std::vector<char> arg_sendbuffer(num_calls * pack_size, 0);
      std::vector<ReturnType> return_receivebuffer(num_calls);
      {
        auto rank_iter = ranks.begin();
        auto index_iter = indices.begin();

        std::size_t send_pos(0);
        std::size_t receive_pos(0);
        while(rank_iter != ranks.end())
        {
          const Index idx = *index_iter++;
          const int rank = *rank_iter++;
          const ArgumentTuple& arg = Intern::nth_of_range(args, idx);

          // NOTE(mmuegge): We need to tag arguments and return types differently,
          // because MPI will match arg-sends and return-receives if the return type
          // and the packed tuple happen to have the same size.

          // Post receive for final result
          result_receives.push_back(_comm.irecv(&return_receivebuffer[receive_pos], 1, rank, tag_return));

          // Pack tuple for sending arg
          int pos = (int)send_pos;
          Intern::pack_tuple(_comm, arg, arg_sendbuffer.data(), num_calls * pack_size, &pos);

          // Post send for sending argument tuple
          arg_sends.push_back(_comm.isend(&arg_sendbuffer[send_pos], pack_size, Dist::dt_packed, rank, tag_argument));

          // Increment buffer positions
          receive_pos++;
          send_pos += pack_size;
        }
      }

      ////////////////////////////////////////////
      // Phase 3 - Receive args and send results
      ////////////////////////////////////////////

      const auto num_args = (std::size_t)std::accumulate(receivecounts.begin(), receivecounts.end(), 0);
      Dist::RequestVector arg_receives;
      std::vector<char> arg_receivebuffer(num_args * pack_size, 0);
      {
        std::size_t pos(0);
        for(std::size_t rank(0); rank < num_ranks; rank++)
        {
          const int receive_count = receivecounts[rank];

          for(int i(0); i < receive_count; i++)
          {
            arg_receives.push_back(
              _comm.irecv(&arg_receivebuffer[pos], pack_size, Dist::dt_packed, (int)rank, tag_argument));
            pos += pack_size;
          }
        }
      }

      Dist::RequestVector result_sends;
      std::vector<ReturnType> return_sendbuffer(num_args);
      {
        std::size_t idx(0);
        Dist::Status status;

        while(arg_receives.wait_any(idx, status))
        {
          ArgumentTuple t;
          int arg_pos = (int)(idx * pack_size);
          Intern::unpack_tuple(_comm, t, arg_receivebuffer.data(), num_args * pack_size, &arg_pos);

          return_sendbuffer[idx] = std::apply(lambda, t);
          result_sends.push_back(_comm.isend(&return_sendbuffer[idx], 1, status.source(), tag_return));
        }
      }

      arg_sends.wait_all();
      arg_receives.wait_all();
      result_sends.wait_all();
      result_receives.wait_all();

      return return_receivebuffer;
    }
  };

  // NOTE(mmuegge): Originally both ScalarRemoteLambda and VectorRemoteLambda were
  // just specialization of RemoteLambda:
  // RemoteLambda<R(Args...)> and RemoteLambda<Range<R>(Args...)>
  // But icx < 2025 is unable to properly choose the Range specialization for lambdas returning
  // e.g. std::vectors.
  // Thus the two cases were split into separate classes, both taking a simple function type parameter.

  /**
   * \brief MPI-abstraction for running vector-valued lambdas on other ranks
   *
   * This is the base template declaration and must not be instantiated.
   * See the function-type specialization for the actual implementation.
   */
  template<typename T_>
  class VectorRemoteLambda;

  /**
   * \brief MPI-abstraction for running vector-valued lambdas on other ranks
   *
   * \tparam Range_ Type of range returned by the lambda function. Must be a contiguous range!
   * \tparam R_ Inner type of ranges returned by the lambda function
   * \tparam Args_ Parameter-pack of argument types of the lambda function
   *
   * This class-template allows calling vector-valued lambda functions on other ranks.
   * We call a lambda function vector-valued if it returns a rang of MPI-transmittable types.
   *
   * \warning
   * Values sent directly from the produced ranges. The \c Range_ type must thus be a contiguous range!
   */
  template<typename R_, typename... Args_>
  class VectorRemoteLambda<R_(Args_...)>
  {
  public:
    /// Range type returned by lambda
    using ReturnType = R_;
    /// Argument tuple type for lambda
    using ArgumentTuple = Intern::NormalizedTuple<Args_...>;

  private:
    /// Tag for messages carrying an argument tuple
    static constexpr int tag_argument = 1;
    /// Tag for messages carrying a range
    static constexpr int tag_return = 2;
    /// Tag for messages carrying the size of a range
    static constexpr int tag_size = 3;
    /// Communicator used by this remote lambda
    Dist::Comm _comm;

  public:
    /// Constructor
    explicit VectorRemoteLambda(const Dist::Comm& comm) : _comm(comm.comm_dup())
    {
    }

    /**
     * \brief Remote lambda call
     *
     * \tparam RankRange_
     * Range of `int`s, must have `begin` and `end` methods.
     * \tparam IndexRange_
     * Range of `Index`, must have `begin` and `end` methods.
     * \tparam ArgsRange_ Range of argument tuples
     * Range of `ArgumentTuples`s, must have `begin` and `end` methods.
     * \tparam Lambda_
     * Type of lambda function. Must be compatible with \c ReturnType and \c ArgumentTuple
     *
     * \param[in] ranks
     * Range of ranks, for each entry a remote call will be made to that rank
     * \param[in] indices
     * Range of indices, must have same size as \c ranks and contain valid indices into the \c args range
     * \param[in] args
     * Range of argument tuples
     * \param[in] lambda
     * Lambda function to call on this rank
     *
     * \returns A vector containing the return values of all calls described by `ranks`, `indices`, and
     * `args`, in the same order.
     *
     * Assume we have ranges `ranks`, `indices`, and `args`.
     * Assume `result = remote_lambda.call(ranks, indices, args, lambda)` for some `remote_lambda` and some `lambda`.
     * Let `n` be the size of ranks. Then for 0 <= i < n:
     * `result[i] = lambda(args[indices[i]])` executed on rank `rank[i]`
     *
     * \note
     * This method works by exchanging argument tuples between ranks, calling the local lambdas, and then exchanging
     * return values back. Returns values are exchanged using two messages each, one containing the size of the range
     * and one containing the range itself. Note that only arguments and return values are exchanged, not the lambdas
     * themselves. This is both easier/possible to implement and allows the lambda to capture from its environment.
     *
     * \note
     * This is a collective operation.
     */
    template<typename RankRange_, typename IndexRange_, typename ArgsRange_, typename Lambda_>
    std::vector<ReturnType>
    call(const RankRange_& ranks, const IndexRange_& indices, const ArgsRange_& args, Lambda_ lambda)
    {
      // Inner type of ranges returned by lambda
      using ValueType = std::decay_t<std::remove_reference_t<decltype(*std::declval<ReturnType>().begin())>>;

      XASSERTM(
        std::distance(ranks.begin(), ranks.end()) == std::distance(indices.begin(), indices.end()),
        "size mismatch between ranks and indices");

      const auto num_ranks = (std::size_t)_comm.size();
      const auto num_calls = (std::size_t)std::distance(ranks.begin(), ranks.end());
      const auto pack_size = (std::size_t)Intern::tuple_pack_size<Args_...>(_comm);

      ////////////////////////////////////////////
      // Phase 1 - Exchange number of calls
      ////////////////////////////////////////////

      // Each rank needs to know how many arguments it is to receive from each rank
      // so that it can post the necessary sends and receives.

      std::vector<int> receivecounts(num_ranks, 0);
      {
        std::vector<int> sendcounts(num_ranks, 0);
        for(auto i : ranks)
        {
          sendcounts[(std::size_t)i] += 1;
        }
        _comm.alltoall(sendcounts.data(), 1, receivecounts.data(), 1);
      }

      ////////////////////////////////////////////////////
      // Phase 2 - Send args and receive sizes of ranges
      ////////////////////////////////////////////////////

      // Sends and receives are posted in the same order as given by the user
      Dist::RequestVector size_receives;
      Dist::RequestVector arg_sends;
      std::vector<std::size_t> range_message_sizes(num_calls, 0);
      std::vector<char> arg_sendbuffer(num_calls * pack_size, 0);
      {
        auto rank_iter = ranks.begin();
        auto index_iter = indices.begin();

        std::size_t send_pos(0);
        std::size_t receive_pos(0);
        while(rank_iter != ranks.end())
        {
          const Index idx = *index_iter++;
          const int rank = *rank_iter++;
          const ArgumentTuple& arg = Intern::nth_of_range(args, idx);

          // Pack tuple for sending arg
          int pos = (int)send_pos;
          Intern::pack_tuple(_comm, arg, arg_sendbuffer.data(), num_calls * pack_size, &pos);

          // Post send for sending argument tuple
          arg_sends.push_back(_comm.isend(&arg_sendbuffer[send_pos], pack_size, Dist::dt_packed, rank, tag_argument));

          // Post receive for size of resulting range
          size_receives.push_back(_comm.irecv(&range_message_sizes[receive_pos], 1, rank, tag_size));

          // Increment buffer positions
          send_pos += pack_size;
          receive_pos++;
        }
      }

      ///////////////////////////////////////////////////
      // Phase 3 - Receive args and send range sizes and ranges
      ///////////////////////////////////////////////////

      const auto num_args = (std::size_t)std::accumulate(receivecounts.begin(), receivecounts.end(), 0);
      Dist::RequestVector arg_receives;
      std::vector<char> arg_receivebuffer(num_args * pack_size, 0);
      {
        std::size_t pos(0);
        for(std::size_t rank(0); rank < num_ranks; rank++)
        {
          const int receive_count = receivecounts[rank];

          for(int i(0); i < receive_count; i++)
          {
            // Post receive for argument
            arg_receives.push_back(
              _comm.irecv(&arg_receivebuffer[pos], pack_size, Dist::dt_packed, (int)rank, tag_argument));
            pos += pack_size;
          }
        }
      }

      Dist::RequestVector size_sends;
      Dist::RequestVector range_sends;
      std::vector<ReturnType> ranges(num_args);
      std::vector<std::size_t> range_sizes(num_args, 0);
      std::vector<std::size_t> range_buffer_sizes(num_ranks, 0);
      {
        std::size_t idx(0);
        Dist::Status status;

        while(arg_receives.wait_any(idx, status))
        {
          int arg_pos = (int)(idx * pack_size);

          ArgumentTuple t;
          Intern::unpack_tuple(_comm, t, arg_receivebuffer.data(), num_args * pack_size, &arg_pos);

          ranges[idx] = std::apply(lambda, t);
          range_sizes[idx] = (std::size_t)std::distance(ranges[idx].begin(), ranges[idx].end());
          range_buffer_sizes[std::size_t(status.source())] += range_sizes[idx];

          // Post send for sizes
          size_sends.push_back(_comm.isend(&range_sizes[idx], 1, status.source(), tag_size));

          // Pointer to beginning of range, or nullptr if range is empty
          // NOTE: This assumes the range is contiguous in memory
          const void* range_pointer = range_sizes[idx] > 0 ? &*ranges[idx].begin() : nullptr;

          // Post send for range
          range_sends.push_back(
            _comm.isend(range_pointer, range_sizes[idx], Dist::autotype<ValueType>(), status.source(), tag_return));
        }
      }

      /////////////////////////////////////////////////////
      // Phase 4 - Exchange buffer sizes for all ranks
      /////////////////////////////////////////////////////

      // At this point each ranks knows how many values it will send to each other rank.
      // Exchange these values so that reach ranks can allocate its range_buffer before receiving individual ranges.

      // TODO: Is there a smarter way than doing _comm.size() reduces?
      std::size_t range_buffer_size(0);
      for(std::size_t rank(0); rank < num_ranks; rank++)
      {
        _comm.reduce(&range_buffer_sizes[rank], &range_buffer_size, 1, Dist::op_sum, (int)rank);
      }

      ////////////////////////////////////////////
      // Phase 5 - Receive range sizes and ranges
      ////////////////////////////////////////////

      size_sends.wait_all();

      Dist::RequestVector range_receives;
      std::vector<ValueType> range_buffer(range_buffer_size);
      std::vector<std::size_t> range_buffer_offsets;
      std::exclusive_scan(
        range_message_sizes.begin(),
        range_message_sizes.end(),
        std::back_inserter(range_buffer_offsets),
        std::size_t(0));
      {
        std::size_t idx(0);
        Dist::Status status;

        // For each message containing a range size, post the corresponding range receive
        while(size_receives.wait_any(idx, status))
        {
          std::size_t offset = range_buffer_offsets[idx];
          std::size_t size = range_message_sizes[idx];
          void* range_pointer = size > 0 ? &range_buffer[offset] : nullptr;
          range_receives.push_back(
            _comm.irecv(range_pointer, size, Dist::autotype<ValueType>(), status.source(), tag_return));
        }
      }

      ////////////////////////////////////////////
      // Phase 6 - Reconstruct ranges
      ////////////////////////////////////////////

      arg_sends.wait_all();
      arg_receives.wait_all();
      range_receives.wait_all();
      range_sends.wait_all();

      std::vector<ReturnType> result;

      // NOTE: We are keeping the same message order throughout all of the above code,
      // so we can just walk linearly over the range_buffer and the reconstructed ranges
      // will have the same order as the input parameters.
      std::size_t pos(0);
      for(std::size_t size : range_message_sizes)
      {
        auto begin = range_buffer.data() + pos;
        auto end = begin + size;
        result.emplace_back(begin, end);
        pos += size;
      }

      return result;
    }
  };
#else // FEAT_HAVE_MPI
  /**
   * \brief MPI-abstraction for running scalar lambdas on other ranks
   *
   * This is the base template declaration and must not be instantiated.
   * See the function-type specialization for the actual implementation.
   */
  template<typename T_>
  class ScalarRemoteLambda;

  /**
   * \brief MPI-abstraction for running scalar lambdas on other ranks
   *
   * \tparam R_ Return type of the lambda function
   * \tparam Args_ Parameter-pack of argument types of the lambda function
   *
   * This class-template allows calling scalar lambda functions on other ranks.
   * We call a lambda function scalar if it returns a single MPI-transmittable type.
   *
   * \see RemoteLambda<Range_<R_>(Args_...)> for the specialization for vector-valued lambdas
   */
  template<typename R_, typename... Args_>
  class ScalarRemoteLambda<R_(Args_...)>
  {
  public:
    /// Return type of the lambda function
    using ReturnType = R_;

    /// Tuple of arguments for the lambda function
    using ArgumentTuple = Intern::NormalizedTuple<Args_...>;

  private:
    /// Tag for messages carrying an argument tuple
    static constexpr int tag_argument = 1;
    /// Tag for messages carrying a return value
    static constexpr int tag_return = 2;
    /// The communicator used by the remote lambda
    Dist::Comm _comm;

  public:
    /// Constructor
    explicit ScalarRemoteLambda(const Dist::Comm& comm) : _comm(comm.comm_dup())
    {
    }

    /**
     * \brief Remote lambda call
     *
     * \tparam RankRange_
     * Range of `int`s, must have `begin` and `end` methods.
     * \tparam IndexRange_
     * Range of `Index`, must have `begin` and `end` methods.
     * \tparam ArgsRange_ Range of argument tuples
     * Range of `ArgumentTuples`s, must have `begin` and `end` methods.
     * \tparam Lambda_
     * Type of lambda function. Must be compatible with \c ReturnType and \c ArgumentTuple
     *
     * \param[in] ranks
     * Range of ranks, for each entry a remote call will be made to that rank
     * \param[in] indices
     * Range of indices, must have same size as \c ranks and contain valid indices into the \c args range
     * \param[in] args
     * Range of argument tuples
     * \param[in] lambda
     * Lambda function to call on this rank
     *
     * \returns A vector containing the return values of all calls described by `ranks`, `indices`, and
     * `args`, in the same order.
     *
     * Assume we have ranges `ranks`, `indices`, and `args`.
     * Assume `result = remote_lambda.call(ranks, indices, args, lambda)` for some `remote_lambda` and some `lambda`.
     * Let `n` be the size of ranks. Then for 0 <= i < n:
     * `result[i] = lambda(args[indices[i]])`
     */
    template<typename RankRange_, typename IndexRange_, typename ArgsRange_, typename Lambda_>
    std::vector<ReturnType>
    call(const RankRange_& ranks, const IndexRange_& indices, const ArgsRange_& args, Lambda_ lambda)
    {
      XASSERTM(
        std::distance(ranks.begin(), ranks.end()) == std::distance(indices.begin(), indices.end()),
        "size mismatch between ranks and indices");

      const auto num_calls = (std::size_t)std::distance(ranks.begin(), ranks.end());

      std::vector<ReturnType> result;
      result.reserve(num_calls);

      auto rank_iter = ranks.begin();
      auto index_iter = indices.begin();
      while(rank_iter != ranks.end())
      {
        const Index idx = *index_iter++;
        const int rank = *rank_iter++;
        const ArgumentTuple& arg = Intern::nth_of_range(args, idx);

        XASSERTM(rank == 0, "invalid rank for non-mpi RemoteLambda");

        result.push_back(std::apply(lambda, arg));
      }

      return result;
    }
  };

  /**
   * \brief MPI-abstraction for running vector-valued lambdas on other ranks
   *
   * This is the base template declaration and must not be instantiated.
   * See the function-type specialization for the actual implementation.
   */
  template<typename T_>
  class VectorRemoteLambda;

  /**
   * \brief MPI-abstraction for running vector-values lambdas on other ranks
   *
   * \tparam Range_ Type of range returned by the lambda function. Must be a contiguous range!
   * \tparam R_ Inner type of ranges returned by the lambda function
   * \tparam Args_ Parameter-pack of argument types of the lambda function
   *
   * This class-template allows calling vector-valued lambda functions on other ranks.
   * We call a lambda function vector-valued if it returns a rang of MPI-transmittable types.
   *
   * \warning
   * Values sent directly from the produced ranges. The \c Range_ type must thus be a contiguous range!
   *
   * \see RemoteLambda<Range_<R_>(Args_...)> for the specialization for scalar-valued lambdas
   */
  template<typename R_, typename... Args_>
  class VectorRemoteLambda<R_(Args_...)>
  {
  public:
    /// Range type returned by lambda
    using ReturnType = R_;
    /// Argument tuple type for lambda
    using ArgumentTuple = Intern::NormalizedTuple<Args_...>;

  private:
    /// Tag for messages carrying an argument tuple
    static constexpr int tag_argument = 1;
    /// Tag for messages carrying a range
    static constexpr int tag_return = 2;
    /// Tag for messages carrying the size of a range
    static constexpr int tag_size = 3;
    /// Communicator used by this remote lambda
    Dist::Comm _comm;

  public:
    /// Constructor
    explicit VectorRemoteLambda(const Dist::Comm& comm) : _comm(comm.comm_dup())
    {
    }

    /**
     * \brief Remote lambda call
     *
     * \tparam RankRange_
     * Range of `int`s, must have `begin` and `end` methods.
     * \tparam IndexRange_
     * Range of `Index`, must have `begin` and `end` methods.
     * \tparam ArgsRange_ Range of argument tuples
     * Range of `ArgumentTuples`s, must have `begin` and `end` methods.
     * \tparam Lambda_
     * Type of lambda function. Must be compatible with \c ReturnType and \c ArgumentTuple
     *
     * \param[in] ranks
     * Range of ranks, for each entry a remote call will be made to that rank
     * \param[in] indices
     * Range of indices, must have same size as \c ranks and contain valid indices into the \c args range
     * \param[in] args
     * Range of argument tuples
     * \param[in] lambda
     * Lambda function to call on this rank
     *
     * \returns A vector containing the return values of all calls described by `ranks`, `indices`, and
     * `args`, in the same order.
     *
     * Assume we have ranges `ranks`, `indices`, and `args`.
     * Assume `result = remote_lambda.call(ranks, indices, args, lambda)` for some `remote_lambda` and some `lambda`.
     * Let `n` be the size of ranks. Then for 0 <= i < n:
     * `result[i] = lambda(args[indices[i]])`
     */
    template<typename RankRange_, typename IndexRange_, typename ArgsRange_, typename Lambda_>
    std::vector<ReturnType>
    call(const RankRange_& ranks, const IndexRange_& indices, const ArgsRange_& args, Lambda_ lambda)
    {
      XASSERTM(
        std::distance(ranks.begin(), ranks.end()) == std::distance(indices.begin(), indices.end()),
        "size mismatch between ranks and indices");

      const auto num_calls = (std::size_t)std::distance(ranks.begin(), ranks.end());

      std::vector<ReturnType> result;
      result.reserve(num_calls);

      auto rank_iter = ranks.begin();
      auto index_iter = indices.begin();
      while(rank_iter != ranks.end())
      {
        const Index idx = *index_iter++;
        const int rank = *rank_iter++;
        const ArgumentTuple& arg = Intern::nth_of_range(args, idx);

        XASSERTM(rank == 0, "invalid rank for non-mpi RemoteLambda");

        result.push_back(std::apply(lambda, arg));
      }

      return result;
    }
  };
#endif
} // namespace FEAT::Global
