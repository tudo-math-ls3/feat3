// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/string.hpp>

#include <array>
#include <numeric>
#include <tuple>
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

    template<typename... Ts_>
    using TypeCache = std::array<const Dist::Datatype*, sizeof...(Ts_)>;

    /// Actual recursive implementation of mpi_datatypes()
    template<std::size_t i = 0, typename... Ts_>
    void _fill_type_cache(TypeCache<Ts_...>& cache)
    {
      if constexpr(i < sizeof...(Ts_))
      {
        using T = std::decay_t<std::remove_reference_t<std::tuple_element_t<i, std::tuple<Ts_...>>>>;
        cache[i] = &Dist::autotype<T>();
        _fill_type_cache<i + 1, Ts_...>(cache);
      }
    }

    /**
     * \brief Retrieve Dist::Datatypes for all types of a parameter pack
     *
     * \tparam Ts_ Types to retrieve Dist::Datatypes for
     *
     * \return An array containing references to the Dist::Datatypes for all requested types
     */
    template<typename... Ts_>
    TypeCache<Ts_...> make_type_cache()
    {
      TypeCache<Ts_...> result;
      _fill_type_cache<0, Ts_...>(result);
      return result;
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
      const TypeCache<Ts_...>& cache,
      void* outbuf,
      std::size_t outsize,
      int* position)
    {
      using T = std::decay_t<std::remove_reference_t<decltype(std::get<i>(tuple))>>;

      const T& t = std::get<i>(tuple);
      comm.pack(&t, 1, *cache[i], outbuf, outsize, position);

      if constexpr(i + 1 != sizeof...(Ts_))
      {
        pack_tuple<i + 1>(comm, tuple, cache, outbuf, outsize, position);
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
      const TypeCache<Ts_...>& cache,
      const void* inbuf,
      std::size_t insize,
      int* position)
    {
      using T = std::decay_t<std::remove_reference_t<decltype(std::get<i>(tuple))>>;

      T& t = std::get<i>(tuple);
      comm.unpack(inbuf, insize, position, &t, 1, *cache[i]);

      if constexpr(i + 1 != sizeof...(Ts_))
      {
        unpack_tuple<i + 1>(comm, tuple, cache, inbuf, insize, position);
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
      // Phase 1 - Exchange number of calls
      ////////////////////////////////////////////

      // NOTE: Arguments are batched. If rank A wants to (in total) perform
      // three lambda calls on rank B, then a single message containing three
      // argument tuples is sent. The number of calls is only relevant to
      // allocate appropriate buffers.

      // Number of incoming lambda calls
      std::vector<int> calls_from(num_ranks, 0);
      // Number of outgoing lambda calls
      std::vector<int> calls_to(num_ranks, 0);
      {
        for(auto i : ranks)
        {
          calls_to[(std::size_t)i] += 1;
        }

        // Sync number of calls across all ranks
        _comm.alltoall(calls_to.data(), 1, calls_from.data(), 1);
      }

      ////////////////////////////////////////////
      // Phase 2 - Create required buffers
      ////////////////////////////////////////////

      // Buffer for sending arguments
      std::vector<std::vector<char>> argument_s_buffers(num_ranks);
      // Buffer for receiving arguments
      std::vector<std::vector<char>> argument_r_buffers(num_ranks);
      // Buffer for sending return values
      std::vector<std::vector<ReturnType>> return_s_buffers(num_ranks);
      // Buffer for receiving return values
      std::vector<std::vector<ReturnType>> return_r_buffers(num_ranks);

      for(Index rank(0); rank < num_ranks; rank++)
      {
        // We send as many arguments to rank r as our local sendcount for rank r
        argument_s_buffers[rank].resize(std::size_t(calls_to[rank]) * pack_size);
        // We receive as many arguments from rank r as our local receivecount for rank r
        argument_r_buffers[rank].resize(std::size_t(calls_from[rank]) * pack_size);

        // We send as many return values to rank r as we received arguments from rank r
        return_s_buffers[rank] = std::vector<ReturnType>(std::size_t(calls_from[rank]));
        // We receive as many return values from rank r as we sent arguments to rank r
        return_r_buffers[rank] = std::vector<ReturnType>(std::size_t(calls_to[rank]));
      }

      ////////////////////////////////////////////
      // Phase 3 - Pack argument buffers
      ////////////////////////////////////////////

      // Packing and unpacking tuples requires access to the Dist::Datatypes of the tuple elements.
      // For non-standard types this requires looking them up in the TypeRegistry, which is behind a mutex.
      // We cache the datatypes once and then re-use them for packing and unpacking each tuple.

      Intern::TypeCache<Args_...> cache = Intern::make_type_cache<Args_...>();

      // argument_buffer_pos[r] is the current offset into argument_r_buffer[r]
      std::vector<int> argument_buffer_pos(num_ranks, 0);
      {
        auto rank_iter = ranks.begin();
        auto index_iter = indices.begin();

        while(rank_iter != ranks.end())
        {
          const Index idx = *index_iter++;
          const Index rank = Index(*rank_iter++);
          const ArgumentTuple& arg = Intern::nth_of_range(args, idx);

          // Pack tuple for sending arg, appends to current buffer
          Intern::pack_tuple(_comm, arg, cache, argument_s_buffers[rank].data(), argument_s_buffers[rank].size(), &argument_buffer_pos[rank]);
        }
      }

      ////////////////////////////////////////////
      // Phase 4 - Send args and receive results
      ////////////////////////////////////////////

      Dist::RequestVector arg_sends;
      Dist::RequestVector result_receives;
      {
        for(Index rank(0); rank < num_ranks; rank++)
        {
          // Send packed argument buffer to intended rank
          arg_sends.push_back(_comm.isend(argument_s_buffers[rank].data(), (std::size_t)argument_buffer_pos[rank], (int)rank, tag_argument));

          // Post receive for results from this rank
          result_receives.push_back(_comm.irecv(return_r_buffers[rank].data(), return_r_buffers[rank].size(), (int)rank, tag_return));
        }
      }

      ////////////////////////////////////////////
      // Phase 5 - Receive args and send results
      ////////////////////////////////////////////

      Dist::RequestVector arg_receives;
      for(Index rank(0); rank < num_ranks; rank++)
      {
        // Post receive for incoming arguments
        arg_receives.push_back(_comm.irecv(argument_r_buffers[rank].data(), argument_r_buffers[rank].size(), (int)rank, tag_argument));
      }

      // For all incoming arguments, handle them in order and write their results into the appropriate send buffer
      Dist::RequestVector result_sends;
      {
        std::size_t idx(0);
        Dist::Status status;

        while(arg_receives.wait_any(idx, status))
        {
          auto rank = Index(status.source());

          int pos(0);
          for(Index arg(0); arg < Index(calls_from[rank]); arg++)
          {
            ArgumentTuple t;
            Intern::unpack_tuple(_comm, t, cache, argument_r_buffers[rank].data(), argument_r_buffers[rank].size(), &pos);

            return_s_buffers[rank][arg] = std::apply(lambda, t);
          }
        }
      }

      // Post sends for return values
      for(Index rank(0); rank < num_ranks; rank++)
      {
        result_sends.push_back(_comm.isend(return_s_buffers[rank].data(), return_s_buffers[rank].size(), (int)rank, tag_return));
      }

      ////////////////////////////////////////////
      // Phase 6 - Merge batches back into final result
      ////////////////////////////////////////////

      arg_sends.wait_all();
      arg_receives.wait_all();
      result_sends.wait_all();
      result_receives.wait_all();

      // We now have the results of all lambda calls, sorted by rank We need to
      // remap these into the order of lambda calls originally given by the
      // user, so that the user knows which return value belongs to which call
      std::vector<ReturnType> result(num_calls);
      std::vector<Index> result_pos(num_ranks, 0);
      {
        auto rank_iter = ranks.begin();

        Index idx(0);
        while(rank_iter != ranks.end())
        {
          const Index rank = Index(*rank_iter++);

          result[idx++] = return_r_buffers[rank][result_pos[rank]++];
        }
      }

      return result;
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

      // NOTE: Arguments are batched. If rank A wants to (in total) perform
      // three lambda calls on rank B, then a single message containing three
      // argument tuples is sent. The number of calls is only relevant to
      // allocate appropriate buffers.

      // Number of incoming lambda calls
      std::vector<int> calls_from(num_ranks, 0);
      // Number of outgoing lambda calls
      std::vector<int> calls_to(num_ranks, 0);
      {
        for(auto i : ranks)
        {
          calls_to[(std::size_t)i] += 1;
        }

        // Sync number of calls across all ranks
        _comm.alltoall(calls_to.data(), 1, calls_from.data(), 1);
      }

      ////////////////////////////////////////////
      // Phase 2 - Create required buffers
      ////////////////////////////////////////////

      std::vector<std::vector<char>> argument_s_buffers(num_ranks);
      std::vector<std::vector<char>> argument_r_buffers(num_ranks);
      std::vector<std::vector<std::size_t>> range_size_s_buffers(num_ranks);
      std::vector<std::vector<std::size_t>> range_size_r_buffers(num_ranks);
      for(Index rank(0); rank < num_ranks; rank++)
      {
        // We send as many arguments to rank r as our local sendcount for rank r
        argument_s_buffers[rank].resize(std::size_t(calls_to[rank]) * pack_size);
        // We receive as many arguments from rank r as our local receivecount for rank r
        argument_r_buffers[rank].resize(std::size_t(calls_from[rank]) * pack_size);

        // We send as many range sizes to rank r as we received arguments from rank r
        range_size_s_buffers[rank].resize(std::size_t(calls_from[rank]));
        // We receive as many range sizes from rank r as we sent argument to rank r
        range_size_r_buffers[rank].resize(std::size_t(calls_to[rank]));
      }

      ////////////////////////////////////////////
      // Phase 3 - Pack argument buffers
      ////////////////////////////////////////////

      // Packing and unpacking tuples requires access to the Dist::Datatypes of the tuple elements.
      // For non-standard types this requires looking them up in the TypeRegistry, which is behind a mutex.
      // We cache the datatypes once and then re-use them for packing and unpacking each tuple.

      Intern::TypeCache<Args_...> cache = Intern::make_type_cache<Args_...>();

      // argument_buffer_pos[r] is the current offset into argument_r_buffer[r]
      std::vector<int> argument_buffer_pos(num_ranks, 0);
      {
        auto rank_iter = ranks.begin();
        auto index_iter = indices.begin();

        while(rank_iter != ranks.end())
        {
          const Index idx = *index_iter++;
          const Index rank = Index(*rank_iter++);
          const ArgumentTuple& arg = Intern::nth_of_range(args, idx);

          // Pack tuple for sending arg
          Intern::pack_tuple(_comm, arg, cache, argument_s_buffers[rank].data(), argument_s_buffers[rank].size(), &argument_buffer_pos[rank]);
        }
      }

      ////////////////////////////////////////////////////
      // Phase 4 - Send args and receive sizes of ranges
      ////////////////////////////////////////////////////

      Dist::RequestVector arg_sends;
      Dist::RequestVector size_receives;
      {
        for(Index rank(0); rank < num_ranks; rank++)
        {
          // Send packed argument buffer to intended rank
          arg_sends.push_back(_comm.isend(argument_s_buffers[rank].data(), (std::size_t)argument_buffer_pos[rank], (int)rank, tag_argument));

          // Post receive for range sizes from this rank
          size_receives.push_back(_comm.irecv(range_size_r_buffers[rank].data(), range_size_r_buffers[rank].size(), (int)rank, tag_size));
        }
      }

      ///////////////////////////////////////////////////
      // Phase 5 - Receive args and send range sizes and ranges
      ///////////////////////////////////////////////////

      Dist::RequestVector arg_receives;
      for(Index rank(0); rank < num_ranks; rank++)
      {
        // Post receive for argument buffer
        arg_receives.push_back(_comm.irecv(argument_r_buffers[rank].data(), argument_r_buffers[rank].size(), (int)rank, tag_argument));
      }

      Dist::RequestVector range_sends;
      Dist::RequestVector size_sends;

      std::vector<std::vector<ValueType>> return_s_buffers(num_ranks);
      {
        std::size_t idx(0);
        Dist::Status status;

        while(arg_receives.wait_any(idx, status))
        {
          auto rank = Index(status.source());

          int pos(0);
          for(Index arg(0); arg < Index(calls_from[rank]); arg++)
          {
            // Unpack argument
            ArgumentTuple t;
            Intern::unpack_tuple(_comm, t, cache, argument_r_buffers[rank].data(), argument_r_buffers[rank].size(), &pos);

            // Create result range
            auto range = std::apply(lambda, t);

            // Record range size so that range can be reconstructed after sending
            range_size_s_buffers[rank][arg] = (std::size_t)std::distance(range.begin(), range.end());

            // Append range values to sendbuffer for that rank
            return_s_buffers[rank].insert(return_s_buffers[rank].end(), range.begin(), range.end());
          }
        }
      }

      for(Index rank(0); rank < num_ranks; rank++)
      {
        // Post sends for range sizes
        size_sends.push_back(_comm.isend(range_size_s_buffers[rank].data(), range_size_s_buffers[rank].size(), (int)rank, tag_size));

        // Post sends for ranges
        range_sends.push_back(_comm.isend(return_s_buffers[rank].data(), return_s_buffers[rank].size(), (int)rank, tag_return));
      }

      ////////////////////////////////////////////
      // Phase 5 - Receive range sizes and ranges
      ////////////////////////////////////////////

      size_sends.wait_all();
      size_receives.wait_all();

      // Allocate return values buffer for each rank
      // Each buffer must be large enough to receive all values of all generated ranges
      std::vector<std::vector<ValueType>> return_r_buffers(num_ranks);
      for(Index rank(0); rank < num_ranks; rank++)
      {
        return_r_buffers[rank].resize(std::accumulate(range_size_r_buffers[rank].begin(), range_size_r_buffers[rank].end(), 0ULL));
      }

      Dist::RequestVector range_receives;
      for(Index rank(0); rank < num_ranks; rank++)
      {
        // Post receives for ranges
        range_receives.push_back(_comm.irecv(return_r_buffers[rank].data(), return_r_buffers[rank].size(), (int)rank, tag_return));
      }

      ////////////////////////////////////////////
      // Phase 6 - Reconstruct ranges
      ////////////////////////////////////////////

      arg_sends.wait_all();
      arg_receives.wait_all();
      range_receives.wait_all();
      range_sends.wait_all();

      // We now have the returned ranges from each rank as one long array, plus the sizes of the individual ranges
      // We now reconstruct the original ranges, by taking subranges according to the range sizes.
      std::vector<std::vector<ReturnType>> ranges(num_ranks);

      for(Index rank(0); rank < num_ranks; rank++)
      {
        Index start(0);

        for(Index arg(0); arg < Index(calls_to[rank]); arg++)
        {
          auto size = range_size_r_buffers[rank][arg];
          auto begin = return_r_buffers[rank].data() + start;
          auto end = return_r_buffers[rank].data() + start + size;
          ranges[rank].emplace_back(begin, end);
          start += size;
        }
      }

      // Same as in the scalar case, we now have the returned ranges sorted by rank
      // We need to reorder them so that they match the order of the original lambda calls,
      // so that the user knows which range belongs to which call.
      std::vector<ReturnType> result(num_calls);
      std::vector<Index> result_pos(num_ranks, 0);
      {
        auto rank_iter = ranks.begin();

        Index idx(0);
        while(rank_iter != ranks.end())
        {
          const Index rank = Index(*rank_iter++);

          result[idx++] = std::move(ranges[rank][result_pos[rank]++]);
        }
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
