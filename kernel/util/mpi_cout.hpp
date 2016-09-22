#pragma once
#ifndef UTIL_GUARD_MPI_COUT_HPP
#define UTIL_GUARD_MPI_COUT_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/comm_base.hpp>

#include <functional>

namespace FEAT
{
  namespace Util
  {
    /// Run cout if condition evaluates to true
    void mpi_cout(FEAT::String string,
    std::function<bool (Index, Index)> func =
    [](Index r, Index /*rs*/) -> bool{ return (r == 0); });

    /**
     * \brief Pads mpi_cout output between two given strings to a given width
     *
     * \param[in] s
     * String on the left of the output.
     *
     * \param[in] t
     * String on the right of the output.
     *
     * \param[in] width
     * Width to pad the line to.
     */
    static inline void mpi_cout_pad_line(const String& s, const String& t, const size_t width = 30)
    {
      mpi_cout(s.pad_back(width, '.') + ": " + t + "\n");
    }

    /**
     * \brief Pads mpi_cout output between a string and a stringify-able type to a given width
     *
     * \param[in] s
     * String on the left of the output.
     *
     * \param[in] t
     * Type that gets written to the right of the output through stringify()
     *
     * \param[in] width
     * Width to pad the line to.
     */
    template<typename T_>
    static inline void mpi_cout_pad_line(const String& s, const T_& t, const size_t width = 30)
    {
      mpi_cout(s.pad_back(width, '.') + ": " + stringify(t) + "\n");
    }

    /**
     * \brief Writes an array, ordered by rank
     *
     * \tparam DT_ Data type, could be anything that can be sent through MPI
     *
     * \param[in] val
     * The data array.
     *
     * \param[in] str
     * String the prepend to the output.
     *
     * \param[in] size
     * The number of entries of val to print.
     *
     * \warning To order the output, it is neccessary to send all information to one rank, which then prints
     * everything. This is a slow communication bottleneck and should ONLY be used for debugging purposes.
     *
     * Sometimes, a series of barriers can be used to get the same effect without first gathering all data in one
     * rank, but this is highly dependent on the MPI implementation. This is the most inefficient, but most portable
     * variant.
     *
     */
    template<typename DT_>
    static inline void mpi_cout_ordered(const DT_* val, const String& str= "", Index size = Index(1))
    {
      Index* other_size(new Index[Util::Comm::size()]);

#ifdef FEAT_HAVE_MPI
      DT_* buffer;
      Index max_size(0);

      // Ignore all comm statuses
      CommStatusIgnore st;
      // This is the rank that receives everything
      if(Util::Comm::rank() == Index(0))
      {
        // First, write out our own stuff
        std::cout << Util::Comm::rank() << " " << str << ": ";

        // Enclose within brackets if we write more than one value
        if(size > Index(1))
          std::cout << "[";

        for(Index i(0); i < size; ++i)
          std::cout << " " << val[i];

        if(size > Index(1))
          std::cout << " ]";

        std::cout << std::endl;

        // Get the sizes of what the other ranks want to send
        for(Index r(1); r < Util::Comm::size(); ++r)
        {
          Util::Comm::recv(&(other_size[r]), Index(1), r, st);
          if(other_size[r] > max_size)
            max_size = other_size[r];
        }

        // Now that we know the size we can allocate the buffer for receiving
        buffer = new DT_[max_size];

        // Receive data from other ranks and print it
        for(Index r(1); r < Util::Comm::size(); ++r)
        {

          Util::Comm::recv(buffer, other_size[r], r, st);

          std::cout << r << " " << str << ": ";

          // Enclose within brackets if we write more than one value
          if(size > Index(1))
            std::cout << "[";

          for(Index i(0); i < other_size[r]; ++i)
          {
            std::cout << " " << buffer[i];
          }

          if(size > Index(1))
            std::cout << " ]";

          std::cout << std::endl;

        }

      }
      // All other ranks
      else
      {
        Comm::send(&size, Index(1), Index(0));

        // Because MPI doesn't care for const, we copy our data to the buffer here and send the buffer
        buffer = new DT_[size];
        for(Index i(0); i < size; ++i)
          buffer[i] = val[i];

        Util::Comm::send(buffer, size, Index(0));

      }
      delete[] buffer;
#else
      // Write out our own stuff
      std::cout << Util::Comm::rank() << " " << str << ": ";

      // Enclose within brackets if we write more than one value
      if(size > Index(1))
        std::cout << "[";

      for(Index i(0); i < size; ++i)
        std::cout << " " << val[i];

      if(size > Index(1))
        std::cout << " ]";

      std::cout << std::endl;
#endif

      delete[] other_size;

    } // void mpi_cout_ordered

  } // namespace Util

} // namespace FEAT


#endif // UTIL_GUARD_MPI_COUT_HPP
