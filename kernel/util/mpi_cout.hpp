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

  } // namespace Util

} // namespace FEAT


#endif // UTIL_GUARD_MPI_COUT_HPP
