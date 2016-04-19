#pragma once
#ifndef UTIL_GUARD_MPI_COUT_HPP
#define UTIL_GUARD_MPI_COUT_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/comm_base.hpp>

#include <functional>

namespace FEAST
{
  namespace Util
  {
    /// Run cout if condition evaluates to true
    void mpi_cout(FEAST::String string,
        std::function<bool (Index, Index)> func =
        [](Index r, Index /*rs*/) -> bool{ return (r == 0); });

  } // namespace Util
} // namespace FEAST


#endif // UTIL_GUARD_MPI_COUT_HPP
