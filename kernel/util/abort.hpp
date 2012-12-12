#pragma once
#ifndef KERNEL_UTIL_ABORT_HPP
#define KERNEL_UTIL_ABORT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/pretty_printer.hpp>

#ifndef SERIAL
#include <mpi.h>
#endif // SERIAL

namespace FEAST
{
  /**
   * \brief Aborts program execution.
   *
   * \note In a parallel simulation this function will not only abort the calling process but the whole universe.
   *
   * \param[in] message
   * A message describing the reason for the program abortion.
   *
   * \param[in] exit_code
   * The exit code which is to be passed to the exit() function.
   *
   * \author Hilmar Wobker
   * \author Peter Zajac
   */
  inline void abort(
    const String& message,
    int exit_code = -1)
  {
    CONTEXT("abort()");

    // flush cout and cerr
    std::cout.flush();
    std::cerr.flush();

    // print error message to logfile and stderr
    PrettyPrinter pp(40, '#');
    pp.add_line_sep();
    if(!message.empty())
    {
      pp.add_line_no_right_delim(message);
      pp.add_line_sep();
    }
    else
    {
      pp.add_line_centered("Aborting program execution...");
    }
    pp.add_line_sep();
    pp.print(std::cerr);
    std::cerr.flush();

#ifndef SERIAL
    // shut down
    int mpi_is_initialised;
    MPI_Initialized(&mpi_is_initialised);
    if (mpi_is_initialised)
    {
      // TODO: Mapping to Feast error codes like in FEAST1? [dom 25.8.2010]
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif // SERIAL

    // abort program execution
    exit(exit_code);
  }
} // namespace FEAST

#endif // KERNEL_UTIL_ABORT_HPP
