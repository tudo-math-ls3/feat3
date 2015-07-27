#pragma once
#ifndef KERNEL_FEAST_HPP
#define KERNEL_FEAST_HPP 1

#include <kernel/base_header.hpp>

namespace FEAST
{
  /**
   * \brief FEAST initialisation
   *
   * This function performs the basic initialisation of the FEAST library.
   *
   * \attention
   * This function should be the first functional called in an application's
   * \p main function.
   *
   * \param[in] argc, argv
   * The argument parameters of the calling \p main function.
   *
   * \param[out] rank
   * The MPI rank of this process. Is always 0 for serial builds.
   *
   * \param[out] nprocs
   * The total number of MPI processes. Is always 0 for serial builds.
   */
  void initialise(int& argc, char**& argv, int& rank, int& nprocs);

  /**
   * \brief FEAST initialisation
   *
   * This function calls the four-argument initialise function and discards the
   * \p rank and \p nprocs values.
   */
  void initialise(int& argc, char**& argv);

  /**
   * \brief FEAST abortion
   *
   * This function terminates this process and, in a MPI-based run, also
   * all other processes belonging to this group.
   */
  void abort();

  /**
   * \brief FEAST finalisation
   *
   * This function finalises the FEAST library.
   *
   * \attention
   * This function should be the last function called in an application's
   * \p main function.
   *
   * \returns
   * An exit code (<c>EXIT_SUCCESS</c>) that can be returned by the \p main function.
   */
  int finalise();
} // namespace FEAST

#endif // KERNEL_FEAST_HPP
