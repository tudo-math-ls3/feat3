#pragma once
#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#include <string>

/**
 * collection of various MPI utilities
 */
class MPIUtils
{
  /* ****************
   * public members *
   ******************/
   public:
    /* ******************
     * member functions *
     ********************/

    /**
     * \brief aborts the program and the MPI universe
     *
     * This function aborts the program and especially blackholes the MPI universe.
     *
     * \param[in] msg
     * message explaining the reason for the abortion
     */
    static void abort(std::string msg)
    {
      std::cerr << msg << " Aborting program..." << std::endl;
      int mpi_is_initialised;
      MPI_Initialized(&mpi_is_initialised);
      if (mpi_is_initialised)
      {
        // TODO: Mapping to Feast error codes like in FEAST1? [dom 25.8.2010]
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      exit(1);
    }


    /**
     * \brief evaluates the return value of MPI routines
     *
     * Some MPI routines return an integer error code. This routine checks if the error code is MPI_SUCCESS.
     * If not, then it aborts the program.
     *
     * \param[in] error_code
     * MPI error code
     *
     * \param[in] mpi_function_name
     * name of the calling routine
     */
    static void validate_mpi_error_code(
      int error_code,
      std::string mpi_function_name)
    {
      if (error_code != MPI_SUCCESS)
      {
        abort("Function " + mpi_function_name + " failed with error code " + StringUtils::stringify(error_code) + ".");
      }
    }
};
#endif //  #ifndef MPI_UTILS_HPP
