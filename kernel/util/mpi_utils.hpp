#pragma once
#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#include <string>
#include <mpi.h>

/// collection of various MPI utilities
class MPIUtils
{

 public:

  /* *****************
  * member variables *
  *******************/
// COMMENT_HILMAR: JUST TEMPORARILY stored in this class until a more suitable place has been found
  /// global buffer used for COMM_WORLD communication
  static char* buffer;
  /// size of buffer #buffer
  static int BUFFERSIZE_BYTES;
  /// current position in buffer #buffer
  static int buffer_pos;
  /// current size of buffer #buffer
  static int received_bytes;
// COMMENT_HILMAR: JUST TEMPORARILY

  /* *****************
  * member functions *
  *******************/

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

  // init message buffer
//  static void init_msg(int id, MPI_Comm& comm, void* buffer, int buffersize_bytes, int pos_bytes)
//  {
//    int mpi_error_code = mpi_pack(id, 1, MPI_INTEGER, buffer, buffersize_bytes, pos_bytes, comm)
//  }
}; // class MPIUtils

// COMMENT_HILMAR: JUST TEMPORARILY
// initialisation of static members
// COMMENT_HILMAR: Use some arbitrary size for the time being. This has to be parameterised somehow...
int MPIUtils::BUFFERSIZE_BYTES = 4194304;
char* MPIUtils::buffer = new char[BUFFERSIZE_BYTES];
int MPIUtils::buffer_pos = 0;
int MPIUtils::received_bytes = 0;
// COMMENT_HILMAR: JUST TEMPORARILY

#endif //  #ifndef MPI_UTILS_HPP
