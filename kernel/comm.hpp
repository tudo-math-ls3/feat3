#pragma once
#ifndef KERNEL_COMM_HPP
#define KERNEL_COMM_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/process.hpp>
#include <kernel/service_ids.hpp>


namespace FEAST
{

  /**
  * \brief class providing basic communication functions
  *
  * \author Hilmar Wobker
  */
  class Comm
  {

  private:

    /// MPI status object needed in MPI_Recv call
    MPI_Status status;

  public:

    /* *****************
    * member variables *
    *******************/
    /// global buffer used for MPI_COMM_WORLD (MCW) communication
    static char* MCW_buffer;
    /// size of buffer Comm::MCW_buffer in bytes
    static int MCW_BUFFERSIZE;
    /// current position in buffer Comm::MCW_buffer
    static int MCW_buffer_pos;
    /// current size of buffer Comm::MCW_buffer
    static int MCW_received_bytes;

    /* *************************
    * constructor & destructor *
    ***************************/

    /* *****************
    * member functions *
    *******************/

    /**
    * \brief init a new MPI_COMM_WORLD message
    *
    * This function resets the MPI_COMM_WORLD buffer #MCW_buffer and writes the message ID to its first position.
    *
    * \param[in] service_id
    * ID of the message, where only values of the enumeration ServiceIDs::service_id are accepted.
    */
    static void init(ServiceIDs::service_id service_id)
    {
      // reset buffer
      Comm::MCW_buffer_pos = 0;
      // write the message id to the buffer
      int mpi_error_code = MPI_Pack(&service_id, 1, MPI_INTEGER, Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }

    /// send the current MPI_COMM_WORLD buffer #MCW_buffer to the master
    static void send()
    {
      // send message
      int mpi_error_code = MPI_Send(Comm::MCW_buffer, Comm::MCW_buffer_pos, MPI_PACKED, Process::rank_master, 0,
                                    MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Send");
    }

  //TODO: use templates for the following functions?

    /**
    * \brief write one integer to the MPI_COMM_WORLD buffer
    *
    * \param[in] msg
    * integer to be written to the buffer
    */
    static void write(int msg)
    {
      // write the integer to the current position of the buffer
      int mpi_error_code = MPI_Pack(&msg, 1, MPI_INTEGER, Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }

    /**
    * \brief read one integer from the MPI_COMM_WORLD buffer
    *
    * \param[out] msg
    * integer read from the buffer
    */
    static void read(int& msg)
    {
      int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, &msg,
                                      1, MPI_INTEGER, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");
    }

    /**
    * \brief write an array of integers to the MPI_COMM_WORLD buffer
    *
    * \param[in] msg
    * integer array to be written to the buffer
    *
    * \param[in] size
    * size of the array
    */
    static void write(int msg[], int const size)
    {
      // write the array to the current position of the buffer
      int mpi_error_code = MPI_Pack(msg, size, MPI_INTEGER, Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }

    /**
    * \brief read an array of integers from the MPI_COMM_WORLD buffer
    *
    * \param[in] size
    * size of the array to be read
    *
    * \param[out] msg
    * array of integers read from the buffer (has to be allocated before function call)
    */
    static void read(int const size, int msg[])
    {
      int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, msg,
                                  size, MPI_INTEGER, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");
    }

    /**
    * \brief write one char to the MPI_COMM_WORLD buffer
    *
    * \param[in] msg
    * integer to be written to the buffer
    */
    static void write(char msg)
    {
      // write the char to the current position of the buffer
      int mpi_error_code = MPI_Pack(&msg, 1, MPI_CHARACTER, Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }

    /**
    * \brief read one char from the MPI_COMM_WORLD buffer
    *
    * \param[out] msg
    * char read from the buffer
    */
    static void read(char& msg)
    {
      int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, &msg,
                                      1, MPI_CHAR, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");
    }

    /**
    * \brief write an array of chars to the MPI_COMM_WORLD buffer
    *
    * \param[in] msg
    * char array to be written to the buffer
    *
    * \param[in] size
    * size of the array
    */
    static void write(char msg[], int const size)
    {
      // write the array to the current position of the buffer
      int mpi_error_code = MPI_Pack(msg, size, MPI_CHAR, Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }

    /**
    * \brief read an array of chars from the MPI_COMM_WORLD buffer
    *
    * \param[in] size
    * size of the array to be read
    *
    * \param[out] msg
    * array of chars read from the buffer (has to be allocated before function call)
    */
    static void read(int const size, char msg[])
    {
      int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, msg,
                                      size, MPI_CHAR, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");
    }

    /**
    * \brief write a string to the MPI_COMM_WORLD buffer
    *
    * \param[in] msg
    * string to be written to the buffer
    */
    static void write(std::string const &msg)
    {
      // Write the string as char array to the current position of the buffer (size + 1, since the null termination
      // symbol is appended). The const of the resulting char array has to be cast away.
      int mpi_error_code = MPI_Pack(const_cast<char *>(msg.c_str()), msg.size()+1, MPI_CHARACTER, Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }

    /**
    * \brief read a string from the MPI_COMM_WORLD buffer
    *
    * \param[in] size
    * size of the string to be read
    *
    * \param[out] msg
    * string read from the buffer
    */
    static void read(int const size, std::string& msg)
    {
      // create char array for storing the sent message (add 1 to the size since the sent message ends with null
      // termination symbol)
      char msg_char[size+1];
      int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, msg_char,
                                      size, MPI_CHAR, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");
      // copy the char array to the string
      msg.assign(msg_char);
    }

  // TODO: implement write(...)/read(...) for single + double

  }; // class Comm

  // COMMENT_HILMAR: JUST TEMPORARILY
  // initialisation of static members
  // COMMENT_HILMAR: Use some arbitrary size for the time being. This has to be parameterised somehow...
  int Comm::MCW_BUFFERSIZE = 4194304;
  char* Comm::MCW_buffer = new char[MCW_BUFFERSIZE];
  int Comm::MCW_buffer_pos = 0;
  int Comm::MCW_received_bytes = 0;
  // COMMENT_HILMAR: JUST TEMPORARILY

} // namespace FEAST

#endif // guard KERNEL_COMM_HPP
