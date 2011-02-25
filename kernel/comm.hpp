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
#include <kernel/util/exception.hpp>
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
    static unsigned int MCW_BUFFERSIZE;
    /// current position in buffer Comm::MCW_buffer
    static int MCW_buffer_pos;
    /// current size of buffer Comm::MCW_buffer
    static int MCW_received_bytes;

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
      CONTEXT("Comm::init()");
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
      CONTEXT("Comm::send()");
      // send message
      int mpi_error_code = MPI_Send(Comm::MCW_buffer, Comm::MCW_buffer_pos, MPI_PACKED, Process::rank_master, 0,
                                    MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Send");
    }


    /**
    * \brief write one item of template type T_ to the MPI_COMM_WORLD buffer
    *
    * \param[in] msg
    * item of template type T_ to be written to the buffer
    */
    template <typename T_>
    static void write(T_ msg)
    {
      CONTEXT("Comm::write()");
      // write the integer to the current position of the buffer
      int mpi_error_code = MPI_Pack(&msg, 1, MPIType<T_>::value(), Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }


    /**
    * \brief overloaded function for writing a string to the MPI_COMM_WORLD buffer
    *
    * Since strings are always sent as char arrays, the string variant cannot be realised via the template version
    * of the write function above.
    *
    * \param[in] msg
    * string to be written to the buffer
    */
    static void write(std::string const &msg)
    {
      CONTEXT("Comm::write()");
      // Write the string as char array to the current position of the buffer (size + 1, since the null termination
      // symbol is appended). The const of the resulting char array has to be cast away.
      int mpi_error_code = MPI_Pack(const_cast<char *>(msg.c_str()), msg.size()+1, MPI_CHARACTER, Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }


    /**
    * \brief read one integer from the MPI_COMM_WORLD buffer
    *
    * \param[out] msg
    * integer read from the buffer
    */
    template <typename T_>
    static void read(T_& msg)
    {
      CONTEXT("Comm::read()");
      int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, &msg,
                                      1, MPIType<T_>::value(), MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");
    }


    /**
    * \brief read a string from the MPI_COMM_WORLD buffer
    *
    * Since strings are always sent as char arrays, one has to know the size in advance. Hence, the string variant
    * cannot be realised via the template version of the read function above.
    *
    * \param[in] size
    * size of the string to be read
    *
    * \param[out] msg
    * string read from the buffer
    */
    static void read(unsigned int const size, std::string& msg)
    {
      CONTEXT("Comm::read()");
      // create char array for storing the sent message (add 1 to the size since the sent message ends with null
      // termination symbol)
      char* msg_char = new char[size+1];
      int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, msg_char,
                                      size, MPI_CHAR, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");
      // copy the char array to the string
      msg.assign(msg_char);
      delete [] msg_char;
    }


    /**
    * \brief write an array of items of template type T_ to the MPI_COMM_WORLD buffer
    *
    * \param[in] msg
    * array of items of template type T_ to be written to the buffer
    *
    * \param[in] size
    * size of the array
    */
    template <typename T_>
    static void write(T_ msg[], int const size)
    {
      CONTEXT("Comm::write()");
      // write the array to the current position of the buffer
      int mpi_error_code = MPI_Pack(msg, size, MPIType<T_>::value(), Comm::MCW_buffer,
                                    Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
    }


    /**
    * \brief read an array of items of template type T_ from the MPI_COMM_WORLD buffer
    *
    * \param[in] size
    * size of the array to be read
    *
    * \param[out] msg
    * array of items of template type T_ read from the buffer (has to be allocated before function call)
    */
    template <typename T_>
    static void read(int const size, T_ msg[])
    {
      CONTEXT("Comm::read()");
      int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, msg,
                                      size, MPIType<T_>::value(), MPI_COMM_WORLD);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");
    }
  }; // class Comm

  // COMMENT_HILMAR: JUST TEMPORARILY
  // initialisation of static members
  // COMMENT_HILMAR: Use some arbitrary size for the time being. This has to be parameterised somehow...
  unsigned int Comm::MCW_BUFFERSIZE = 4194304;
  char* Comm::MCW_buffer = new char[MCW_BUFFERSIZE];
  int Comm::MCW_buffer_pos = 0;
  int Comm::MCW_received_bytes = 0;
  // COMMENT_HILMAR: JUST TEMPORARILY

} // namespace FEAST

#endif // guard KERNEL_COMM_HPP
