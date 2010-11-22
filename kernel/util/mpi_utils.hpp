#pragma once
#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

// includes, system
#include <stdlib.h>
#include <string>
#include <mpi.h>

// includes, Feast
#include <kernel/util/string_utils.hpp>

namespace FEAST
{

  /// collection of various MPI utilities
  class MPIUtils
  {

  public:

    /* *****************
    * member functions *
    *******************/

    /**
    * \brief aborts the program and the MPI universe
    *
    * This function aborts the program and especially blackholes the MPI universe.
    *
    * \param[in] msg
    * message explaining the reason for the abortion (default "")
    */
    static void abort(const std::string& msg = "")
    {
      // flush cout and cerr
      std::cout.flush();
      std::cerr.flush();
      // print error message to logfile and stderr
      PrettyPrinter pp(40, '#');
      pp.add_line_sep();
      if(msg.size() > 0)
      {
        pp.add_line_no_right_delim(msg);
        pp.add_line_sep();
      }
      else
      {
        pp.add_line_centered("Aborting the MPI universe...");
      }
      pp.add_line_sep();
      pp.print(std::cerr);
      std::cerr.flush();
      // shut down
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

  }; // class MPIUtils


  /**
  * \brief template class for translating built-in types into MPI types (e.g. int to MPI_INTEGER)
  *
  * \author Dirk Ribbrock
  * \author Hilmar Wobker
  */
  template <typename T_>
  class MPIType
  {
  };

  /// template specialisation: bool to MPI_LOGICAL
  template <>
  class MPIType<bool>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_LOGICAL;
    }
  };

  /// template specialisation: char to MPI_CHAR
  template <>
  class MPIType<char>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_CHAR;
    }
  };

  /// template specialisation: short int to MPI_SHORT
  template <>
  class MPIType<short int>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_SHORT;
    }
  };

  /// template specialisation: unsigned short int to MPI_UNSIGNED_SHORT
  template <>
  class MPIType<unsigned short int>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_UNSIGNED_SHORT;
    }
  };

  /// template specialisation: int to MPI_INTEGER
  template <>
  class MPIType<int>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_INTEGER;
    }
  };

  /// template specialisation: unsigned int to MPI_UNSIGNED
  template <>
  class MPIType<unsigned int>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_UNSIGNED;
    }
  };

  /// template specialisation: long to MPI_LONG
  template <>
  class MPIType<long>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_LONG;
    }
  };

  /// template specialisation: unsigned long to MPI_UNSIGNED_LONG
  template <>
  class MPIType<unsigned long>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_UNSIGNED_LONG;
    }
  };

  /// template specialisation: long long to MPI_LONG_LONG
  template <>
  class MPIType<long long>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_LONG_LONG;
    }
  };

  /// template specialisation: unsigned long long to MPI_UNSIGNED_LONG_LONG
  template <>
  class MPIType<unsigned long long>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_UNSIGNED_LONG_LONG;
    }
  };

  /// template specialisation: float to MPI_FLOAT
  template <>
  class MPIType<float>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_FLOAT;
    }
  };

  /// template specialisation: double to MPI_DOUBLE
  template <>
  class MPIType<double>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_DOUBLE;
    }
  };

  /// template specialisation: long double to MPI_LONG_DOUBLE
  template <>
  class MPIType<long double>
  {
  public:
    static inline MPI_Datatype value()
    {
      return MPI_LONG_DOUBLE;
    }
  };
} // namespace FEAST

#endif //  #ifndef MPI_UTILS_HPP
