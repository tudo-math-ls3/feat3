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

//wird im code dann aufgerufen per MPIType<DT_>::value()
//mit DT_ halt grad der DataType nachdem du deine Funktion templatisierst hast.
//
//Hm fuer ohne internet einfach mal ausgeschnitten:
//
//  template <typename DT_>
//  class MPIType
//  {
//  };
//
//  template <>
//  class MPIType<bool>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_LOGICAL;
//// Or MPI_C_BOOL?
//    }
//  };
//
//  template <>
//  class MPIType<char>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_CHAR;
//    }
//  };
//
//  template <>
//  class MPIType<short int>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_SHORT;
//    }
//  };
//
//  template <>
//  class MPIType< short int>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI__SHORT;
//    }
//  };
//
//  template <>
//  class MPIType<int>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_INTEGER;
//    }
//  };
//
//  template <>
//  class MPIType< int>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_;
//    }
//  };
//
//  template <>
//  class MPIType<long>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_LONG;
//    }
//  };
//
//  template <>
//  class MPIType< long>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI__LONG;
//    }
//  };
//
//  template <>
//  class MPIType<long long>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_LONG_LONG;
//    }
//  };
//
//  template <>
//  class MPIType< long long>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI__LONG_LONG;
//    }
//  };
//
//  template <>
//  class MPIType<float>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_FLOAT;
//    }
//  };
//
//  template <>
//  class MPIType<double>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_DOUBLE;
//    }
//  };
//
//  template <>
//  class MPIType<long double>
//  {
//  public:
//    static inline MPI_Datatype value()
//    {
//      return MPI_LONG_DOUBLE;
//    }
//  };
} // namespace FEAST

#endif //  #ifndef MPI_UTILS_HPP
