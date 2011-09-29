/* GENERAL_REMARK_BY_HILMAR:
 * For the case a more general 'MPI functionality class' is written (see my comment in comm.hpp), then maybe some of
 * the functionality implemented here will move to that class.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_UTIL_MPI_UTILS_HPP
#define KERNEL_UTIL_MPI_UTILS_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/string_utils.hpp>
#include <kernel/util/pretty_printer.hpp>
#include <kernel/util/exception.hpp>

// includes, system
#include <stdlib.h>
#include <string>
#include <mpi.h>

/**
* \file collection of various MPI utilities
*
* \note Following the FEAST naming convention, the functions defined here should actually be called mpi_abort,
* mpi_init, etc. (the distinguishing part of the name has to be appended to the common part). But to avoid confusion
* with the original MPI functions, we deliberately break this rule here.
*
* \author Hilmar Wobker
* \author Dirk Ribbrock
*/

namespace FEAST
{
  /**
  * \brief aborts the program and the MPI universe
  *
  * This function aborts the program and especially blackholes the MPI universe.
  *
  * \param[in] msg
  * message explaining the reason for the abortion (default "")
  */
  void abort_mpi(const std::string& msg = "")
  {
    CONTEXT("abort_mpi()");
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
  void validate_error_code_mpi(
    int error_code,
    std::string mpi_function_name)
  {
    CONTEXT("validate_error_code_mpi()");
    if (error_code != MPI_SUCCESS)
    {
      abort_mpi("Function " + mpi_function_name + " failed with error code " + stringify(error_code) + ".");
    }
  }


  /**
  * \brief initialises MPI
  *
  * This function only calls MPI_Init(...). It has to be called before anything else happens.
  *
  * \note If FEAST's test system is used, it calls MPI_Init(). To prevent calling the function twice,
  * we inquire if MPI is already initialised.
  */
  void init_mpi()
  {
    CONTEXT("init_mpi()");

    // inquire if MPI is already finalised
    int mpi_is_finalised;
    MPI_Finalized(&mpi_is_finalised);
    if(mpi_is_finalised)
    {
      throw InternalError("MPI cannot be initialised, since it has already been finalised.");
    }

    // init MPI
    int mpi_is_initialised;
    MPI_Initialized(&mpi_is_initialised);
    if(!mpi_is_initialised)
    {
      int mpi_error_code = MPI_Init(NULL, NULL);
      validate_error_code_mpi(mpi_error_code, "MPI_Init");
    }
  }


  /**
  * \brief finalises MPI
  *
  * This function only calls MPI_Finalize(...). It has to be called before anything else happens.
  *
  * \note If FEAST's test system is used, it needs to call MPI_Finalize(). To prevent calling the function twice,
  * we inquire if MPI is already finalised.
  */
  void finalise_mpi()
  {
    CONTEXT("finalise_mpi()");
    // inquire if MPI is initialised
    int mpi_is_initialised;
    MPI_Initialized(&mpi_is_initialised);
    if(!mpi_is_initialised)
    {
      throw InternalError("MPI cannot be finalised, since it has not been initialised.");
    }
    // finalise MPI
    int mpi_is_finalised;
    MPI_Finalized(&mpi_is_finalised);
    if(!mpi_is_finalised)
    {
      int mpi_error_code = MPI_Finalize();
      validate_error_code_mpi(mpi_error_code, "MPI_Finalize");
    }
  }


  /**
  * \brief initialises MPI
  *
  * This function only calls MPI_Init(...). It has to be called before anything else happens.
  *
  * \param[in] argc
  * argument count passed to the main() method
  *
  * \param[in] argv
  * arguments passed to the main() method
  */
  void init_mpi(
    int argc,
    char* argv[])
  {
    CONTEXT("init_mpi()");
    // init MPI
    int mpi_error_code = MPI_Init(&argc, &argv);
    validate_error_code_mpi(mpi_error_code, "MPI_Init");
  }


  /**
  * \brief template class for translating built-in types into MPI types (e.g. int to MPI_INTEGER)
  *
  * Empty class definition, to be specialised w.r.t. type T_.
  *
  * \tparam T_
  * type to be translated into MPI type
  *
  * \author Dirk Ribbrock
  * \author Hilmar Wobker
  */
  template<typename T_>
  class MPIType
  {
  };

  /// template specialisation: bool to MPI_LOGICAL
  template<>
  class MPIType<bool>
  {
  public:
    /// calls an MPI_Abort() since there is no logical/boolean MPI datatype.
    static inline MPI_Datatype value()
    {
      // call abort_mpi() (error handler is not available here)
      abort_mpi("Don't try to translate bool into an MPI datatype! There is no such type!");
      return MPI_LOGICAL;
    }
  };

  /// template specialisation: char to MPI_CHARACTER
  template<>
  class MPIType<char>
  {
  public:
    /// translates char to MPI_CHARACTER
    static inline MPI_Datatype value()
    {
      return MPI_CHARACTER;
    }
  };

  /// template specialisation: short int to MPI_SHORT
  template<>
  class MPIType<short int>
  {
  public:
    /// translates short int to MPI_SHORT
    static inline MPI_Datatype value()
    {
      return MPI_SHORT;
    }
  };

  /// template specialisation: unsigned short int to MPI_UNSIGNED_SHORT
  template<>
  class MPIType<unsigned short int>
  {
  public:
    /// translates unsigned short int to MPI_UNSIGNED_SHORT
    static inline MPI_Datatype value()
    {
      return MPI_UNSIGNED_SHORT;
    }
  };

  /// template specialisation: int to MPI_INTEGER
  template<>
  class MPIType<int>
  {
  public:
    /// translates int to MPI_INTEGER
    static inline MPI_Datatype value()
    {
      return MPI_INTEGER;
    }
  };

  /// template specialisation: unsigned int to MPI_UNSIGNED
  template<>
  class MPIType<unsigned int>
  {
  public:
    /// translates unsigned int to MPI_UNSIGNED
    static inline MPI_Datatype value()
    {
      return MPI_UNSIGNED;
    }
  };

  /// template specialisation: long to MPI_LONG
  template<>
  class MPIType<long>
  {
  public:
    /// translates long to MPI_LONG
    static inline MPI_Datatype value()
    {
      return MPI_LONG;
    }
  };

  /// template specialisation: unsigned long to MPI_UNSIGNED_LONG
  template<>
  class MPIType<unsigned long>
  {
  public:
    /// translates unsigned long to MPI_UNSIGNED_LONG
    static inline MPI_Datatype value()
    {
      return MPI_UNSIGNED_LONG;
    }
  };

  /// template specialisation: long long to MPI_LONG_LONG
  template<>
  class MPIType<long long>
  {
  public:
    /// translates long long to MPI_LONG_LONG
    static inline MPI_Datatype value()
    {
      return MPI_LONG_LONG;
    }
  };

  /// template specialisation: unsigned long long to MPI_UNSIGNED_LONG_LONG
  template<>
  class MPIType<unsigned long long>
  {
  public:
    /// translates unsigned long long to MPI_UNSIGNED_LONG_LONG
    static inline MPI_Datatype value()
    {
      return MPI_UNSIGNED_LONG_LONG;
    }
  };

  /// template specialisation: float to MPI_FLOAT
  template<>
  class MPIType<float>
  {
  public:
    /// translates float to MPI_FLOAT
    static inline MPI_Datatype value()
    {
      return MPI_FLOAT;
    }
  };

  /// template specialisation: double to MPI_DOUBLE
  template<>
  class MPIType<double>
  {
  public:
    /// translates double to MPI_DOUBLE
    static inline MPI_Datatype value()
    {
      return MPI_DOUBLE;
    }
  };

  /// template specialisation: long double to MPI_LONG_DOUBLE
  template<>
  class MPIType<long double>
  {
  public:
    /// translates long double to MPI_LONG_DOUBLE
    static inline MPI_Datatype value()
    {
      return MPI_LONG_DOUBLE;
    }
  };
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_UTIL_MPI_UTILS_HPP
