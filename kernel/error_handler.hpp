#pragma once
#ifndef UTIL_ERROR_HANDLER_HHP
/// Header guard
#define UTIL_ERROR_HANDLER_HHP 1

// includes, system

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/logger.hpp>
#include <kernel/util/exception.hpp>

namespace FEAST
{

  /**
  * \brief class providing some functions for handling warnings, errors and exceptions
  *
  * \todo This class is not finished and not tested yet!
  *
  * \author Hilmar Wobker
  */
  class ErrorHandler
  {
  private:

  public:

    /**
    * \brief severity of the error
    */
    enum severity
    {
      /// a critical error triggers abortion of the program
      CRITICAL,
      /// a non-critical error only triggers output
      NON_CRITICAL
    };

/*
COMMENT_HILMAR: corresponding to the function exception_occured(...), there should be something like the following
function:
    /// function reacting to an error
    static void exception_occured(int error_code, severity sev)
    {
      // create pretty printed error message with the prefix EXCEPTION which can be grepped for
      PrettyPrinter pp(40, '#', "ERROR ");
      pp.add_line_sep();
      if(sev == CRITICAL)
      {
        pp.add_line_centered("Error occured on process " + StringUtils::stringify(Process::rank) + "!");
      }
      else
      {
        pp.add_line_centered("Warning occured on process " + StringUtils::stringify(Process::rank) + "!");
      }
      pp.add_line_sep();
      pp.add_line_no_right_delim(...); // here we need some mapping of error_code to corresponding error message
      pp.add_line_sep();
      pp.print(Logger::file);

      // If the error is critical, then write the error to std::cerr and abort the program.
      if(sev == CRITICAL)
      {
        pp.print(std::cerr);
        MPIUtils::abort();
      }
    }

Der error_code wäre dann der return value einer Funktion. error_code = 0 heisst typischerweise "kein Fehler".
Die Frage ist, wie man ein geschicktes Mapping zwischen error_code und error_message hinbekommt.
Einfach ein statisches String-Array mit error_code = array index?
*/

    /// function reacting to an exception
    static void exception_occured(Exception const& e)
    {
      // create pretty printed error message with the prefix EXCEPTION which can be grepped for
      PrettyPrinter pp(40, '#', "EXCEPTION ");
      pp.add_line_sep();
      pp.add_line_centered("Exception occured on process " + StringUtils::stringify(Process::rank) + "!");
      pp.add_line_sep();
      pp.add_line_no_right_delim(e.message());
      pp.add_line_sep();
      pp.print(Logger::file);

      // An exception always leads to program abortion via MPI_Abort(). Its message is written to stderr. If the error
      // occurs simultaneously on several processes, the error message appears several times on stderr. There is no
      // way to prevent this.
      pp.print(std::cerr);
      MPIUtils::abort();
    }
  };
} // namespace FEAST

#endif //UTIL_ERROR_HANDLER_HPP
