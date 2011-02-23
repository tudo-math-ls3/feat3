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
  * \brief class providing some functions for handling warnings, exceptions and errors
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

// COMMENT_HILMAR: wenn wir beschlieﬂen, dass eine geworfene Exception automatisch zum Programm-Abbruch fuehren soll,
// dann wird das severity-Zeuch hier raugeschmissen.

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
