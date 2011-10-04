/* GENERAL_REMARK_BY_HILMAR:
 * This class is meant to provide auxiliary functionality for reacting on errors, exceptions and similar. "React" means:
 * Processing status objects (which do not exist yet; see issue 00028), performing corresponding file/screen output,
 * eventually aborting the program. Until now, I only realised some sort of "exception handler" (which simply
 * performs pretty printing of the exception message and aborts the program). Maybe "ErrorHandler" is not the best
 * name for the class, since it should also process warnings and exceptions... maybe you find a better name.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_ERROR_HANDLER_HHP
#define KERNEL_ERROR_HANDLER_HHP 1

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/logger.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/pretty_printer.hpp>
#include <kernel/util/exception.hpp>

// includes, system

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
    /// function reacting to error/warnings
    static void error_occured(StatusObject status)
    {
      // create pretty printed error message with the prefix Error which can be grepped for
      PrettyPrinter pp(40, '#', "ERROR ");
      pp.add_line_sep();
      if(status.severity == CRITICAL)
      {
        pp.add_line_centered("Error occured on process " + stringify(Process::rank) + "!");
      }
      else
      {
        pp.add_line_centered("Warning occured on process " + stringify(Process::rank) + "!");
      }
      pp.add_line_sep();
      pp.add_line_no_right_delim(status.message);
      pp.add_line_sep();
      pp.print(Logger::file);

      // If the error is critical, then write the error to std::cerr and abort the program.
      if(sev == CRITICAL)
      {
        pp.print(std::cerr);
        abort_mpi();
      }
    }
*/

    /// function reacting to an exception
    static void exception_occured(Exception const& e)
    {
      // create pretty printed error message with the prefix EXCEPTION which can be grepped for
      PrettyPrinter pp(40, '#', "EXCEPTION ");
      pp.add_line_sep();
#ifdef PARALLEL
      pp.add_line_centered("Exception occured on process " + stringify(Process::rank) + "!");
#endif
      pp.add_line_sep();
      pp.add_line_no_right_delim(e.message());
      pp.add_line_sep();
#ifndef FEAST_NO_CONTEXT
      pp.add_line("Backtrace:");
      pp.add_line(e.backtrace("\nEXCEPTION # "));
      pp.add_line_sep();
#endif
      //pp.print(Logger::file);
      Logger::log(pp.block());

      // An exception always leads to program abortion via MPI_Abort(). Its message is written to stderr. If the error
      // occurs simultaneously on several processes, the error message appears several times on stderr. There is no
      // way to prevent this.
      pp.print(std::cerr);
      abort("Unhandled exception");
    }
  };
} // namespace FEAST

#endif // KERNEL_UTIL_ERROR_HANDLER_HPP
