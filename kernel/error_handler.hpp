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

    /// function reacting to an exception
    static void exception_occured(Exception const& e, severity sev)
    {
      // create pretty printed error message with the prefix EXCEPTION which can be grepped for
      PrettyPrinter pp(40, '#', "EXCEPTION ");
      pp.add_line_sep();
      if(sev == CRITICAL)
      {
        pp.add_line_centered("CRITICAL exception occured");
      }
      else
      {
        pp.add_line_centered("NON-CRITICAL exception occured");
      }
      pp.add_line_centered("on process " + StringUtils::stringify(Process::rank) + "!");
      pp.add_line_sep();
      pp.add_line_no_right_delim(e.message());
      pp.add_line_sep();
      pp.print(Logger::file);

      // If the error is critical, it is written to stderr and the program is aborted. If the error occurs
      // simultaneously on several processes, the error message appears several times on stderr. There is no way to
      // prevent this.
      if(sev == CRITICAL)
      {
        pp.print(std::cerr);
        MPIUtils::abort("");
      }
    }
  };
} // namespace FEAST

#endif //UTIL_ERROR_HANDLER_HPP
