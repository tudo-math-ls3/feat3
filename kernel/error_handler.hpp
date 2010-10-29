#pragma once
#ifndef UTIL_ERROR_HANDLER_HHP
/// Header guard
#define UTIL_ERROR_HANDLER_HHP 1

// includes, system

// includes, Feast
#include <kernel/util/string_utils.hpp>
#include <kernel/base_header.hpp>
#include <kernel/logger.hpp>

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
    /// function reacting to an exception
    void exception_occured(Exception e, bool critical, bool send_to_master) const
    {
      PrettyPrinter pp(40, '#');
      pp.add_line_sep();
      if(critical)
      {
        pp.add_line_centered("CRITICAL exception occured!");
        // ... further stuff
      }
      else
      {
        pp.add_line_centered("NON-CRITICAL exception occured!");
        // ... further stuff
      }
      pp.add_line_sep();
      // ... further stuff
      pp.print(Logger::file);
    }
  };

} // namespace FEAST

#endif //UTIL_ERROR_HANDLER_HPP
