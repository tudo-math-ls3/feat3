/* GENERAL_REMARK_BY_HILMAR:
 * We did not really decide yet, how we deal with global constants... so this file was intended to be the preliminary
 * place where to collect global constants. Either continue using it, or move the constants already defined here to a
 * more appropriate place.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_CONSTANTS_HPP
#define KERNEL_CONSTANTS_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#include <limits> // for numeric_limits

namespace FEAST
{
  /// providing global constants
  class Constants
  {

  public:

    /// global constant: max counting number for an index (if a given index evaluates to this value then it is not valid)
    static const index_glob_t MAX_INDEX;

    /// global constant: max counting number for a number (if a given number evaluates to this value then it is not valid)
    static const index_glob_t MAX_NUMBER;
  };

  /// global constant: max counting number for an index (if a given index evaluates to this value then it is not valid)
  const index_glob_t Constants::MAX_INDEX(std::numeric_limits<index_glob_t>::max());

  /// global constant: max counting number for a number (if a given number evaluates to this value then it is not valid)
  const index_glob_t Constants::MAX_NUMBER(std::numeric_limits<index_glob_t>::max());

} // namespace FEAST


#endif // KERNEL_CONSTANTS_HPP
