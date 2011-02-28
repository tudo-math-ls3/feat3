#pragma once
#ifndef KERNEL_CONSTANTS_HPP
#define KERNEL_CONSTANTS_HPP 1

// includes, system
#include <limits> // for numeric_limits

// includes, FEAST
#include <kernel/base_header.hpp>

/**
 * \brief FEAST namespace.
 */
namespace FEAST
{

  /// providing global constants
  class Constants
  {

  public:

    /// global constant: max counting number for an index (if a given index evaluates to this value then it is not valid)
    static const index_t_glob MAX_INDEX;

    /// global constant: max counting number for a number (if a given number evaluates to this value then it is not valid)
    static const index_t_glob MAX_NUMBER;
  };

  /// global constant: max counting number for an index (if a given index evaluates to this value then it is not valid)
  const index_t_glob Constants::MAX_INDEX(std::numeric_limits<index_t_glob>::max());

  /// global constant: max counting number for a number (if a given number evaluates to this value then it is not valid)
  const index_t_glob Constants::MAX_NUMBER(std::numeric_limits<index_t_glob>::max());

} // namespace FEAST


#endif // #define KERNEL_CONSTANTS_HPP
