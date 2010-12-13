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

  class Constants
  {

  public:

    /// global constant: max counting number for an index (if a given index evaluates to this value then it is not valid)
    static const global_index_t MAX_INDEX;

    /// global constant: max counting number for a number (if a given number evaluates to this value then it is not valid)
    static const global_index_t MAX_NUMBER;
  };

  /// global constant: max counting number for an index (if a given index evaluates to this value then it is not valid)
  const global_index_t Constants::MAX_INDEX(std::numeric_limits<global_index_t>::max());

  /// global constant: max counting number for a number (if a given number evaluates to this value then it is not valid)
  const global_index_t Constants::MAX_NUMBER(std::numeric_limits<global_index_t>::max());

} // namespace FEAST


#endif // #define KERNEL_CONSTANTS_HPP
