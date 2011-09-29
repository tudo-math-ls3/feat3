
// includes, FEAST
#include <kernel/constants.hpp>

// includes, system
#include <limits> // for numeric_limits

namespace FEAST
{
  /// global constant: max counting number for an index (if a given index evaluates to this value then it is not valid)
  const index_glob_t Constants::MAX_INDEX(std::numeric_limits<index_glob_t>::max());

  /// global constant: max counting number for a number (if a given number evaluates to this value then it is not valid)
  const index_glob_t Constants::MAX_NUMBER(std::numeric_limits<index_glob_t>::max());
} // namespace FEAST