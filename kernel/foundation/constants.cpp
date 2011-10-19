
// includes, FEAST
#include <kernel/foundation/constants.hpp>

// includes, system
#include <limits> // for numeric_limits

namespace FEAST
{
  /// global constant: max counting number for an index (if a given index evaluates to this value then it is not valid)
  const Index Constants::MAX_INDEX(std::numeric_limits<Index>::max());

  /// global constant: max counting number for a number (if a given number evaluates to this value then it is not valid)
  const Index Constants::MAX_NUMBER(std::numeric_limits<Index>::max());
} // namespace FEAST