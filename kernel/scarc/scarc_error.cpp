#include <kernel/scarc/scarc_error.hpp>
#include <string>

using namespace FEAST;
using namespace ScaRC;

ScaRCError::ScaRCError(const std::string & message) throw () :
    Exception(message)
{
}
