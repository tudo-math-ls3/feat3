#include <kernel/foundation/attribute_error.hpp>
#include <string>

using namespace FEAST;
using namespace Foundation;

AttributeError::AttributeError(const std::string & message_in) :
    Exception(message_in)
{
}

AttributeTypeMismatch::AttributeTypeMismatch() :
    AttributeError("Type mismatch in attribute")
{
}
