#include <kernel/foundation/attribute_error.hpp>
#include <string>

using namespace FEAST;
using namespace Foundation;

AttributeError::AttributeError(const std::string & message_in) throw () :
    Exception(message_in)
{
}

AttributeTypeMismatch::AttributeTypeMismatch() throw () :
    AttributeError("Type mismatch in attribute")
{
}
