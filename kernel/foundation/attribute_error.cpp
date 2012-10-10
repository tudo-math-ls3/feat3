#include <kernel/foundation/attribute_error.hpp>
#include <string>

using namespace FEAST;
using namespace Foundation;

AttributeError::AttributeError(const std::string & message) throw () :
    Exception(message)
{
}

AttributeTypeMismatch::AttributeTypeMismatch() throw () :
    AttributeError("Type mismatch in attribute")
{
}
