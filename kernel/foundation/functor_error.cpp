#include <kernel/foundation/functor_error.hpp>
#include <string>

using namespace FEAST;

FunctorError::FunctorError(const std::string & message) throw () :
    Exception(message)
{
}
