#include <kernel/foundation/mesh_error.hpp>
#include <string>

using namespace FEAST;
using namespace Foundation;

MeshError::MeshError(const std::string & message_in) :
    Exception(message_in)
{
}

MeshInternalIndexOutOfBounds::MeshInternalIndexOutOfBounds(Index index, Index max_index) :
    MeshError("Internal index '" + stringify(index) + "' exceeds max index '"
            + stringify(max_index) + "'")
{
}
