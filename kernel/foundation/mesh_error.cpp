#include <kernel/foundation/mesh_error.hpp>
#include <string>

using namespace FEAST;
using namespace Foundation;

MeshError::MeshError(const std::string & message) throw () :
    Exception(message)
{
}

MeshInternalIndexOutOfBounds::MeshInternalIndexOutOfBounds(Index index, Index max_index) throw () :
    MeshError("Internal index '" + stringify(index) + "' exceeds max index '"
            + stringify(max_index) + "'")
{
}
