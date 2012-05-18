#include <kernel/foundation/mesh_error.hpp>
#include <string>

using namespace FEAST;

MeshError::MeshError(const std::string & message) throw () :
    Exception(message)
{
}

MeshInternalIndexOutOfBounds::MeshInternalIndexOutOfBounds(unsigned long index, unsigned long max_index) throw () :
    MeshError("Internal index '" + stringify(index) + "' exceeds max index '"
            + stringify(max_index) + "'")
{
}
