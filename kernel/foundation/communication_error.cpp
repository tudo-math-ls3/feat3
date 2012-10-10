#include <kernel/foundation/communication_error.hpp>
#include <string>

using namespace FEAST;
using namespace Foundation;

CommunicationError::CommunicationError(const std::string & message) throw () :
    Exception(message)
{
}

CommunicationHaloOverlapMismatch::CommunicationHaloOverlapMismatch(Index comm_overlap, Index halo_overlap) throw () :
    CommunicationError("Comm overlap '" + stringify(comm_overlap) + "' does not match halo overlap '"
            + stringify(halo_overlap) + "'")
{
}
