#include <kernel/foundation/communication_error.hpp>
#include <string>

using namespace FEAST;
using namespace Foundation;

CommunicationError::CommunicationError(const std::string & message_in) :
    Exception(message_in)
{
}

CommunicationHaloOverlapMismatch::CommunicationHaloOverlapMismatch(Index comm_overlap, Index halo_overlap) :
    CommunicationError("Comm overlap '" + stringify(comm_overlap) + "' does not match halo overlap '"
            + stringify(halo_overlap) + "'")
{
}
