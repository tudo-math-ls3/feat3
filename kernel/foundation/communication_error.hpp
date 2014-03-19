#pragma once
#ifndef FEM_GUARD_COMMUNICATION_ERROR_HH
#define FEM_GUARD_COMMUNICATION_ERROR_HH 1

#include <kernel/util/exception.hpp>

#include <string>

namespace FEAST
{
  namespace Foundation
  {
    class CommunicationError :
        public Exception
    {
        public:
            CommunicationError(const std::string & message_in) throw ();
    };

    class CommunicationHaloOverlapMismatch :
        public CommunicationError
    {
        public:
            CommunicationHaloOverlapMismatch(Index comm_overlap, Index halo_overlap) throw ();
    };
  }
}

#endif
