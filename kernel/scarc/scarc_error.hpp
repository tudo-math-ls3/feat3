#pragma once
#ifndef SCARC_GUARD_SCARC_ERROR_HH
#define SCARC_GUARD_SCARC_ERROR_HH 1

#include <kernel/util/exception.hpp>
#include <string>

namespace FEAST
{
  namespace ScaRC
  {
    class ScaRCError :
        public Exception
    {
        public:
            ScaRCError(const std::string & message_in);
    };

  }
}

#endif
