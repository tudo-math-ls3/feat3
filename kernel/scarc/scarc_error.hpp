#pragma once
#ifndef FEM_GUARD_SCARC_ERROR_HH
#define FEM_GUARD_SCARC_ERROR_HH 1

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
            ScaRCError(const std::string & message) throw ();
    };

  }
}

#endif
