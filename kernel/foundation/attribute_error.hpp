#pragma once
#ifndef FEM_GUARD_ATTRIBUTE_ERROR_HH
#define FEM_GUARD_ATTRIBUTE_ERROR_HH 1

#include <kernel/util/exception.hpp>

#include <string>

namespace FEAST
{
  namespace Foundation
  {
    class AttributeError :
        public Exception
    {
        public:
            AttributeError(const std::string & message_in) throw ();
    };

    class AttributeTypeMismatch :
        public AttributeError
    {
        public:
            AttributeTypeMismatch() throw ();
    };
  }
}

#endif
