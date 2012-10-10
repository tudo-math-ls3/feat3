#pragma once
#ifndef FEM_GUARD_FUNCTOR_ERROR_HH
#define FEM_GUARD_FUNCTOR_ERROR_HH 1

#include <kernel/util/exception.hpp>

#include <string>

namespace FEAST
{
  namespace Foundation
  {
    class FunctorError :
        public Exception
    {
        public:
            FunctorError(const std::string & message) throw ();
    };
  }
}

#endif
