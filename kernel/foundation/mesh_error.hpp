#pragma once
#ifndef KERNEL_FOUNDATION_MESH_ERROR_HH
#define KERNEL_FOUNDATION_MESH_ERROR_HH 1

#include <kernel/util/exception.hpp>
#include <string>

namespace FEAST
{
  namespace Foundation
  {
    class MeshError :
      public Exception
    {
      public:
        MeshError(const std::string & message_in) throw ();
    };

    class MeshInternalIndexOutOfBounds :
      public MeshError
    {
      public:
        MeshInternalIndexOutOfBounds(Index index, Index max_index) throw ();
    };
  }
}

#endif
