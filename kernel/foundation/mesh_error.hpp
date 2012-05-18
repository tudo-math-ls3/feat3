#pragma once
#ifndef KERNEL_FOUNDATION_MESH_ERROR_HH
#define KERNEL_FOUNDATION_MESH_ERROR_HH 1

#include <kernel/util/exception.hpp>
#include <string>

namespace FEAST
{
  class MeshError :
    public Exception
  {
    public:
      MeshError(const std::string & message) throw ();
  };

  class MeshInternalIndexOutOfBounds :
    public MeshError
  {
    public:
      MeshInternalIndexOutOfBounds(unsigned long index, unsigned long max_index) throw ();
  };
}

#endif
