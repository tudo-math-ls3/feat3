#pragma once
#ifndef UTIL_GUARD_ENVIRONMENT_HPP
#define UTIL_GUARD_ENVIRONMENT_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/util/comm_base.hpp>

#ifndef SERIAL
#include<mpi.h>
#endif

namespace FEAST
{
  class Environment
  {
    public:
      static Index reserve_tag()
      {
#ifndef SERIAL
        void* max_tag_loc;
        int flag(0);
        MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &max_tag_loc, &flag);
        int max_tag(*(int*)max_tag_loc);
        if(!flag)
          throw InternalError("Environment is unable to retrieve TAG_UB!");
        if(Index(max_tag) < MIN_TAG_UB)
          throw InternalError("Environment gets too small value for TAG_UB!");
        Index result(tags_reserved_low % Index(max_tag) == 0 ? tags_reserved_low : Index(0));
#else
        Index result(0);
#endif
        ++tags_reserved_low;
        return result;
      }

    private:
      static const Index MIN_TAG_UB;
      static Index tags_reserved_low;
  };

}

#endif
