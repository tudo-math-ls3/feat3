#include <kernel/util/environment.hpp>
#include <string>

using namespace FEAST;

#ifndef SERIAL
    const Index Environment::MIN_TAG_UB = Index(32767); //given by MPI 3 standard
#else
    const Index Environment::MIN_TAG_UB = Index(0); //no op should use tags
#endif
    Index Environment::tags_reserved_low = Index(0);
