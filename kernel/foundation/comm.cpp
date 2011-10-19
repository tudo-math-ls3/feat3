#include <kernel/foundation/comm.hpp>
#ifdef PARALLEL

namespace FEAST
{
  // COMMENT_HILMAR: JUST TEMPORARILY
  // initialisation of static members
  // COMMENT_HILMAR: Use some arbitrary size for the time being. This has to be parameterised somehow...
  unsigned int Comm::MCW_BUFFERSIZE = 4194304;
  char* Comm::MCW_buffer = new char[MCW_BUFFERSIZE];
  int Comm::MCW_buffer_pos = 0;
  int Comm::MCW_received_bytes = 0;
  // COMMENT_HILMAR: JUST TEMPORARILY

} // namespace FEAST

#endif // PARALLEL
