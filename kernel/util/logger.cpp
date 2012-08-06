#include <kernel/util/logger.hpp>

namespace FEAST
{
  std::ofstream Logger::_stream[Logger::max_files];

#ifdef PARALLEL
  Logger::MasterSender* Logger::_master_sender = nullptr;
#endif // PARALLEL

} // namespace FEAST