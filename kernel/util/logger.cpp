#include <kernel/util/logger.hpp>

namespace FEAST
{
  std::ofstream Logger::_stream[Logger::max_files];

#ifndef SERIAL
  Logger::MasterSender* Logger::_master_sender = nullptr;
#endif // SERIAL

} // namespace FEAST
