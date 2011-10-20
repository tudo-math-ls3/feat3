#include <kernel/util/logger.hpp>

namespace FEAST
{
#ifdef OLD_LOGGER
  std::string Logger::file_name;
  std::string Logger::file_base_name;
  std::ofstream Logger::file;
#else

  std::ofstream Logger::_stream[Logger::max_files];

#ifdef PARALLEL
  Logger::MasterSender* Logger::_master_sender = nullptr;
#endif // PARALLEL

#endif
} // namespace FEAST