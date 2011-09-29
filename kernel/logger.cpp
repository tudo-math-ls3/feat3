#include <kernel/logger.hpp>

namespace FEAST
{
  std::string Logger::file_name;
  std::string Logger::file_base_name;
  std::ofstream Logger::file;
} // namespace FEAST