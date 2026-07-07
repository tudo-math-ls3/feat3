#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/dist_file_io.hpp>

#include <ctime>

namespace Gendie
{
  enum class PrintLevel : std::uint8_t
  {
    debug = 0,
    verbose,
    info,
    warning,
    error,
    off
  };

  namespace Intern
  {
    inline const char* get_print_char(PrintLevel lvl)
    {
      switch(lvl)
      {
        case PrintLevel::debug:
          return "Debug";
        case PrintLevel::verbose:
          return "Verbose";
        case PrintLevel::info:
          return "Info";
        case PrintLevel::warning:
          return "Warning";
        case PrintLevel::error:
          return "Error";
        case PrintLevel::off:
          return "Off";
      }
      return "Unknown";
    }
  }

  std::ostream& operator<<(std::ostream& os, PrintLevel lvl)
  {
    return os << Intern::get_print_char(lvl);
  }

  class Logger
  {
  public:
    const FEAT::Dist::Comm& comm;
    mutable FEAT::String log;
    PrintLevel std_out_lvl = PrintLevel::info;
    PrintLevel log_lvl = PrintLevel::verbose;

    bool print_time, always_flush;


    Logger(const FEAT::Dist::Comm& comm_, PrintLevel out_lvl_ = PrintLevel::info, PrintLevel log_lvl_ = PrintLevel::verbose, bool print_time_ = true, bool always_flush_ = false) :
      comm(comm_),
      std_out_lvl(out_lvl_),
      log_lvl(log_lvl_),
      print_time(print_time_),
      always_flush(always_flush_)
    {}

    ~Logger() = default;

    Logger(const Logger&) = delete;
    Logger(Logger&&) = delete;
    Logger& operator=(const Logger&) = delete;
    Logger& operator=(Logger&&) = delete;

    FEAT::String get_utc_time() const
    {
      std::time_t timer = (std::time(nullptr));
      char timeString[std::size("yyyy-mm-ddThh:mm:ssZ")];
      std::strftime(std::data(timeString), std::size(timeString),
                  "%FT%TZ", std::localtime(&timer));
      return {timeString};

    }

    void print(const FEAT::String& msg, PrintLevel lvl = PrintLevel::info) const
    {
      FEAT::String pre_line = FEAT::String(Intern::get_print_char(lvl)).pad_back(6, '.');
      if(print_time)
      {
        pre_line += "[" + get_utc_time() + "]";
      }

      pre_line += ": ";

      if(lvl >= std_out_lvl)
      {
        for(const auto& line : msg.split_by_charset("\n"))
        {
          comm.print(pre_line + line);
        }
      }
      if(lvl >= log_lvl)
      {
        for(const auto& line : msg.split_by_charset("\n"))
        {
          log += "\n" + pre_line + line;
        }
      }

      if(always_flush)
        this->flush_print();
    }

    void flush_print() const
    {
      comm.print_flush();
      // todo: write out log to file...
    }

  };
}