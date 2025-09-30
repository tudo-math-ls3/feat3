#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/dist_file_io.hpp>

#include <ctime>

namespace Gendie
{
  enum print_level
  {
    debug = 0,
    verbose = 1,
    info = 2,
    warning = 3,
    error = 4,
    off = 5
  };

  namespace Intern
  {
    inline const char* get_print_char(print_level lvl)
    {
      switch(lvl)
      {
        case print_level::debug:
          return "Debug";
        case print_level::verbose:
          return "Verbose";
        case print_level::info:
          return "Info";
        case print_level::warning:
          return "Warning";
        case print_level::error:
          return "Error";
        case print_level::off:
          return "Off";
      }
      return "Unknown";
    }
  }

  class Logger
  {
  public:
    const FEAT::Dist::Comm& comm;
    mutable FEAT::String log;
    print_level std_out_lvl = print_level::info;
    print_level log_lvl = print_level::verbose;

    bool print_time;


    Logger(const FEAT::Dist::Comm& comm_, print_level out_lvl_ = print_level::info, print_level log_lvl_ = print_level::verbose, bool print_time_ = true) :
      comm(comm_),
      std_out_lvl(out_lvl_),
      log_lvl(log_lvl_),

      print_time(print_time_)
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

    void print(const FEAT::String& msg, print_level lvl = print_level::info) const
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
    }

    void flush_print() const
    {
      comm.print_flush();
      // todo: write out log to file...
    }

  };
}