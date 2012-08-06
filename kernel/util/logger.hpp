#pragma once
#ifndef KERNEL_UTIL_LOGGER_HPP
#define KERNEL_UTIL_LOGGER_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/abort.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <iostream>
#include <fstream>

namespace FEAST
{
  /**
   * \brief Logger class implementation
   *
   * \todo detailed documentation
   * \author Peter Zajac
   */
  class Logger
  {
  public:
    /**
     * \brief Logger channel enumeration
     *
     * This enumeration describes the channels to which a message can be logged. Different values in this enumeration
     * can be combined by using the bit-wise \c OR ( \c | ) and \c AND ( \c & ) operators.
     *
     * \internal \b Internal: \n
     * Each channel encodes a 32 bit mask, where the low-order 16 bits are reserved for local channels, whereas
     * the high-order 16 bits represent channels on the master process.\n
     * Each master channel (e.g. \c master_file_0), can be translated into the corresponding local channel (e.g.
     * \c local_file_0) by simply performing a right bit shift by 16 bits and vice versa.
     */
    enum Channel
    {
      /**
       * \brief Null channel.
       */
      none                    = 0x00000000,

      /**
       * \brief Default local channel.
       *
       * \li In a parallel build, this channel is equivalent to \link local_file_0\endlink.
       * \li In a serial build, this channel is an OR-ed combination of \link local_file_0\endlink and
       *     \link local_standard\endlink.
       */
#ifdef PARALLEL
      local                   = 0x00000001, // = local_file_0
#else
      local                   = 0x00008001, // = local_file_0 | local_standard
#endif

      /**
       * \brief Default master channel.
       * \li In a parallel build, this channel is an OR-ed combination of \link master_file_0\endlink and
       *     \link master_standard\endlink.
       * \li In a serial build, this channel is (currently) equivalent to \link none\endlink.
       */
#ifdef PARALLEL
      master                  = 0x80010000, // = master_file_0 | master_standard
#else
      master                  = 0,          // = none; maybe change that...
#endif

      /**
       * \brief Local standard output channel.
       *
       * This channel represents the standard output stream (i.e. <c>std::cout</c>) of the calling process.
       */
      local_standard          = 0x00008000,

      /**
       * \brief Local log file channels.
       *
       * Each log file channel represents a separate log file, which can be targeted for output by setting
       * the \p channels parameter of Logger::log.\n
       * Local log files are created <em>per process</em>, i.e. each process in a parallel simulation has
       * its own set of local log files.
       */
      local_file_0            = 0x00000001,
      /// see \link local_file_0\endlink
      local_file_1            = 0x00000002,
      /// see \link local_file_0\endlink
      local_file_2            = 0x00000004,
      /// see \link local_file_0\endlink
      local_file_3            = 0x00000008,

      /**
       * \brief Master standard output channel.
       *
       * This channel represents the standard output stream  (i.e. <c>std::cout</c>) of the master process.
       */
      master_standard         = 0x80000000,

      /**
       * \brief Master log file channels.
       */
      master_file_0           = 0x00010000,
      /// see \link master_file_0\endlink
      master_file_1           = 0x00020000,
      /// see \link master_file_0\endlink
      master_file_2           = 0x00040000,
      /// see \link master_file_0\endlink
      master_file_3           = 0x00080000,

      /**
       * \brief Screen output channel.
       *
       * \li In a parallel build, this channel refers to the master standard output stream.\n
       * \li In a serial build, this channel refers to the local standard output stream.
       */
#ifdef PARALLEL
      screen                  = master_standard,
#else
      screen                  = local_standard,
#endif
      /// \cond internal
      /** bit mask for local channels */
      local_mask              = 0x0000FFFF,
      /** bit mask for master channels */
      master_mask             = 0xFFFF0000
      /// \endcond
    }; // enum Channel

    enum
    {
      /**
       * \brief Maximum number of log files.
       */
      max_files = 4
    };

#ifdef PARALLEL
    /**
     * \brief Master sender interface
     *
     * This abstract class acts as an interface for the master sender in a parallel build, which implements
     * the functionality of sending log messages to the master process.
     *
     * \author Peter Zajac
     */
    class MasterSender
    {
    public:
      /**
       * \brief virtual destructor
       */
      virtual ~MasterSender()
      {
      }

      /**
       * \brief Sends a message to the master.
       *
       * \param[in] message
       * The message to be sent.
       *
       * \param[in] channels
       * The channels to be logged.
       */
      virtual void send(
        const String& message,
        Channel channels) = 0;
    }; // class MasterSender
#endif // PARALLEL

  private:

    /// local log file streams
    static std::ofstream _stream[max_files];

#ifdef PARALLEL
    /// master sender object pointer
    static MasterSender* _master_sender;
#endif // PARALLEL

  public:
    /**
     * \brief Opens a local log file.
     *
     * This function opens a local log file, which can then be targeted by the <c>log_channel_<index></c> channel
     * for logging.
     *
     * \param[in] file_name
     * The file name of the log file.
     *
     * \param[in] index
     * The index of the log file that is to be opened. Must be in range {0, ..., Logger::max_files - 1}.
     */
    static void open(
      const String& file_name,
      int index = 0)
    {
      CONTEXT("Logger::open()");

      // check index range
      ASSERT((index >= 0) && (index < max_files), "Index out of range");

      // ensure that the specified channel isn't already open
      if(_stream[index].is_open())
      {
        abort("Log file stream is already open!");
      }

      // try to open the desired log file stream
      _stream[index].open(file_name.c_str());
      if(_stream[index].fail())
      {
        abort("Failed to open log file '" + file_name + "'");
      }
    }

    /**
     * \brief Closes a local log file.
     *
     * \param[in] index
     * The index of the local log file to be closed.
     */
    static void close(int index)
    {
      CONTEXT("Logger::close()");
      ASSERT((index >= 0) && (index < max_files), "Index out of range");
      if(!_stream[index].is_open())
      {
        abort("Log file is not open (anymore)");
      }
      _stream[index].close();
    }

    /**
     * \brief Closes all currently open local log files.
     */
    static void close_all()
    {
      CONTEXT("Logger::close_all()");

      // loop over all streams and close the open ones.
      for(int i = 0; i < max_files; ++i)
      {
        if(_stream[i].is_open())
        {
          _stream[i].close();
        }
      }
    }

#ifdef PARALLEL
    /**
     * \brief Sets the MasterSender object.
     *
     * \param[in] master_sender
     * A pointer to a class implementing the MasterSender interface.
     */
    static void set_master_sender(MasterSender* master_sender)
    {
      _master_sender = master_sender;
    }
#endif // PARALLEL

    /**
     * \brief Returns a log file channel.
     *
     * This is an auxiliary function which returns the channel for a local or master log file.
     *
     * \param[in] index
     * The index of the log file channel to be returned.
     *
     * \param[in] master
     * If \c true, then the returned channel will be a master channel, otherwise a local channel.
     *
     * \returns
     * A Logger::Channel for the desired log file.
     *
     * Examples:
     * \li <c>file_channel(0,false)</c> returns Logger::local_file_0
     * \li <c>file_channel(2,true)</c> returns Logger::master_file_2
     */
    static Channel file_channel(
      int index,
      bool master = false)
    {
      CONTEXT("Logger::file_channel()");
      ASSERT((index >= 0) && (index < max_files), "Index out of range");

      // encode channel
      return (Channel)(1 << (index + (master ? 16 : 0)));
    }

    /**
     * \brief Builds a filename.
     *
     * This auxiliary function can be used to generate unique file names for log files in a parallel
     * build.
     *
     * \param[in] base_name
     * The base name of the file name.
     *
     * \param[in] rank
     * The rank of the calling process.
     *
     * \param[in] num_ranks
     * The total number of ranks.
     *
     * \returns
     * The file name of the form <c>\<base-name\>_\<rank\>.log</c> as a String.
     */
    static String build_name(
      const String& base_name,
      unsigned int rank = 0,
      unsigned int num_ranks = 0)
    {
      String file_name(base_name);

      // count number of decimal digits
      unsigned int num_digits(0);
      unsigned int i = num_ranks;
      while(i != 0)
      {
        ++num_digits;
        i /= 10u;
      }
      num_digits = std::max(3u, num_digits);
      // print rank suffix
      char* tmp = new char[num_digits+3u];
      sprintf(tmp, "_%0*i", num_digits, rank);
      file_name.append(tmp);
      delete [] tmp;

      file_name.append(".log");
      return file_name;
    }

    /**
     * \brief Logs a message.
     *
     * This function writes a message string to the specified output channels.
     *
     * \param[in] message
     * The string message that is to be logged.
     *
     * \param[in] channels
     * An OR-ed combination of log channels to which the message should be written.
     */
    static void log(
      const String& message,
      Channel channels = local)
    {
      CONTEXT("Logger::log()");

      // trivial call?
      if(message.empty() || (channels == none))
        return;

      // log to local stdout?
      if((channels & local_standard) != none)
      {
        std::cout << message;
      }

      // log to local files?
      for(int i(0); i < max_files; ++i)
      {
        if((channels & file_channel(i, false)) != none)
        {
          if(_stream[i].is_open())
          {
            _stream[i] << message;
          }
          // Note: The lack of an else-block spitting out error messages about closed streams is intended.
        }
      }

      /// \todo implement sending to master
#ifdef PARALLEL
      // check whether we need to send something to the master, otherwise return here
      if((channels & master_mask) == none)
        return;

      if(_master_sender != nullptr)
        _master_sender->send(message, channels);
#endif // PARALLEL
    }

    /**
     * \brief Flushes all local channels.
     *
     * This function flushes the local standard output stream as well as all currently open local log file streams.
     *
     * \todo Maybe this function should also send a message to the master telling him to flush his channels, too?
     */
    static void flush()
    {
      CONTEXT("Logger::flush()");

      // flush local stdout
      std::cout.flush();

      // flush local files
      for(int i(0); i < max_files; ++i)
      {
        if(_stream[i].is_open())
        {
          _stream[i].flush();
        }
      }
    }
  }; // class Logger

  /// \cond nodoxy
  inline Logger::Channel operator|(Logger::Channel a, Logger::Channel b)
  {
    // The following type-casts should be safe, as the C++ standard describes in [dcl.enum] that an enumeration
    // has to use the int type as an underlying type.
    return (Logger::Channel)(((int)a) | ((int)b));
  }

  inline Logger::Channel operator&(Logger::Channel a, Logger::Channel b)
  {
    return (Logger::Channel)(((int)a) & ((int)b));
  }
  /// \endcond
} // namespace FEAST

#endif // KERNEL_UTIL_LOGGER_HPP
