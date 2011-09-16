/* GENERAL_REMARK_BY_HILMAR:
 * If you decide to kick out the master process, you can probably also remove this class.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_SERVICE_IDS_HPP
#define KERNEL_SERVICE_IDS_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>

// includes, system

namespace FEAST
{
  /**
  * \brief class defining an enumeration of service IDs for the master service loop
  *
  * \author Hilmar Wobker
  */
  class ServiceIDs
  {

  public:

    /// enumeration of service IDs for the master service loop
    enum service_id
    {
      /// receive single log message (see class Logger)
      LOG_RECEIVE,
      /// receive array of log messages (see class Logger)
      LOG_RECEIVE_ARRAY,
      /// finish the master's service loop
      MASTER_FINISH_SERVICE
    };

  }; // class ServiceIDs
} // namespace FEAST

#endif // KERNEL_SERVICE_IDS_HPP
