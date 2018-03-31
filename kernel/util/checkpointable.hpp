// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_CHECKPOINTABLE_HPP
#define KERNEL_UTIL_CHECKPOINTABLE_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>


namespace FEAT
{
  /**
   * \brief Checkpoint interface class.
   *
   * This abstract class provides the necessary interface for all checkpointable objects.
   *
   */
  class Checkpointable
  {
    public:
      /**
       * \brief Calculate size
       *
       * Calculate size of complete object as it is stored in the checkpoint
       *
       * \return size of object
       *
       */
      virtual uint64_t get_checkpoint_size() = 0;

      /**
       * \brief Extract object from checkpoint
       *
       * Restores complete object with all its contents from checkpoint data
       *
       * \param[out] data object as bytestrem
       *
       */
      virtual void restore_from_checkpoint_data(std::vector<char> & data) = 0;

      /**
       * \brief Collect checkpoint data from object
       *
       * Adds the condensed complete object with all its contents to the end of the checkpoint buffer
       *
       * \param[out] data object as bytestream
       *
       */
      virtual void set_checkpoint_data(std::vector<char>& data) = 0;

    protected:
      virtual ~Checkpointable(){}
  }; // class Checkpointable

} // namespace FEAT
#endif // KERNEL_UTIL_CHECKPOINTABLE_HPP
