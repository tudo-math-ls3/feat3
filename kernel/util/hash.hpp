// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <cstddef>
#include <cstdint>

namespace FEAT
{
  class Hash
  {
  public:
    /**
     * \brief Computes the CRC-32 checksum of a buffer
     *
     * The polynomial for the CRC table is 0x04C11DB7.
     *
     * \param[in] n
     * The length of the input buffer in bytes.
     *
     * \param[in] buffer
     * A \transient pointer to the buffer that the checksum is to be computed for
     *
     * \returns
     * The CRC-32 checksum of the input buffer.
     */
    static std::uint32_t crc32(std::size_t n, const void* buffer);

    /**
     * \brief Updates a CRC-32 checksum of a buffer
     *
     * The polynomial for the CRC table is 0x04C11DB7.
     *
     * \param[in] n
     * The length of the input buffer in bytes.
     *
     * \param[in] buffer
     * A \transient pointer to the buffer that the checksum is to be updated for
     *
     * \returns
     * The updated raw CRC-32 checksum of the input buffer.
     */
    static std::uint32_t crc32_update(std::size_t n, const void* buffer, std::uint32_t crc_in);
  };
} // namespace FEAT
